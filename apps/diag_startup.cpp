// Start-up diagnostics for the Nektar++ flow map.
//
// Usage: diag_startup mesh.xml session.xml   (IC: ideally a converged EQ)
//
//  1. our time loop with StartMode::Standard reproduces DoSolve();
//  2. each start mode is deterministic;
//  3. invariance of the IC under k = 1, 2, 4, 8 multiples of T and the
//     semigroup defect |Phi_2T(u) - Phi_T(Phi_T(u))|, per start mode.
// For a steady state u of the regular scheme, ConstantHistory gives
// Phi_kT(u) = u for every k; Standard does not (first step is order 1).
#include <nekecs/nektar/flowmap.hpp>
#include <nekecs/space.hpp>

#include <iomanip>
#include <iostream>

using namespace nekecs;
using namespace nekecs::nektar;

int main(int argc, char *argv[])
{
    auto session = Nektar::LibUtilities::SessionReader::CreateInstance(argc, argv);
    auto graph   = Nektar::SpatialDomains::MeshGraphIO::Read(session);
    int N;
    session->LoadParameter("ECS_Steps", N, 100);

    VCSFlowMap fm(session, graph);
    const NektarSpace S = fm.space();
    const FieldVector u = fm.state();
    const double un     = norm(S, u);
    std::cout << std::scientific << std::setprecision(3);
    auto diff = [&](const FieldVector &a, const FieldVector &b) {
        FieldVector d = clone(S, a);
        S.axpy(-1.0, b, d);
        return norm(S, d) / un;
    };
    FieldVector a = S.zeros_like(u), b = S.zeros_like(u), c = S.zeros_like(u);

    // 1. loop vs DoSolve
    fm.set_use_dosolve(true);
    fm.advance(u, N, a);
    fm.set_use_dosolve(false);
    fm.set_start_mode(StartMode::Standard);
    fm.advance(u, N, b);
    std::cout << "[diag] loop(Standard) vs DoSolve: " << diff(a, b) << std::endl;

    for (auto mode : {StartMode::Standard, StartMode::ConstantHistory})
    {
        const char *name = mode == StartMode::Standard ? "Standard       "
                                                       : "ConstantHistory";
        fm.set_start_mode(mode);
        // 2. determinism
        fm.advance(u, N, a);
        fm.advance(u, N, b);
        std::cout << "[diag] " << name << "  repeat-call diff = " << diff(a, b)
                  << std::endl;
        // 3. invariance and semigroup defect
        for (int k : {1, 2, 4, 8})
        {
            fm.advance(u, k * N, c);
            std::cout << "[diag] " << name << "  |Phi_" << k
                      << "T(u) - u|/|u| = " << diff(c, u) << std::endl;
        }
        fm.advance(u, N, a);
        fm.advance(a, N, b);
        fm.advance(u, 2 * N, c);
        std::cout << "[diag] " << name << "  semigroup defect = " << diff(c, b)
                  << std::endl;
    }
    // 4. Two-level (History) state: exact semigroup for an UNSTEADY state.
    //    x0 = (w, w) with w = u + noise, so the flow is far from steady.
    {
        FieldVector w     = fm.random_field(3);
        S.scale(0.05 * un / norm(S, w), w);
        S.axpy(1.0, u, w);
        FieldVector x0 = fm.with_history(w);
        const NektarSpace S2 = fm.space();
        auto d2 = [&](const FieldVector &p, const FieldVector &q) {
            FieldVector d = clone(S2, p);
            S2.axpy(-1.0, q, d);
            return norm(S2, d) / norm(S2, q);
        };
        FieldVector A = S2.zeros_like(x0), B = S2.zeros_like(x0),
                    C = S2.zeros_like(x0);
        fm.advance(x0, N, A);     // N steps in one call
        fm.advance(x0, N / 2, B); // N/2 + N/2 in two calls
        fm.advance(B, N - N / 2, C);
        std::cout << "[diag] History   |Phi_T(x) - Phi_T/2(Phi_T/2(x))| / |.| = "
                  << d2(C, A) << "  (unsteady: |Phi_T(x) - x|/|x| = "
                  << d2(A, x0) << ")" << std::endl;
        fm.advance(x0, N, B);
        std::cout << "[diag] History   repeat-call diff = " << d2(B, A) << std::endl;

        // 5. dt change and back leaves no trace.
        const double dt = fm.timestep();
        fm.set_timestep(1.02 * dt);
        fm.advance(x0, N, B);
        std::cout << "[diag] History   dt x 1.02: |Phi - Phi(dt)|/|.| = " << d2(B, A)
                  << " (should be O(1e-2) different)" << std::endl;
        fm.set_timestep(dt);
        fm.advance(x0, N, B);
        std::cout << "[diag] History   dt restored: repeat diff = " << d2(B, A)
                  << std::endl;

        // 6. ConstantHistory start == History start with (u, u).
        fm.set_start_mode(StartMode::ConstantHistory);
        FieldVector c1 = S.zeros_like(u);
        fm.advance(w, N, c1);
        std::cout << "[diag] ConstantHistory(w) vs History(w, w): "
                  << diff(c1, fm.level0(A)) << std::endl;
        for (int n : {0, 1, 2})
        {
            fm.advance(w, n, c1);
            fm.advance(x0, n, B);
            std::cout << "[diag]   after " << n << " steps: CH vs History level0 = "
                      << diff(c1, fm.level0(B)) << std::endl;
        }
        // Semigroup at the step level: 1 + 1 vs 2.
        fm.advance(x0, 2, A);
        fm.advance(x0, 1, B);
        fm.advance(B, 1, C);
        std::cout << "[diag]   History 1+1 vs 2 steps: " << d2(C, A) << std::endl;
    }

    session->Finalise();
    return 0;
}
