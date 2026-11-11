// Checks on the Nektar++ backend that the Newton solver relies on.
//
// Usage: test_flowmap mesh.xml session.xml   (ECS_Steps from the session)
//
//  1. Flow map is deterministic (repeat calls agree to round-off).
//  2. Inner product is symmetric and positive.
//  3. Finite-difference Jacobian is consistent: D(eps) = (Phi(u+eps v) -
//     Phi(u)) / eps converges at first order in eps until round-off, so the
//     default eps (~1e-8 relative) gives J v to ~1e-6 or better.
//  4. Projection: set_state is idempotent on states produced by the solver,
//     so Krylov vectors (differences of solver outputs) are not distorted.
#include <nekecs/nektar/flowmap.hpp>
#include <nekecs/space.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace nekecs;
using namespace nekecs::nektar;

static int failures = 0;
#define CHECK(cond, msg)                                                       \
    do                                                                         \
    {                                                                          \
        std::cout << ((cond) ? "ok   " : "FAIL ") << msg << std::endl;         \
        if (!(cond))                                                           \
            ++failures;                                                        \
    } while (0)

int main(int argc, char *argv[])
{
    auto session = Nektar::LibUtilities::SessionReader::CreateInstance(argc, argv);
    auto graph   = Nektar::SpatialDomains::MeshGraphIO::Read(session);
    int nsteps;
    session->LoadParameter("ECS_Steps", nsteps, 20);

    VCSFlowMap fm(session, graph);
    const NektarSpace S = fm.space();
    const FieldVector u = fm.state();
    const double un     = norm(S, u);
    std::cout << std::scientific << std::setprecision(3);

    // 1. determinism
    FieldVector a = S.zeros_like(u), b = S.zeros_like(u);
    fm.advance(u, nsteps, a);
    fm.advance(u, nsteps, b);
    FieldVector d = clone(S, a);
    S.axpy(-1.0, b, d);
    CHECK(norm(S, d) <= 1e-13 * un,
          "flow map deterministic: |Phi(u)-Phi(u)'|/|u| = " << norm(S, d) / un);

    // 2. inner product
    FieldVector r = clone(S, a);
    S.axpy(-1.0, u, r); // residual direction Phi(u) - u
    const double rn = norm(S, r);
    CHECK(std::abs(S.dot(u, r) - S.dot(r, u)) <= 1e-14 * un * rn,
          "dot symmetric");
    CHECK(S.dot(u, u) > 0.0 && rn > 0.0, "dot positive (|u| = " << un
                                             << ", |Phi(u)-u| = " << rn << ")");

    // 3. FD Jacobian convergence along v = r / |r|
    FieldVector v = clone(S, r);
    S.scale(1.0 / rn, v);
    FieldVector prev;
    double prev_diff = -1.0, first_ratio = 0.0;
    for (int k = 2; k <= 8; ++k)
    {
        const double eps = un * std::pow(10.0, -k);
        FieldVector up = clone(S, u);
        S.axpy(eps, v, up);
        FieldVector D = S.zeros_like(u);
        fm.advance(up, nsteps, D);
        S.axpy(-1.0, a, D);
        S.scale(1.0 / eps, D);
        if (!prev.u.empty())
        {
            FieldVector dd = clone(S, D);
            S.axpy(-1.0, prev, dd);
            const double diff = norm(S, dd) / norm(S, D);
            std::cout << "     eps/|u| = 1e-" << k
                      << "  |D(eps) - D(10 eps)| / |D| = " << diff
                      << "  |Jv| = " << norm(S, D) << std::endl;
            if (prev_diff > 0.0 && first_ratio == 0.0)
            {
                first_ratio = prev_diff / diff;
            }
            prev_diff = diff;
        }
        prev = D;
    }
    CHECK(first_ratio > 5.0,
          "FD Jacobian converges at first order (ratio per decade = "
              << first_ratio << ")");
    CHECK(prev_diff < 1e-5,
          "FD Jacobian at eps ~1e-8 accurate to " << prev_diff);

    // 4. projection idempotent on solver output
    FieldVector a2 = S.zeros_like(u);
    fm.advance(a, 0, a2); // zero steps: set_state + read back
    FieldVector pd = clone(S, a2);
    S.axpy(-1.0, a, pd);
    CHECK(norm(S, pd) <= 1e-12 * un,
          "state projection idempotent (" << norm(S, pd) / un << ")");

    std::cout << (failures ? "test_flowmap: FAILED" : "test_flowmap: all checks passed")
              << " (" << fm.calls() << " flow-map calls)" << std::endl;
    session->Finalise();
    return failures ? 1 : 0;
}
