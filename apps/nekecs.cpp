// nekecs: exact coherent structures of the incompressible Navier-Stokes
// equations with Newton-Krylov-hookstep around the Nektar++ flow map.
//
//   equilibrium (ECS_Period <= 0):  F(u) = Phi_T(u) - u,  T = ECS_Steps * TimeStep
//   periodic orbit (ECS_Period > 0):
//       x = (u_n, u_{n-1}, T)   (two-level state: exact IMEX-2 history)
//       F(x) = (Phi_T(u_n, u_{n-1}) - (u_n, u_{n-1}),  <b, u_n - u_ref>)
//       with dt = T / ECS_Steps, b the unit flow direction at u_ref
//       (Poincare hyperplane; moved to each new point in continuation).
//
// With SOLVERINFO ECS_ContParam = <name>, the solution is continued in that
// parameter by pseudo-arclength continuation (core/continuation.hpp).
//
// Usage: nekecs mesh.xml session.xml
//
// Session PARAMETERS (defaults in brackets):
//   ECS_Period     initial period guess; > 0 selects a periodic orbit [0]
//   ECS_Steps      time steps per flow-map call (per period for POs) [100]
//   ECS_Tol        stop when ||F|| / ||x|| <= ECS_Tol       [1e-10]
//   ECS_MaxNewton  Newton iterations                        [30]
//   ECS_MaxKrylov  GMRES dimension per Newton step          [60]
//   ECS_GmresTol   GMRES relative tolerance                 [1e-3]
//   ECS_Delta0     initial trust radius / ||x0|| (<=0: none) [-1]
//   ECS_Verify     extra multiples of T for the final check [4]
//   ECS_FDEps      finite-difference step for J v, relative  [1e-7]
//   ECS_Start      EQ integrator start per call (0 standard, 1 constant
//                  history)                                  [1]
//   ECS_ProbeX, ECS_ProbeY  probe point for diagnostics      [2, 0]
//
// Linear stability / Floquet analysis after the solve (ECS_MaxNewton = 0 to
// only analyse the initial condition):
//   ECS_Eigs       leading eigenvalues to compute            [0 = off]
//   ECS_EigMaxDim  Arnoldi Krylov dimension limit            [80]
//   ECS_EigTol     relative Ritz residual                    [1e-6]
//   ECS_EigSteps   EQ: time steps per Jacobian application   [ECS_Steps]
//   ECS_EigStart   EQ: 0 standard start-up, 1 constant history,
//                  2 two-level history state (exact, no start-up bias) [2]
// POs always use the two-level state and one period. Reports multipliers mu
// and lambda = log(mu)/T (frequency modulo 2 pi / T); writes
// <session>_eig<i>_{re,im}.fld (level u_n).
//
// Continuation (SOLVERINFO ECS_ContParam: "Re", "Kinvis", or any session
// parameter that appears directly in boundary-condition expressions, e.g. an
// inflow angle; name as written in the session file):
//   ECS_ContEnd    stop when the parameter passes this value       [required]
//   ECS_ContDmu    first (natural-parameter) step, signed           [required]
//   ECS_ContDs     arclength step, in parameter units               [|ECS_ContDmu|]
//   ECS_ContDsMax  largest arclength step, parameter units          [4 ECS_ContDs]
//   ECS_ContSteps  maximum arclength steps                          [50]
//   ECS_ContWeight relative weight of the parameter in the arclength [1]
//                  (1: a 1% change of the parameter counts like a 1%
//                  change of ||x||)
// With ECS_Eigs > 0 the leading eigenvalue (EQ) or the largest non-trivial
// Floquet multiplier (PO) is computed at every point and crossings of the
// stability boundary are reported (EQ: Hopf / steady bifurcation; PO:
// |mu| = 1). Output: <session>_branch.dat (one line per point),
// <session>_c<NNN>.fld (state at each point, level u_n).
//
// Run with OPENBLAS_NUM_THREADS=1: threaded OpenBLAS on the small elemental
// matrices is ~8x slower (thread contention). Parallelise with MPI instead.
//
// Output (<session> = basename of the last XML file on the command line):
// <session>_eq.fld or <session>_po.fld, <session>_ecs_latest.fld (last
// accepted iterate), log on stdout.
#include <nekecs/arnoldi.hpp>
#include <nekecs/continuation.hpp>
#include <nekecs/nektar/flowmap.hpp>
#include <nekecs/newton.hpp>

#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace nekecs;
using namespace nekecs::nektar;

namespace
{

struct Params
{
    int nsteps, max_newton, max_krylov, verify, neig, eig_dim, start,
        eig_steps, eig_start, cont_steps;
    NekDouble tol, gmres_tol, delta0, fd_eps, eig_tol, period, probe_x,
        probe_y, cont_end, cont_dmu, cont_ds, cont_dsmax, cont_weight;
    std::string cont_param;
};

double rel_diff(const NektarSpace &S, const FieldVector &a, const FieldVector &b)
{
    FieldVector d = clone(S, a);
    S.axpy(-1.0, b, d);
    return norm(S, d) / norm(S, b);
}

struct Spectrum
{
    std::vector<std::complex<double>> mu;  // multipliers of D Phi_T
    std::vector<double> residual;
    double T = 0.0;
    bool converged = false;
};

// Leading eigenvalues of D Phi at `base` (any number of levels), Phi = nsteps
// steps at the current dt. Optionally writes eigenvectors (level 0).
Spectrum stability(VCSFlowMap &fm, const FieldVector &base, int nsteps,
                   const Params &P, const std::string &name, bool root,
                   bool write_vectors)
{
    const NektarSpace S = fm.space();
    const double T      = nsteps * fm.timestep();
    const double bnorm  = norm(S, base);

    FieldVector phi = S.zeros_like(base);
    fm.advance(base, nsteps, phi);
    auto J = [&](const FieldVector &v, FieldVector &Jv) {
        const double eps = P.fd_eps * (1.0 + bnorm) / norm(S, v);
        FieldVector up   = clone(S, base);
        S.axpy(eps, v, up);
        fm.advance(up, nsteps, Jv);
        S.axpy(-1.0, phi, Jv);
        S.scale(1.0 / eps, Jv);
    };
    // Admissible start vector: one Jacobian application to noise.
    FieldVector noise = fm.random_field(1);
    if (base.u.size() > noise.u.size())
    {
        noise = fm.with_history(noise);
    }
    noise.extra.clear();
    S.scale(1.0 / norm(S, noise), noise);
    FieldVector v0 = S.zeros_like(base);
    J(noise, v0);

    ArnoldiOptions aopt;
    aopt.nev     = P.neig;
    aopt.max_dim = P.eig_dim;
    aopt.tol     = P.eig_tol;
    aopt.log     = root ? &std::cout : nullptr;
    const int calls0 = fm.calls();
    const auto E     = arnoldi(S, J, v0, aopt);

    if (root)
    {
        std::cout << "[nekecs] stability: "
                  << (E.converged ? "converged" : "NOT converged")
                  << ", Krylov dim = " << E.k << ", flow-map calls = "
                  << fm.calls() - calls0 << ", T = " << T << ", levels = "
                  << base.u.size() / fm.nvel() << "\n"
                  << "[nekecs]   i            mu            |mu|    "
                     "lambda = log(mu)/T        residual"
                  << std::endl;
    }
    Spectrum out;
    out.T         = T;
    out.converged = E.converged;
    for (std::size_t i = 0; i < E.values.size(); ++i)
    {
        const auto mu  = E.values[i];
        const auto lam = std::log(mu) / T;
        out.mu.push_back(mu);
        out.residual.push_back(E.residuals[i]);
        if (root)
        {
            std::cout << std::scientific << std::setprecision(5)
                      << "[nekecs]   " << i << "  " << std::showpos << mu.real()
                      << " " << mu.imag() << "i  " << std::noshowpos
                      << std::abs(mu) << "  " << std::showpos << lam.real()
                      << " " << lam.imag() << "i  " << std::noshowpos
                      << std::setprecision(1) << E.residuals[i] << std::endl;
        }
        if (write_vectors)
        {
            fm.write_perturbation(fm.level0(E.re[i]),
                                  name + "_eig" + std::to_string(i) + "_re.fld");
            fm.write_perturbation(fm.level0(E.im[i]),
                                  name + "_eig" + std::to_string(i) + "_im.fld");
        }
    }
    return out;
}

NewtonOptions<NektarSpace> newton_options(const Params &P, double xnorm,
                                          bool root)
{
    NewtonOptions<NektarSpace> opt;
    opt.tol        = P.tol * xnorm;
    opt.max_newton = P.max_newton;
    opt.max_krylov = P.max_krylov;
    opt.gmres_rtol = P.gmres_tol;
    opt.fd_eps     = P.fd_eps;
    opt.delta0     = P.delta0 > 0.0 ? P.delta0 * xnorm : -1.0;
    opt.log        = root ? &std::cout : nullptr;
    return opt;
}

// One ECS problem (EQ or PO) on the flow map: everything the single solve
// and the continuation driver need.
struct Problem
{
    bool po = false;
    FieldVector x0;
    // F(x) at the flow map's current parameter values.
    std::function<void(const FieldVector &, FieldVector &)> F;
    // J v at x (F(x) given); null: finite differences on F.
    std::function<void(const FieldVector &, const FieldVector &,
                       const FieldVector &, FieldVector &)>
        jac;
    // Invalidate cached Jacobian data (call when a parameter changes).
    std::function<void()> invalidate = [] {};
    // Move the phase condition to the solution x (PO); no-op for EQ.
    std::function<void(const FieldVector &)> rebase = [](const FieldVector &) {};
    // |Phi_{verify T}(x) - x| / |x| at a converged x.
    std::function<double(const FieldVector &)> verify;
    // Spectrum at a converged x.
    std::function<Spectrum(const FieldVector &, bool write_vectors)> spectrum;
    double period(const FieldVector &x) const
    {
        return po ? x.extra[0] : 0.0;
    }
};

Problem make_equilibrium(VCSFlowMap &fm, const Params &P,
                         const std::string &name, bool root)
{
    const StartMode start =
        P.start == 0 ? StartMode::Standard : StartMode::ConstantHistory;
    fm.set_start_mode(start);
    const NektarSpace S = fm.space();
    Problem pb;
    pb.x0 = fm.state();
    pb.F  = [&fm, S, P](const FieldVector &x, FieldVector &Fx) {
        fm.advance(x, P.nsteps, Fx);
        S.axpy(-1.0, x, Fx);
    };
    pb.verify = [&fm, S, P](const FieldVector &x) {
        FieldVector y = S.zeros_like(x);
        fm.advance(x, P.verify * P.nsteps, y);
        return rel_diff(S, y, x);
    };
    pb.spectrum = [&fm, P, name, root, start](const FieldVector &x, bool wv) {
        if (P.eig_start == 2)
        {
            return stability(fm, fm.with_history(x), P.eig_steps, P, name,
                             root, wv);
        }
        fm.set_start_mode(P.eig_start == 0 ? StartMode::Standard
                                           : StartMode::ConstantHistory);
        auto s = stability(fm, x, P.eig_steps, P, name, root, wv);
        fm.set_start_mode(start);
        return s;
    };
    return pb;
}

// Periodic-orbit problem. Shared state (phase reference, Jacobian cache)
// lives in a heap struct captured by the closures.
Problem make_periodic_orbit(VCSFlowMap &fm, const Params &P,
                            const std::string &name, bool root)
{
    struct State
    {
        FieldVector uref, b;          // phase hyperplane
        FieldVector cache_x, dPhidT;  // period column of the Jacobian
        double cache_T = -1.0;
        bool cache_ok  = false;
    };
    auto st             = std::make_shared<State>();
    const NektarSpace S = fm.space();
    const int N         = P.nsteps;
    auto set_period     = [&fm, N](double T) { fm.set_timestep(T / N); };
    auto levels         = [S](const FieldVector &x) {
        FieldVector s = clone(S, x);
        s.extra.clear();
        return s;
    };
    auto flow = [&fm, N, set_period, levels](const FieldVector &x, double T,
                                             FieldVector &out) {
        set_period(T);
        fm.advance(levels(x), N, out);
    };

    Problem pb;
    pb.po       = true;
    pb.x0       = fm.with_history(fm.state());
    pb.x0.extra = {P.period};

    // Hyperplane through u_n of x, normal to d u_n / dt (one regular step).
    pb.rebase = [&fm, S, st, set_period, levels](const FieldVector &x) {
        set_period(x.extra[0]);
        const FieldVector s = levels(x);
        FieldVector y1      = S.zeros_like(s);
        fm.advance(s, 1, y1);
        st->uref = fm.level0(s);
        st->b    = fm.level0(y1);
        S.axpy(-1.0, st->uref, st->b);
        S.scale(1.0 / norm(S, st->b), st->b);
    };
    pb.rebase(pb.x0);
    auto phase = [&fm, S, st](const FieldVector &x) {
        FieldVector d = fm.level0(x);
        S.axpy(-1.0, st->uref, d);
        return S.dot(st->b, d);
    };

    pb.F = [S, flow, levels, phase](const FieldVector &x, FieldVector &Fx) {
        FieldVector s = levels(x);
        FieldVector y = S.zeros_like(s);
        flow(x, x.extra[0], y);
        S.axpy(-1.0, s, y);
        Fx.u     = y.u;
        Fx.extra = {phase(x)};
    };

    // Jacobian: FD in the state at fixed T (no time-step change, no new
    // Helmholtz factorisation per GMRES iteration) plus the period column
    // dPhi/dT, computed once per Newton iterate by FD in T.
    pb.invalidate = [st] { st->cache_ok = false; };
    pb.jac = [&fm, S, P, N, st, set_period, levels, flow](
                 const FieldVector &x, const FieldVector &Fx,
                 const FieldVector &v, FieldVector &Jv) {
        const double T = x.extra[0];
        FieldVector s  = levels(x);
        FieldVector phi = levels(Fx); // Phi_T(x) = F_state(x) + x
        S.axpy(1.0, s, phi);
        if (!st->cache_ok || st->cache_T != T || rel_diff(S, s, st->cache_x) != 0.0)
        {
            const double dT = P.fd_eps * T;
            FieldVector yT  = S.zeros_like(s);
            flow(x, T + dT, yT);
            S.axpy(-1.0, phi, yT);
            S.scale(1.0 / dT, yT);
            st->dPhidT   = yT;
            st->cache_x  = s;
            st->cache_T  = T;
            st->cache_ok = true;
        }
        FieldVector vs  = levels(v);
        const double vn = norm(S, vs);
        FieldVector Ju  = S.zeros_like(s);
        if (vn > 0.0)
        {
            const double eps = P.fd_eps * (1.0 + norm(S, s)) / vn;
            FieldVector up   = clone(S, s);
            S.axpy(eps, vs, up);
            FieldVector yp = S.zeros_like(s);
            set_period(T);
            fm.advance(up, N, yp);
            S.axpy(-1.0, phi, yp);
            S.scale(1.0 / eps, yp);
            Ju = yp;
            S.axpy(-1.0, vs, Ju);
        }
        S.axpy(v.extra[0], st->dPhidT, Ju);
        Jv.u     = Ju.u;
        Jv.extra = {S.dot(st->b, fm.level0(v))};
    };

    // Closure after `verify` periods in ONE call (exact two-level map).
    pb.verify = [&fm, S, P, N, set_period, levels](const FieldVector &x) {
        FieldVector s = levels(x), y = S.zeros_like(levels(x));
        set_period(x.extra[0]);
        fm.advance(s, P.verify * N, y);
        return rel_diff(S, y, s);
    };
    pb.spectrum = [&fm, P, N, name, root, set_period, levels](
                      const FieldVector &x, bool wv) {
        set_period(x.extra[0]);
        return stability(fm, levels(x), N, P, name, root, wv); // Floquet
    };
    return pb;
}

int run_single(VCSFlowMap &fm, Problem &pb, const Params &P,
               const std::string &name, bool root)
{
    const NektarSpace S = fm.space();
    const double xnorm  = norm(S, pb.x0);
    auto opt            = newton_options(P, xnorm, root);
    opt.on_accept       = [&](const FieldVector &x, double r, int it) {
        fm.write(fm.level0(x), name + "_ecs_latest.fld");
        if (root)
        {
            std::cout << "[nekecs] it " << it << "  |F|/|x| = "
                      << std::scientific << std::setprecision(3)
                      << r / norm(S, x);
            if (pb.po)
            {
                std::cout << "  T = " << std::setprecision(8) << x.extra[0];
            }
            std::cout << "  flow-map calls = " << fm.calls() << std::endl;
        }
        return false;
    };

    if (root)
    {
        FieldVector F0 = S.zeros_like(pb.x0);
        pb.F(pb.x0, F0);
        std::cout << std::scientific << std::setprecision(3) << "[nekecs] "
                  << (pb.po ? "periodic-orbit" : "equilibrium")
                  << " search: nvel = " << fm.nvel() << ", steps/call = "
                  << P.nsteps << ", dt = " << fm.timestep();
        if (pb.po)
        {
            std::cout << ", levels = 2, T0 = " << P.period;
        }
        else
        {
            std::cout << ", T = " << P.nsteps * fm.timestep() << ", start = "
                      << (P.start == 0 ? "standard" : "constant-history");
        }
        std::cout << ", |x0| = " << xnorm
                  << ", |F(x0)|/|x0| = " << norm(S, F0) / xnorm << std::endl;
    }

    std::function<void(const FieldVector &, const FieldVector &,
                       const FieldVector &, FieldVector &)>
        J = pb.jac;
    const auto R       = newton_hookstep(S, pb.F, pb.x0, opt, J);
    const double check = pb.verify(R.x);
    const std::string out = name + (pb.po ? "_po.fld" : "_eq.fld");
    fm.write(fm.level0(R.x), out);

    if (root)
    {
        std::cout << std::scientific << std::setprecision(3)
                  << "[nekecs] status = " << to_string(R.status)
                  << ", newton its = " << R.iterations
                  << ", rejected = " << R.rejected_steps
                  << ", flow-map calls = " << fm.calls() << "\n"
                  << "[nekecs] |F|/|x| = "
                  << R.residual_history.back() / norm(S, R.x) << ", |Phi_"
                  << P.verify << "T(x) - x|/|x| = " << check;
        if (pb.po)
        {
            const double T = R.x.extra[0];
            std::cout << "\n[nekecs] T = " << std::setprecision(10) << T
                      << ", St = f D / U = " << 1.0 / T
                      << ", dt = " << T / P.nsteps;
        }
        else
        {
            std::cout << ", |u - u0|/|u0| = " << rel_diff(S, R.x, pb.x0);
        }
        std::cout << "\n[nekecs] wrote " << out << std::endl;
        if (!R.converged() && P.neig > 0)
        {
            std::cout << "[nekecs] warning: stability of a non-converged state"
                      << std::endl;
        }
    }
    if (P.neig > 0)
    {
        pb.spectrum(R.x, true);
    }
    return R.converged() ? 0 : 1;
}

// Stability indicator of a point: EQ: largest Re(lambda), lambda =
// log(mu)/T; PO: largest |mu| - 1 over the non-trivial multipliers (the
// multiplier closest to 1 is the time-translation mode and is skipped).
struct Indicator
{
    double value = 0.0;          // > 0 unstable
    std::complex<double> mode;   // lambda (EQ) or mu (PO) of the leading mode
};

Indicator indicator(const Problem &pb, const Spectrum &s)
{
    Indicator ind;
    ind.value = -std::numeric_limits<double>::infinity();
    std::size_t trivial = s.mu.size();
    if (pb.po)
    {
        double best = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < s.mu.size(); ++i)
        {
            if (std::abs(s.mu[i] - 1.0) < best)
            {
                best    = std::abs(s.mu[i] - 1.0);
                trivial = i;
            }
        }
    }
    for (std::size_t i = 0; i < s.mu.size(); ++i)
    {
        if (i == trivial)
        {
            continue;
        }
        const double v = pb.po ? std::abs(s.mu[i]) - 1.0
                               : std::log(std::abs(s.mu[i])) / s.T;
        if (v > ind.value)
        {
            ind.value = v;
            ind.mode  = pb.po ? s.mu[i] : std::log(s.mu[i]) / s.T;
        }
    }
    return ind;
}

int run_continuation(VCSFlowMap &fm, Problem &pb, const Params &P,
                     const std::string &name, bool root)
{
    const NektarSpace S  = fm.space();
    const std::string pn = P.cont_param;
    const double mu0     = fm.parameter(pn);
    const double xnorm   = norm(S, pb.x0);

    ContinuationOptions<NektarSpace> co;
    co.newton = newton_options(P, xnorm, root);
    // Arclength metric: a relative change of the parameter weighs like the
    // same relative change of the state (times ECS_ContWeight).
    const double mscale = std::abs(mu0) > 0.0 ? std::abs(mu0) : 1.0;
    co.mu_weight        = P.cont_weight * (xnorm / mscale) * (xnorm / mscale);
    const double sw     = std::sqrt(co.mu_weight); // parameter units -> arclength
    co.ds               = P.cont_ds * sw;
    co.ds_max           = P.cont_dsmax * sw;
    co.ds_min           = 1e-3 * co.ds;
    co.max_steps        = P.cont_steps;
    co.fd_eps_mu        = P.fd_eps;
    const double range  = std::abs(P.cont_end - mu0);
    if (P.cont_dmu > 0.0)
    {
        co.mu_max = P.cont_end;
        co.mu_min = mu0 - range; // room to come back round a fold
    }
    else
    {
        co.mu_min = P.cont_end;
        co.mu_max = mu0 + range;
    }
    co.log       = root ? &std::cout : nullptr;
    co.on_newton = [&](const FieldVector &x, double mu, double r, int it) {
        fm.write(fm.level0(x), name + "_ecs_latest.fld");
        if (root)
        {
            std::cout << "[nekecs] it " << it << "  " << pn << " = "
                      << std::setprecision(8) << mu << "  |H|/|x| = "
                      << std::scientific << std::setprecision(3)
                      << r / norm(S, x) << "  flow-map calls = " << fm.calls()
                      << std::endl;
        }
    };

    double cur_mu = std::numeric_limits<double>::quiet_NaN();
    auto set_mu   = [&](double mu) {
        if (mu != cur_mu)
        {
            fm.set_parameter(pn, mu);
            pb.invalidate();
            cur_mu = mu;
        }
    };
    std::function<void(const FieldVector &, double, FieldVector &)> G =
        [&](const FieldVector &x, double mu, FieldVector &Gx) {
            set_mu(mu);
            pb.F(x, Gx);
        };
    std::function<void(const FieldVector &, double, const FieldVector &,
                       const FieldVector &, FieldVector &)>
        jx = nullptr;
    if (pb.jac)
    {
        jx = [&](const FieldVector &x, double mu, const FieldVector &Gx,
                 const FieldVector &v, FieldVector &Jv) {
            set_mu(mu);
            pb.jac(x, Gx, v, Jv);
        };
    }

    std::ofstream branch(name + "_branch.dat");
    if (root)
    {
        branch << "# nekecs continuation in " << pn << " ("
               << (pb.po ? "periodic orbit" : "equilibrium") << ")\n"
               << "# point " << pn << " arclength dmu/ds newton_its |F|/|x| "
               << "verify |u| u_probe v_probe T St"
               << " stab_indicator mode_re mode_im   (probe at ("
               << P.probe_x << ", " << P.probe_y << "); indicator: EQ max "
               << "Re(lambda), PO max |mu|-1 non-trivial; mode: EQ lambda, "
               << "PO mu)\n";
    }

    bool have_prev = false;
    Indicator prev_ind;
    double prev_mu = 0.0;
    std::vector<std::string> events;
    auto on_point = [&](const ContinuationPoint<NektarSpace> &p) {
        set_mu(p.mu);
        const double check = pb.verify(p.x);
        const double T     = pb.period(p.x);
        char fname[64];
        std::snprintf(fname, sizeof(fname), "_c%03d.fld", p.index);
        fm.write(fm.level0(p.x), name + fname);

        Indicator ind;
        bool have_ind = false;
        if (P.neig > 0)
        {
            const Spectrum s = pb.spectrum(p.x, false);
            ind              = indicator(pb, s);
            have_ind         = true;
            if (have_prev && (ind.value > 0.0) != (prev_ind.value > 0.0))
            {
                // Linear interpolation of the indicator between the points.
                const double a  = prev_ind.value, b = ind.value;
                const double mc = prev_mu + (p.mu - prev_mu) * a / (a - b);
                std::ostringstream os;
                os << std::setprecision(8) << "stability change between "
                   << pn << " = " << prev_mu << " and " << p.mu
                   << ", interpolated " << pn << "_c = " << mc;
                if (pb.po)
                {
                    os << ", multiplier " << ind.mode
                       << (std::abs(ind.mode.imag()) > 1e-6
                               ? " (Neimark-Sacker / torus)"
                               : (ind.mode.real() < 0 ? " (period doubling)"
                                                      : " (fold / pitchfork)"));
                }
                else
                {
                    os << ", lambda " << ind.mode
                       << (std::abs(ind.mode.imag()) > 1e-6
                               ? " (Hopf, omega = " +
                                     std::to_string(std::abs(ind.mode.imag())) + ")"
                               : " (steady: fold / pitchfork)");
                }
                events.push_back(os.str());
                if (root)
                {
                    std::cout << "[nekecs] *** " << os.str() << std::endl;
                }
            }
            prev_ind  = ind;
            have_prev = true;
        }
        prev_mu = p.mu;

        if (root)
        {
            const FieldVector u0 = fm.level0(p.x);
            branch << std::defaultfloat << std::setprecision(10) << p.index
                   << " " << p.mu << " "
                   << std::scientific << std::setprecision(6) << p.arclength / sw
                   << " " << p.dmu_ds * sw << " " << p.newton_its << " "
                   << std::setprecision(3) << p.residual / norm(S, p.x) << " "
                   << check << " " << std::setprecision(10) << norm(S, u0)
                   << " " << fm.probe(u0, 0, P.probe_x, P.probe_y) << " "
                   << fm.probe(u0, 1, P.probe_x, P.probe_y) << " " << T << " "
                   << (T > 0 ? 1.0 / T : 0.0);
            if (have_ind)
            {
                branch << " " << ind.value << " " << ind.mode.real() << " "
                       << ind.mode.imag();
            }
            branch << std::endl;
            std::cout << "[nekecs] point " << p.index << ": " << pn << " = "
                      << std::setprecision(8) << std::defaultfloat << p.mu
                      << std::scientific << std::setprecision(3)
                      << ", |F|/|x| = " << p.residual / norm(S, p.x)
                      << ", verify = " << check;
            if (pb.po)
            {
                std::cout << ", T = " << std::setprecision(8) << T;
            }
            if (have_ind)
            {
                std::cout << ", indicator = " << std::setprecision(4)
                          << ind.value;
            }
            std::cout << ", flow-map calls = " << fm.calls() << ", wrote "
                      << name << fname << std::endl;
        }
        pb.rebase(p.x); // next phase condition through this point
        return true;
    };

    if (root)
    {
        std::cout << "[nekecs] continuation in " << pn << " from " << mu0
                  << " towards " << P.cont_end << ": first step "
                  << P.cont_dmu << ", ds = " << P.cont_ds << " (max "
                  << P.cont_dsmax << ") in parameter units, weight w = "
                  << co.mu_weight << std::endl;
    }
    const auto R = continuation<NektarSpace>(S, G, pb.x0, mu0, P.cont_dmu, co,
                                             jx, on_point);
    if (root)
    {
        std::cout << "[nekecs] continuation " << to_string(R.status) << ": "
                  << R.mu.size() << " points, " << pn << " from "
                  << std::defaultfloat << std::setprecision(8) << R.mu.front()
                  << " to " << R.mu.back() << ", failed correctors = "
                  << R.failed_correctors << ", flow-map calls = " << fm.calls()
                  << "\n[nekecs] wrote " << name << "_branch.dat" << std::endl;
        for (const auto &e : events)
        {
            std::cout << "[nekecs] summary: " << e << std::endl;
        }
    }
    const bool ok = R.status == ContinuationStatus::ReachedBound ||
                    R.status == ContinuationStatus::MaxSteps;
    return ok ? 0 : 1;
}

} // namespace

int main(int argc, char *argv[])
{
    auto session = Nektar::LibUtilities::SessionReader::CreateInstance(argc, argv);
    auto graph   = Nektar::SpatialDomains::MeshGraphIO::Read(session);
    const bool root = session->GetComm()->GetRank() == 0;

    Params P;
    session->LoadParameter("ECS_Steps", P.nsteps, 100);
    session->LoadParameter("ECS_Tol", P.tol, 1e-10);
    session->LoadParameter("ECS_MaxNewton", P.max_newton, 30);
    session->LoadParameter("ECS_MaxKrylov", P.max_krylov, 60);
    session->LoadParameter("ECS_GmresTol", P.gmres_tol, 1e-3);
    session->LoadParameter("ECS_Delta0", P.delta0, -1.0);
    session->LoadParameter("ECS_Verify", P.verify, 4);
    session->LoadParameter("ECS_FDEps", P.fd_eps, 1e-7);
    session->LoadParameter("ECS_Eigs", P.neig, 0);
    session->LoadParameter("ECS_EigMaxDim", P.eig_dim, 80);
    session->LoadParameter("ECS_EigTol", P.eig_tol, 1e-6);
    session->LoadParameter("ECS_Start", P.start, 1);
    session->LoadParameter("ECS_EigSteps", P.eig_steps, P.nsteps);
    session->LoadParameter("ECS_EigStart", P.eig_start, 2);
    session->LoadParameter("ECS_Period", P.period, 0.0);
    session->LoadParameter("ECS_ProbeX", P.probe_x, 2.0);
    session->LoadParameter("ECS_ProbeY", P.probe_y, 0.0);

    P.cont_param = session->DefinesSolverInfo("ECS_ContParam")
                       ? session->GetSolverInfo("ECS_ContParam")
                       : "";
    if (!P.cont_param.empty())
    {
        if (!session->DefinesParameter("ECS_ContEnd") ||
            !session->DefinesParameter("ECS_ContDmu"))
        {
            throw std::runtime_error(
                "nekecs: continuation needs ECS_ContEnd and ECS_ContDmu");
        }
        session->LoadParameter("ECS_ContEnd", P.cont_end);
        session->LoadParameter("ECS_ContDmu", P.cont_dmu);
        session->LoadParameter("ECS_ContDs", P.cont_ds, std::abs(P.cont_dmu));
        session->LoadParameter("ECS_ContDsMax", P.cont_dsmax, 4.0 * P.cont_ds);
        session->LoadParameter("ECS_ContSteps", P.cont_steps, 50);
        session->LoadParameter("ECS_ContWeight", P.cont_weight, 1.0);
    }

    VCSFlowMap fm(session, graph);
    // Name outputs after the last XML file (the conditions file), not the
    // mesh: e.g. session.xml -> session_eq.fld.
    std::string name = session->GetFilenames().back();
    name = name.substr(name.find_last_of('/') + 1);
    name = name.substr(0, name.find('.'));

    Problem pb = P.period > 0.0 ? make_periodic_orbit(fm, P, name, root)
                                : make_equilibrium(fm, P, name, root);
    const int status = P.cont_param.empty()
                           ? run_single(fm, pb, P, name, root)
                           : run_continuation(fm, pb, P, name, root);
    session->Finalise();
    return status;
}
