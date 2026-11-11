// Pseudo-arclength continuation of solutions of G(x; mu) = 0 in a scalar
// parameter mu (Reynolds number, angle of attack, ...).
//
// G: X x R -> X is the square problem solved by newton_hookstep at fixed mu
// (equilibrium: Phi_T(u; mu) - u; periodic orbit: x = (u, T) with a phase
// condition). Continuation works on z = (x, mu) in the augmented space
// X x R with inner product <z1, z2> = <x1, x2> + w mu1 mu2 (w = mu_weight,
// so the arclength mixes state and parameter with a chosen scale).
//
// Steps
//   0. Newton at fixed mu0 from the initial guess.
//   1. Natural-parameter step: Newton at fixed mu0 + dmu0, predictor x0.
//   2. Repeat: secant tangent t = (z1 - z0) / |z1 - z0|, predictor
//      z = z1 + ds t, corrector Newton on the augmented system
//          H(z) = ( G(x; mu),  <t, z - z1> - ds ),
//      which passes through folds (mu is an unknown, the arclength is not).
//      The step ds adapts to the Newton iteration count; a failed corrector
//      halves ds and retries from z1.
//
// Jacobian of H: J v = (J_x v_x + v_mu dG/dmu, <t, v>). The column dG/dmu
// is computed by finite differences ONCE per Newton iterate, so every
// Krylov vector is applied at fixed mu (for Nektar++ this keeps the viscous
// Helmholtz matrices; a mu change refactorises them). J_x v is the user's
// operator if given, else finite differences at fixed mu.
#pragma once

#include "newton.hpp"
#include "space.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace nekecs
{

template <class Space> struct Augmented
{
    typename Space::Vector x;
    double mu = 0.0;
};

// X x R with <z1, z2> = <x1, x2> + w mu1 mu2.
template <class Space> class AugmentedSpace
{
public:
    using Vector = Augmented<Space>;

    AugmentedSpace(const Space &S, double w) : m_S(S), m_w(w)
    {
    }

    Vector zeros_like(const Vector &z) const
    {
        return Vector{m_S.zeros_like(z.x), 0.0};
    }
    void copy(const Vector &src, Vector &dst) const
    {
        m_S.copy(src.x, dst.x);
        dst.mu = src.mu;
    }
    double dot(const Vector &a, const Vector &b) const
    {
        return m_S.dot(a.x, b.x) + m_w * a.mu * b.mu;
    }
    void axpy(double a, const Vector &x, Vector &y) const
    {
        m_S.axpy(a, x.x, y.x);
        y.mu += a * x.mu;
    }
    void scale(double a, Vector &x) const
    {
        m_S.scale(a, x.x);
        x.mu *= a;
    }
    const Space &base() const
    {
        return m_S;
    }

private:
    Space m_S;
    double m_w;
};

template <class Space> struct ContinuationOptions
{
    using Vector = typename Space::Vector;

    // Newton settings for every solve (tol, Krylov dimension, trust region,
    // fd_eps for J_x). on_accept is ignored; use on_newton below.
    NewtonOptions<Space> newton;

    double mu_weight = 1.0;  // w in the augmented inner product
    double ds        = 1.0;  // initial arclength step (augmented norm)
    double ds_min    = 1e-6;
    double ds_max    = std::numeric_limits<double>::infinity();
    double grow      = 1.5;  // ds *= grow if the corrector took <= easy_its
    double shrink    = 0.5;  // ds *= shrink if it took >= hard_its, or failed
    int easy_its     = 2;
    int hard_its     = 5;
    int max_retries  = 6;    // corrector failures per step before giving up
    int max_steps    = 100;  // arclength steps after the two starting points
    double mu_min    = -std::numeric_limits<double>::infinity();
    double mu_max    = std::numeric_limits<double>::infinity();
    double fd_eps_mu = 1e-7; // dG/dmu step: h = fd_eps_mu (1 + |mu|)

    // After each accepted Newton step: (x, mu, |H|, iteration).
    std::function<void(const Vector &, double, double, int)> on_newton;

    std::ostream *log = nullptr;
};

template <class Space> struct ContinuationPoint
{
    typename Space::Vector x;
    double mu        = 0.0;
    int index        = 0;    // 0, 1: starting points; then arclength steps
    int newton_its   = 0;
    int retries      = 0;    // failed correctors before this point
    double ds        = 0.0;  // arclength step that produced it (0, 1: none)
    double arclength = 0.0;  // accumulated |z_i - z_{i-1}|
    double dmu_ds    = 0.0;  // secant d mu / ds into this point (fold: sign change)
    double residual  = 0.0;  // |G(x; mu)|
};

enum class ContinuationStatus
{
    ReachedBound,   // mu left [mu_min, mu_max]
    MaxSteps,
    StoppedByUser,  // on_point returned false
    StepTooSmall,   // corrector kept failing down to ds_min
    StartFailed,    // Newton at mu0 or mu0 + dmu0 did not converge
};

inline std::string to_string(ContinuationStatus s)
{
    switch (s)
    {
        case ContinuationStatus::ReachedBound:
            return "reached-bound";
        case ContinuationStatus::MaxSteps:
            return "max-steps";
        case ContinuationStatus::StoppedByUser:
            return "stopped-by-user";
        case ContinuationStatus::StepTooSmall:
            return "step-too-small";
        case ContinuationStatus::StartFailed:
            return "start-failed";
    }
    return "?";
}

template <class Space> struct ContinuationResult
{
    ContinuationStatus status = ContinuationStatus::MaxSteps;
    std::vector<double> mu;      // parameter at each accepted point
    std::vector<int> newton_its; // corrector iterations at each point
    int failed_correctors = 0;
    ContinuationPoint<Space> last;
};

// G(x, mu, Gx): Gx pre-allocated. jac_x(x, mu, Gx, v, Jv) optional.
// on_point(point) is called for every accepted point (including the two
// starting points); return false to stop. It may change the definition of G
// for later points (e.g. move a phase-condition reference to point.x).
template <class Space>
ContinuationResult<Space> continuation(
    const Space &S,
    std::function<void(const typename Space::Vector &, double,
                       typename Space::Vector &)>
        G,
    const typename Space::Vector &x0, double mu0, double dmu0,
    const ContinuationOptions<Space> &opt,
    std::function<void(const typename Space::Vector &, double,
                       const typename Space::Vector &,
                       const typename Space::Vector &,
                       typename Space::Vector &)>
        jac_x = nullptr,
    std::function<bool(const ContinuationPoint<Space> &)> on_point = nullptr)
{
    using X  = typename Space::Vector;
    using AS = AugmentedSpace<Space>;
    using Z  = typename AS::Vector;
    const AS A(S, opt.mu_weight);

    ContinuationResult<Space> out;
    auto logline = [&](const std::string &s) {
        if (opt.log)
        {
            *opt.log << s << std::endl;
        }
    };
    auto fmt = [](double v, int p = 6) {
        std::ostringstream os;
        os << std::scientific << std::setprecision(p) << v;
        return os.str();
    };

    // ---- Newton at fixed mu (starting points) ------------------------------
    auto solve_fixed = [&](const X &guess, double mu) {
        NewtonOptions<Space> o = opt.newton;
        o.on_accept            = [&](const X &x, double r, int it) {
            if (opt.on_newton)
            {
                opt.on_newton(x, mu, r, it);
            }
            return false;
        };
        std::function<void(const X &, const X &, const X &, X &)> J = nullptr;
        if (jac_x)
        {
            J = [&](const X &x, const X &Gx, const X &v, X &Jv) {
                jac_x(x, mu, Gx, v, Jv);
            };
        }
        return newton_hookstep(
            S, [&](const X &x, X &Gx) { G(x, mu, Gx); }, guess, o, J);
    };

    ContinuationPoint<Space> p0, p1;
    {
        logline("[cont] point 0: Newton at mu = " + fmt(mu0));
        const auto R = solve_fixed(x0, mu0);
        if (!R.converged())
        {
            logline("[cont] Newton at mu0 failed: " + to_string(R.status));
            out.status = ContinuationStatus::StartFailed;
            return out;
        }
        p0.x          = clone(S, R.x);
        p0.mu         = mu0;
        p0.index      = 0;
        p0.newton_its = R.iterations;
        p0.residual   = R.residual_history.back();
    }
    out.mu.push_back(p0.mu);
    out.newton_its.push_back(p0.newton_its);
    out.last = p0;
    if (on_point && !on_point(p0))
    {
        out.status = ContinuationStatus::StoppedByUser;
        return out;
    }

    {
        double dmu = dmu0;
        bool ok    = false;
        for (int tries = 0; tries <= opt.max_retries && !ok; ++tries)
        {
            logline("[cont] point 1: Newton at mu = " + fmt(mu0 + dmu));
            const auto R = solve_fixed(p0.x, mu0 + dmu);
            if (R.converged())
            {
                p1.x          = clone(S, R.x);
                p1.mu         = mu0 + dmu;
                p1.newton_its = R.iterations;
                p1.retries    = tries;
                p1.residual   = R.residual_history.back();
                ok            = true;
            }
            else
            {
                ++out.failed_correctors;
                dmu *= opt.shrink;
            }
        }
        if (!ok)
        {
            logline("[cont] natural-parameter step failed");
            out.status = ContinuationStatus::StartFailed;
            return out;
        }
        X d = clone(S, p1.x);
        S.axpy(-1.0, p0.x, d);
        const double len =
            std::sqrt(S.dot(d, d) + opt.mu_weight * (p1.mu - p0.mu) * (p1.mu - p0.mu));
        p1.index     = 1;
        p1.arclength = len;
        p1.dmu_ds    = (p1.mu - p0.mu) / len;
    }
    out.mu.push_back(p1.mu);
    out.newton_its.push_back(p1.newton_its);
    out.last = p1;
    if (on_point && !on_point(p1))
    {
        out.status = ContinuationStatus::StoppedByUser;
        return out;
    }

    // ---- pseudo-arclength steps ---------------------------------------------
    double ds = opt.ds;
    for (int step = 1; step <= opt.max_steps; ++step)
    {
        if (p1.mu < opt.mu_min || p1.mu > opt.mu_max)
        {
            out.status = ContinuationStatus::ReachedBound;
            return out;
        }

        // Secant tangent through the last two points.
        Z z0{clone(S, p0.x), p0.mu}, z1{clone(S, p1.x), p1.mu};
        Z t = clone(A, z1);
        A.axpy(-1.0, z0, t);
        A.scale(1.0 / norm(A, t), t);

        bool accepted = false;
        int retries   = 0;
        for (; retries <= opt.max_retries; ++retries)
        {
            Z zp = clone(A, z1);
            A.axpy(ds, t, zp);

            auto H = [&](const Z &z, Z &Hz) {
                G(z.x, z.mu, Hz.x);
                Z d = clone(A, z);
                A.axpy(-1.0, z1, d);
                Hz.mu = A.dot(t, d) - ds;
            };

            // dG/dmu at the current iterate, cached per iterate.
            Z cache_z;
            bool have_cache = false;
            X dGdmu;
            auto J = [&](const Z &z, const Z &Hz, const Z &v, Z &Jv) {
                bool same = have_cache && cache_z.mu == z.mu;
                if (same)
                {
                    X d = clone(S, z.x);
                    S.axpy(-1.0, cache_z.x, d);
                    same = S.dot(d, d) == 0.0;
                }
                if (!same)
                {
                    const double h = opt.fd_eps_mu * (1.0 + std::abs(z.mu));
                    dGdmu          = S.zeros_like(z.x);
                    G(z.x, z.mu + h, dGdmu);
                    S.axpy(-1.0, Hz.x, dGdmu);
                    S.scale(1.0 / h, dGdmu);
                    cache_z    = clone(A, z);
                    have_cache = true;
                }
                // J_x v_x at fixed mu.
                if (jac_x)
                {
                    jac_x(z.x, z.mu, Hz.x, v.x, Jv.x);
                }
                else
                {
                    const double vn = norm(S, v.x);
                    if (vn == 0.0)
                    {
                        S.scale(0.0, Jv.x);
                    }
                    else
                    {
                        const double eps =
                            opt.newton.fd_eps * (1.0 + norm(S, z.x)) / vn;
                        X xp = clone(S, z.x);
                        S.axpy(eps, v.x, xp);
                        G(xp, z.mu, Jv.x);
                        S.axpy(-1.0, Hz.x, Jv.x);
                        S.scale(1.0 / eps, Jv.x);
                    }
                }
                S.axpy(v.mu, dGdmu, Jv.x);
                Jv.mu = A.dot(t, v);
            };

            NewtonOptions<AS> o;
            o.max_newton     = opt.newton.max_newton;
            o.tol            = opt.newton.tol;
            o.max_krylov     = opt.newton.max_krylov;
            o.gmres_rtol     = opt.newton.gmres_rtol;
            o.delta0         = opt.newton.delta0;
            o.delta_min      = opt.newton.delta_min;
            o.delta_max      = opt.newton.delta_max;
            o.max_hook_tries = opt.newton.max_hook_tries;
            o.fd_eps         = opt.newton.fd_eps;
            o.log            = opt.newton.log;
            o.on_accept      = [&](const Z &z, double r, int it) {
                if (opt.on_newton)
                {
                    opt.on_newton(z.x, z.mu, r, it);
                }
                return false;
            };

            logline("[cont] step " + std::to_string(step) + ": ds = " +
                    fmt(ds, 3) + ", predictor mu = " + fmt(zp.mu));
            const auto R = newton_hookstep(A, H, zp, o, J);
            if (R.converged())
            {
                ContinuationPoint<Space> p;
                p.x          = clone(S, R.x.x);
                p.mu         = R.x.mu;
                p.index      = p1.index + 1;
                p.newton_its = R.iterations;
                p.retries    = retries;
                p.ds         = ds;
                p.residual   = R.residual_history.back();
                Z d          = clone(A, R.x);
                A.axpy(-1.0, z1, d);
                const double len = norm(A, d);
                p.arclength      = p1.arclength + len;
                p.dmu_ds         = (p.mu - p1.mu) / len;
                if (p.dmu_ds * p1.dmu_ds < 0.0)
                {
                    logline("[cont] fold between mu = " + fmt(p1.mu) +
                            " and " + fmt(p.mu));
                }
                p0 = std::move(p1);
                p1 = std::move(p);
                accepted = true;

                if (R.iterations <= opt.easy_its)
                {
                    ds = std::min(ds * opt.grow, opt.ds_max);
                }
                else if (R.iterations >= opt.hard_its)
                {
                    ds = std::max(ds * opt.shrink, opt.ds_min);
                }
                break;
            }
            ++out.failed_correctors;
            logline("[cont] corrector failed (" + to_string(R.status) +
                    "), halving ds");
            ds *= opt.shrink;
            if (ds < opt.ds_min)
            {
                break;
            }
        }
        if (!accepted)
        {
            out.status = ContinuationStatus::StepTooSmall;
            return out;
        }

        out.mu.push_back(p1.mu);
        out.newton_its.push_back(p1.newton_its);
        out.last = p1;
        logline("[cont] point " + std::to_string(p1.index) + ": mu = " +
                fmt(p1.mu) + ", newton its = " + std::to_string(p1.newton_its) +
                ", dmu/ds = " + fmt(p1.dmu_ds, 3));
        if (on_point && !on_point(p1))
        {
            out.status = ContinuationStatus::StoppedByUser;
            return out;
        }
    }
    out.status = (p1.mu < opt.mu_min || p1.mu > opt.mu_max)
                     ? ContinuationStatus::ReachedBound
                     : ContinuationStatus::MaxSteps;
    return out;
}

} // namespace nekecs
