// Newton-Krylov-hookstep solver for F(x) = 0, F: X -> X (square system).
//
// The backend packs everything into X and F: for a periodic orbit,
// x = (u, T) and F = (Phi_T(u) - u, phase condition). The core only sees
// vectors, an inner product, and F.
//
// Each Newton iteration:
//   1. GMRES on J dx = -F(x), J v by finite differences (or a user operator).
//   2. Hookstep in the Krylov space with trust radius delta.
//   3. Evaluate F(x + dx). Compare the actual decrease of ||F||^2 with the
//      decrease predicted by the linear model, rho = actual / predicted.
//      Reject (rho < accept_ratio or ||F|| grew): shrink delta, re-solve the
//      hookstep in the SAME Krylov space (no new GMRES), retry.
//      Accept: shrink delta if rho < shrink_ratio; grow it if rho > grow_ratio
//      and the step was constrained by the trust region.
#pragma once

#include "gmres.hpp"
#include "hookstep.hpp"
#include "space.hpp"

#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace nekecs
{

template <class Space> struct NewtonOptions
{
    using Vector = typename Space::Vector;

    int max_newton    = 50;
    double tol        = 1e-10; // stop when ||F|| <= tol (absolute; scale F in the backend)
    int max_krylov    = 100;   // GMRES dimension (no restarts)
    double gmres_rtol = 1e-3;  // inexact Newton: ||J dx + F|| <= gmres_rtol ||F||

    // Trust region. delta0 <= 0 means "no initial constraint": the first
    // step is the full Newton step and delta adapts from there.
    double delta0      = -1.0;
    double delta_min   = 1e-12;
    double delta_max   = std::numeric_limits<double>::infinity();
    int max_hook_tries = 12; // hookstep retries per Newton iteration
    double accept_ratio = 0.01;
    double shrink_ratio = 0.25;
    double grow_ratio   = 0.75;
    double shrink       = 0.5;
    double grow         = 2.0;

    // Finite-difference Jacobian: J v ~ (F(x + eps v) - F(x)) / eps,
    // eps = fd_eps * (1 + ||x||) / ||v||.
    double fd_eps = 1.4901161193847656e-08; // sqrt(machine epsilon)

    // Called after each accepted step (e.g. to write a checkpoint or move a
    // phase-condition reference point). Arguments: x, ||F||, iteration.
    // Return true if the callback changed the definition of F (e.g. moved a
    // phase-condition reference), so F(x) is re-evaluated.
    std::function<bool(const Vector &, double, int)> on_accept;

    std::ostream *log = nullptr;
};

enum class NewtonStatus
{
    Converged,
    MaxIterations,
    TrustRegionCollapsed,
};

inline std::string to_string(NewtonStatus s)
{
    switch (s)
    {
        case NewtonStatus::Converged:
            return "converged";
        case NewtonStatus::MaxIterations:
            return "max-iterations";
        case NewtonStatus::TrustRegionCollapsed:
            return "trust-region-collapsed";
    }
    return "?";
}

template <class Space> struct NewtonResult
{
    typename Space::Vector x;
    NewtonStatus status = NewtonStatus::MaxIterations;
    int iterations      = 0;
    int residual_evals  = 0; // calls to F (includes Jacobian FD evaluations)
    int rejected_steps  = 0; // hooksteps rejected and retried with smaller delta
    std::vector<double> residual_history; // ||F|| at each accepted iterate
    bool converged() const
    {
        return status == NewtonStatus::Converged;
    }
};

// F: callable (const Vector &x, Vector &Fx), Fx pre-allocated.
// Jac (optional): callable (const Vector &x, const Vector &Fx,
//                           const Vector &v, Vector &Jv).
template <class Space, class Residual>
NewtonResult<Space> newton_hookstep(
    const Space &S, Residual &&F, const typename Space::Vector &x0,
    const NewtonOptions<Space> &opt,
    std::function<void(const typename Space::Vector &,
                       const typename Space::Vector &,
                       const typename Space::Vector &,
                       typename Space::Vector &)>
        jacobian = nullptr)
{
    using Vector = typename Space::Vector;
    NewtonResult<Space> out;
    out.x = clone(S, x0);

    auto evalF = [&](const Vector &x, Vector &Fx) {
        F(x, Fx);
        ++out.residual_evals;
    };

    Vector Fx = S.zeros_like(x0);
    evalF(out.x, Fx);
    double r = norm(S, Fx);
    out.residual_history.push_back(r);

    double delta = opt.delta0 > 0.0 ? opt.delta0
                                    : std::numeric_limits<double>::infinity();

    auto logline = [&](const std::string &s) {
        if (opt.log)
        {
            *opt.log << s << std::endl;
        }
    };
    auto fmt = [](double v) {
        std::ostringstream os;
        os << std::scientific << std::setprecision(3) << v;
        return os.str();
    };

    logline("[newton] it 0  |F| = " + fmt(r));

    for (int it = 1; it <= opt.max_newton; ++it)
    {
        if (r <= opt.tol)
        {
            out.status     = NewtonStatus::Converged;
            out.iterations = it - 1;
            return out;
        }

        // J v operator at the current iterate.
        const double xnorm = norm(S, out.x);
        auto Jop = [&](const Vector &v, Vector &Jv) {
            if (jacobian)
            {
                jacobian(out.x, Fx, v, Jv);
                return;
            }
            const double vnorm = norm(S, v);
            if (vnorm == 0.0)
            {
                S.scale(0.0, Jv);
                return;
            }
            const double eps = opt.fd_eps * (1.0 + xnorm) / vnorm;
            Vector xp        = clone(S, out.x);
            S.axpy(eps, v, xp);
            evalF(xp, Jv);
            S.axpy(-1.0, Fx, Jv);
            S.scale(1.0 / eps, Jv);
        };

        Vector b = clone(S, Fx);
        S.scale(-1.0, b);
        const auto K = gmres(S, Jop, b, opt.max_krylov, opt.gmres_rtol);
        const HookstepSolver hook(K.H, K.beta);

        logline("[newton] it " + std::to_string(it) + "  gmres k = " +
                std::to_string(K.k) + "  lin.res = " +
                fmt(K.beta > 0 ? K.residual / K.beta : 0.0) +
                "  |dx_newton| = " + fmt(hook.newton_step_norm()));

        bool accepted = false;
        for (int tries = 0; tries < opt.max_hook_tries; ++tries)
        {
            const HookstepResult h = hook.solve(delta);
            Vector xnew            = clone(S, out.x);
            S.axpy(1.0, expand(S, K, h.y), xnew);
            Vector Fnew = S.zeros_like(x0);
            evalF(xnew, Fnew);
            const double rnew = norm(S, Fnew);

            const double pred = r * r - h.predicted_residual *
                                            h.predicted_residual;
            const double act = r * r - rnew * rnew;
            const double rho = pred > 0.0 ? act / pred : -1.0;

            logline("[newton]   try " + std::to_string(tries) +
                    "  delta = " + fmt(delta) + "  |dx| = " + fmt(h.ynorm) +
                    (h.constrained ? " (hook)" : " (newton)") +
                    "  |F_new| = " + fmt(rnew) + "  rho = " + fmt(rho));

            if (std::isfinite(rnew) && rnew < r && rho >= opt.accept_ratio)
            {
                // Trust-region update for the next iteration.
                if (rho < opt.shrink_ratio)
                {
                    delta = opt.shrink * h.ynorm;
                }
                else if (rho > opt.grow_ratio && h.constrained)
                {
                    delta = std::min(opt.grow * delta, opt.delta_max);
                }
                out.x = std::move(xnew);
                Fx    = std::move(Fnew);
                r     = rnew;
                out.residual_history.push_back(r);
                accepted = true;
                break;
            }

            ++out.rejected_steps;
            delta = opt.shrink * std::min(delta, h.ynorm);
            if (delta < opt.delta_min)
            {
                break;
            }
        }

        if (!accepted)
        {
            logline("[newton] trust region collapsed (delta = " +
                    fmt(delta) + ")");
            out.status     = NewtonStatus::TrustRegionCollapsed;
            out.iterations = it;
            return out;
        }

        logline("[newton] it " + std::to_string(it) + "  |F| = " + fmt(r));
        if (opt.on_accept && opt.on_accept(out.x, r, it))
        {
            evalF(out.x, Fx);
            r = norm(S, Fx);
            out.residual_history.back() = r;
        }
        out.iterations = it;
    }

    out.status = r <= opt.tol ? NewtonStatus::Converged
                              : NewtonStatus::MaxIterations;
    return out;
}

} // namespace nekecs
