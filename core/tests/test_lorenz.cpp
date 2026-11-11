// End-to-end Newton-Krylov-hookstep on the Lorenz system, set up exactly as
// the Navier-Stokes problem will be: the residual is built from a fixed-step
// flow map Phi_T with dt = T / N (N fixed, so Phi is smooth in T), Jacobian
// actions by finite differences.
//
//   equilibrium:     F(x) = Phi_T(x) - x,                         T fixed
//   periodic orbit:  F(x, T) = (Phi_T(x) - x, f(x0) . (x - x0)),  T unknown
//
// The periodic-orbit initial guess comes from a recurrence search on a
// trajectory, as it would for a turbulent flow.
#include "check.hpp"

#include <nekecs/newton.hpp>
#include <nekecs/space.hpp>

#include <array>
#include <limits>

using namespace nekecs;
using Vec = std::vector<double>;

#include "lorenz.hpp"

using namespace lorenz;


static void test_equilibrium()
{
    StdVectorSpace S;
    const double T = 0.5;
    const int N    = 500;
    auto F = [&](const Vec &x, Vec &Fx) {
        auto y = flow({x[0], x[1], x[2]}, T, N);
        for (int i = 0; i < 3; ++i) Fx[i] = y[i] - x[i];
    };

    NewtonOptions<StdVectorSpace> opt;
    opt.tol = 1e-11;
    const auto R = newton_hookstep(S, F, Vec{7.0, 9.5, 25.0}, opt);

    const double c = std::sqrt(beta * (rho - 1.0));
    CHECK(R.converged(), "Lorenz EQ: converged in " << R.iterations
                                                    << " Newton iterations");
    CHECK_NEAR(R.x[0], c, 1e-8, "Lorenz EQ x = sqrt(beta(rho-1))");
    CHECK_NEAR(R.x[1], c, 1e-8, "Lorenz EQ y = sqrt(beta(rho-1))");
    CHECK_NEAR(R.x[2], rho - 1.0, 1e-8, "Lorenz EQ z = rho - 1");
}


// expectedT < 0: any periodic orbit is acceptable (only convergence and
// closure are checked).
static void test_periodic_orbit(const char *name, double Tmin, double Tmax,
                                double expectedT, bool verbose)
{
    const Guess g = recurrence_guess(Tmin, Tmax);
    std::cout << "     " << name << " recurrence guess: T = " << g.T
              << ", |x(T) - x(0)| = " << g.recurrence << std::endl;

    StdVectorSpace S;
    const int N            = 4000;
    const auto x0          = g.x;
    const auto f0          = f(x0);
    auto F = [&](const Vec &x, Vec &Fx) {
        auto y = flow({x[0], x[1], x[2]}, x[3], N);
        for (int i = 0; i < 3; ++i) Fx[i] = y[i] - x[i];
        // Poincare section through the initial guess, normal to the flow.
        Fx[3] = f0[0] * (x[0] - x0[0]) + f0[1] * (x[1] - x0[1]) +
                f0[2] * (x[2] - x0[2]);
    };

    NewtonOptions<StdVectorSpace> opt;
    opt.tol        = 1e-10;
    opt.max_newton = 40;
    opt.delta0     = 0.1 * g.recurrence + 1e-3; // cautious start: guess is rough
    opt.log        = verbose ? &std::cout : nullptr;
    const auto R = newton_hookstep(S, F, Vec{x0[0], x0[1], x0[2], g.T}, opt);

    CHECK(R.converged(), name << ": " << to_string(R.status) << " after "
                                       << R.iterations << " iterations, "
                                       << R.residual_evals << " F evals");
    std::cout << "     " << name << " converged T = " << R.x[3] << std::endl;
    if (expectedT > 0.0)
    {
        CHECK_NEAR(R.x[3], expectedT, 2e-6, name << " period T");
    }

    // Independent check: integrate the converged orbit with a finer step.
    const auto y = flow({R.x[0], R.x[1], R.x[2]}, R.x[3], 4 * N);
    const double closure = dist(y, {R.x[0], R.x[1], R.x[2]});
    CHECK(closure < 1e-8, name << " closes under 4x finer integration ("
                              << closure << ")");

    // Quadratic convergence in the final iterations.
    const auto &h = R.residual_history;
    if (h.size() >= 3)
    {
        const double a = h[h.size() - 3], b = h[h.size() - 2];
        CHECK(b < 1e-2 * a || b < 1e-8,
              "fast terminal convergence (" << a << " -> " << b << ")");
    }
}

// Guesses where the full Newton step overshoots (delta0 = inf, so the first
// step is unconstrained): the hookstep must reject, shrink and retry.
//
//  - "far" guess (|F0| ~ 30): no convergence is expected from a local method,
//    but ||F|| must decrease monotonically and rejections must occur.
//  - "perturbed" AB guess: must converge to the AB orbit.
static nekecs::NewtonResult<StdVectorSpace>
solve_from(const std::array<double, 3> &x0, double T0, int max_newton)
{
    StdVectorSpace S;
    const int N   = 4000;
    const auto f0 = f(x0);
    auto F = [&](const Vec &x, Vec &Fx) {
        auto y = flow({x[0], x[1], x[2]}, x[3], N);
        for (int i = 0; i < 3; ++i) Fx[i] = y[i] - x[i];
        Fx[3] = f0[0] * (x[0] - x0[0]) + f0[1] * (x[1] - x0[1]) +
                f0[2] * (x[2] - x0[2]);
    };
    NewtonOptions<StdVectorSpace> opt;
    opt.tol        = 1e-10;
    opt.max_newton = max_newton;
    return newton_hookstep(S, F, Vec{x0[0], x0[1], x0[2], T0}, opt);
}

static void test_bad_guesses()
{
    const auto far = solve_from({-13.0, -18.0, 28.5}, 1.7, 30);
    bool monotone  = true;
    for (std::size_t i = 1; i < far.residual_history.size(); ++i)
    {
        monotone &= far.residual_history[i] < far.residual_history[i - 1];
    }
    std::cout << "     far guess: |F| " << far.residual_history.front()
              << " -> " << far.residual_history.back() << " in "
              << far.iterations << " its, " << far.rejected_steps
              << " rejected steps" << std::endl;
    CHECK(monotone, "far guess: |F| decreases monotonically");
    CHECK(far.rejected_steps > 0, "far guess: steps rejected and retried");
    CHECK(far.residual_history.back() < 0.05 * far.residual_history.front(),
          "far guess: |F| reduced by > 20x");

    // Start from the AB orbit perturbed by O(1) in state and 3% in period.
    const Guess g = recurrence_guess(1.3, 1.8);
    const auto R  = solve_from({g.x[0] + 0.8, g.x[1] - 0.8, g.x[2] + 0.8},
                               1.05 * g.T, 60);
    std::cout << "     perturbed AB guess: |F0| = "
              << R.residual_history.front() << ", " << to_string(R.status)
              << ", T = " << R.x[3] << ", " << R.iterations << " its, "
              << R.rejected_steps << " rejected" << std::endl;
    CHECK(R.converged(), "perturbed AB guess converges");
    CHECK_NEAR(R.x[3], 1.5586522, 2e-6, "perturbed AB guess finds AB orbit");
}

int main()
{
    test_equilibrium();
    // Shortest orbit "AB" (one loop around each wing): T = 1.5586522
    // (Viswanath 2003).
    test_periodic_orbit("PO AB", 1.3, 1.8, 1.5586522, true);
    // Longer orbits from rougher recurrences exercise the trust region.
    test_periodic_orbit("PO T~2.3", 2.1, 2.5, -1.0, true);
    test_periodic_orbit("PO T~3.0", 2.8, 3.3, -1.0, false);
    test_bad_guesses();
    return test::report("test_lorenz");
}
