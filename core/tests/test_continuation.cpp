// Pseudo-arclength continuation, set up like the Navier-Stokes problems:
// residuals built from fixed-step RK4 flow maps, Jacobians by finite
// differences.
//
//  1. Saddle-node: x' = mu - x^2, y' = x - y. Equilibria x = y = +-sqrt(mu)
//     meet at the fold mu = 0. Starting on the upper branch at mu = 1 with
//     decreasing mu, continuation must go round the fold (natural-parameter
//     continuation cannot) and come back up the lower branch.
//  2. Lorenz periodic orbit AB continued in rho from 28 to 32: x = (u, T),
//     phase condition through the previous point (moved by on_point).
#include "check.hpp"

#include <nekecs/continuation.hpp>
#include <nekecs/newton.hpp>

#include <algorithm>
#include <array>
#include <cmath>

using namespace nekecs;
using Vec = std::vector<double>;

#include "lorenz.hpp"

namespace
{

template <class Rhs>
std::array<double, 3> rk4(Rhs f, std::array<double, 3> u, double T, int N)
{
    const double h = T / N;
    for (int n = 0; n < N; ++n)
    {
        std::array<double, 3> t;
        auto k1 = f(u);
        for (int i = 0; i < 3; ++i) t[i] = u[i] + 0.5 * h * k1[i];
        auto k2 = f(t);
        for (int i = 0; i < 3; ++i) t[i] = u[i] + 0.5 * h * k2[i];
        auto k3 = f(t);
        for (int i = 0; i < 3; ++i) t[i] = u[i] + h * k3[i];
        auto k4 = f(t);
        for (int i = 0; i < 3; ++i)
            u[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return u;
}

std::array<double, 3> lorenz_f(const std::array<double, 3> &u, double r)
{
    using lorenz::beta;
    using lorenz::sigma;
    return {sigma * (u[1] - u[0]), u[0] * (r - u[2]) - u[1],
            u[0] * u[1] - beta * u[2]};
}

void test_fold()
{
    StdVectorSpace S;
    const double T = 1.0;
    const int N    = 200;
    auto G = [&](const Vec &x, double mu, Vec &Gx) {
        auto f = [mu](const std::array<double, 3> &u) {
            return std::array<double, 3>{mu - u[0] * u[0], u[0] - u[1], 0.0};
        };
        auto y = rk4(f, {x[0], x[1], 0.0}, T, N);
        Gx[0]  = y[0] - x[0];
        Gx[1]  = y[1] - x[1];
    };

    ContinuationOptions<StdVectorSpace> opt;
    opt.newton.tol = 1e-12;
    opt.ds         = 0.1;
    opt.ds_max     = 0.2;
    opt.mu_max     = 1.0;
    opt.max_steps  = 80;

    std::vector<std::array<double, 3>> pts; // (x, y, mu)
    const auto R = continuation<StdVectorSpace>(
        S, G, Vec{0.9, 1.1}, 1.0, -0.1, opt, nullptr,
        [&](const ContinuationPoint<StdVectorSpace> &p) {
            pts.push_back({p.x[0], p.x[1], p.mu});
            return true;
        });

    double on_branch = 0.0, mu_lo = 1e9;
    int folds        = 0;
    for (std::size_t i = 0; i < pts.size(); ++i)
    {
        const auto &p = pts[i];
        on_branch = std::max({on_branch, std::abs(p[0] * p[0] - p[2]),
                              std::abs(p[1] - p[0])});
        mu_lo     = std::min(mu_lo, p[2]);
        if (i >= 2 && (p[2] - pts[i - 1][2]) * (pts[i - 1][2] - pts[i - 2][2]) < 0)
        {
            ++folds;
        }
    }
    const auto &e = pts.back();
    std::cout << "     fold: " << pts.size() << " points, "
              << to_string(R.status) << ", min mu = " << mu_lo
              << ", end (x, mu) = (" << e[0] << ", " << e[2] << "), "
              << R.failed_correctors << " failed correctors" << std::endl;
    CHECK(R.status == ContinuationStatus::ReachedBound,
          "fold: continuation reaches mu_max on the far side");
    CHECK(on_branch < 1e-9, "fold: every point is an equilibrium (" << on_branch << ")");
    CHECK(folds == 1, "fold: exactly one fold traversed");
    CHECK(mu_lo >= 0.0 && mu_lo < 0.02, "fold: passes close to mu = 0 (" << mu_lo << ")");
    CHECK(e[0] < 0.0 && e[2] > 1.0, "fold: ends on the lower branch x = -sqrt(mu)");
}

void test_lorenz_po()
{
    // AB orbit at rho = 28 from a recurrence (as in test_lorenz).
    const lorenz::Guess g = lorenz::recurrence_guess(1.3, 1.8);
    StdVectorSpace S;
    const int N = 4000;
    std::array<double, 3> xref = g.x;
    std::array<double, 3> fref = lorenz_f(xref, 28.0);

    auto G = [&](const Vec &x, double r, Vec &Gx) {
        auto f = [r](const std::array<double, 3> &u) { return lorenz_f(u, r); };
        auto y = rk4(f, {x[0], x[1], x[2]}, x[3], N);
        for (int i = 0; i < 3; ++i) Gx[i] = y[i] - x[i];
        Gx[3] = fref[0] * (x[0] - xref[0]) + fref[1] * (x[1] - xref[1]) +
                fref[2] * (x[2] - xref[2]);
    };

    ContinuationOptions<StdVectorSpace> opt;
    opt.newton.tol        = 1e-10;
    opt.newton.max_newton = 40;
    opt.ds                = 1.0;
    opt.ds_max            = 2.0;
    opt.mu_max            = 32.0;
    opt.max_steps         = 40;

    std::vector<std::array<double, 2>> rT;
    std::vector<Vec> xs;
    double worst_closure = 0.0;
    const auto R = continuation<StdVectorSpace>(
        S, G, Vec{g.x[0], g.x[1], g.x[2], g.T}, 28.0, 0.25, opt, nullptr,
        [&](const ContinuationPoint<StdVectorSpace> &p) {
            rT.push_back({p.mu, p.x[3]});
            xs.push_back(p.x);
            // Independent check: 4x finer integration closes the orbit.
            const double r = p.mu;
            auto f = [r](const std::array<double, 3> &u) { return lorenz_f(u, r); };
            const std::array<double, 3> u{p.x[0], p.x[1], p.x[2]};
            const auto y = rk4(f, u, p.x[3], 4 * N);
            worst_closure = std::max(worst_closure, lorenz::dist(y, u));
            // Next section through this point, normal to the flow.
            xref = u;
            fref = lorenz_f(u, r);
            return true;
        });

    std::cout << "     Lorenz AB: " << rT.size() << " points, "
              << to_string(R.status) << ", failed correctors "
              << R.failed_correctors << std::endl;
    for (const auto &p : rT)
    {
        std::cout << "       rho = " << p[0] << "  T = " << p[1] << std::endl;
    }
    CHECK(R.status == ContinuationStatus::ReachedBound, "Lorenz AB: reaches rho = 32");
    CHECK_NEAR(rT.front()[1], 1.5586522, 2e-6, "Lorenz AB: T at rho = 28");
    CHECK(worst_closure < 1e-8, "Lorenz AB: every point closes under 4x finer "
                                "integration (" << worst_closure << ")");
    bool monotone = true;
    for (std::size_t i = 1; i < rT.size(); ++i)
    {
        monotone &= rT[i][0] > rT[i - 1][0];
    }
    CHECK(monotone, "Lorenz AB: rho increases monotonically (no spurious fold)");

    // Cross-check the last point with a plain fixed-rho Newton from the
    // previous point (natural-parameter step, no arclength constraint).
    const double rend = rT.back()[0];
    const Vec &xp     = xs[xs.size() - 2];
    xref              = {xp[0], xp[1], xp[2]};
    fref              = lorenz_f(xref, rT[rT.size() - 2][0]);
    NewtonOptions<StdVectorSpace> no;
    no.tol        = 1e-10;
    no.max_newton = 40;
    const auto Rn = newton_hookstep(
        S, [&](const Vec &x, Vec &Gx) { G(x, rend, Gx); }, xp, no);
    CHECK(Rn.converged(), "Lorenz AB: direct Newton at rho_end converges");
    CHECK_NEAR(Rn.x[3], rT.back()[1], 1e-8,
               "Lorenz AB: period at rho_end matches a direct Newton solve");
}

} // namespace

int main()
{
    test_fold();
    test_lorenz_po();
    return test::report("test_continuation");
}
