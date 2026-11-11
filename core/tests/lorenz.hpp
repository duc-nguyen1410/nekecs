// Lorenz system helpers shared by the core tests.
#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace lorenz
{
constexpr double sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;

inline std::array<double, 3> f(const std::array<double, 3> &u)
{
    return {sigma * (u[1] - u[0]), u[0] * (rho - u[2]) - u[1],
            u[0] * u[1] - beta * u[2]};
}

// Classical RK4 with N fixed steps of size T/N.
inline std::array<double, 3> flow(std::array<double, 3> u, double T, int N)
{
    const double h = T / N;
    for (int n = 0; n < N; ++n)
    {
        auto k1 = f(u);
        std::array<double, 3> t;
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

// Tangent-linear flow: D Phi_T(u) v, by RK4 on the augmented system
// (u, v)' = (f(u), Df(u) v) with the same N steps as `flow`. Exact Jacobian
// of the discrete map up to round-off (no finite-difference error).
inline std::array<double, 3> flow_tangent(std::array<double, 3> u,
                                          std::array<double, 3> v, double T,
                                          int N)
{
    using A3 = std::array<double, 3>;
    auto Jv = [](const A3 &x, const A3 &w) -> A3 {
        return {sigma * (w[1] - w[0]),
                (rho - x[2]) * w[0] - w[1] - x[0] * w[2],
                x[1] * w[0] + x[0] * w[1] - beta * w[2]};
    };
    const double h = T / N;
    for (int n = 0; n < N; ++n)
    {
        A3 tu, tv;
        auto k1 = f(u);
        auto l1 = Jv(u, v);
        for (int i = 0; i < 3; ++i) { tu[i] = u[i] + 0.5 * h * k1[i]; tv[i] = v[i] + 0.5 * h * l1[i]; }
        auto k2 = f(tu);
        auto l2 = Jv(tu, tv);
        for (int i = 0; i < 3; ++i) { tu[i] = u[i] + 0.5 * h * k2[i]; tv[i] = v[i] + 0.5 * h * l2[i]; }
        auto k3 = f(tu);
        auto l3 = Jv(tu, tv);
        for (int i = 0; i < 3; ++i) { tu[i] = u[i] + h * k3[i]; tv[i] = v[i] + h * l3[i]; }
        auto k4 = f(tu);
        auto l4 = Jv(tu, tv);
        for (int i = 0; i < 3; ++i)
        {
            u[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            v[i] += h / 6.0 * (l1[i] + 2.0 * l2[i] + 2.0 * l3[i] + l4[i]);
        }
    }
    return v;
}

inline double dist(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) +
                     (a[1] - b[1]) * (a[1] - b[1]) +
                     (a[2] - b[2]) * (a[2] - b[2]));
}

struct Guess
{
    std::array<double, 3> x;
    double T;
    double recurrence;
};

// Close recurrence |x(t+T) - x(t)| with T in [Tmin, Tmax] on the attractor.
inline Guess recurrence_guess(double Tmin, double Tmax)
{
    const double dt = 0.01;
    std::array<double, 3> u{1.0, 1.0, 1.0};
    u = flow(u, 20.0, 4000); // settle onto the attractor
    std::vector<std::array<double, 3>> traj;
    for (int n = 0; n < 6000; ++n)
    {
        traj.push_back(u);
        u = flow(u, dt, 2);
    }
    Guess best{traj[0], 0.0, std::numeric_limits<double>::infinity()};
    const int lo = int(Tmin / dt), hi = int(Tmax / dt);
    for (std::size_t i = 0; i + hi < traj.size(); ++i)
    {
        for (int l = lo; l <= hi; ++l)
        {
            const double d = dist(traj[i], traj[i + l]);
            if (d < best.recurrence)
            {
                best = {traj[i], l * dt, d};
            }
        }
    }
    return best;
}

} // namespace lorenz
