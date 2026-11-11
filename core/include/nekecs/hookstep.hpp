// Hookstep (Viswanath 2007, after Dennis & Schnabel): the trust-region
// subproblem restricted to the GMRES Krylov space,
//
//     min_y || beta e_0 - H y ||   subject to   ||y|| <= delta,
//
// with H the (k+1) x k Hessenberg matrix. Since V_{k+1} has orthonormal
// columns, ||V_k y|| = ||y|| and ||F + J dx|| ~ ||beta e_0 - H y|| for
// dx = V_k y, so this is the full-space trust-region step restricted to the
// Krylov space.
//
// With H = U diag(s) W^T and p = beta U^T e_0, writing y = W z:
//     unconstrained:  z_i = p_i / s_i
//     constrained:    z_i = s_i p_i / (s_i^2 + mu),  mu > 0 chosen so ||z|| = delta
// ||z(mu)|| decreases monotonically in mu, so mu is found by bracketing and
// bisection. The SVD is computed once per Krylov space; each new delta (after
// a rejected step) costs O(k^2).
#pragma once

#include "dense.hpp"

#include <cmath>
#include <vector>

namespace nekecs
{

struct HookstepResult
{
    std::vector<double> y;           // Krylov coefficients of the step
    double ynorm              = 0.0; // ||y|| = ||dx||
    double mu                 = 0.0; // Levenberg-Marquardt parameter (0 = Newton)
    bool constrained          = false;
    double predicted_residual = 0.0; // || beta e_0 - H y ||
};

class HookstepSolver
{
public:
    // rank_tol: singular values below rank_tol * s_max are treated as zero
    // (the unconstrained step is then the minimum-norm least-squares step).
    HookstepSolver(const Matrix &H, double beta, double rank_tol = 1e-13)
        : m_H(H), m_svd(svd(H)), m_beta(beta)
    {
        const int k = H.cols();
        m_p.assign(k, 0.0);
        for (int i = 0; i < k; ++i)
        {
            m_p[i] = beta * m_svd.U(0, i);
        }
        const double smax = k > 0 ? m_svd.s[0] : 0.0;
        m_rank            = 0;
        for (int i = 0; i < k; ++i)
        {
            if (m_svd.s[i] > rank_tol * smax)
            {
                ++m_rank;
            }
        }
    }

    const std::vector<double> &singular_values() const
    {
        return m_svd.s;
    }

    // Norm of the unconstrained (Newton / least-squares) step.
    double newton_step_norm() const
    {
        return znorm(0.0);
    }

    HookstepResult solve(double delta) const
    {
        HookstepResult r;
        double mu = 0.0;
        if (!(znorm(0.0) <= delta))
        {
            r.constrained = true;
            double lo = 0.0, hi = 1.0;
            const double s0 = m_svd.s.empty() ? 1.0 : m_svd.s[0];
            hi = std::max(1e-300, s0 * s0 * 1e-8);
            while (znorm(hi) > delta && hi < 1e300)
            {
                hi *= 4.0;
            }
            for (int it = 0; it < 200; ++it)
            {
                const double mid = 0.5 * (lo + hi);
                if (znorm(mid) > delta)
                {
                    lo = mid;
                }
                else
                {
                    hi = mid;
                }
                if (hi - lo <= 1e-14 * hi)
                {
                    break;
                }
            }
            mu = hi; // ||z(hi)|| <= delta: never exceed the trust region
        }
        r.mu = mu;

        const std::vector<double> z = zvec(mu);
        const int k                 = static_cast<int>(z.size());
        r.y.assign(k, 0.0);
        for (int i = 0; i < k; ++i)
        {
            for (int j = 0; j < k; ++j)
            {
                r.y[i] += m_svd.W(i, j) * z[j];
            }
        }
        r.ynorm = norm2(z);

        // Evaluate ||beta e_0 - H y|| directly: the SVD form
        // beta^2 - ||p||^2 + sum (p_i - s_i z_i)^2 cancels catastrophically
        // once the linear residual is small.
        std::vector<double> res = m_H.multiply(r.y);
        for (auto &v : res)
        {
            v = -v;
        }
        res[0] += m_beta;
        r.predicted_residual = norm2(res);
        return r;
    }

private:
    std::vector<double> zvec(double mu) const
    {
        const int k = static_cast<int>(m_p.size());
        std::vector<double> z(k, 0.0);
        for (int i = 0; i < m_rank; ++i)
        {
            const double s = m_svd.s[i];
            z[i]           = mu == 0.0 ? m_p[i] / s : s * m_p[i] / (s * s + mu);
        }
        return z;
    }
    double znorm(double mu) const
    {
        return norm2(zvec(mu));
    }

    Matrix m_H;
    SVD m_svd;
    double m_beta;
    std::vector<double> m_p;
    int m_rank = 0;
};

} // namespace nekecs
