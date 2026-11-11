// Unrestarted GMRES for A x = b with x0 = 0, keeping the Arnoldi basis and the
// Hessenberg matrix so the caller can re-solve the projected problem (the
// hookstep) without new operator applications.
//
// Arnoldi relation on exit, with k = result.k:
//     A V_k = V_{k+1} H,   H is (k+1) x k,   V_{k+1}^T V_{k+1} = I,
//     b = beta v_0.
// Orthogonalisation is modified Gram-Schmidt applied twice ("twice is
// enough"): operator applications dominate the cost, and loss of
// orthogonality would corrupt the hookstep, which assumes ||V y|| = ||y||.
#pragma once

#include "dense.hpp"
#include "space.hpp"

#include <cmath>
#include <vector>

namespace nekecs
{

template <class Space> struct KrylovResult
{
    using Vector = typename Space::Vector;

    std::vector<Vector> basis; // k+1 vectors (k if the last was not needed)
    Matrix H;                  // (k+1) x k Hessenberg matrix
    double beta = 0.0;         // ||b||
    int k = 0;                 // Krylov dimension used
    std::vector<double> y;     // least-squares coefficients: x = V_k y
    double residual = 0.0;     // ||b - A x|| estimated from the projection
    bool converged = false;    // residual <= rtol * beta
    bool breakdown = false;    // invariant subspace found (exact solve)
};

// x = V_k y for the first k basis vectors.
template <class Space>
typename Space::Vector expand(const Space &S, const KrylovResult<Space> &K,
                              const std::vector<double> &y)
{
    auto x = S.zeros_like(K.basis.at(0));
    for (std::size_t i = 0; i < y.size(); ++i)
    {
        S.axpy(y[i], K.basis[i], x);
    }
    return x;
}

// Op: callable (const Vector &v, Vector &out), out <- A v. `out` is
// pre-allocated with zeros_like.
template <class Space, class Op>
KrylovResult<Space> gmres(const Space &S, const Op &A,
                          const typename Space::Vector &b, int max_iter,
                          double rtol)
{
    using Vector = typename Space::Vector;
    KrylovResult<Space> R;
    R.beta = norm(S, b);
    if (R.beta == 0.0)
    {
        R.basis.push_back(S.zeros_like(b));
        R.H = Matrix(1, 0);
        R.converged = true;
        return R;
    }

    const int m = max_iter;
    Matrix H(m + 1, m);
    Matrix Rt(m + 1, m); // H with Givens rotations applied (upper triangular)
    std::vector<double> cs(m, 0.0), sn(m, 0.0), g(m + 1, 0.0);
    g[0] = R.beta;

    Vector v0 = clone(S, b);
    S.scale(1.0 / R.beta, v0);
    R.basis.push_back(std::move(v0));

    int k = 0;
    double resid = R.beta;
    for (int j = 0; j < m; ++j)
    {
        Vector w = S.zeros_like(b);
        A(R.basis[j], w);

        for (int pass = 0; pass < 2; ++pass)
        {
            for (int i = 0; i <= j; ++i)
            {
                const double h = S.dot(R.basis[i], w);
                H(i, j) += h;
                S.axpy(-h, R.basis[i], w);
            }
        }
        const double hnext = norm(S, w);
        H(j + 1, j) = hnext;
        k = j + 1;

        // Apply previous rotations to the new column, then a new rotation to
        // zero H(j+1, j).
        for (int i = 0; i <= j + 1; ++i)
        {
            Rt(i, j) = H(i, j);
        }
        for (int i = 0; i < j; ++i)
        {
            const double a = Rt(i, j), c = Rt(i + 1, j);
            Rt(i, j)     = cs[i] * a + sn[i] * c;
            Rt(i + 1, j) = -sn[i] * a + cs[i] * c;
        }
        const double a = Rt(j, j), c = Rt(j + 1, j);
        const double r = std::hypot(a, c);
        cs[j] = r > 0.0 ? a / r : 1.0;
        sn[j] = r > 0.0 ? c / r : 0.0;
        Rt(j, j)     = r;
        Rt(j + 1, j) = 0.0;
        g[j + 1] = -sn[j] * g[j];
        g[j]     = cs[j] * g[j];
        resid    = std::abs(g[j + 1]);

        // Happy breakdown: the Krylov space is invariant, solve is exact.
        const double hscale = std::abs(H(j, j)) + hnext;
        if (hnext <= 1e-14 * (hscale > 0.0 ? hscale : 1.0))
        {
            R.breakdown = true;
            break;
        }

        S.scale(1.0 / hnext, w);
        R.basis.push_back(std::move(w));

        if (resid <= rtol * R.beta)
        {
            break;
        }
    }

    // Back substitution on the k x k triangle.
    std::vector<double> y(k, 0.0);
    for (int i = k - 1; i >= 0; --i)
    {
        double s = g[i];
        for (int l = i + 1; l < k; ++l)
        {
            s -= Rt(i, l) * y[l];
        }
        y[i] = Rt(i, i) != 0.0 ? s / Rt(i, i) : 0.0;
    }

    R.k         = k;
    R.H         = H.block(k + 1, k);
    R.y         = std::move(y);
    R.residual  = resid;
    R.converged = resid <= rtol * R.beta || R.breakdown;
    return R;
}

} // namespace nekecs
