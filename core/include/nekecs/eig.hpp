// Eigenvalues and eigenvectors of a small real upper-Hessenberg matrix
// (the Arnoldi projection), using complex arithmetic throughout:
//
//  - eigenvalues: shifted QR iteration on the complex Hessenberg matrix with
//    Wilkinson shifts and Givens rotations, deflating from the bottom;
//  - eigenvectors: inverse iteration with a complex LU (partial pivoting).
//
// Complex arithmetic avoids the bookkeeping of real double-shift (Francis)
// steps; for the sizes used here (k <= a few hundred) the cost is negligible
// next to one flow-map evaluation.
#pragma once

#include "dense.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace nekecs
{

using cplx = std::complex<double>;

class CMatrix
{
public:
    CMatrix() = default;
    CMatrix(int n, int m) : m_rows(n), m_cols(m), m_data(std::size_t(n) * m)
    {
    }
    explicit CMatrix(const Matrix &A) : CMatrix(A.rows(), A.cols())
    {
        for (int i = 0; i < m_rows; ++i)
        {
            for (int j = 0; j < m_cols; ++j)
            {
                (*this)(i, j) = A(i, j);
            }
        }
    }
    int rows() const
    {
        return m_rows;
    }
    int cols() const
    {
        return m_cols;
    }
    cplx &operator()(int i, int j)
    {
        return m_data[std::size_t(i) * m_cols + j];
    }
    cplx operator()(int i, int j) const
    {
        return m_data[std::size_t(i) * m_cols + j];
    }

private:
    int m_rows = 0, m_cols = 0;
    std::vector<cplx> m_data;
};

// All eigenvalues of the leading n x n block of an upper-Hessenberg matrix.
inline std::vector<cplx> hessenberg_eigenvalues(const Matrix &Hin, int n,
                                                int max_iter_per_eig = 100)
{
    CMatrix H(n, n);
    double scale = 0.0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            H(i, j) = Hin(i, j);
            scale   = std::max(scale, std::abs(Hin(i, j)));
        }
    }
    std::vector<cplx> eig(n);
    if (n == 0)
    {
        return eig;
    }
    if (scale == 0.0)
    {
        return eig; // zero matrix
    }

    const double eps = 1e-15;
    int hi           = n - 1;
    int iter         = 0;
    while (hi >= 0)
    {
        if (hi == 0)
        {
            eig[0] = H(0, 0);
            break;
        }
        // Find the start of the active unreduced block [lo, hi].
        int lo = hi;
        while (lo > 0)
        {
            const double s = std::abs(H(lo - 1, lo - 1)) + std::abs(H(lo, lo));
            if (std::abs(H(lo, lo - 1)) <= eps * (s > 0.0 ? s : scale))
            {
                H(lo, lo - 1) = 0.0;
                break;
            }
            --lo;
        }
        if (lo == hi)
        {
            eig[hi] = H(hi, hi);
            --hi;
            iter = 0;
            continue;
        }
        if (++iter > max_iter_per_eig)
        {
            throw std::runtime_error("hessenberg_eigenvalues: no convergence");
        }

        // Wilkinson shift from the trailing 2x2 block; exceptional shift
        // every 10 iterations to break cycles.
        const cplx a = H(hi - 1, hi - 1), b = H(hi - 1, hi), c = H(hi, hi - 1),
                   d = H(hi, hi);
        cplx mu;
        if (iter % 10 == 0)
        {
            mu = d + cplx(std::abs(c), 0.0) * 0.75;
        }
        else
        {
            const cplx tr = 0.5 * (a + d);
            const cplx disc = std::sqrt(0.25 * (a - d) * (a - d) + b * c);
            const cplx m1 = tr + disc, m2 = tr - disc;
            mu = std::abs(m1 - d) < std::abs(m2 - d) ? m1 : m2;
        }

        // QR step on H[lo..hi] - mu I via Givens rotations: H <- R Q + mu I.
        const int m = hi - lo + 1;
        std::vector<cplx> cs(m - 1);
        std::vector<cplx> sn(m - 1);
        for (int i = lo; i <= hi; ++i)
        {
            H(i, i) -= mu;
        }
        for (int k = lo; k < hi; ++k)
        {
            const cplx x = H(k, k), y = H(k + 1, k);
            const double r = std::sqrt(std::norm(x) + std::norm(y));
            cplx c_ = 1.0, s_ = 0.0;
            if (r > 0.0)
            {
                c_ = x / r;
                s_ = y / r;
            }
            cs[k - lo] = c_;
            sn[k - lo] = s_;
            // Apply G^H = [[conj(c), conj(s)], [-s, c]] to rows k, k+1:
            // maps (x, y) to (r, 0).
            for (int j = k; j < n; ++j)
            {
                const cplx t1 = H(k, j), t2 = H(k + 1, j);
                H(k, j)     = std::conj(c_) * t1 + std::conj(s_) * t2;
                H(k + 1, j) = -s_ * t1 + c_ * t2;
            }
        }
        for (int k = lo; k < hi; ++k)
        {
            const cplx c_ = cs[k - lo], s_ = sn[k - lo];
            // Right-multiply columns k, k+1 by G = (G^H)^H =
            // [[c, -conj(s)], [s, conj(c)]]:
            //   col_k <- c col_k + s col_{k+1},
            //   col_{k+1} <- -conj(s) col_k + conj(c) col_{k+1}.
            const int imax = std::min(hi, k + 2);
            for (int i = 0; i <= imax; ++i)
            {
                const cplx t1 = H(i, k), t2 = H(i, k + 1);
                H(i, k)     = t1 * c_ + t2 * s_;
                H(i, k + 1) = -t1 * std::conj(s_) + t2 * std::conj(c_);
            }
        }
        for (int i = lo; i <= hi; ++i)
        {
            H(i, i) += mu;
        }
    }
    return eig;
}

// Unit eigenvector y of the leading n x n block of H for eigenvalue lambda,
// by inverse iteration.
inline std::vector<cplx> hessenberg_eigenvector(const Matrix &Hin, int n,
                                                cplx lambda)
{
    double scale = 0.0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            scale = std::max(scale, std::abs(Hin(i, j)));
        }
    }
    // Tiny perturbation keeps (H - lambda I) numerically invertible.
    const cplx shift = lambda + cplx(1e-13 * (scale > 0 ? scale : 1.0), 0.0);

    CMatrix A(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A(i, j) = Hin(i, j);
        }
        A(i, i) -= shift;
    }
    // LU with partial pivoting.
    std::vector<int> piv(n);
    for (int k = 0; k < n; ++k)
    {
        int p = k;
        for (int i = k + 1; i < n; ++i)
        {
            if (std::abs(A(i, k)) > std::abs(A(p, k)))
            {
                p = i;
            }
        }
        piv[k] = p;
        if (p != k)
        {
            for (int j = 0; j < n; ++j)
            {
                std::swap(A(k, j), A(p, j));
            }
        }
        if (std::abs(A(k, k)) == 0.0)
        {
            A(k, k) = 1e-300;
        }
        for (int i = k + 1; i < n; ++i)
        {
            A(i, k) /= A(k, k);
            for (int j = k + 1; j < n; ++j)
            {
                A(i, j) -= A(i, k) * A(k, j);
            }
        }
    }
    auto solve = [&](std::vector<cplx> b) {
        for (int k = 0; k < n; ++k)
        {
            std::swap(b[k], b[piv[k]]);
            for (int i = k + 1; i < n; ++i)
            {
                b[i] -= A(i, k) * b[k];
            }
        }
        for (int i = n - 1; i >= 0; --i)
        {
            for (int j = i + 1; j < n; ++j)
            {
                b[i] -= A(i, j) * b[j];
            }
            b[i] /= A(i, i);
        }
        return b;
    };
    auto normalise = [](std::vector<cplx> &y) {
        double s = 0.0;
        for (auto &v : y)
        {
            s += std::norm(v);
        }
        s = std::sqrt(s);
        for (auto &v : y)
        {
            v /= s;
        }
    };

    std::vector<cplx> y(n);
    for (int i = 0; i < n; ++i)
    {
        y[i] = cplx(1.0 + 0.1 * i, 0.3 - 0.01 * i); // generic start
    }
    normalise(y);
    for (int it = 0; it < 3; ++it)
    {
        y = solve(y);
        normalise(y);
    }
    return y;
}

} // namespace nekecs
