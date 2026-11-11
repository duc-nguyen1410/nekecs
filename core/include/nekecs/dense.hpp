// Small dense linear algebra for the Krylov-subspace problems (sizes ~10-500).
// Self-contained on purpose: no LAPACK dependency in core.
#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace nekecs
{

// Row-major dense matrix.
class Matrix
{
public:
    Matrix() = default;
    Matrix(int rows, int cols, double value = 0.0)
        : m_rows(rows), m_cols(cols),
          m_data(static_cast<std::size_t>(rows) * cols, value)
    {
    }

    int rows() const
    {
        return m_rows;
    }
    int cols() const
    {
        return m_cols;
    }

    double &operator()(int i, int j)
    {
        assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
        return m_data[static_cast<std::size_t>(i) * m_cols + j];
    }
    double operator()(int i, int j) const
    {
        assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
        return m_data[static_cast<std::size_t>(i) * m_cols + j];
    }

    // Leading rows x cols block (copy).
    Matrix block(int rows, int cols) const
    {
        Matrix B(rows, cols);
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                B(i, j) = (*this)(i, j);
            }
        }
        return B;
    }

    std::vector<double> multiply(const std::vector<double> &x) const
    {
        assert(static_cast<int>(x.size()) == m_cols);
        std::vector<double> y(m_rows, 0.0);
        for (int i = 0; i < m_rows; ++i)
        {
            for (int j = 0; j < m_cols; ++j)
            {
                y[i] += (*this)(i, j) * x[j];
            }
        }
        return y;
    }

private:
    int m_rows = 0;
    int m_cols = 0;
    std::vector<double> m_data;
};

inline double norm2(const std::vector<double> &x)
{
    double s = 0.0;
    for (double v : x)
    {
        s += v * v;
    }
    return std::sqrt(s);
}

// Thin SVD A = U diag(s) W^T of an m x n matrix with m >= n, via one-sided
// Jacobi rotations (accurate for small singular values, which matters for the
// hookstep). U is m x n, s has n entries sorted descending, W is n x n.
struct SVD
{
    Matrix U;
    std::vector<double> s;
    Matrix W;
};

inline SVD svd(const Matrix &A, int max_sweeps = 60, double tol = 1e-15)
{
    const int m = A.rows();
    const int n = A.cols();
    if (m < n)
    {
        throw std::invalid_argument("nekecs::svd requires rows >= cols");
    }

    Matrix U = A;      // columns rotated into orthogonal directions
    Matrix W(n, n);    // accumulated right rotations
    for (int i = 0; i < n; ++i)
    {
        W(i, i) = 1.0;
    }

    for (int sweep = 0; sweep < max_sweeps; ++sweep)
    {
        bool rotated = false;
        for (int p = 0; p < n - 1; ++p)
        {
            for (int q = p + 1; q < n; ++q)
            {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (int i = 0; i < m; ++i)
                {
                    alpha += U(i, p) * U(i, p);
                    beta += U(i, q) * U(i, q);
                    gamma += U(i, p) * U(i, q);
                }
                if (std::abs(gamma) <= tol * std::sqrt(alpha * beta) ||
                    gamma == 0.0)
                {
                    continue;
                }
                rotated = true;
                const double zeta = (beta - alpha) / (2.0 * gamma);
                const double t = std::copysign(1.0, zeta) /
                                 (std::abs(zeta) + std::sqrt(1.0 + zeta * zeta));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double sn = c * t;
                for (int i = 0; i < m; ++i)
                {
                    const double up = U(i, p), uq = U(i, q);
                    U(i, p) = c * up - sn * uq;
                    U(i, q) = sn * up + c * uq;
                }
                for (int i = 0; i < n; ++i)
                {
                    const double wp = W(i, p), wq = W(i, q);
                    W(i, p) = c * wp - sn * wq;
                    W(i, q) = sn * wp + c * wq;
                }
            }
        }
        if (!rotated)
        {
            break;
        }
    }

    std::vector<double> s(n);
    for (int j = 0; j < n; ++j)
    {
        double nrm = 0.0;
        for (int i = 0; i < m; ++i)
        {
            nrm += U(i, j) * U(i, j);
        }
        s[j] = std::sqrt(nrm);
    }

    // Sort descending and normalise U's columns.
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return s[a] > s[b]; });

    SVD out{Matrix(m, n), std::vector<double>(n), Matrix(n, n)};
    for (int k = 0; k < n; ++k)
    {
        const int j = order[k];
        out.s[k] = s[j];
        for (int i = 0; i < m; ++i)
        {
            out.U(i, k) = s[j] > 0.0 ? U(i, j) / s[j] : 0.0;
        }
        for (int i = 0; i < n; ++i)
        {
            out.W(i, k) = W(i, j);
        }
    }
    return out;
}

} // namespace nekecs
