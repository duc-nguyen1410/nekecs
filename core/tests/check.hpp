// Minimal test helpers (no framework dependency).
#pragma once

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include <nekecs/dense.hpp>

namespace test
{

inline int &failures()
{
    static int n = 0;
    return n;
}

#define CHECK(cond, msg)                                                       \
    do                                                                         \
    {                                                                          \
        if (!(cond))                                                           \
        {                                                                      \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << ": " << msg \
                      << "  [" #cond "]" << std::endl;                         \
            ++test::failures();                                                \
        }                                                                      \
        else                                                                   \
        {                                                                      \
            std::cout << "ok   " << msg << std::endl;                          \
        }                                                                      \
    } while (0)

#define CHECK_NEAR(a, b, tol, msg)                                             \
    CHECK(std::abs((a) - (b)) <= (tol),                                        \
          msg << " (" << (a) << " vs " << (b) << ", tol " << (tol) << ")")

inline nekecs::Matrix random_matrix(int m, int n, std::mt19937 &rng)
{
    std::normal_distribution<double> N(0.0, 1.0);
    nekecs::Matrix A(m, n);
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A(i, j) = N(rng);
        }
    }
    return A;
}

inline std::vector<double> random_vector(int n, std::mt19937 &rng)
{
    std::normal_distribution<double> N(0.0, 1.0);
    std::vector<double> v(n);
    for (auto &x : v)
    {
        x = N(rng);
    }
    return v;
}

// Dense solve A x = b by Gaussian elimination with partial pivoting.
inline std::vector<double> dense_solve(nekecs::Matrix A, std::vector<double> b)
{
    const int n = A.rows();
    for (int c = 0; c < n; ++c)
    {
        int p = c;
        for (int i = c + 1; i < n; ++i)
        {
            if (std::abs(A(i, c)) > std::abs(A(p, c)))
            {
                p = i;
            }
        }
        for (int j = 0; j < n; ++j)
        {
            std::swap(A(c, j), A(p, j));
        }
        std::swap(b[c], b[p]);
        for (int i = c + 1; i < n; ++i)
        {
            const double f = A(i, c) / A(c, c);
            for (int j = c; j < n; ++j)
            {
                A(i, j) -= f * A(c, j);
            }
            b[i] -= f * b[c];
        }
    }
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i)
    {
        double s = b[i];
        for (int j = i + 1; j < n; ++j)
        {
            s -= A(i, j) * x[j];
        }
        x[i] = s / A(i, i);
    }
    return x;
}

inline int report(const char *name)
{
    if (failures() == 0)
    {
        std::cout << name << ": all checks passed" << std::endl;
    }
    else
    {
        std::cout << name << ": " << failures() << " check(s) FAILED"
                  << std::endl;
    }
    return failures() == 0 ? 0 : 1;
}

} // namespace test
