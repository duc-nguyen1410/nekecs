// Vector-space abstraction used by every algorithm in nekecs/core.
//
// The algorithms never touch vector internals. They are templated on a
// "Space" type that provides:
//
//   using Vector = ...;                                   // any type
//   Vector zeros_like(const Vector &x) const;             // NEW storage, zero-filled
//   void   copy(const Vector &src, Vector &dst) const;    // dst <- src (same shape)
//   double dot(const Vector &a, const Vector &b) const;   // inner product
//   void   axpy(double a, const Vector &x, Vector &y) const; // y <- y + a x
//   void   scale(double a, Vector &x) const;              // x <- a x
//
// Notes for backends:
//  - zeros_like must allocate: the algorithms keep Krylov bases of many vectors.
//    (Nektar++ Array<> copies are shallow, so a backend may not rely on the
//    Vector copy constructor making a deep copy; the core never does.)
//  - dot defines the geometry: GMRES orthogonality, trust-region radius and the
//    Newton residual norm all use it. For spectral/hp fields it must be the
//    quadrature-weighted L2 inner product, reduced over MPI ranks.
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace nekecs
{

template <class Space>
double norm(const Space &S, const typename Space::Vector &x)
{
    return std::sqrt(S.dot(x, x));
}

template <class Space>
typename Space::Vector clone(const Space &S, const typename Space::Vector &x)
{
    auto y = S.zeros_like(x);
    S.copy(x, y);
    return y;
}

// Plain Euclidean space over std::vector<double>; used by the tests and for
// small ODE problems.
struct StdVectorSpace
{
    using Vector = std::vector<double>;

    Vector zeros_like(const Vector &x) const
    {
        return Vector(x.size(), 0.0);
    }
    void copy(const Vector &src, Vector &dst) const
    {
        dst = src;
    }
    double dot(const Vector &a, const Vector &b) const
    {
        double s = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i)
        {
            s += a[i] * b[i];
        }
        return s;
    }
    void axpy(double a, const Vector &x, Vector &y) const
    {
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            y[i] += a * x[i];
        }
    }
    void scale(double a, Vector &x) const
    {
        for (auto &v : x)
        {
            v *= a;
        }
    }
};

} // namespace nekecs
