// Arnoldi eigensolver for the leading eigenvalues (largest magnitude) of a
// linear operator given only by its action, e.g. the flow-map Jacobian
// D Phi_T or D(sigma Phi_T) evaluated by finite differences.
//
// Unrestarted: the Krylov space grows until the `nev` leading Ritz pairs
// converge or `max_dim` is reached. For flow-map Jacobians the leading
// multipliers separate quickly from the strongly damped bulk, so this is
// usually enough; the basis costs max_dim vectors of memory.
//
// Ritz pair (theta, x = V_k y), with H_k y = theta y and ||y|| = 1:
//     ||A x - theta x|| = |h_{k+1,k}| |y_{k-1}|        (exact in exact arithmetic)
// Converged when that residual <= tol * |theta|.
//
// For a flow map over time T, multipliers mu relate to growth rates by
// lambda = log(mu) / T; the imaginary part is only defined modulo 2 pi / T.
#pragma once

#include "eig.hpp"
#include "space.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <iomanip>
#include <ostream>
#include <vector>

namespace nekecs
{

template <class Space> struct EigenResult
{
    using Vector = typename Space::Vector;

    std::vector<cplx> values;   // leading eigenvalues, |.| descending
    std::vector<double> residuals; // ||A x - theta x|| / |theta|
    std::vector<Vector> re, im; // Ritz vectors: x = re + i im, ||x|| = 1
    int k          = 0;         // Krylov dimension used
    bool converged = false;
};

struct ArnoldiOptions
{
    int nev         = 4;    // eigenvalues wanted
    int max_dim     = 100;  // Krylov dimension limit
    double tol      = 1e-6; // relative Ritz residual
    int check_every = 5;    // test convergence every this many iterations
    std::ostream *log = nullptr;
};

// Op: callable (const Vector &v, Vector &Av).
template <class Space, class Op>
EigenResult<Space> arnoldi(const Space &S, const Op &A,
                           const typename Space::Vector &v0,
                           const ArnoldiOptions &opt)
{
    using Vector = typename Space::Vector;
    EigenResult<Space> R;

    std::vector<Vector> V;
    Vector v = clone(S, v0);
    S.scale(1.0 / norm(S, v), v);
    V.push_back(std::move(v));
    Matrix H(opt.max_dim + 1, opt.max_dim);

    struct Ritz
    {
        cplx theta;
        double res;
        std::vector<cplx> y;
    };
    auto ritz = [&](int k) {
        std::vector<cplx> th = hessenberg_eigenvalues(H, k);
        std::sort(th.begin(), th.end(), [](cplx a, cplx b) {
            return std::abs(a) > std::abs(b);
        });
        std::vector<Ritz> out;
        const int want = std::min<int>(opt.nev, k);
        for (int i = 0; i < want; ++i)
        {
            Ritz r{th[i], 0.0, hessenberg_eigenvector(H, k, th[i])};
            r.res = std::abs(H(k, k - 1)) * std::abs(r.y[k - 1]) /
                    std::max(std::abs(th[i]), 1e-300);
            out.push_back(std::move(r));
        }
        return out;
    };

    std::vector<Ritz> pairs;
    int k = 0;
    for (int j = 0; j < opt.max_dim; ++j)
    {
        Vector w = S.zeros_like(V[0]);
        A(V[j], w);
        for (int pass = 0; pass < 2; ++pass)
        {
            for (int i = 0; i <= j; ++i)
            {
                const double h = S.dot(V[i], w);
                H(i, j) += h;
                S.axpy(-h, V[i], w);
            }
        }
        const double hn = norm(S, w);
        H(j + 1, j)     = hn;
        k               = j + 1;
        const bool breakdown = hn <= 1e-14 * std::abs(H(j, j)) || hn == 0.0;
        if (!breakdown)
        {
            S.scale(1.0 / hn, w);
            V.push_back(std::move(w));
        }

        if (breakdown || k % opt.check_every == 0 || k == opt.max_dim)
        {
            pairs = ritz(k);
            bool ok = static_cast<int>(pairs.size()) == std::min(opt.nev, k);
            for (const auto &p : pairs)
            {
                ok &= p.res <= opt.tol || breakdown;
            }
            if (opt.log)
            {
                *opt.log << "[arnoldi] k = " << k;
                for (const auto &p : pairs)
                {
                    *opt.log << std::scientific << std::setprecision(4)
                             << "  (" << p.theta.real() << ", "
                             << p.theta.imag() << ") res " << std::setprecision(1)
                             << p.res;
                }
                *opt.log << std::endl;
            }
            if (ok && k >= opt.nev)
            {
                R.converged = true;
                break;
            }
            if (breakdown)
            {
                break;
            }
        }
    }

    // Ritz vectors x = V_k y.
    R.k = k;
    for (const auto &p : pairs)
    {
        Vector xr = S.zeros_like(V[0]), xi = S.zeros_like(V[0]);
        for (int i = 0; i < k; ++i)
        {
            S.axpy(p.y[i].real(), V[i], xr);
            S.axpy(p.y[i].imag(), V[i], xi);
        }
        R.values.push_back(p.theta);
        R.residuals.push_back(p.res);
        R.re.push_back(std::move(xr));
        R.im.push_back(std::move(xi));
    }
    return R;
}

} // namespace nekecs
