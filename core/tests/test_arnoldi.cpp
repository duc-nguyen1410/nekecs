// Dense Hessenberg eigensolver and Arnoldi, against known spectra:
//  1. random Hessenberg matrices: H y = theta y residuals, trace and det;
//  2. Arnoldi on A = Q D Q^T with prescribed leading eigenvalues;
//  3. Lorenz equilibrium C+: multipliers exp(lambda T) of the flow-map
//     Jacobian (FD) vs the roots of the characteristic cubic;
//  4. Lorenz AB periodic orbit: Floquet multipliers contain 1 (time
//     translation) and their product is exp(-(sigma + 1 + beta) T)
//     (Liouville: the flow contracts volume at a constant rate).
#include "check.hpp"
#include "lorenz.hpp"

#include <nekecs/arnoldi.hpp>
#include <nekecs/newton.hpp>
#include <nekecs/space.hpp>

#include <algorithm>

using namespace nekecs;
using Vec = std::vector<double>;

static void test_dense()
{
    std::mt19937 rng(7);
    double worst = 0.0, trace_err = 0.0, det_err = 0.0;
    for (int trial = 0; trial < 20; ++trial)
    {
        const int n = 3 + trial; // 3 .. 22
        Matrix R = test::random_matrix(n, n, rng);
        Matrix H(n + 1, n);
        for (int i = 0; i < n; ++i)
            for (int j = std::max(0, i - 1); j < n; ++j)
                H(i, j) = R(i, j);

        const auto th = hessenberg_eigenvalues(H, n);
        cplx tr = 0.0, prod = 1.0;
        for (auto t : th)
        {
            tr += t;
            prod *= t;
            const auto y = hessenberg_eigenvector(H, n, t);
            double r = 0.0;
            for (int i = 0; i < n; ++i)
            {
                cplx s = -t * y[i];
                for (int j = 0; j < n; ++j) s += H(i, j) * y[j];
                r += std::norm(s);
            }
            worst = std::max(worst, std::sqrt(r));
        }
        double trH = 0.0;
        for (int i = 0; i < n; ++i) trH += H(i, i);
        // det via the same LU-free route: product of eigenvalues vs dense det
        Matrix A = H.block(n, n);
        double det = 1.0;
        {
            Matrix M = A;
            for (int c = 0; c < n; ++c)
            {
                int p = c;
                for (int i = c + 1; i < n; ++i)
                    if (std::abs(M(i, c)) > std::abs(M(p, c))) p = i;
                if (p != c)
                {
                    for (int j = 0; j < n; ++j) std::swap(M(c, j), M(p, j));
                    det = -det;
                }
                det *= M(c, c);
                for (int i = c + 1; i < n; ++i)
                {
                    const double f = M(i, c) / M(c, c);
                    for (int j = c; j < n; ++j) M(i, j) -= f * M(c, j);
                }
            }
        }
        trace_err = std::max(trace_err, std::abs(tr - trH) / (1.0 + std::abs(trH)));
        det_err   = std::max(det_err, std::abs(prod - det) / (1.0 + std::abs(det)));
    }
    // Inverse iteration near clustered eigenvalues: ~1e-9 on these matrices.
    CHECK(worst < 1e-8, "dense eig: max ||H y - theta y|| = " << worst);
    CHECK(trace_err < 1e-12, "dense eig: sum of eigenvalues = trace (" << trace_err << ")");
    CHECK(det_err < 1e-10, "dense eig: product of eigenvalues = det (" << det_err << ")");
}

static void test_arnoldi_matrix()
{
    // A = Q D Q^T, D block diagonal: 0.9 +- 0.3i, 0.8, 0.5, rest in [-0.3, 0.3].
    std::mt19937 rng(11);
    const int n = 80;
    Matrix D(n, n);
    D(0, 0) = 0.9;  D(0, 1) = 0.3;
    D(1, 0) = -0.3; D(1, 1) = 0.9;
    D(2, 2) = 0.8;
    D(3, 3) = 0.5;
    std::uniform_real_distribution<double> U(-0.3, 0.3);
    for (int i = 4; i < n; ++i) D(i, i) = U(rng);
    // Random orthogonal Q by Gram-Schmidt.
    Matrix Q = test::random_matrix(n, n, rng);
    for (int j = 0; j < n; ++j)
    {
        for (int pass = 0; pass < 2; ++pass)
            for (int k = 0; k < j; ++k)
            {
                double d = 0.0;
                for (int i = 0; i < n; ++i) d += Q(i, k) * Q(i, j);
                for (int i = 0; i < n; ++i) Q(i, j) -= d * Q(i, k);
            }
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += Q(i, j) * Q(i, j);
        s = std::sqrt(s);
        for (int i = 0; i < n; ++i) Q(i, j) /= s;
    }
    Matrix A(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        {
            double s = 0.0;
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b) s += Q(i, a) * D(a, b) * Q(j, b);
            A(i, j) = s;
        }

    StdVectorSpace S;
    auto op = [&](const Vec &v, Vec &out) { out = A.multiply(v); };
    ArnoldiOptions opt;
    opt.nev = 4;
    opt.tol = 1e-10;
    opt.max_dim = 80;
    const auto R = arnoldi(S, op, test::random_vector(n, rng), opt);

    CHECK(R.converged, "Arnoldi (80x80) converged at k = " << R.k);
    const cplx expect[4] = {{0.9, 0.3}, {0.9, -0.3}, {0.8, 0.0}, {0.5, 0.0}};
    double err = 0.0;
    for (int i = 0; i < 4 && i < (int)R.values.size(); ++i)
    {
        double best = 1e300;
        for (auto e : expect) best = std::min(best, std::abs(R.values[i] - e));
        err = std::max(err, best);
    }
    CHECK(R.values.size() == 4 && err < 1e-9,
          "Arnoldi finds 0.9+-0.3i, 0.8, 0.5 (max err " << err << ")");

    // Ritz vector check: ||A x - theta x|| for the leading pair.
    if (!R.values.empty())
    {
        const auto &xr = R.re[0], &xi = R.im[0];
        const auto Ar = A.multiply(xr), Ai = A.multiply(xi);
        const cplx t = R.values[0];
        double r = 0.0;
        for (int i = 0; i < n; ++i)
        {
            const cplx lhs(Ar[i], Ai[i]), x(xr[i], xi[i]);
            r += std::norm(lhs - t * x);
        }
        CHECK(std::sqrt(r) < 1e-9, "Ritz vector residual " << std::sqrt(r));
    }
}

// Jacobian of the Lorenz flow map at x, exact via the tangent-linear flow.
static auto flow_jacobian(const std::array<double, 3> &x, double T, int N)
{
    return [x, T, N](const Vec &v, Vec &out) {
        auto w = lorenz::flow_tangent(x, {v[0], v[1], v[2]}, T, N);
        out.assign(w.begin(), w.end());
    };
}

static void test_lorenz_equilibrium()
{
    using namespace lorenz;
    // Roots of lambda^3 + a lambda^2 + b lambda + c at C+.
    const double a = sigma + beta + 1.0, b = beta * (sigma + rho),
                 c = 2.0 * sigma * beta * (rho - 1.0);
    double r = -14.0;
    for (int i = 0; i < 50; ++i)
        r -= (((r + a) * r + b) * r + c) / ((3.0 * r + 2.0 * a) * r + b);
    const double p = a + r, q = b + r * p; // (lambda - r)(lambda^2 + p lambda + q)
    const cplx l1(-0.5 * p, 0.5 * std::sqrt(4.0 * q - p * p));
    std::cout << "     C+ eigenvalues: " << r << ", " << l1.real() << " +- "
              << l1.imag() << "i" << std::endl;

    const double T = 0.5;
    const int N    = 2000;
    const double e = std::sqrt(beta * (rho - 1.0));
    StdVectorSpace S;
    ArnoldiOptions opt;
    opt.nev = 3;
    opt.tol = 1e-8;
    opt.max_dim = 3;
    opt.check_every = 1;
    const auto R = arnoldi(S, flow_jacobian({e, e, rho - 1.0}, T, N), Vec{1.0, 0.3, -0.2}, opt);

    const cplx m1 = std::exp(l1 * T), m2 = std::conj(m1), m3 = std::exp(r * T);
    double err = 0.0;
    for (auto mu : R.values)
    {
        err = std::max(err, std::min({std::abs(mu - m1), std::abs(mu - m2),
                                      std::abs(mu - m3)}));
    }
    CHECK(R.values.size() == 3 && err < 1e-6,
          "Lorenz C+ multipliers exp(lambda T) (max err " << err << ")");
    if (!R.values.empty())
    {
        // log(mu)/T fixes the frequency only modulo 2 pi / T (here
        // omega T = 5.1 > pi, so the principal branch aliases).
        const cplx lam = std::log(R.values[0]) / T;
        const double w = 2.0 * M_PI / T;
        auto alias = [w](double a, double b) {
            const double d = std::remainder(a - b, w);
            return std::abs(d);
        };
        const bool ok = std::abs(lam.real() - l1.real()) < 1e-9 &&
                        (alias(lam.imag(), l1.imag()) < 1e-8 ||
                         alias(lam.imag(), -l1.imag()) < 1e-8);
        CHECK(ok, "leading growth rate log(mu)/T = " << lam
                                                     << " (= C+ pair mod 2pi/T)");
    }
}

static void test_lorenz_floquet()
{
    using namespace lorenz;
    // Converge the AB orbit first (as in test_lorenz).
    const Guess g = recurrence_guess(1.3, 1.8);
    const int N   = 4000;
    const auto x0 = g.x;
    const auto f0 = f(x0);
    auto F = [&](const Vec &x, Vec &Fx) {
        auto y = flow({x[0], x[1], x[2]}, x[3], N);
        for (int i = 0; i < 3; ++i) Fx[i] = y[i] - x[i];
        Fx[3] = f0[0] * (x[0] - x0[0]) + f0[1] * (x[1] - x0[1]) + f0[2] * (x[2] - x0[2]);
    };
    StdVectorSpace S;
    NewtonOptions<StdVectorSpace> nopt;
    nopt.tol = 1e-11;
    const auto P = newton_hookstep(S, F, Vec{x0[0], x0[1], x0[2], g.T}, nopt);
    CHECK(P.converged(), "AB orbit converged (T = " << P.x[3] << ")");

    ArnoldiOptions opt;
    opt.nev = 3;
    opt.tol = 1e-8;
    opt.max_dim = 3;
    opt.check_every = 1;
    const double T = P.x[3];
    const auto R = arnoldi(S, flow_jacobian({P.x[0], P.x[1], P.x[2]}, T, N),
                           Vec{0.2, 1.0, 0.5}, opt);
    cplx prod = 1.0;
    double dist1 = 1e300;
    std::cout << "     AB Floquet multipliers:";
    for (auto mu : R.values)
    {
        std::cout << " " << mu;
        prod *= mu;
        dist1 = std::min(dist1, std::abs(mu - 1.0));
    }
    std::cout << std::endl;
    const double liouville = std::exp(-(sigma + 1.0 + beta) * T);
    CHECK(R.values.size() == 3 && dist1 < 1e-6,
          "one Floquet multiplier = 1 (|mu - 1| = " << dist1 << ")");
    // The smallest multiplier (~1e-10) has absolute accuracy ~eps * 4.7,
    // i.e. ~1e-5 relative, which bounds the product's accuracy.
    CHECK(std::abs(prod - liouville) < 1e-4 * liouville,
          "product of multipliers = exp(-(sigma+1+beta)T) = " << liouville
              << " (got " << prod.real() << ")");
}

int main()
{
    test_dense();
    test_arnoldi_matrix();
    test_lorenz_equilibrium();
    test_lorenz_floquet();
    return test::report("test_arnoldi");
}
