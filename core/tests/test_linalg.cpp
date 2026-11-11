// SVD, GMRES and hookstep against dense reference computations.
#include "check.hpp"

#include <nekecs/gmres.hpp>
#include <nekecs/hookstep.hpp>
#include <nekecs/space.hpp>

#include <algorithm>
#include <limits>

using namespace nekecs;

static void test_svd()
{
    std::mt19937 rng(1);
    const int m = 9, n = 6;
    Matrix A = test::random_matrix(m, n, rng);
    // Make it ill-conditioned: scale columns over 10 orders of magnitude.
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < m; ++i)
        {
            A(i, j) *= std::pow(10.0, -2.0 * j);
        }
    }
    const SVD d = svd(A);

    double recon = 0.0, orthU = 0.0, orthW = 0.0;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double s = 0.0;
            for (int k = 0; k < n; ++k)
            {
                s += d.U(i, k) * d.s[k] * d.W(j, k);
            }
            recon = std::max(recon, std::abs(s - A(i, j)));
        }
    }
    for (int a = 0; a < n; ++a)
    {
        for (int b = 0; b < n; ++b)
        {
            double u = 0.0, w = 0.0;
            for (int i = 0; i < m; ++i)
            {
                u += d.U(i, a) * d.U(i, b);
            }
            for (int i = 0; i < n; ++i)
            {
                w += d.W(i, a) * d.W(i, b);
            }
            orthU = std::max(orthU, std::abs(u - (a == b)));
            orthW = std::max(orthW, std::abs(w - (a == b)));
        }
    }
    bool sorted = std::is_sorted(d.s.rbegin(), d.s.rend());
    CHECK(recon < 1e-13, "svd reconstructs A (max err " << recon << ")");
    CHECK(orthU < 1e-12, "svd U orthonormal (" << orthU << ")");
    CHECK(orthW < 1e-12, "svd W orthogonal (" << orthW << ")");
    CHECK(sorted, "svd singular values sorted descending");
    CHECK(d.s[n - 1] > 0.0 && d.s[0] / d.s[n - 1] > 1e8,
          "svd resolves condition number " << d.s[0] / d.s[n - 1]);
}

static void test_gmres()
{
    std::mt19937 rng(2);
    const int n = 40;
    Matrix A = test::random_matrix(n, n, rng);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A(i, j) = 0.3 * A(i, j) / std::sqrt(double(n)) + (i == j ? 1.0 : 0.0);
        }
    }
    const auto b = test::random_vector(n, rng);
    const auto xref = test::dense_solve(A, b);

    StdVectorSpace S;
    auto op = [&](const std::vector<double> &v, std::vector<double> &out) {
        out = A.multiply(v);
    };

    const auto K = gmres(S, op, b, n, 1e-12);
    const auto x = expand(S, K, K.y);
    double err = 0.0;
    for (int i = 0; i < n; ++i)
    {
        err = std::max(err, std::abs(x[i] - xref[i]));
    }
    CHECK(K.converged, "gmres converged in k = " << K.k << " iterations");
    CHECK(err < 1e-10, "gmres matches dense solve (max err " << err << ")");

    // True residual agrees with the projected estimate.
    auto Ax = A.multiply(x);
    double rtrue = 0.0;
    for (int i = 0; i < n; ++i)
    {
        rtrue += (b[i] - Ax[i]) * (b[i] - Ax[i]);
    }
    rtrue = std::sqrt(rtrue);
    CHECK(std::abs(rtrue - K.residual) < 1e-10 * K.beta,
          "gmres residual estimate " << K.residual << " ~ true " << rtrue);

    // Arnoldi relation A V_k = V_{k+1} H and orthonormal basis, on a
    // truncated run (k = 12) like a Newton step would use.
    const auto K2 = gmres(S, op, b, 12, 0.0);
    double arn = 0.0, orth = 0.0;
    for (int j = 0; j < K2.k; ++j)
    {
        auto lhs = A.multiply(K2.basis[j]);
        for (int i = 0; i <= K2.k; ++i)
        {
            S.axpy(-K2.H(i, j), K2.basis[i], lhs);
        }
        arn = std::max(arn, norm(S, lhs));
    }
    for (int a = 0; a <= K2.k; ++a)
    {
        for (int c = 0; c <= K2.k; ++c)
        {
            orth = std::max(orth, std::abs(S.dot(K2.basis[a], K2.basis[c]) -
                                           (a == c)));
        }
    }
    CHECK(K2.k == 12 && K2.H.rows() == 13 && K2.H.cols() == 12,
          "truncated gmres returns 13 x 12 Hessenberg");
    CHECK(arn < 1e-12, "Arnoldi relation holds (" << arn << ")");
    CHECK(orth < 1e-13, "Krylov basis orthonormal (" << orth << ")");

    // Hookstep with no constraint reproduces the GMRES least-squares step.
    const HookstepSolver hook(K2.H, K2.beta);
    const auto h = hook.solve(std::numeric_limits<double>::infinity());
    double dy = 0.0;
    for (int i = 0; i < K2.k; ++i)
    {
        dy = std::max(dy, std::abs(h.y[i] - K2.y[i]));
    }
    CHECK(!h.constrained && dy < 1e-10,
          "unconstrained hookstep == GMRES solution (" << dy << ")");
    CHECK(std::abs(h.predicted_residual - K2.residual) < 1e-10 * K2.beta,
          "hookstep predicted residual == GMRES residual");
}

static void test_hookstep()
{
    std::mt19937 rng(3);
    const int k = 10;
    Matrix H(k + 1, k);
    Matrix R = test::random_matrix(k + 1, k, rng);
    for (int i = 0; i <= k; ++i)
    {
        for (int j = std::max(0, i - 1); j < k; ++j)
        {
            H(i, j) = R(i, j);
        }
    }
    const double beta = 2.5;
    const HookstepSolver hook(H, beta);
    const double full = hook.newton_step_norm();

    auto direct_residual = [&](const std::vector<double> &y) {
        auto Hy = H.multiply(y);
        double s = 0.0;
        for (int i = 0; i <= k; ++i)
        {
            const double r = (i == 0 ? beta : 0.0) - Hy[i];
            s += r * r;
        }
        return std::sqrt(s);
    };

    double prev_res = std::numeric_limits<double>::infinity();
    bool monotone = true, on_boundary = true, pred_ok = true;
    for (double frac : {0.01, 0.1, 0.3, 0.6, 0.9})
    {
        const auto h = hook.solve(frac * full);
        on_boundary &= h.constrained &&
                       std::abs(norm2(h.y) - frac * full) <= 1e-8 * full;
        pred_ok &= std::abs(direct_residual(h.y) - h.predicted_residual) <
                   1e-12 * beta;
        monotone &= h.predicted_residual <= prev_res;
        prev_res = h.predicted_residual;
    }
    CHECK(on_boundary, "constrained hookstep lies on ||y|| = delta");
    CHECK(pred_ok, "predicted residual equals ||beta e0 - H y||");
    CHECK(monotone, "model residual decreases as delta grows");

    const auto h = hook.solve(2.0 * full);
    CHECK(!h.constrained && std::abs(norm2(h.y) - full) < 1e-12 * full,
          "delta above Newton step gives Newton step");
}

int main()
{
    test_svd();
    test_gmres();
    test_hookstep();
    return test::report("test_linalg");
}
