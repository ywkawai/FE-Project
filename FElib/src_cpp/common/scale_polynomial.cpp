#include "scale_polynomial.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

// LAPACK's symmetric tridiagonal eigensolver, as linked via SCALE_MATHLIB_LIBS. 
extern "C" {
void dstev_(const char* jobz, const int* n, double* d, double* e, double* z, const int* ldz,
            double* work, int* info);
}

namespace FElib {
namespace common {
namespace polynomial {

namespace {

void RequireOrderAtLeast(int order, int minOrder, const char* funcName)
{
    if (order < minOrder) {
        throw std::invalid_argument(std::string(funcName) + ": order must be >= " +
                                     std::to_string(minOrder));
    }
}

template <typename View>
void RequireExtent1(const char* funcName, const View& v, std::size_t n0)
{
    if (v.extent(0) != n0) {
        throw std::invalid_argument(std::string(funcName) + ": array shape mismatch");
    }
}

template <typename View>
void RequireExtent2(const char* funcName, const View& v, std::size_t n0, std::size_t n1)
{
    if (v.extent(0) != n0 || v.extent(1) != n1) {
        throw std::invalid_argument(std::string(funcName) + ": array shape mismatch");
    }
}

//! Golub-Welsch quadrature points: the N+1 eigenvalues of the Jacobi matrix for parameters (alpha, beta). 
void GenJacobiGaussQuadraturePts(int alpha, int beta, int N, ArrayView<double, 1> x)
{
    if (N < 0) {
        throw std::invalid_argument("GenJacobiGaussQuadraturePts: N must be non-negative");
    }

    const auto npts = static_cast<std::size_t>(N) + 1;
    if (x.size() != npts) {
        throw std::invalid_argument("GenJacobiGaussQuadraturePts: output array size mismatch");
    }

    const double a = static_cast<double>(alpha);
    const double b = static_cast<double>(beta);

    if (N == 0) {
        x(0) = -(a - b) / (a + b + 2);
        return;
    }

    std::vector<double> h1(npts);
    for (std::size_t i = 0; i < npts; ++i) {
        h1[i] = 2.0 * static_cast<double>(i) + a + b;
    }

    std::vector<double> d(npts);
    for (std::size_t i = 0; i < npts; ++i) {
        d[i] = -(a * a - b * b) / (h1[i] * (h1[i] + 2.0));
    }

    std::vector<double> e(npts - 1);
    for (std::size_t i = 1; i < npts; ++i) {
        const double ii = static_cast<double>(i);
        const double hi = h1[i - 1];
        e[i - 1] = 2.0 / (hi + 2.0) *
            std::sqrt(ii * (ii + a + b) * (ii + a) * (ii + b) / ((hi + 1.0) * (hi + 3.0)));
    }

    if (alpha + beta == 0) d[0] = 0.0;

    const int n = N + 1;
    int ldz = 1;
    int info = 0;
    std::vector<double> work(static_cast<std::size_t>(std::max(1, 2 * n - 2)));    
    dstev_("N", &n, d.data(), e.data(), /*z=*/nullptr, &ldz, work.data(), &info);
    if (info != 0) {
        throw std::runtime_error("GenJacobiGaussQuadraturePts: dstev failed with info=" + std::to_string(info));
    }

    for (std::size_t i = 0; i < npts; ++i) {
        x(i) = d[i];
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// Legendre polynomials
// ---------------------------------------------------------------------------

void GenLegendrePoly(int order, ArrayView<const double, 1> x, ArrayView<double, 2> P)
{
    RequireOrderAtLeast(order, 0, "GenLegendrePoly");
    const std::size_t nx = x.extent(0);
    RequireExtent2("GenLegendrePoly", P, nx, static_cast<std::size_t>(order) + 1);

    for (std::size_t i = 0; i < nx; ++i) P(i, 0) = 1.0;
    if (order == 0) return;

    for (std::size_t i = 0; i < nx; ++i) P(i, 1) = x(i);
    for (int degree = 2; degree <= order; ++degree) {
        const double d = static_cast<double>(degree);
        const auto k = static_cast<std::size_t>(degree);
        for (std::size_t i = 0; i < nx; ++i) {
            P(i, k) = ((2.0 * d - 1.0) * x(i) * P(i, k - 1) - (d - 1.0) * P(i, k - 2)) / d;
        }
    }
}

std::vector<double> GenLegendrePoly(int order, ArrayView<const double, 1> x)
{
    RequireOrderAtLeast(order, 0, "GenLegendrePoly");
    const std::size_t nx = x.extent(0);
    const auto ncol = static_cast<std::size_t>(order) + 1;
    std::vector<double> buf(nx * ncol);
    GenLegendrePoly(order, x, ArrayView<double, 2>(buf.data(), nx, ncol));
    return buf;
}

void GenDLegendrePoly(int order, ArrayView<const double, 1> x, ArrayView<const double, 2> P,
                       ArrayView<double, 2> gradP)
{
    RequireOrderAtLeast(order, 0, "GenDLegendrePoly");
    const std::size_t nx = x.extent(0);
    const auto ncol = static_cast<std::size_t>(order) + 1;
    RequireExtent2("GenDLegendrePoly", P, nx, ncol);
    RequireExtent2("GenDLegendrePoly", gradP, nx, ncol);

    for (std::size_t i = 0; i < nx; ++i) {
        const double xi = x(i);

        gradP(i, 0) = 0.0;
        if (order >= 1) gradP(i, 1) = 1.0;

        for (int degree = 2; degree <= order; ++degree) {
            const auto k = static_cast<std::size_t>(degree);
            gradP(i, k) = 2.0 * xi * gradP(i, k - 1) - gradP(i, k - 2) + P(i, k - 1);
        }
    }
}

std::vector<double> GenDLegendrePoly(int order, ArrayView<const double, 1> x, ArrayView<const double, 2> P)
{
    RequireOrderAtLeast(order, 0, "GenDLegendrePoly");
    const std::size_t nx = x.extent(0);
    const auto ncol = static_cast<std::size_t>(order) + 1;
    std::vector<double> buf(nx * ncol);
    GenDLegendrePoly(order, x, P, ArrayView<double, 2>(buf.data(), nx, ncol));
    return buf;
}

// ---------------------------------------------------------------------------
// Gauss-Lobatto-Legendre (GLL) points and weights
// ---------------------------------------------------------------------------

void GenGaussLobattoPt(int order, ArrayView<double, 1> pts)
{
    RequireOrderAtLeast(order, 1, "GenGaussLobattoPt");
    RequireExtent1("GenGaussLobattoPt", pts, static_cast<std::size_t>(order) + 1);

    pts(0) = -1.0;
    pts(static_cast<std::size_t>(order)) = 1.0;
    if (order == 1) return;

    ArrayView<double, 1> interior(pts.data() + 1, static_cast<std::size_t>(order - 1));
    GenJacobiGaussQuadraturePts(1, 1, order - 2, interior);
}

std::vector<double> GenGaussLobattoPt(int order)
{
    RequireOrderAtLeast(order, 1, "GenGaussLobattoPt");
    std::vector<double> buf(static_cast<std::size_t>(order) + 1);
    GenGaussLobattoPt(order, ArrayView<double, 1>(buf.data(), buf.size()));
    return buf;
}

void GenGaussLobattoPtIntWeight(int order, ArrayView<double, 1> weights)
{
    RequireOrderAtLeast(order, 1, "GenGaussLobattoPtIntWeight");
    const auto npts = static_cast<std::size_t>(order) + 1;
    RequireExtent1("GenGaussLobattoPtIntWeight", weights, npts);

    const std::vector<double> lglPts = GenGaussLobattoPt(order);
    const std::vector<double> P = GenLegendrePoly(order, ArrayView<const double, 1>(lglPts.data(), npts));
    const ArrayView<const double, 2> Pview(P.data(), npts, npts);

    const double denomCoeff = static_cast<double>(order * (order + 1));
    for (std::size_t i = 0; i < npts; ++i) {
        const double p = Pview(i, static_cast<std::size_t>(order));
        weights(i) = 2.0 / (denomCoeff * p * p);
    }
}

std::vector<double> GenGaussLobattoPtIntWeight(int order)
{
    RequireOrderAtLeast(order, 1, "GenGaussLobattoPtIntWeight");
    std::vector<double> buf(static_cast<std::size_t>(order) + 1);
    GenGaussLobattoPtIntWeight(order, ArrayView<double, 1>(buf.data(), buf.size()));
    return buf;
}

// ---------------------------------------------------------------------------
// Gauss-Legendre (GL) points and weights
// ---------------------------------------------------------------------------

void GenGaussLegendrePt(int order, ArrayView<double, 1> pts)
{
    RequireOrderAtLeast(order, 1, "GenGaussLegendrePt");
    RequireExtent1("GenGaussLegendrePt", pts, static_cast<std::size_t>(order));

    GenJacobiGaussQuadraturePts(0, 0, order - 1, pts);
}

std::vector<double> GenGaussLegendrePt(int order)
{
    RequireOrderAtLeast(order, 1, "GenGaussLegendrePt");
    std::vector<double> buf(static_cast<std::size_t>(order));
    GenGaussLegendrePt(order, ArrayView<double, 1>(buf.data(), buf.size()));
    return buf;
}

void GenGaussLegendrePtIntWeight(int order, ArrayView<double, 1> weights)
{
    RequireOrderAtLeast(order, 1, "GenGaussLegendrePtIntWeight");
    const auto npts = static_cast<std::size_t>(order);
    RequireExtent1("GenGaussLegendrePtIntWeight", weights, npts);

    const std::vector<double> glPts = GenGaussLegendrePt(order);
    ArrayView<const double, 1> glPtsView(glPts.data(), npts);
    const std::vector<double> P = GenLegendrePoly(order, glPtsView);
    ArrayView<const double, 2> Pview(P.data(), npts, static_cast<std::size_t>(order) + 1);
    const std::vector<double> dP = GenDLegendrePoly(order, glPtsView, Pview);
    ArrayView<const double, 2> dPview(dP.data(), npts, static_cast<std::size_t>(order) + 1);

    for (std::size_t i = 0; i < npts; ++i) {
        const double xi = glPtsView(i);
        const double dp = dPview(i, static_cast<std::size_t>(order));
        weights(i) = 2.0 / ((1.0 - xi * xi) * dp * dp);
    }
}

std::vector<double> GenGaussLegendrePtIntWeight(int order)
{
    RequireOrderAtLeast(order, 1, "GenGaussLegendrePtIntWeight");
    std::vector<double> buf(static_cast<std::size_t>(order));
    GenGaussLegendrePtIntWeight(order, ArrayView<double, 1>(buf.data(), buf.size()));
    return buf;
}

// ---------------------------------------------------------------------------
// Lagrange basis functions on GLL points
// ---------------------------------------------------------------------------

void GenLagrangePoly(int order, ArrayView<const double, 1> xLgl, ArrayView<const double, 1> x,
                      ArrayView<double, 2> L)
{
    // order >= 1 is required: the formula below divides by order*(order+1),
    // which is identically 0 for order == 0.
    RequireOrderAtLeast(order, 1, "GenLagrangePoly");
    const auto ncol = static_cast<std::size_t>(order) + 1;
    const std::size_t nx = x.extent(0);
    RequireExtent1("GenLagrangePoly", xLgl, ncol);
    RequireExtent2("GenLagrangePoly", L, nx, ncol);

    const std::vector<double> pLgl = GenLegendrePoly(order, xLgl);
    ArrayView<const double, 2> pLglView(pLgl.data(), ncol, ncol);

    const std::vector<double> p = GenLegendrePoly(order, x);
    ArrayView<const double, 2> pView(p.data(), nx, ncol);
    const std::vector<double> pr = GenDLegendrePoly(order, x, pView);
    ArrayView<const double, 2> prView(pr.data(), nx, ncol);

    const double denomCoeff = static_cast<double>(order) * static_cast<double>(order + 1);

    for (std::size_t n = 0; n < ncol; ++n) {
        for (std::size_t i = 0; i < nx; ++i) {
            if (std::abs(x(i) - xLgl(n)) < 1e-15) {
                L(i, n) = 1.0;
            } else {
                L(i, n) = (x(i) - 1.0) * (x(i) + 1.0) * prView(i, static_cast<std::size_t>(order)) /
                          (denomCoeff * pLglView(n, static_cast<std::size_t>(order)) * (x(i) - xLgl(n)));
            }
        }
    }
}

std::vector<double> GenLagrangePoly(int order, ArrayView<const double, 1> xLgl,
                                     ArrayView<const double, 1> x)
{
    RequireOrderAtLeast(order, 1, "GenLagrangePoly");
    const std::size_t nx = x.extent(0);
    const auto ncol = static_cast<std::size_t>(order) + 1;
    std::vector<double> buf(nx * ncol);
    GenLagrangePoly(order, xLgl, x, ArrayView<double, 2>(buf.data(), nx, ncol));
    return buf;
}

void GenDLagrangePoly_lglpt(int order, ArrayView<const double, 1> xLgl, ArrayView<double, 2> Lr)
{
    RequireOrderAtLeast(order, 0, "GenDLagrangePoly_lglpt");
    const auto ncol = static_cast<std::size_t>(order) + 1;
    RequireExtent1("GenDLagrangePoly_lglpt", xLgl, ncol);
    RequireExtent2("GenDLagrangePoly_lglpt", Lr, ncol, ncol);

    const std::vector<double> p = GenLegendrePoly(order, xLgl);
    ArrayView<const double, 2> pView(p.data(), ncol, ncol);
    const auto lastCol = static_cast<std::size_t>(order);

    for (std::size_t n = 0; n < ncol; ++n) {
        double lrnn = 0.0;
        for (std::size_t k = 0; k < ncol; ++k) {
            double v;
            if (k == 0 && n == 0) {
                v = -0.25 * static_cast<double>(order) * static_cast<double>(order + 1);
            } else if (k == lastCol && n == lastCol) {
                v = 0.25 * static_cast<double>(order) * static_cast<double>(order + 1);
            } else if (k == n) {
                v = 0.0;
            } else {
                v = pView(n, lastCol) / (pView(k, lastCol) * (xLgl(n) - xLgl(k)));
            }
            Lr(k, n) = v;
            if (k != n) lrnn += v;
        }
        Lr(n, n) = -lrnn;
    }
}

std::vector<double> GenDLagrangePoly_lglpt(int order, ArrayView<const double, 1> xLgl)
{
    RequireOrderAtLeast(order, 0, "GenDLagrangePoly_lglpt");
    const auto ncol = static_cast<std::size_t>(order) + 1;
    std::vector<double> buf(ncol * ncol);
    GenDLagrangePoly_lglpt(order, xLgl, ArrayView<double, 2>(buf.data(), ncol, ncol));
    return buf;
}

}  // namespace polynomial
}  // namespace common
}  // namespace FElib
