// Unit test of scale_polynomial

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "cpp_test_harness.hpp"
#include "scale_polynomial.hpp"

using FElib::common::polynomial::GenDLagrangePoly_lglpt;
using FElib::common::polynomial::GenDLegendrePoly;
using FElib::common::polynomial::GenGaussLegendrePt;
using FElib::common::polynomial::GenGaussLegendrePtIntWeight;
using FElib::common::polynomial::GenGaussLobattoPt;
using FElib::common::polynomial::GenGaussLobattoPtIntWeight;
using FElib::common::polynomial::GenLagrangePoly;
using FElib::common::polynomial::GenLegendrePoly;
using FElib::utility::ArrayView;
using namespace FElib::test;

namespace {

constexpr double kCheckEps = 5.0e-15;

// -- Closed-form references, copied from FElib/test/FE/polynomial/test_polynomial.f90 --

void TestClosedFormReferences()
{
    std::cout << "-- TestClosedFormReferences\n";

    // n = 1
    {
        std::vector<double> lglPtsAns = {-1.0, 1.0};
        std::vector<double> lglIntwAns = {1.0, 1.0};
        std::vector<double> glPtsAns = {0.0};
        std::vector<double> glIntwAns = {2.0};

        Check(MaxAbsDiff(GenGaussLobattoPt(1), lglPtsAns) <= kCheckEps, "n=1: GLL points");
        Check(MaxAbsDiff(GenGaussLobattoPtIntWeight(1), lglIntwAns) <= kCheckEps, "n=1: GLL weights");
        Check(MaxAbsDiff(GenGaussLegendrePt(1), glPtsAns) <= kCheckEps, "n=1: GL points");
        Check(MaxAbsDiff(GenGaussLegendrePtIntWeight(1), glIntwAns) <= kCheckEps, "n=1: GL weights");
    }

    // n = 2
    {
        const double glN2 = std::sqrt(1.0 / 3.0);
        std::vector<double> lglPtsAns = {-1.0, 0.0, 1.0};
        std::vector<double> lglIntwAns = {1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0};
        std::vector<double> glPtsAns = {-glN2, glN2};
        std::vector<double> glIntwAns = {1.0, 1.0};

        Check(MaxAbsDiff(GenGaussLobattoPt(2), lglPtsAns) <= kCheckEps, "n=2: GLL points");
        Check(MaxAbsDiff(GenGaussLobattoPtIntWeight(2), lglIntwAns) <= kCheckEps, "n=2: GLL weights");
        Check(MaxAbsDiff(GenGaussLegendrePt(2), glPtsAns) <= kCheckEps, "n=2: GL points");
        Check(MaxAbsDiff(GenGaussLegendrePtIntWeight(2), glIntwAns) <= kCheckEps, "n=2: GL weights");
    }

    // n = 3
    {
        const double xN3 = std::sqrt(1.0 / 5.0);
        const double glN3 = std::sqrt(3.0 / 5.0);
        std::vector<double> lglPtsAns = {-1.0, -xN3, xN3, 1.0};
        std::vector<double> lglIntwAns = {1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0};
        std::vector<double> glPtsAns = {-glN3, 0.0, glN3};
        std::vector<double> glIntwAns = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

        Check(MaxAbsDiff(GenGaussLobattoPt(3), lglPtsAns) <= kCheckEps, "n=3: GLL points");
        Check(MaxAbsDiff(GenGaussLobattoPtIntWeight(3), lglIntwAns) <= kCheckEps, "n=3: GLL weights");
        Check(MaxAbsDiff(GenGaussLegendrePt(3), glPtsAns) <= kCheckEps, "n=3: GL points");
        Check(MaxAbsDiff(GenGaussLegendrePtIntWeight(3), glIntwAns) <= kCheckEps, "n=3: GL weights");
    }

    // n = 4
    {
        const double xN4 = std::sqrt(21.0) / 7.0;
        const double gl1N4 = std::sqrt(3.0 + 2.0 * std::sqrt(6.0 / 5.0)) / std::sqrt(7.0);
        const double gl2N4 = std::sqrt(3.0 - 2.0 * std::sqrt(6.0 / 5.0)) / std::sqrt(7.0);
        const double intw1N4 = (18.0 + std::sqrt(30.0)) / 36.0;
        const double intw2N4 = (18.0 - std::sqrt(30.0)) / 36.0;

        std::vector<double> lglPtsAns = {-1.0, -xN4, 0.0, xN4, 1.0};
        std::vector<double> lglIntwAns = {1.0 / 10.0, 49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0, 1.0 / 10.0};
        std::vector<double> glPtsAns = {-gl1N4, -gl2N4, gl2N4, gl1N4};
        std::vector<double> glIntwAns = {intw2N4, intw1N4, intw1N4, intw2N4};

        Check(MaxAbsDiff(GenGaussLobattoPt(4), lglPtsAns) <= kCheckEps, "n=4: GLL points");
        Check(MaxAbsDiff(GenGaussLobattoPtIntWeight(4), lglIntwAns) <= kCheckEps, "n=4: GLL weights");
        Check(MaxAbsDiff(GenGaussLegendrePt(4), glPtsAns) <= kCheckEps, "n=4: GL points");
        Check(MaxAbsDiff(GenGaussLegendrePtIntWeight(4), glIntwAns) <= kCheckEps, "n=4: GL weights");
    }
}

// -- Legendre polynomial sanity checks (P_0=1, P_1(x)=x, endpoint values) --

void TestLegendreSanity()
{
    std::cout << "-- TestLegendreSanity\n";

    const int order = 6;
    std::vector<double> x = {-1.0, -0.5, 0.0, 0.3, 0.7, 1.0};
    ArrayView<const double, 1> xView(x.data(), x.size());

    const std::vector<double> P = GenLegendrePoly(order, xView);
    ArrayView<const double, 2> Pview(P.data(), x.size(), static_cast<std::size_t>(order) + 1);

    bool p0IsOne = true, p1IsX = true;
    for (std::size_t i = 0; i < x.size(); ++i) {
        p0IsOne = p0IsOne && std::abs(Pview(i, 0) - 1.0) < 1e-14;
        p1IsX = p1IsX && std::abs(Pview(i, 1) - x[i]) < 1e-14;
    }
    Check(p0IsOne, "P_0(x) == 1 for all x");
    Check(p1IsX, "P_1(x) == x for all x");

    // Endpoint identities: P_n(1) == 1, P_n(-1) == (-1)^n, for all n.
    std::vector<double> endpoints = {-1.0, 1.0};
    const std::vector<double> Pend = GenLegendrePoly(order, ArrayView<const double, 1>(endpoints.data(), 2));
    ArrayView<const double, 2> Pendview(Pend.data(), 2, static_cast<std::size_t>(order) + 1);
    bool endpointsOk = true;
    for (int n = 0; n <= order; ++n) {
        const double expectedAtMinus1 = (n % 2 == 0) ? 1.0 : -1.0;
        endpointsOk = endpointsOk && std::abs(Pendview(1, static_cast<std::size_t>(n)) - 1.0) < 1e-13;
        endpointsOk = endpointsOk &&
                      std::abs(Pendview(0, static_cast<std::size_t>(n)) - expectedAtMinus1) < 1e-13;
    }
    Check(endpointsOk, "P_n(1) == 1 and P_n(-1) == (-1)^n for n=0..order");
}

// -- Lagrange basis: Kronecker-delta property at the GLL points --

void TestLagrangeKroneckerDelta()
{
    std::cout << "-- TestLagrangeKroneckerDelta\n";

    for (int order = 1; order <= 6; ++order) {
        const std::vector<double> xLgl = GenGaussLobattoPt(order);
        ArrayView<const double, 1> xLglView(xLgl.data(), xLgl.size());

        const std::vector<double> L = GenLagrangePoly(order, xLglView, xLglView);
        ArrayView<const double, 2> Lview(L.data(), xLgl.size(), xLgl.size());

        bool ok = true;
        for (std::size_t n = 0; n < xLgl.size(); ++n) {
            for (std::size_t i = 0; i < xLgl.size(); ++i) {
                const double expected = (i == n) ? 1.0 : 0.0;
                ok = ok && std::abs(Lview(i, n) - expected) < 1e-12;
            }
        }
        Check(ok, "order=" + std::to_string(order) + ": Lagrange basis is Kronecker-delta at GLL points");
    }
}

// -- Differentiation matrix identities --

void TestDifferentiationMatrixIdentities()
{
    std::cout << "-- TestDifferentiationMatrixIdentities\n";

    for (int order = 1; order <= 6; ++order) {
        const std::vector<double> xLgl = GenGaussLobattoPt(order);
        ArrayView<const double, 1> xLglView(xLgl.data(), xLgl.size());
        const std::vector<double> Lr = GenDLagrangePoly_lglpt(order, xLglView);
        ArrayView<const double, 2> LrView(Lr.data(), xLgl.size(), xLgl.size());

        // Lr(k, n) = derivative of the k-th basis function at point n (row = basis,
        // column = point), so for a fixed evaluation point n we sum over the basis
        // index k (rows). d/dx of the constant function 1 must be 0 at every point.
        bool constDerivIsZero = true;
        for (std::size_t n = 0; n < xLgl.size(); ++n) {
            double s = 0.0;
            for (std::size_t k = 0; k < xLgl.size(); ++k) s += LrView(k, n);
            constDerivIsZero = constDerivIsZero && std::abs(s) < 1e-12;
        }
        Check(constDerivIsZero,
              "order=" + std::to_string(order) + ": differentiation matrix annihilates constants");

        // d/dx of f(x) = x must be 1 at every point.
        bool linearDerivIsOne = true;
        for (std::size_t n = 0; n < xLgl.size(); ++n) {
            double s = 0.0;
            for (std::size_t k = 0; k < xLgl.size(); ++k) s += LrView(k, n) * xLglView(k);
            linearDerivIsOne = linearDerivIsOne && std::abs(s - 1.0) < 1e-12;
        }
        Check(linearDerivIsOne,
              "order=" + std::to_string(order) + ": differentiation matrix reproduces d/dx(x)=1");
    }
}

// -- Gauss-Legendre quadrature exactness, at orders with no simple closed form --

void TestGaussLegendreExactness()
{
    std::cout << "-- TestGaussLegendreExactness\n";

    for (int order : {5, 6, 7, 8}) {
        const std::vector<double> pts = GenGaussLegendrePt(order);
        const std::vector<double> w = GenGaussLegendrePtIntWeight(order);

        // An order-point Gauss-Legendre rule integrates polynomials up to
        // degree 2*order-1 exactly. Integral of x^m over [-1,1] is
        // 0 for odd m, 2/(m+1) for even m.
        bool exact = true;
        for (int m = 0; m <= 2 * order - 1; ++m) {
            double quad = 0.0;
            for (std::size_t i = 0; i < pts.size(); ++i) quad += w[i] * std::pow(pts[i], m);
            const double exactVal = (m % 2 == 0) ? 2.0 / static_cast<double>(m + 1) : 0.0;
            exact = exact && std::abs(quad - exactVal) < 1e-11;
        }
        Check(exact, "order=" + std::to_string(order) +
                         ": Gauss-Legendre rule is exact for polynomials up to degree 2*order-1");
    }
}

// -- Dual API: core (out-param) and convenience (returning) overloads must agree --

void TestDualApiConsistency()
{
    std::cout << "-- TestDualApiConsistency\n";

    const int order = 5;
    std::vector<double> x = {-0.9, -0.3, 0.1, 0.6, 1.0};
    ArrayView<const double, 1> xView(x.data(), x.size());

    const std::vector<double> viaWrapper = GenLegendrePoly(order, xView);
    std::vector<double> viaCore(x.size() * (static_cast<std::size_t>(order) + 1));
    GenLegendrePoly(order, xView, ArrayView<double, 2>(viaCore.data(), x.size(),
                                                        static_cast<std::size_t>(order) + 1));
    Check(viaWrapper == viaCore, "GenLegendrePoly: core and convenience overloads agree exactly");

    const std::vector<double> lglWrapper = GenGaussLobattoPt(order);
    std::vector<double> lglCore(static_cast<std::size_t>(order) + 1);
    GenGaussLobattoPt(order, ArrayView<double, 1>(lglCore.data(), lglCore.size()));
    Check(lglWrapper == lglCore, "GenGaussLobattoPt: core and convenience overloads agree exactly");
}

// -- Validation: inputs the Fortran source mishandles are now rejected --

void TestValidation()
{
    std::cout << "-- TestValidation\n";

    std::vector<double> x = {0.0, 0.5};
    ArrayView<const double, 1> xView(x.data(), x.size());
    std::vector<double> P(x.size() * 3);

    CheckThrows([&] { GenLegendrePoly(-1, xView); }, "GenLegendrePoly(order=-1, ...) throws");
    CheckThrows([&] { GenGaussLobattoPt(0); }, "GenGaussLobattoPt(order=0) throws");
    CheckThrows([&] { GenGaussLobattoPtIntWeight(0); }, "GenGaussLobattoPtIntWeight(order=0) throws");
    CheckThrows([&] { GenGaussLegendrePt(0); }, "GenGaussLegendrePt(order=0) throws");
    CheckThrows([&] { GenGaussLegendrePtIntWeight(0); }, "GenGaussLegendrePtIntWeight(order=0) throws");

    std::vector<double> xLgl2 = {-1.0, 1.0};  // order=1 GLL points
    ArrayView<const double, 1> xLgl2View(xLgl2.data(), 2);
    CheckThrows([&] { GenLagrangePoly(0, xLgl2View, xView); }, "GenLagrangePoly(order=0, ...) throws");

    // Shape mismatch: P sized for the wrong order.
    CheckThrows(
        [&] {
            GenLegendrePoly(2, xView, ArrayView<double, 2>(P.data(), x.size(), 2 /* should be 3 */));
        },
        "GenLegendrePoly with a mis-shaped output array throws");

    // order == 0 is valid (not rejected) for the plain Legendre/D-Lagrange-matrix
    // routines, since their formulas do not degenerate at order == 0.
    std::vector<double> onePt = {0.3};
    ArrayView<const double, 1> onePtView(onePt.data(), 1);
    bool orderZeroOk = true;
    try {
        auto P0 = GenLegendrePoly(0, onePtView);
        orderZeroOk = (P0.size() == 1) && std::abs(P0[0] - 1.0) < 1e-15;
    } catch (const std::exception&) {
        orderZeroOk = false;
    }
    Check(orderZeroOk, "GenLegendrePoly(order=0, ...) is accepted and returns P_0(x)=1");
}

}  // namespace

int main()
{
    return RunTestMain("test_polynomial_cpp", [] {
        TestClosedFormReferences();
        TestLegendreSanity();
        TestLagrangeKroneckerDelta();
        TestDifferentiationMatrixIdentities();
        TestGaussLegendreExactness();
        TestDualApiConsistency();
        TestValidation();
    });
}
