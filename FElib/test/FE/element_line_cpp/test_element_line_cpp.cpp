// Unit test of scale_element_base / scale_element_line (LineElement).
//
// There is no narrow Fortran unit test for scale_element_base/scale_element_line
// alone (only large integration tests like test_mesh_linedom1d.f90, which
// build a LineElement but never inspect its matrices), so this test's
// reference values are derived independently rather than ported from
// Fortran:
//
//  - For p=1 and p=2, the mass/differentiation/lift matrices have simple
//    exact rational closed forms, hand-derived from the construction
//    formulas and cross-checked against the standard consistent-mass-matrix
//    formula and the well-known Legendre-Gauss-Lobatto differentiation
//    matrix (see e.g. Canuto et al., "Spectral Methods", or Kopriva,
//    "Implementing Spectral Methods for PDEs").
//  - For p=1..3, the differentiation matrix is additionally checked against
//    the classical closed-form GLL differentiation matrix formula (Costa &
//    Don, 2000), and the mass matrix against numerical integration of the
//    Lagrange basis products via Gauss-Legendre quadrature -- both
//    implemented independently below, i.e. via a different computational
//    path than LineElement's own (GenDLagrangePoly_lglpt-based / V*V^T-based)
//    construction.

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "cpp_test_harness.hpp"
#include "element/scale_element_line.hpp"
#include "scale_linalgebra.hpp"
#include "scale_polynomial.hpp"

using FElib::common::linalgebra::Inv;
using FElib::common::polynomial::GenGaussLegendrePt;
using FElib::common::polynomial::GenGaussLegendrePtIntWeight;
using FElib::common::polynomial::GenLegendrePoly;
using FElib::element::LineElement;
using FElib::utility::ArrayView;
using namespace FElib::test;

namespace {

constexpr double kExactEps = 1.0e-13;
constexpr double kQuadEps = 1.0e-11;

// -- Independent reference computations (deliberately NOT reusing
//    LineElement's own construction path) -----------------------------------

// L_i(x), the i-th Lagrange basis function built from `nodes`, at x.
double LagrangeBasisAt(const std::vector<double>& nodes, std::size_t i, double x)
{
    double v = 1.0;
    for (std::size_t k = 0; k < nodes.size(); ++k) {
        if (k == i) continue;
        v *= (x - nodes[k]) / (nodes[i] - nodes[k]);
    }
    return v;
}

// Classical closed-form Legendre-Gauss-Lobatto differentiation matrix
// (Costa & Don, 2000): D(i,j) = [P_N(x_i)/P_N(x_j)] / (x_i - x_j) for i!=j,
// D(0,0) = -N(N+1)/4, D(N,N) = N(N+1)/4, D(i,i) = 0 for interior i.
std::vector<double> CostaDonGllDiffMatrix(const std::vector<double>& nodes, int order)
{
    const std::size_t n = nodes.size();
    ArrayView<const double, 1> nodes_view(nodes.data(), n);
    const std::vector<double> P = GenLegendrePoly(order, nodes_view);
    ArrayView<const double, 2> Pview(P.data(), n, static_cast<std::size_t>(order) + 1);

    std::vector<double> D(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            D[i * n + j] = (Pview(i, static_cast<std::size_t>(order)) / Pview(j, static_cast<std::size_t>(order))) *
                          (1.0 / (nodes[i] - nodes[j]));
        }
    }
    D[0 * n + 0] = -static_cast<double>(order) * static_cast<double>(order + 1) / 4.0;
    D[(n - 1) * n + (n - 1)] = static_cast<double>(order) * static_cast<double>(order + 1) / 4.0;
    return D;
}

// Exact mass matrix M(i,j) = integral_{-1}^{1} L_i(x) L_j(x) dx, evaluated by
// Gauss-Legendre quadrature with enough points to integrate the degree-2*order
// integrand exactly (an order+1-point rule is already exact; order+2 is used
// for a small margin).
std::vector<double> QuadratureMassMatrix(const std::vector<double>& nodes, int order)
{
    const std::size_t n = nodes.size();
    const std::vector<double> qpts = GenGaussLegendrePt(order + 2);
    const std::vector<double> qw = GenGaussLegendrePtIntWeight(order + 2);

    std::vector<double> M(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (std::size_t q = 0; q < qpts.size(); ++q) {
                s += qw[q] * LagrangeBasisAt(nodes, i, qpts[q]) * LagrangeBasisAt(nodes, j, qpts[q]);
            }
            M[i * n + j] = s;
        }
    }
    return M;
}

double MaxAbsMatDiff(const std::vector<double>& a, ArrayView<const double, 2> b)
{
    double m = 0.0;
    const std::size_t n0 = b.extent(0), n1 = b.extent(1);
    for (std::size_t i = 0; i < n0; ++i) {
        for (std::size_t j = 0; j < n1; ++j) {
            m = std::max(m, std::abs(a[i * n1 + j] - b(i, j)));
        }
    }
    return m;
}

// -- i) p=1..3: mass/lift/differentiation matrices match analytically
//    pre-derived reference values ---------------------------------------

void TestAnalyticPropertiesForOrder(int order)
{
    LineElement elem;
    elem.Init(order, /*lumped_mass_mat_flag=*/false);

    const std::size_t np = static_cast<std::size_t>(order) + 1;
    const std::vector<double> nodes(elem.x1().data(), elem.x1().data() + np);

    const std::string tag = "p=" + std::to_string(order) + ": ";

    // Differentiation matrix, vs. the classical closed-form GLL formula.
    const std::vector<double> D_ref = CostaDonGllDiffMatrix(nodes, order);
    Check(MaxAbsMatDiff(D_ref, elem.Dx1()) < kQuadEps, tag + "Dx1 matches the closed-form GLL differentiation matrix");

    // Mass matrix, vs. Gauss-Legendre quadrature of the Lagrange basis products.
    const std::vector<double> M_ref = QuadratureMassMatrix(nodes, order);
    Check(MaxAbsMatDiff(M_ref, elem.M()) < kQuadEps, tag + "M matches the quadrature-integrated mass matrix");

    // Lift matrix: Lift = inv(M) * Emat, where Emat just selects the two
    // boundary nodes (Nfp=1, Nfaces=2 for every line element), so Lift's
    // two columns must equal inv(M_ref)'s column 0 and column (np-1).
    const std::vector<double> invM_ref = Inv(ArrayView<const double, 2>(M_ref.data(), np, np));
    ArrayView<const double, 2> invM_ref_view(invM_ref.data(), np, np);
    double lift_err = 0.0;
    for (std::size_t i = 0; i < np; ++i) {
        lift_err = std::max(lift_err, std::abs(elem.Lift()(i, 0) - invM_ref_view(i, 0)));
        lift_err = std::max(lift_err, std::abs(elem.Lift()(i, 1) - invM_ref_view(i, np - 1)));
    }
    Check(lift_err < kQuadEps, tag + "Lift matches inv(quadrature mass matrix) restricted to the boundary nodes");
}

// -- Exact rational reference values for p=1 and p=2, hand-derived from the
//    construction formulas (see the file-level comment). These duplicate
//    (with tighter tolerance and exact literals) part of what
//    TestAnalyticPropertiesForOrder already checks numerically. --------------

void TestExactValuesP1()
{
    LineElement elem;
    elem.Init(1, /*lumped_mass_mat_flag=*/false);

    // clang-format off
    const std::vector<double> M_ans  = {2.0/3.0, 1.0/3.0,
                                        1.0/3.0, 2.0/3.0};
    const std::vector<double> Dx1_ans = {-0.5, 0.5,
                                         -0.5, 0.5};
    const std::vector<double> Lift_ans = {2.0, -1.0,
                                          -1.0,  2.0};
    // clang-format on

    Check(MaxAbsMatDiff(M_ans, elem.M()) < kExactEps, "p=1: M matches the exact closed form");
    Check(MaxAbsMatDiff(Dx1_ans, elem.Dx1()) < kExactEps, "p=1: Dx1 matches the exact closed form");
    Check(MaxAbsMatDiff(Lift_ans, elem.Lift()) < kExactEps, "p=1: Lift matches the exact closed form");
}

void TestExactValuesP2()
{
    LineElement elem;
    elem.Init(2, /*lumped_mass_mat_flag=*/false);

    // clang-format off
    const std::vector<double> M_ans = {
        4.0/15.0,  2.0/15.0, -1.0/15.0,
        2.0/15.0, 16.0/15.0,  2.0/15.0,
       -1.0/15.0,  2.0/15.0,  4.0/15.0,
    };
    const std::vector<double> Dx1_ans = {
        -1.5,  2.0, -0.5,
        -0.5,  0.0,  0.5,
         0.5, -2.0,  1.5,
    };
    const std::vector<double> Lift_ans = {
        4.5,   1.5,
       -0.75, -0.75,
        1.5,   4.5,
    };
    // clang-format on

    Check(MaxAbsMatDiff(M_ans, elem.M()) < kExactEps, "p=2: M matches the exact closed form");
    Check(MaxAbsMatDiff(Dx1_ans, elem.Dx1()) < kExactEps, "p=2: Dx1 matches the exact closed form");
    Check(MaxAbsMatDiff(Lift_ans, elem.Lift()) < kExactEps, "p=2: Lift matches the exact closed form");
}

// -- Basic structural checks not covered above -------------------------------

void TestBasicProperties()
{
    for (int order = 1; order <= 3; ++order) {
        LineElement elem;
        elem.Init(order, /*lumped_mass_mat_flag=*/false);

        Check(elem.Np() == order + 1, "p=" + std::to_string(order) + ": Np == order+1");
        Check(elem.Nv() == 2, "p=" + std::to_string(order) + ": Nv == 2");
        Check(elem.Nfaces() == 2, "p=" + std::to_string(order) + ": Nfaces == 2");
        Check(elem.NfpTot() == 2, "p=" + std::to_string(order) + ": NfpTot == 2");
        Check(!elem.IsLumpedMatrix(), "p=" + std::to_string(order) + ": IsLumpedMatrix() reflects the Init() flag");
        Check(elem.Fmask()(0, 0) == 0 && elem.Fmask()(0, 1) == elem.Np() - 1,
              "p=" + std::to_string(order) + ": Fmask selects the two boundary nodes");

        // invV * V must be the identity (V is invertible by construction).
        const auto V = elem.V();
        const auto invV = elem.InvV();
        double max_err = 0.0;
        for (int i = 0; i < elem.Np(); ++i) {
            for (int j = 0; j < elem.Np(); ++j) {
                double s = 0.0;
                for (int k = 0; k < elem.Np(); ++k) s += invV(static_cast<std::size_t>(i), static_cast<std::size_t>(k)) * V(static_cast<std::size_t>(k), static_cast<std::size_t>(j));
                max_err = std::max(max_err, std::abs(s - (i == j ? 1.0 : 0.0)));
            }
        }
        Check(max_err < kExactEps, "p=" + std::to_string(order) + ": invV * V == I");
    }
}

void TestLumpedMassMatrixIsDiagonalAndMatchesIntWeights()
{
    for (int order = 1; order <= 3; ++order) {
        LineElement elem;
        elem.Init(order, /*lumped_mass_mat_flag=*/true);

        const auto M = elem.M();
        const auto w = elem.IntWeightLgl();
        double max_err = 0.0;
        for (int i = 0; i < elem.Np(); ++i) {
            for (int j = 0; j < elem.Np(); ++j) {
                const double expect = (i == j) ? w(static_cast<std::size_t>(i)) : 0.0;
                max_err = std::max(max_err, std::abs(M(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) - expect));
            }
        }
        Check(max_err < kExactEps,
              "p=" + std::to_string(order) + " (lumped): M is diagonal with the LGL integration weights");
    }
}

void TestGenIntGaussLegendreIntrpMat()
{
    LineElement elem;
    elem.Init(3, /*lumped_mass_mat_flag=*/false);

    const auto result = elem.GenIntGaussLegendreIntrpMat(4);
    Check(result.points.size() == 4, "GenIntGaussLegendreIntrpMat: points has the requested length");
    Check(result.weights.size() == 4, "GenIntGaussLegendreIntrpMat: weights has the requested length");
    Check(result.matrix.size() == 4 * static_cast<std::size_t>(elem.Np()),
          "GenIntGaussLegendreIntrpMat: matrix has the expected shape");

    // The interpolation matrix maps nodal values at the element's LGL nodes
    // onto the requested Gauss-Legendre points; interpolating the constant
    // function 1 must reproduce 1 at every target point (partition of unity).
    ArrayView<const double, 2> mat(result.matrix.data(), 4, static_cast<std::size_t>(elem.Np()));
    std::vector<double> ones(static_cast<std::size_t>(elem.Np()), 1.0);
    double max_err = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        double s = 0.0;
        for (std::size_t k = 0; k < ones.size(); ++k) s += mat(i, k) * ones[k];
        max_err = std::max(max_err, std::abs(s - 1.0));
    }
    Check(max_err < kQuadEps, "GenIntGaussLegendreIntrpMat: reproduces the constant function 1 exactly");
}

}  // namespace

int main()
{
    return RunTestMain("test_element_line_cpp", [] {
        TestExactValuesP1();
        TestExactValuesP2();
        TestAnalyticPropertiesForOrder(1);
        TestAnalyticPropertiesForOrder(2);
        TestAnalyticPropertiesForOrder(3);
        TestBasicProperties();
        TestLumpedMassMatrixIsDiagonalAndMatchesIntWeights();
        TestGenIntGaussLegendreIntrpMat();
    });
}
