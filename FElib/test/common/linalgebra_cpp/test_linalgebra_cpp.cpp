// Unit test of scale_linalgebra (SolveLinEq only -- see PORTING_STATUS.md
// for why linalgebra_inv/LU/SolveLinEq_BndMat/SolveLinEq_GMRES are not
// ported yet). The 4x4 system/answer below is copied from
// FElib/test/common/linalgebra/test_linalgebra.f90's test_SolveLinEq so the
// two versions are checked against the same reference values.

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "cpp_test_harness.hpp"
#include "scale_linalgebra.hpp"

using FElib::common::linalgebra::SolveLinEq;
using FElib::utility::ArrayView;
using namespace FElib::test;

namespace {

constexpr double kCheckEps = 1.0e-12;

// -- Fixture, copied from FElib/test/common/linalgebra/test_linalgebra.f90 --

constexpr int kN = 4;
// clang-format off
const std::vector<double> kA = {
    -0.23, 2.54, -3.66,  0.0,
    -6.98, 2.46, -2.73, -2.13,
     0.0,  2.56,  2.46,  4.07,
     0.0,  0.0,  -4.78, -3.82,
};
// clang-format on
const std::vector<double> kB = {4.42, 27.13, -6.14, 10.50};
const std::vector<double> kXAns = {-2.0, 3.0, 1.0, -4.0};

void TestSolveLinEqVector()
{
    ArrayView<const double, 2> A(kA.data(), kN, kN);
    ArrayView<const double, 1> b(kB.data(), kN);

    std::vector<double> x(kN);
    SolveLinEq(A, b, ArrayView<double, 1>(x.data(), kN));
    Check(MaxAbsDiff(x, kXAns) < kCheckEps, "SolveLinEq(vector, out-arg core) matches the known answer");

    const auto x2 = SolveLinEq(A, b);
    Check(MaxAbsDiff(x2, kXAns) < kCheckEps, "SolveLinEq(vector, convenience wrapper) matches the known answer");
}

void TestSolveLinEqMatrix()
{
    ArrayView<const double, 2> A(kA.data(), kN, kN);

    // Two RHS columns: b and 2*b. Since Ax=b is linear, the solution for the
    // second column must be exactly 2*XAns -- this is a self-contained check
    // (no need to borrow the dead A_batch/B_batch fixture from the Fortran
    // test, which doesn't actually fit the "one A, several b's" b2D shape).
    std::vector<double> B(kN * 2);
    for (int i = 0; i < kN; ++i) {
        B[static_cast<std::size_t>(i) * 2 + 0] = kB[i];
        B[static_cast<std::size_t>(i) * 2 + 1] = 2.0 * kB[i];
    }
    ArrayView<const double, 2> Bview(B.data(), kN, 2);

    std::vector<double> X(kN * 2);
    SolveLinEq(A, Bview, ArrayView<double, 2>(X.data(), kN, 2));

    double max_err = 0.0;
    for (int i = 0; i < kN; ++i) {
        max_err = std::max(max_err, std::abs(X[static_cast<std::size_t>(i) * 2 + 0] - kXAns[i]));
        max_err = std::max(max_err, std::abs(X[static_cast<std::size_t>(i) * 2 + 1] - 2.0 * kXAns[i]));
    }
    Check(max_err < kCheckEps, "SolveLinEq(matrix RHS, out-arg core): both columns match the known answer");

    const auto X2 = SolveLinEq(A, Bview);
    Check(MaxAbsDiff(X, X2) < kCheckEps, "SolveLinEq(matrix RHS, convenience wrapper) matches the out-arg core");
}

void TestValidation()
{
    // Non-square A.
    const std::vector<double> a_nonsquare(kN * (kN + 1), 1.0);
    ArrayView<const double, 2> A_nonsquare(a_nonsquare.data(), kN, kN + 1);
    ArrayView<const double, 1> b(kB.data(), kN);
    std::vector<double> x(kN);
    CheckThrows([&] { SolveLinEq(A_nonsquare, b, ArrayView<double, 1>(x.data(), kN)); },
                "SolveLinEq with non-square A throws");

    // b/x length mismatch.
    ArrayView<const double, 2> A(kA.data(), kN, kN);
    std::vector<double> x_wrong(kN + 1);
    CheckThrows([&] { SolveLinEq(A, b, ArrayView<double, 1>(x_wrong.data(), kN + 1)); },
                "SolveLinEq with mismatched b/x length throws");

    // Singular A.
    const std::vector<double> a_singular(kN * kN, 0.0);  // all-zero matrix
    ArrayView<const double, 2> A_singular(a_singular.data(), kN, kN);
    CheckThrows([&] { SolveLinEq(A_singular, b, ArrayView<double, 1>(x.data(), kN)); },
                "SolveLinEq with a singular A throws");
}

}  // namespace

int main()
{
    return RunTestMain("test_linalgebra_cpp", [] {
        TestSolveLinEqVector();
        TestSolveLinEqMatrix();
        TestValidation();
    });
}
