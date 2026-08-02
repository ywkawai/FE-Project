// Unit test for the C++ port of scale_sparsemat
// (FElib/src_cpp/common/scale_sparsemat.{hpp,cpp}).
//
// This mirrors FElib/test/common/sparsemat/test_sparsemat.f90 (same 5x5
// matrix, same GetVal/SpMV checks) and additionally covers behavior that
// is specific to, or intentionally different from, the Fortran version --
// see FElib/src_cpp/common/scale_sparsemat.hpp for the list of differences.

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "scale_sparsemat.hpp"

using FElib::common::MatMat;
using FElib::common::MatVec;
using FElib::common::MatVec2;
using FElib::common::ParseStorageFormat;
using FElib::common::SparseMat;
using FElib::common::StorageFormat;
using FElib::common::ToString;
using FElib::utility::ArrayView;

namespace {

int g_checks = 0;
int g_failures = 0;

void Check(bool cond, const std::string& what)
{
    ++g_checks;
    if (cond) {
        std::cout << "  ok: " << what << "\n";
    } else {
        std::cerr << "  FAILED: " << what << "\n";
        ++g_failures;
    }
}

template <typename Fn>
void CheckThrows(Fn&& fn, const std::string& what)
{
    bool threw = false;
    try {
        fn();
    } catch (const std::exception&) {
        threw = true;
    }
    Check(threw, what + " (expected an exception)");
}

// Same 5x5 matrix as FElib/test/common/sparsemat/test_sparsemat.f90.
constexpr std::size_t kN = 5;

std::vector<double> BuildDenseMatrix()
{
    // clang-format off
    return {
        1, 3, 0, 0, 0,
        1, 2, 5, 0, 0,
        4, 1, 3, 0, 0,
        0, 3, 7, 4, 0,
        1, 0, 0, 0, 5,
    };
    // clang-format on
}

constexpr double kEps = 1e-12;

void TestConstructionConsistency(StorageFormat format)
{
    std::cout << "-- TestConstructionConsistency(" << ToString(format) << ")\n";

    auto data = BuildDenseMatrix();
    ArrayView<const double, 2> A(data.data(), kN, kN);

    SparseMat mat;
    mat.Init(A, kEps, format);

    Check(mat.Rows() == kN, "Rows() == N");
    Check(mat.Cols() == kN, "Cols() == N");
    Check(mat.GetStorageFormat() == format, "GetStorageFormat() matches requested format");

    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t j = 0; j < kN; ++j) {
            const double expect = A(i, j);
            const auto v = mat.GetVal(i, j);
            if (std::abs(expect) > kEps) {
                Check(v.has_value() && std::abs(*v - expect) < kEps,
                      "GetVal(" + std::to_string(i) + "," + std::to_string(j) + ") matches dense value");
            } else {
                Check(!v.has_value(),
                      "GetVal(" + std::to_string(i) + "," + std::to_string(j) + ") reports not-stored");
            }
        }
    }
}

void TestMatVec(StorageFormat format)
{
    std::cout << "-- TestMatVec(" << ToString(format) << ")\n";

    auto data = BuildDenseMatrix();
    ArrayView<const double, 2> A(data.data(), kN, kN);
    SparseMat mat;
    mat.Init(A, kEps, format);

    std::vector<double> x(kN, 1.0), b(kN, 0.0), bAns(kN, 0.0);
    for (std::size_t i = 0; i < kN; ++i) {
        double s = 0;
        for (std::size_t j = 0; j < kN; ++j) s += A(i, j) * x[j];
        bAns[i] = s;
    }

    MatVec(mat, ArrayView<const double, 1>(x.data(), kN), ArrayView<double, 1>(b.data(), kN));

    bool ok = true;
    for (std::size_t i = 0; i < kN; ++i) ok = ok && std::abs(b[i] - bAns[i]) < kEps;
    Check(ok, "MatVec result matches dense matmul");
}

void TestMatVec2(StorageFormat format)
{
    std::cout << "-- TestMatVec2(" << ToString(format) << ")\n";

    auto data = BuildDenseMatrix();
    ArrayView<const double, 2> A(data.data(), kN, kN);
    SparseMat mat;
    mat.Init(A, kEps, format);

    std::vector<double> b1(kN, 2.0), b2(kN, 3.0), c(kN, 0.0), cAns(kN, 0.0);
    for (std::size_t i = 0; i < kN; ++i) {
        double s = 0;
        for (std::size_t j = 0; j < kN; ++j) s += A(i, j) * b1[j] * b2[j];
        cAns[i] = s;
    }

    MatVec2(mat, ArrayView<const double, 1>(b1.data(), kN), ArrayView<const double, 1>(b2.data(), kN),
            ArrayView<double, 1>(c.data(), kN));

    bool ok = true;
    for (std::size_t i = 0; i < kN; ++i) ok = ok && std::abs(c[i] - cAns[i]) < kEps;
    Check(ok, "MatVec2 result matches dense c=A*(b1.*b2)");
}

void TestMatMat(StorageFormat format)
{
    std::cout << "-- TestMatMat(" << ToString(format) << ")\n";

    auto data = BuildDenseMatrix();
    ArrayView<const double, 2> A(data.data(), kN, kN);
    SparseMat mat;
    mat.Init(A, kEps, format);

    constexpr std::size_t kNQ = 2;
    std::vector<double> B(kN * kNQ), C(kN * kNQ, 0.0), CAns(kN * kNQ, 0.0);
    for (std::size_t j = 0; j < kN; ++j) {
        B[j * kNQ + 0] = 1.0;
        B[j * kNQ + 1] = static_cast<double>(j + 1);
    }
    for (std::size_t i = 0; i < kN; ++i) {
        for (std::size_t q = 0; q < kNQ; ++q) {
            double s = 0;
            for (std::size_t j = 0; j < kN; ++j) s += A(i, j) * B[j * kNQ + q];
            CAns[i * kNQ + q] = s;
        }
    }

    MatMat(mat, ArrayView<const double, 2>(B.data(), kN, kNQ), ArrayView<double, 2>(C.data(), kN, kNQ));

    bool ok = true;
    for (std::size_t idx = 0; idx < kN * kNQ; ++idx) ok = ok && std::abs(C[idx] - CAns[idx]) < kEps;
    Check(ok, "MatMat (batched, NQ=2) result matches dense reference");
}

void TestReplaceVal(StorageFormat format)
{
    std::cout << "-- TestReplaceVal(" << ToString(format) << ")\n";

    auto data = BuildDenseMatrix();
    ArrayView<const double, 2> A(data.data(), kN, kN);
    SparseMat mat;
    mat.Init(A, kEps, format);

    // A(0,1) == 3 is stored -> ReplaceVal must succeed and update it.
    Check(mat.ReplaceVal(0, 1, 99.0), "ReplaceVal on a stored entry returns true");
    auto v = mat.GetVal(0, 1);
    Check(v.has_value() && *v == 99.0, "GetVal reflects the value written by ReplaceVal");

    // A(0,2) == 0 is not stored -> ReplaceVal must no-op and return false.
    Check(!mat.ReplaceVal(0, 2, 42.0), "ReplaceVal on a non-stored entry returns false");
    Check(!mat.GetVal(0, 2).has_value(), "GetVal still reports not-stored after a failed ReplaceVal");
}

void TestOutOfRangeIndex(StorageFormat format)
{
    std::cout << "-- TestOutOfRangeIndex(" << ToString(format) << ")\n";

    auto data = BuildDenseMatrix();
    ArrayView<const double, 2> A(data.data(), kN, kN);
    SparseMat mat;
    mat.Init(A, kEps, format);

    CheckThrows([&] { mat.GetVal(kN, 0); }, "GetVal(Rows(), 0) throws std::out_of_range");
    CheckThrows([&] { mat.GetVal(0, kN); }, "GetVal(0, Cols()) throws std::out_of_range");
    CheckThrows([&] { mat.ReplaceVal(kN, 0, 1.0); }, "ReplaceVal(Rows(), 0, ...) throws std::out_of_range");
}

void TestShapeValidation()
{
    std::cout << "-- TestShapeValidation\n";

    SparseMat mat;
    ArrayView<const double, 2> empty(nullptr, 0, 0);
    CheckThrows([&] { mat.Init(empty); }, "Init on an empty (0x0) matrix throws");

    auto data = BuildDenseMatrix();
    ArrayView<const double, 2> A(data.data(), kN, kN);
    CheckThrows([&] { mat.Init(A, -1.0); }, "Init with a negative eps throws");

    // Mismatched vector length in MatVec.
    SparseMat csr;
    csr.Init(A, kEps, StorageFormat::CSR);
    std::vector<double> wrongSize(kN - 1, 0.0), c(kN, 0.0);
    CheckThrows(
        [&] {
            MatVec(csr, ArrayView<const double, 1>(wrongSize.data(), wrongSize.size()),
                   ArrayView<double, 1>(c.data(), kN));
        },
        "MatVec with a mismatched vector length throws");
}

void TestInvalidStorageFormatString()
{
    std::cout << "-- TestInvalidStorageFormatString\n";

    CheckThrows([&] { ParseStorageFormat("XYZ"); }, "ParseStorageFormat(\"XYZ\") throws");
    Check(ParseStorageFormat("CSR") == StorageFormat::CSR, "ParseStorageFormat(\"CSR\") == StorageFormat::CSR");
    Check(ParseStorageFormat("ELL") == StorageFormat::ELL, "ParseStorageFormat(\"ELL\") == StorageFormat::ELL");
}

// Documents the threshold asymmetry inherited from the Fortran source:
// CSR treats a value as nonzero only when abs(v) > eps (strict), while ELL
// uses abs(v) >= eps (inclusive). This test intentionally exercises a value
// exactly equal to eps, so CSR and ELL disagree on whether it is stored.
// This is NOT a bug in the C++ port -- it is a deliberately preserved
// asymmetry from FElib/src/common/scale_sparsemat.F90 (see the class
// comment in scale_sparsemat.hpp).
void TestEpsBoundaryAsymmetry()
{
    std::cout << "-- TestEpsBoundaryAsymmetry\n";

    constexpr double kBoundaryEps = 1e-6;
    std::vector<double> data(kN * kN, 0.0);
    data[0 * kN + 1] = kBoundaryEps;  // exactly at the threshold

    ArrayView<const double, 2> A(data.data(), kN, kN);

    SparseMat csr, ell;
    csr.Init(A, kBoundaryEps, StorageFormat::CSR);
    ell.Init(A, kBoundaryEps, StorageFormat::ELL);

    Check(!csr.GetVal(0, 1).has_value(),
          "CSR excludes a value exactly == eps (strict '>' in the Fortran source)");
    Check(ell.GetVal(0, 1).has_value() && *ell.GetVal(0, 1) == kBoundaryEps,
          "ELL includes a value exactly == eps (inclusive '>=' in the Fortran source)");
}

// The Fortran ELL padding-fixup loop (scale_sparsemat.F90) starts at flat
// index 2 and never touches flat index 1, which corresponds to (row=0,
// slot=0). If row 0 has no stored nonzeros at all, that slot's column index
// is left at its sentinel value, and a subsequent SpMV could read out of
// bounds. The C++ port assigns every padding slot's sentinel independently
// (not by copying a neighboring flat index), so this case must not crash.
void TestEllAllZeroFirstRow()
{
    std::cout << "-- TestEllAllZeroFirstRow\n";

    std::vector<double> data(kN * kN, 0.0);
    for (std::size_t j = 1; j < kN; ++j) {
        data[1 * kN + j] = static_cast<double>(j);  // row 0 stays all-zero
    }
    ArrayView<const double, 2> A(data.data(), kN, kN);

    SparseMat ell;
    ell.Init(A, kEps, StorageFormat::ELL);

    bool row0Empty = true;
    for (std::size_t j = 0; j < kN; ++j) {
        row0Empty = row0Empty && !ell.GetVal(0, j).has_value();
    }
    Check(row0Empty, "GetVal reports no stored entries for an all-zero first row");

    std::vector<double> x(kN, 1.0), b(kN, 0.0);
    MatVec(ell, ArrayView<const double, 1>(x.data(), kN), ArrayView<double, 1>(b.data(), kN));
    Check(b[0] == 0.0, "MatVec does not crash and yields 0 for the all-zero first row");
}

}  // namespace

int main()
{
    std::cout << "- Start test_sparsemat_cpp ..\n";

    TestConstructionConsistency(StorageFormat::CSR);
    TestConstructionConsistency(StorageFormat::ELL);
    TestMatVec(StorageFormat::CSR);
    TestMatVec(StorageFormat::ELL);
    TestMatVec2(StorageFormat::CSR);
    TestMatVec2(StorageFormat::ELL);
    TestMatMat(StorageFormat::CSR);
    TestMatMat(StorageFormat::ELL);
    TestReplaceVal(StorageFormat::CSR);
    TestReplaceVal(StorageFormat::ELL);
    TestOutOfRangeIndex(StorageFormat::CSR);
    TestOutOfRangeIndex(StorageFormat::ELL);
    TestShapeValidation();
    TestInvalidStorageFormatString();
    TestEpsBoundaryAsymmetry();
    TestEllAllZeroFirstRow();

    std::cout << g_checks << " checks run, " << g_failures << " failed.\n";
    if (g_failures == 0) {
        std::cout << "test_sparsemat_cpp has been succeeded!\n";
        return 0;
    }
    std::cerr << "test_sparsemat_cpp FAILED.\n";
    return 1;
}
