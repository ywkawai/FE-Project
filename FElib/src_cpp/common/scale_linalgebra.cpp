#include "scale_linalgebra.hpp"

#include <stdexcept>
#include <string>

namespace FElib {
namespace common {
namespace linalgebra {

namespace {

// LAPACK's LU factorization/solve, as linked via SCALE_MATHLIB_LIBS.
extern "C" {
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a, const int* lda,
             const int* ipiv, double* b, const int* ldb, int* info);
}

}  // namespace

void SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 2> B, ArrayView<Real, 2> X)
{
    if (A.extent(0) != A.extent(1)) {
        throw std::invalid_argument("linalgebra::SolveLinEq: A must be square");
    }
    const std::size_t n = A.extent(0);
    const std::size_t nrhs = B.extent(1);
    if (B.extent(0) != n) {
        throw std::invalid_argument("linalgebra::SolveLinEq: B's row count must match A's size");
    }
    if (X.extent(0) != n || X.extent(1) != nrhs) {
        throw std::invalid_argument("linalgebra::SolveLinEq: X's shape must match B's shape");
    }

    // A, B, X all use the row-major ArrayView convention (A(i,j), not
    // LAPACK's native column-major layout), so this builds column-major
    // scratch copies for the LAPACK call rather than reinterpreting the
    // row-major buffers directly (which would silently solve a transposed
    // system instead).
    const int n_i = static_cast<int>(n);
    const int nrhs_i = static_cast<int>(nrhs);

    std::vector<Real> a_col(n * n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            a_col[j * n + i] = A(i, j);
        }
    }

    std::vector<int> ipiv(n);
    int info = 0;
    dgetrf_(&n_i, &n_i, a_col.data(), &n_i, ipiv.data(), &info);
    if (info != 0) {
        throw std::runtime_error("linalgebra::SolveLinEq: matrix is singular (dgetrf info=" +
                                  std::to_string(info) + ")");
    }

    std::vector<Real> b_col(n * nrhs);
    for (std::size_t j = 0; j < nrhs; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            b_col[j * n + i] = B(i, j);
        }
    }

    dgetrs_("N", &n_i, &nrhs_i, a_col.data(), &n_i, ipiv.data(), b_col.data(), &n_i, &info);
    if (info != 0) {
        throw std::runtime_error("linalgebra::SolveLinEq: dgetrs failed (info=" + std::to_string(info) + ")");
    }

    for (std::size_t j = 0; j < nrhs; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            X(i, j) = b_col[j * n + i];
        }
    }
}

std::vector<Real> SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 2> B)
{
    std::vector<Real> x(B.extent(0) * B.extent(1));
    SolveLinEq(A, B, ArrayView<Real, 2>(x.data(), B.extent(0), B.extent(1)));
    return x;
}

void SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 1> b, ArrayView<Real, 1> x)
{
    if (x.extent(0) != b.extent(0)) {
        throw std::invalid_argument("linalgebra::SolveLinEq: x's length must match b's length");
    }
    // A length-n vector is the same memory layout as an (n, 1) row-major
    // matrix, so this just reinterprets b/x as single-column matrices and
    // delegates to the matrix form above -- no copy, no duplicated formula.
    SolveLinEq(A, ArrayView<const Real, 2>(b.data(), b.extent(0), 1), ArrayView<Real, 2>(x.data(), x.extent(0), 1));
}

std::vector<Real> SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 1> b)
{
    std::vector<Real> x(b.extent(0));
    SolveLinEq(A, b, ArrayView<Real, 1>(x.data(), x.size()));
    return x;
}

}  // namespace linalgebra
}  // namespace common
}  // namespace FElib
