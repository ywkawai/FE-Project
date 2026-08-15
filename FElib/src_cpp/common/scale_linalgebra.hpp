/*
 * @par Description
 * Linear algebra utilities. Only `SolveLinEq` is ported so far -- the part
 * of FElib/src/common/scale_linalgebra.F90 needed by
 * scale_timeint_rk_butcher_tab.cpp's ShuOsher2Butcher (which mirrors the
 * Fortran source's own dependency: `ShuOsher2Butcher` calls
 * `linalgebra_SolveLinEq`). `linalgebra_inv`, `linalgebra_LU`,
 * `linalgebra_SolveLinEq_BndMat`, and `linalgebra_SolveLinEq_GMRES` are not
 * ported yet; see PORTING_STATUS.md for the remaining scope and open
 * questions (an apparent bug in the Fortran ILU(0) preconditioner, a
 * SparseMat API extension GMRES would need, etc.).
 *
 * @author Yuta Kawai, Team SCALE
 */
#pragma once

#include <vector>

#include "utility/array_view.hpp"

namespace FElib {
namespace common {
namespace linalgebra {

using utility::ArrayView;
using Real = double;

//! Solve the linear equation A X = B for X, where A is n x n and B/X are n x nrhs.
//! @throws std::invalid_argument if A is not square, or B's/X's shapes are inconsistent with A or each other.
//! @throws std::runtime_error if A is singular (LAPACK dgetrf/dgetrs reports a failure).
void SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 2> B, ArrayView<Real, 2> X);

//! Convenience overload: allocates and returns X, shape (B.extent(0), B.extent(1)) row-major.
std::vector<Real> SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 2> B);

//! Solve the linear equation A x = b for x, where A is n x n and b/x are length-n vectors.
//! @throws std::invalid_argument if A is not square, or b's/x's length is inconsistent with A or each other.
//! @throws std::runtime_error if A is singular.
void SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 1> b, ArrayView<Real, 1> x);

//! Convenience overload: allocates and returns x, length b.extent(0).
std::vector<Real> SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 1> b);

}  // namespace linalgebra
}  // namespace common
}  // namespace FElib
