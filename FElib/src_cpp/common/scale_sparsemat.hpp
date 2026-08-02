/* 
 * @par A class to treat sparse matrices and the associated operations.
 *
 * @author Yuta Kawai, Xuanzhengbo Ren, and Team SCALE
 */
#pragma once

#include <cstddef>
#include <iosfwd>
#include <limits>
#include <optional>
#include <string_view>
#include <vector>

#include "utility/array_view.hpp"

namespace FElib {
namespace common {

using utility::ArrayView;

//! Storage format for SparseMat.
enum class StorageFormat {
    CSR,
    ELL,
};

//! Convert a format name ("CSR"/"ELL") to StorageFormat.
//! @throws std::invalid_argument on an unknown name.
StorageFormat ParseStorageFormat(std::string_view name);

//! Convert a StorageFormat to its canonical name ("CSR"/"ELL").
std::string_view ToString(StorageFormat format);

//! Provisional default threshold for treating a dense matrix entry as zero.
//! The Fortran version (FElib/src/common/scale_sparsemat.F90) defaults to `CONST_EPS * 500`, where CONST_EPS is `epsilon(0.0_RP)` set at run time by scale_const (scalelib/src/common/scale_const.F90). 
inline constexpr double kProvisionalDefaultEps = 500.0 * std::numeric_limits<double>::epsilon();

class SparseMat {
public:
    using Real = double;
    using DenseView = ArrayView<const Real, 2>;

    SparseMat() = default;

    //! Build the compressed representation of a dense Rows() x Cols() matrix.
    //! @param mat    dense input; must have mat.extent(0) > 0 and mat.extent(1) > 0.
    //! @param eps    threshold below/at which an entry is treated as zero
    //!               (see class comment for the CSR/ELL threshold difference).
    //! @param format CSR (default) or ELL.
    //! @throws std::invalid_argument if mat is empty or eps < 0.
    void Init(DenseView mat,
              Real eps = kProvisionalDefaultEps,
              StorageFormat format = StorageFormat::CSR);

    //! Print shape and dense-reconstructed contents to `os` (debug aid).
    void Print(std::ostream& os) const;

    //! Look up A(i,j). Returns std::nullopt if the entry was not stored
    //! @throws std::out_of_range if i >= Rows() or j >= Cols().
    std::optional<Real> GetVal(std::size_t i, std::size_t j) const;

    //! Replace A(i,j) with v. Returns false (no-op) if the entry is not stored.
    //! @throws std::out_of_range if i >= Rows() or j >= Cols().
    bool ReplaceVal(std::size_t i, std::size_t j, Real v);

    StorageFormat GetStorageFormat() const noexcept { return format_; }
    std::size_t Rows() const noexcept { return M_; }
    std::size_t Cols() const noexcept { return N_; }

    //! Physical length of the internal value/column-index buffers.
    //! For CSR this equals the actual number of stored nonzeros.
    //! For ELL this equals Rows() * ColSize(), i.e. it INCLUDES the padding slots introduced to give every row the same slot count. 
    std::size_t StorageSize() const noexcept { return val_.size(); }

    //! Number of compressed columns per row in ELL format. 0 for CSR.
    std::size_t ColSize() const noexcept { return colSize_; }

    //! Number of entries in the CSR row-pointer array (Rows()+1). 0 for ELL.
    std::size_t RowPtrSize() const noexcept { return rowPtr_.size(); }

    friend void MatVec(const SparseMat& A, ArrayView<const Real, 1> b, ArrayView<Real, 1> c);
    friend void MatVec2(const SparseMat& A, ArrayView<const Real, 1> b1,
                         ArrayView<const Real, 1> b2, ArrayView<Real, 1> c);
    friend void MatMat(const SparseMat& A, ArrayView<const Real, 2> b, ArrayView<Real, 2> c);

private:
    // Sentinel stored in colIdx_ for unused ELL padding slots. 
    // Never a valid column index (columns are >= 0), so `colIdx_[idx] == kNoEntry` is an unambiguous "no entry here" check.
    static constexpr int kNoEntry = -1;

    void InitCSR(DenseView mat, Real eps);
    void InitELL(DenseView mat, Real eps);

    std::size_t M_ = 0;
    std::size_t N_ = 0;
    StorageFormat format_ = StorageFormat::CSR;

    std::vector<Real> val_;
    std::vector<int> colIdx_;           // 0-based column index; kNoEntry marks ELL padding
    std::vector<std::size_t> rowPtr_;   // CSR only, size Rows()+1 when present
    std::size_t colSize_ = 0;           // ELL only
};

//! c = A * b
void MatVec(const SparseMat& A, ArrayView<const SparseMat::Real, 1> b,
            ArrayView<SparseMat::Real, 1> c);

//! c = A * (b1 elementwise* b2)
void MatVec2(const SparseMat& A, ArrayView<const SparseMat::Real, 1> b1, ArrayView<const SparseMat::Real, 1> b2, 
            ArrayView<SparseMat::Real, 1> c);

//! Batched matrix-vector product: c(:,q) = A * b(:,q) for every column q.
//! b has shape (Cols(), NQ); c has shape (Rows(), NQ). NQ is the fast (inner) index. 
//! This matches the memory layout of the Fortran matmul2's b(NQ,N)/c(NQ,M) arguments once the two Fortran index positions are swapped for row-major ArrayView access (see FElib/src_cpp design notes).
void MatMat(const SparseMat& A, ArrayView<const SparseMat::Real, 2> b,
            ArrayView<SparseMat::Real, 2> c);

}  // namespace common
}  // namespace FElib
