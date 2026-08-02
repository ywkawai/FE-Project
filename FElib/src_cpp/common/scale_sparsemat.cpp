#include "scale_sparsemat.hpp"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>

namespace FElib {
namespace common {

StorageFormat ParseStorageFormat(std::string_view name)
{
    if (name == "CSR") return StorageFormat::CSR;
    if (name == "ELL") return StorageFormat::ELL;
    throw std::invalid_argument(
        "ParseStorageFormat: unknown storage format \"" + std::string(name) + "\"");
}

std::string_view ToString(StorageFormat format)
{
    switch (format) {
    case StorageFormat::CSR: return "CSR";
    case StorageFormat::ELL: return "ELL";
    }
    throw std::invalid_argument("ToString(StorageFormat): invalid enum value");
}

namespace {

void CheckIndex(std::size_t i, std::size_t j, std::size_t rows, std::size_t cols)
{
    if (i >= rows || j >= cols) {
        throw std::out_of_range("SparseMat: index out of range");
    }
}

}  // namespace

void SparseMat::Init(DenseView mat, Real eps, StorageFormat format)
{
    if (mat.extent(0) == 0 || mat.extent(1) == 0) {
        throw std::invalid_argument("SparseMat::Init: mat must be non-empty");
    }
    if (eps < 0.0) {
        throw std::invalid_argument("SparseMat::Init: eps must be >= 0");
    }

    format_ = format;
    switch (format_) {
    case StorageFormat::CSR: InitCSR(mat, eps); return;
    case StorageFormat::ELL: InitELL(mat, eps); return;
    }
    throw std::invalid_argument("SparseMat::Init: invalid storage format");
}

void SparseMat::InitCSR(DenseView mat, Real eps)
{
    M_ = mat.extent(0);
    N_ = mat.extent(1);
    colSize_ = 0;

    val_.clear();
    colIdx_.clear();
    rowPtr_.assign(M_ + 1, 0);

    for (std::size_t i = 0; i < M_; ++i) {
        for (std::size_t j = 0; j < N_; ++j) {
            // Strict '>' matches the CSR branch of the Fortran source.
            if (std::abs(mat(i, j)) > eps) {
                val_.push_back(mat(i, j));
                colIdx_.push_back(static_cast<int>(j));
            }
        }
        rowPtr_[i + 1] = val_.size();
    }
}

void SparseMat::InitELL(DenseView mat, Real eps)
{
    M_ = mat.extent(0);
    N_ = mat.extent(1);
    rowPtr_.clear();

    std::vector<std::size_t> rowCount(M_, 0);
    std::size_t maxCount = 0;
    for (std::size_t i = 0; i < M_; ++i) {
        for (std::size_t j = 0; j < N_; ++j) {
            // '>=' matches the ELL branch of the Fortran source; note this
            // differs from CSR's strict '>' there (kept as-is, see header).
            if (std::abs(mat(i, j)) >= eps) {
                ++rowCount[i];
            }
        }
        maxCount = std::max(maxCount, rowCount[i]);
    }

    colSize_ = maxCount;
    val_.assign(M_ * colSize_, Real{0});
    colIdx_.assign(M_ * colSize_, kNoEntry);

    for (std::size_t i = 0; i < M_; ++i) {
        std::size_t k = 0;
        for (std::size_t j = 0; j < N_ && k < colSize_; ++j) {
            if (std::abs(mat(i, j)) >= eps) {
                const std::size_t idx = i * colSize_ + k;
                val_[idx] = mat(i, j);
                colIdx_[idx] = static_cast<int>(j);
                ++k;
            }
        }
        // Remaining slots [k, colSize_) keep val=0 / colIdx=kNoEntry from the
        // assign() calls above -- unlike the Fortran version, these are never
        // reassigned to a borrowed column index from a neighboring slot.
    }
}

void SparseMat::Print(std::ostream& os) const
{
    os << "-- SparseMat (" << ToString(format_) << ") --\n";
    os << "shape: " << M_ << " x " << N_ << "\n";
    os << "storage size: " << val_.size() << "\n";
    for (std::size_t i = 0; i < M_; ++i) {
        os << "row " << i << ": ";
        for (std::size_t j = 0; j < N_; ++j) {
            const auto v = GetVal(i, j);
            os << (v ? *v : Real{0}) << ' ';
        }
        os << "\n";
    }
}

std::optional<SparseMat::Real> SparseMat::GetVal(std::size_t i, std::size_t j) const
{
    CheckIndex(i, j, M_, N_);

    switch (format_) {
    case StorageFormat::CSR:
        for (std::size_t n = rowPtr_[i]; n < rowPtr_[i + 1]; ++n) {
            if (static_cast<std::size_t>(colIdx_[n]) == j) {
                return val_[n];
            }
        }
        return std::nullopt;

    case StorageFormat::ELL:
        for (std::size_t k = 0; k < colSize_; ++k) {
            const std::size_t idx = i * colSize_ + k;
            if (colIdx_[idx] == kNoEntry) {
                continue;
            }
            if (static_cast<std::size_t>(colIdx_[idx]) == j) {
                return val_[idx];
            }
        }
        return std::nullopt;
    }
    return std::nullopt;
}

bool SparseMat::ReplaceVal(std::size_t i, std::size_t j, Real v)
{
    CheckIndex(i, j, M_, N_);

    switch (format_) {
    case StorageFormat::CSR:
        for (std::size_t n = rowPtr_[i]; n < rowPtr_[i + 1]; ++n) {
            if (static_cast<std::size_t>(colIdx_[n]) == j) {
                val_[n] = v;
                return true;
            }
        }
        return false;

    case StorageFormat::ELL:
        for (std::size_t k = 0; k < colSize_; ++k) {
            const std::size_t idx = i * colSize_ + k;
            if (colIdx_[idx] == kNoEntry) {
                continue;
            }
            if (static_cast<std::size_t>(colIdx_[idx]) == j) {
                val_[idx] = v;
                return true;
            }
        }
        return false;
    }
    return false;
}

void MatVec(const SparseMat& A, ArrayView<const SparseMat::Real, 1> b,
            ArrayView<SparseMat::Real, 1> c)
{
    if (b.extent(0) != A.N_ || c.extent(0) != A.M_) {
        throw std::invalid_argument("MatVec: shape mismatch");
    }

    switch (A.format_) {
    case StorageFormat::CSR:
        for (std::size_t p = 0; p < A.M_; ++p) {
            SparseMat::Real s = 0;
            for (std::size_t n = A.rowPtr_[p]; n < A.rowPtr_[p + 1]; ++n) {
                s += A.val_[n] * b(static_cast<std::size_t>(A.colIdx_[n]));
            }
            c(p) = s;
        }
        return;

    case StorageFormat::ELL:
        for (std::size_t i = 0; i < A.M_; ++i) {
            SparseMat::Real s = 0;
            for (std::size_t k = 0; k < A.colSize_; ++k) {
                const std::size_t idx = i * A.colSize_ + k;
                if (A.colIdx_[idx] == SparseMat::kNoEntry) continue;
                s += A.val_[idx] * b(static_cast<std::size_t>(A.colIdx_[idx]));
            }
            c(i) = s;
        }
        return;
    }
}

void MatVec2(const SparseMat& A, ArrayView<const SparseMat::Real, 1> b1,
             ArrayView<const SparseMat::Real, 1> b2, ArrayView<SparseMat::Real, 1> c)
{
    if (b1.extent(0) != A.N_ || b2.extent(0) != A.N_ || c.extent(0) != A.M_) {
        throw std::invalid_argument("MatVec2: shape mismatch");
    }

    switch (A.format_) {
    case StorageFormat::CSR:
        for (std::size_t p = 0; p < A.M_; ++p) {
            SparseMat::Real s = 0;
            for (std::size_t n = A.rowPtr_[p]; n < A.rowPtr_[p + 1]; ++n) {
                const auto col = static_cast<std::size_t>(A.colIdx_[n]);
                s += A.val_[n] * b1(col) * b2(col);
            }
            c(p) = s;
        }
        return;

    case StorageFormat::ELL:
        for (std::size_t i = 0; i < A.M_; ++i) {
            SparseMat::Real s = 0;
            for (std::size_t k = 0; k < A.colSize_; ++k) {
                const std::size_t idx = i * A.colSize_ + k;
                if (A.colIdx_[idx] == SparseMat::kNoEntry) continue;
                const auto col = static_cast<std::size_t>(A.colIdx_[idx]);
                s += A.val_[idx] * b1(col) * b2(col);
            }
            c(i) = s;
        }
        return;
    }
}

void MatMat(const SparseMat& A, ArrayView<const SparseMat::Real, 2> b,
            ArrayView<SparseMat::Real, 2> c)
{
    if (b.extent(0) != A.N_ || c.extent(0) != A.M_ || b.extent(1) != c.extent(1)) {
        throw std::invalid_argument("MatMat: shape mismatch");
    }
    const std::size_t NQ = b.extent(1);

    switch (A.format_) {
    case StorageFormat::CSR:
        for (std::size_t p = 0; p < A.M_; ++p) {
            for (std::size_t q = 0; q < NQ; ++q) c(p, q) = 0;
            for (std::size_t n = A.rowPtr_[p]; n < A.rowPtr_[p + 1]; ++n) {
                const auto col = static_cast<std::size_t>(A.colIdx_[n]);
                for (std::size_t q = 0; q < NQ; ++q) {
                    c(p, q) += A.val_[n] * b(col, q);
                }
            }
        }
        return;

    case StorageFormat::ELL:
        for (std::size_t i = 0; i < A.M_; ++i) {
            for (std::size_t q = 0; q < NQ; ++q) c(i, q) = 0;
            for (std::size_t k = 0; k < A.colSize_; ++k) {
                const std::size_t idx = i * A.colSize_ + k;
                if (A.colIdx_[idx] == SparseMat::kNoEntry) continue;
                const auto col = static_cast<std::size_t>(A.colIdx_[idx]);
                for (std::size_t q = 0; q < NQ; ++q) {
                    c(i, q) += A.val_[idx] * b(col, q);
                }
            }
        }
        return;
    }
}

}  // namespace common
}  // namespace FElib
