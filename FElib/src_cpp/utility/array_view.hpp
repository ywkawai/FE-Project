#pragma once

#include <cstddef>
#include <stdexcept>

namespace FElib {
namespace utility {

template <typename T, int Rank>
class ArrayView;

template <typename T>
class ArrayView<T, 1> {
public:
    using size_type = std::size_t;
    using value_type = T;

    ArrayView() = default;

    ArrayView(T* data, size_type n0)
        : data_(data), n0_(n0)
    {
        if (data == nullptr && n0 != 0) {
            throw std::invalid_argument("ArrayView<1>: null data");
        }
    }

    T& operator()(size_type i)
    {
        return data_[i];
    }

    const T& operator()(size_type i) const
    {
        return data_[i];
    }

    size_type extent(size_type d) const
    {
        if (d != 0) {
            throw std::out_of_range("ArrayView<1>::extent");
        }
        return n0_;
    }

    size_type size() const noexcept
    {
        return n0_;
    }

    T* data() noexcept
    {
        return data_;
    }

    const T* data() const noexcept
    {
        return data_;
    }

private:
    T* data_ = nullptr;
    size_type n0_ = 0;
};

template <typename T>
class ArrayView<T, 2> {
public:
    using size_type = std::size_t;
    using value_type = T;

    ArrayView() = default;

    ArrayView(T* data, size_type n0, size_type n1)
        : data_(data), n0_(n0), n1_(n1)
    {
        if (data == nullptr && n0 * n1 != 0) {
            throw std::invalid_argument("ArrayView<2>: null data");
        }
    }

    // Row-major layout: (i,j) -> i*n1 + j
    T& operator()(size_type i, size_type j)
    {
        return data_[i * n1_ + j];
    }

    constexpr T& operator()(size_type i, size_type j) const noexcept
    {
        return data_[i * n1_ + j];
    }

    size_type extent(size_type d) const
    {
        if (d == 0) return n0_;
        if (d == 1) return n1_;
        throw std::out_of_range("ArrayView<2>::extent");
    }

    size_type size() const noexcept
    {
        return n0_ * n1_;
    }

    T* data() noexcept
    {
        return data_;
    }

    const T* data() const noexcept
    {
        return data_;
    }

private:
    T* data_ = nullptr;
    size_type n0_ = 0;
    size_type n1_ = 0;
};

template <typename T>
class ArrayView<T, 3> {
public:
    using size_type = std::size_t;
    using value_type = T;

    ArrayView() = default;

    ArrayView(T* data, size_type n0, size_type n1, size_type n2)
        : data_(data), n0_(n0), n1_(n1), n2_(n2)
    {
        if (data == nullptr && n0 * n1 * n2 != 0) {
            throw std::invalid_argument("ArrayView<3>: null data");
        }
    }

    // Row-major layout: (i,j,k) -> (i*n1 + j)*n2 + k
    T& operator()(size_type i, size_type j, size_type k)
    {
        return data_[(i * n1_ + j) * n2_ + k];
    }

    const T& operator()(size_type i, size_type j, size_type k) const
    {
        return data_[(i * n1_ + j) * n2_ + k];
    }

    size_type extent(size_type d) const
    {
        if (d == 0) return n0_;
        if (d == 1) return n1_;
        if (d == 2) return n2_;
        throw std::out_of_range("ArrayView<3>::extent");
    }

    size_type size() const noexcept
    {
        return n0_ * n1_ * n2_;
    }

    T* data() noexcept
    {
        return data_;
    }

    const T* data() const noexcept
    {
        return data_;
    }

private:
    T* data_ = nullptr;
    size_type n0_ = 0;
    size_type n1_ = 0;
    size_type n2_ = 0;
};

}  // namespace utility
}  // namespace FElib