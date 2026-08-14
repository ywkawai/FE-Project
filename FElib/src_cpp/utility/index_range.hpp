/*
 * @par Description
 * A dimension-generic replacement for the erb-based d=1..3 loop-nest duplication used on the Fortran side (e.g. FElib/src/common/scale_timeint_rk.F90.erb).
 * BoxRange<Rank> describes a half-open index box; ForEachIndex visits every point in it, handing the caller's lambda one index per dimension. 
 * Only the (trivial) loop nest is specialized per Rank here -- the numerical formula that a caller passes in as the lambda body is written exactly once,
 * regardless of Rank.
 *
 * @author Yuta Kawai, Team SCALE
 */
#pragma once

#include <array>
#include <cstddef>

namespace FElib {
namespace utility {

//! Half-open box [lo, hi) in each of Rank dimensions, 0-based.
//! Rank is `int` (not `std::size_t`) to match FElib::utility::ArrayView and FElib::common::TimeIntRK,
//! so a single Rank can be deduced consistently across all three in generic code.
template <int Rank>
struct BoxRange {
    std::array<std::size_t, Rank> lo{};
    std::array<std::size_t, Rank> hi{};
};

// EXPERIMENTAL: parallelized with OpenMP when built with -fopenmp (silently ignored otherwise). 
// This is safe for every kernel in scale_timeint_rk.cpp
// because each kernel's per-point work only reads/writes memory addressed by that point's own flat index. 
// If ForEachIndex is ever reused for a kernel that does NOT have this per-point-independent property, 
// that call site would need its own reduction/critical handling. 
//
// Collapse depth deliberately mirrors scale_timeint_rk.F90.erb's OpenMP convention, NOT its OpenACC one: 
// every `!$omp do`/`!$omp parallel do` in that file collapses only `d-1` dimensions 
// (`!$acc loop`/`!$acc parallel loop` collapses the full `d`). 
// The innermost (fastest-varying, unit-stride via ArrayView's row-major layout) dimension is deliberately left OUT of the OpenMP collapse 
// so each thread's share of the work is a contiguous run the compiler can still auto-vectorize internally; 
// OpenACC's full collapse instead maps every point to its own GPU thread. Rank=1 has only one dimension, 
// so there is nothing to leave uncollapsed for it, matching the erb file's d=1 case (which also just uses a plain `!$omp do`).

// EXPERIMENTAL: the threshold for switching to OpenMP is a provisional constant, not a runtime-tunable parameter.
constexpr std::size_t kOmpMinPoints = 10000;

template <typename F>
void ForEachIndex(const BoxRange<1>& box, F&& f)
{
    const std::size_t n0 = box.hi[0] - box.lo[0];

#pragma omp parallel for if (n0 >= kOmpMinPoints)
    for (std::size_t i0 = box.lo[0]; i0 < box.hi[0]; ++i0) {
        f(i0);
    }
}

template <typename F>
void ForEachIndex(const BoxRange<2>& box, F&& f)
{
    const std::size_t n0 = box.hi[0] - box.lo[0];
    const std::size_t n1 = box.hi[1] - box.lo[1];

#pragma omp parallel for if (n0 * n1 >= kOmpMinPoints)
    for (std::size_t i0 = box.lo[0]; i0 < box.hi[0]; ++i0) {
        for (std::size_t i1 = box.lo[1]; i1 < box.hi[1]; ++i1) {
            f(i0, i1);
        }
    }
}

template <typename F>
void ForEachIndex(const BoxRange<3>& box, F&& f)
{
    const std::size_t n0 = box.hi[0] - box.lo[0];
    const std::size_t n1 = box.hi[1] - box.lo[1];
    const std::size_t n2 = box.hi[2] - box.lo[2];

#pragma omp parallel for collapse(2) if (n0 * n1 * n2 >= kOmpMinPoints)
    for (std::size_t i0 = box.lo[0]; i0 < box.hi[0]; ++i0) {
        for (std::size_t i1 = box.lo[1]; i1 < box.hi[1]; ++i1) {
            for (std::size_t i2 = box.lo[2]; i2 < box.hi[2]; ++i2) {
                f(i0, i1, i2);
            }
        }
    }
}

}  // namespace utility
}  // namespace FElib
