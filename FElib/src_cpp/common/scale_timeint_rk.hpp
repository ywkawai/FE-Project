/*
 * @par Description
 * Runge-Kutta time integrator. Provides the explicit RK / diagonally
 * implicit-explicit (IMEX) RK schemes listed in scale_timeint_rk_butcher_tab.hpp.
 *
 * The Fortran source (FElib/src/common/scale_timeint_rk.F90.erb) uses an erb
 * template to generate near-identical 1D/2D/3D subroutines, because Fortran
 * arrays must have a compile-time-fixed rank. Here that duplication is
 * replaced by a single class template on Rank (the number of spatial
 * dimensions of the fields being advanced), instantiated for Rank = 1, 2, 3;
 * the update formulas are written once and reused for all three ranks via
 * FElib::utility::ForEachIndex (utility/index_range.hpp), which is the
 * dimension-generic stand-in for the erb-generated loop nests.
 *
 * @author Yuta Kawai, Team SCALE
 */
#pragma once

#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

#include "utility/array_view.hpp"
#include "utility/index_range.hpp"

namespace FElib {
namespace common {

using utility::ArrayView;
using utility::BoxRange;

template <int Rank>
class TimeIntRK {
public:
    using Real = double;

    TimeIntRK() = default;

    //! @param rk_scheme_name see scale_timeint_rk_butcher_tab.hpp for supported names.
    //! @param size_each_var  extents of every field this integrator will advance
    //!                       (all fields passed to Advance()/StoreVar0()/etc. must share this shape).
    //! @throws std::invalid_argument if rk_scheme_name is not a supported scheme.
    void Init(std::string_view rk_scheme_name, Real dt, int var_num,
              const std::array<std::size_t, Rank>& size_each_var);

    Real GetDeltaTime() const noexcept { return dt_; }

    //! Diagonal coefficient of the implicit part at a stage, i.e. dt * a_im(nowstage, nowstage).
    //! @param nowstage 0-based stage index.
    Real GetImplicitDiagFac(int nowstage) const { return CoefAIm(nowstage, nowstage) * dt_; }

    int Nstage() const noexcept { return nstage_; }
    int TendBufSize() const noexcept { return tend_buf_size_; }
    bool ImexFlag() const noexcept { return imex_flag_; }
    bool LowStorageFlag() const noexcept { return low_storage_flag_; }

    //! 0-based tendency-buffer slot to write to/read from at a given (0-based) stage.
    int TendBufIndex(int nowstage) const { return tend_buf_indmap_[nowstage]; }

    //! Reference to the caller-owned explicit tendency buffer slot for (varID, stage, idx).
    //! The caller writes the explicit tendency here before calling Advance() for that stage.
    //! @throws std::out_of_range if varID is out of [0, var_num).
    Real& TendBufEx(int varID, int stage, std::array<std::size_t, Rank> idx);

    //! Same as TendBufEx(), for the implicit part (only meaningful if ImexFlag()).
    //! @throws std::out_of_range if varID is out of [0, var_num).
    Real& TendBufIm(int varID, int stage, std::array<std::size_t, Rank> idx);

    //! Advance q (a single variable, identified by varID) over box by one RK stage.
    void Advance(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box);

    //! Advance several variables at once. For low-storage schemes this pairs
    //! variables two at a time (mirroring the Fortran source's hand-fused
    //! rk_advance_low_storage*D_var2), which is a loop-fusion optimization,
    //! not a change of formula.
    void AdvanceVarlist(int nowstage, const std::vector<ArrayView<Real, Rank>>& vars,
                        const std::vector<int>& varIDs, const BoxRange<Rank>& box);

    //! Advance a density-weighted tracer variable (q holds a mixing ratio; the
    //! RK update is applied to q*density, then divided back out).
    //! @throws std::logic_error if ImexFlag() (not supported, mirrors the Fortran source).
    void AdvanceTracerVar(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box,
                          ArrayView<const Real, Rank> ddens, ArrayView<const Real, Rank> ddens0,
                          ArrayView<const Real, Rank> dens_hyd);

    //! Store q as the time-level-n state for varID (used by IMEX schemes).
    void StoreVar0(ArrayView<const Real, Rank> q, int varID, const BoxRange<Rank>& box);

    //! Apply the diagonal implicit correction after the implicit tendency for
    //! the current stage has been written via TendBufIm(). No-op if !ImexFlag().
    void StoreImplicit(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box);

private:
    Real CoefAEx(int i, int j) const { return coef_a_ex_[static_cast<std::size_t>(i) * nstage_ + j]; }
    Real CoefBEx(int i) const { return coef_b_ex_[i]; }
    Real CoefCEx(int i) const { return coef_c_ex_[i]; }
    Real CoefSigEx(int i, int j) const { return coef_sig_ex_[static_cast<std::size_t>(i) * nstage_ + j]; }
    Real CoefGamEx(int i, int j) const { return coef_gam_ex_[static_cast<std::size_t>(i) * nstage_ + j]; }
    Real CoefAIm(int i, int j) const { return coef_a_im_[static_cast<std::size_t>(i) * nstage_ + j]; }
    Real CoefBIm(int i) const { return coef_b_im_[i]; }

    std::size_t FlatFieldIndexFromArray(const std::array<std::size_t, Rank>& idx) const;

    template <typename... Idx>
    std::size_t FlatFieldIndex(Idx... idx) const
    {
        static_assert(sizeof...(Idx) == static_cast<std::size_t>(Rank), "index count must equal Rank");
        return FlatFieldIndexFromArray(std::array<std::size_t, Rank>{static_cast<std::size_t>(idx)...});
    }

    Real& Var0Flat(std::size_t flat, int varID) { return var0_[flat * static_cast<std::size_t>(var_num_) + varID]; }
    Real& VarTmpFlat(std::size_t flat, int varID) { return var_tmp_[flat * static_cast<std::size_t>(var_num_) + varID]; }
    Real& TendBufExFlat(std::size_t flat, int varID, int stage)
    {
        return tend_buf_ex_[(flat * static_cast<std::size_t>(var_num_) + varID) * tend_buf_size_ + stage];
    }
    Real& TendBufImFlat(std::size_t flat, int varID, int stage)
    {
        return tend_buf_im_[(flat * static_cast<std::size_t>(var_num_) + varID) * tend_buf_size_ + stage];
    }

    void AdvanceLowStorage(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box);
    void AdvanceLowStoragePair(int nowstage, ArrayView<Real, Rank> q1, ArrayView<Real, Rank> q2, int varID1,
                              int varID2, const BoxRange<Rank>& box);
    void AdvanceGeneral(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box);

    void AdvanceTracerVarLowStorage(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box,
                                    ArrayView<const Real, Rank> ddens, ArrayView<const Real, Rank> ddens0,
                                    ArrayView<const Real, Rank> dens_hyd);
    void AdvanceTracerVarGeneral(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box,
                                 ArrayView<const Real, Rank> ddens, ArrayView<const Real, Rank> ddens0,
                                 ArrayView<const Real, Rank> dens_hyd);

    Real dt_ = 0.0;
    int nstage_ = 0;
    int tend_buf_size_ = 0;
    bool low_storage_flag_ = false;
    bool imex_flag_ = false;
    int var_num_ = 0;
    std::array<std::size_t, Rank> extents_{};
    std::size_t field_size_ = 0;

    std::vector<Real> coef_a_ex_, coef_b_ex_, coef_c_ex_;
    std::vector<Real> coef_sig_ex_, coef_gam_ex_;
    std::vector<Real> coef_a_im_, coef_b_im_, coef_c_im_;
    std::vector<int> tend_buf_indmap_;

    std::vector<Real> var0_;         //< shape (field_size, var_num)
    std::vector<Real> var_tmp_;      //< shape (field_size, var_num)
    std::vector<Real> tend_buf_ex_;  //< shape (field_size, var_num, tend_buf_size)
    std::vector<Real> tend_buf_im_;  //< shape (field_size, var_num, tend_buf_size), empty if !imex_flag_
};

}  // namespace common
}  // namespace FElib
