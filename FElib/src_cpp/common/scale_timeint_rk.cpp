#include "scale_timeint_rk.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "scale_timeint_rk_butcher_tab.hpp"

namespace FElib {
namespace common {

using utility::ForEachIndex;

namespace {

// Provisional threshold for skipping a genuinely-near-zero Shu-Osher
// contribution in the low-storage update (see PORTING_STATUS.md §2-5: the
// Fortran source uses scale_const::CONST_EPS, a value only available at
// SCALE-DG runtime, so this uses the same provisional-eps convention as
// scale_sparsemat.hpp's kProvisionalDefaultEps).
constexpr double kProvisionalEps = 500.0 * std::numeric_limits<double>::epsilon();

}  // namespace

// ---------------------------------------------------------------------------
// Init / small accessors
// ---------------------------------------------------------------------------

template <int Rank>
std::size_t TimeIntRK<Rank>::FlatFieldIndexFromArray(const std::array<std::size_t, Rank>& idx) const
{
    std::size_t flat = 0;
    for (std::size_t d = 0; d < static_cast<std::size_t>(Rank); ++d) {
        flat = flat * extents_[d] + idx[d];
    }
    return flat;
}

template <int Rank>
void TimeIntRK<Rank>::Init(std::string_view rk_scheme_name, Real dt, int var_num,
                          const std::array<std::size_t, Rank>& size_each_var)
{
    dt_ = dt;
    var_num_ = var_num;
    extents_ = size_each_var;
    field_size_ = 1;
    for (std::size_t d = 0; d < static_cast<std::size_t>(Rank); ++d) {
        field_size_ *= extents_[d];
    }

    const auto info = timeint_rk_butcher_tab::GetInfo(rk_scheme_name);
    nstage_ = info.nstage;
    tend_buf_size_ = info.tend_buf_size;
    low_storage_flag_ = info.low_storage_flag;
    imex_flag_ = info.imex_flag;

    const auto table = timeint_rk_butcher_tab::Get(rk_scheme_name);
    coef_a_ex_ = table.coef_a_ex;
    coef_b_ex_ = table.coef_b_ex;
    coef_c_ex_ = table.coef_c_ex;
    coef_sig_ex_ = table.coef_sig_ex;
    coef_gam_ex_ = table.coef_gam_ex;
    coef_a_im_ = table.coef_a_im;
    coef_b_im_ = table.coef_b_im;
    coef_c_im_ = table.coef_c_im;
    tend_buf_indmap_ = table.tend_buf_indmap;

    const std::size_t var_num_sz = static_cast<std::size_t>(var_num_);
    var0_.assign(field_size_ * var_num_sz, 0.0);
    var_tmp_.assign(field_size_ * var_num_sz, 0.0);
    tend_buf_ex_.assign(field_size_ * var_num_sz * static_cast<std::size_t>(tend_buf_size_), 0.0);
    tend_buf_im_.assign(imex_flag_ ? field_size_ * var_num_sz * static_cast<std::size_t>(tend_buf_size_) : 0, 0.0);
}

template <int Rank>
typename TimeIntRK<Rank>::Real& TimeIntRK<Rank>::TendBufEx(int varID, int stage, std::array<std::size_t, Rank> idx)
{
    if (varID < 0 || varID >= var_num_) {
        throw std::out_of_range("TimeIntRK::TendBufEx: varID out of range");
    }
    return TendBufExFlat(FlatFieldIndexFromArray(idx), varID, stage);
}

template <int Rank>
typename TimeIntRK<Rank>::Real& TimeIntRK<Rank>::TendBufIm(int varID, int stage, std::array<std::size_t, Rank> idx)
{
    if (varID < 0 || varID >= var_num_) {
        throw std::out_of_range("TimeIntRK::TendBufIm: varID out of range");
    }
    return TendBufImFlat(FlatFieldIndexFromArray(idx), varID, stage);
}

// ---------------------------------------------------------------------------
// Advance (single variable)
// ---------------------------------------------------------------------------

template <int Rank>
void TimeIntRK<Rank>::Advance(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box)
{
    if (low_storage_flag_) {
        AdvanceLowStorage(nowstage, q, varID, box);
    } else {
        AdvanceGeneral(nowstage, q, varID, box);
    }
}

template <int Rank>
void TimeIntRK<Rank>::AdvanceLowStorage(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box)
{
    const Real sig_ss = CoefSigEx(nowstage + 1, nowstage);
    const Real sig_Ns = CoefSigEx(nstage_, nowstage);
    const Real one_minus_sig_ss = 1.0 - sig_ss;
    const Real gam_ss = dt_ * CoefGamEx(nowstage + 1, nowstage);
    const Real gam_Ns = dt_ * CoefGamEx(nstage_, nowstage);

    if (nowstage == nstage_ - 1) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            q(idx...) = VarTmpFlat(flat, varID) + sig_ss * q(idx...) + gam_ss * TendBufExFlat(flat, varID, 0);
        });
        return;
    }

    if (nowstage == 0) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            Var0Flat(flat, varID) = q(idx...);
            VarTmpFlat(flat, varID) = 0.0;
        });
    }

    if (std::abs(sig_Ns) > kProvisionalEps || std::abs(gam_Ns) > kProvisionalEps) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            VarTmpFlat(flat, varID) += sig_Ns * q(idx...) + gam_Ns * TendBufExFlat(flat, varID, 0);
        });
    }

    ForEachIndex(box, [&](auto... idx) {
        const std::size_t flat = FlatFieldIndex(idx...);
        q(idx...) = one_minus_sig_ss * Var0Flat(flat, varID) + sig_ss * q(idx...) + gam_ss * TendBufExFlat(flat, varID, 0);
    });
}

template <int Rank>
void TimeIntRK<Rank>::AdvanceLowStoragePair(int nowstage, ArrayView<Real, Rank> q1, ArrayView<Real, Rank> q2,
                                           int varID1, int varID2, const BoxRange<Rank>& box)
{
    const Real sig_ss = CoefSigEx(nowstage + 1, nowstage);
    const Real sig_Ns = CoefSigEx(nstage_, nowstage);
    const Real one_minus_sig_ss = 1.0 - sig_ss;
    const Real gam_ss = dt_ * CoefGamEx(nowstage + 1, nowstage);
    const Real gam_Ns = dt_ * CoefGamEx(nstage_, nowstage);

    if (nowstage == nstage_ - 1) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            q1(idx...) = VarTmpFlat(flat, varID1) + sig_ss * q1(idx...) + gam_ss * TendBufExFlat(flat, varID1, 0);
            q2(idx...) = VarTmpFlat(flat, varID2) + sig_ss * q2(idx...) + gam_ss * TendBufExFlat(flat, varID2, 0);
        });
        return;
    }

    if (nowstage == 0) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            Var0Flat(flat, varID1) = q1(idx...);
            Var0Flat(flat, varID2) = q2(idx...);
            VarTmpFlat(flat, varID1) = 0.0;
            VarTmpFlat(flat, varID2) = 0.0;
        });
    }

    if (std::abs(sig_Ns) > kProvisionalEps || std::abs(gam_Ns) > kProvisionalEps) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            VarTmpFlat(flat, varID1) += sig_Ns * q1(idx...) + gam_Ns * TendBufExFlat(flat, varID1, 0);
            VarTmpFlat(flat, varID2) += sig_Ns * q2(idx...) + gam_Ns * TendBufExFlat(flat, varID2, 0);
        });
    }

    ForEachIndex(box, [&](auto... idx) {
        const std::size_t flat = FlatFieldIndex(idx...);
        q1(idx...) = one_minus_sig_ss * Var0Flat(flat, varID1) + sig_ss * q1(idx...) + gam_ss * TendBufExFlat(flat, varID1, 0);
        q2(idx...) = one_minus_sig_ss * Var0Flat(flat, varID2) + sig_ss * q2(idx...) + gam_ss * TendBufExFlat(flat, varID2, 0);
    });
}

template <int Rank>
void TimeIntRK<Rank>::AdvanceGeneral(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box)
{
    const int tintbuf_ind = tend_buf_indmap_[nowstage];

    if (nowstage == 0 && !imex_flag_) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            Var0Flat(flat, varID) = q(idx...);
            VarTmpFlat(flat, varID) = q(idx...);
        });
    }

    if (nowstage == nstage_ - 1) {
        if (imex_flag_) {
            const Real coef_b_ex_dt = CoefBEx(nowstage) * dt_;
            const Real coef_b_im_dt = CoefBIm(nowstage) * dt_;
            ForEachIndex(box, [&](auto... idx) {
                const std::size_t flat = FlatFieldIndex(idx...);
                q(idx...) = VarTmpFlat(flat, varID) + coef_b_ex_dt * TendBufExFlat(flat, varID, tintbuf_ind) +
                          coef_b_im_dt * TendBufImFlat(flat, varID, tintbuf_ind);
            });
        } else {
            const Real coef_b_ex_dt = CoefBEx(nowstage) * dt_;
            ForEachIndex(box, [&](auto... idx) {
                const std::size_t flat = FlatFieldIndex(idx...);
                q(idx...) = VarTmpFlat(flat, varID) + coef_b_ex_dt * TendBufExFlat(flat, varID, tintbuf_ind);
            });
        }
        return;
    }

    if (imex_flag_) {
        const Real coef_b_ex_dt = CoefBEx(nowstage) * dt_;
        const Real coef_b_im_dt = CoefBIm(nowstage) * dt_;
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            q(idx...) = Var0Flat(flat, varID);
            VarTmpFlat(flat, varID) += coef_b_ex_dt * TendBufExFlat(flat, varID, tintbuf_ind) +
                                     coef_b_im_dt * TendBufImFlat(flat, varID, tintbuf_ind);
        });
    } else {
        const Real coef_b_ex_dt = CoefBEx(nowstage) * dt_;
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            q(idx...) = Var0Flat(flat, varID);
            VarTmpFlat(flat, varID) += coef_b_ex_dt * TendBufExFlat(flat, varID, tintbuf_ind);
        });
    }

    if (tend_buf_size_ == 1 && !imex_flag_) {
        const Real coef_a_ex_dt = dt_ * CoefAEx(nowstage + 1, nowstage);
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            q(idx...) = Var0Flat(flat, varID) + coef_a_ex_dt * TendBufExFlat(flat, varID, 0);
        });
    } else if (!imex_flag_) {
        for (int s = 0; s <= nowstage; ++s) {
            const Real coef_a_ex_dt = dt_ * CoefAEx(nowstage + 1, s);
            ForEachIndex(box, [&](auto... idx) {
                const std::size_t flat = FlatFieldIndex(idx...);
                q(idx...) += coef_a_ex_dt * TendBufExFlat(flat, varID, s);
            });
        }
    } else {  // IMEX
        for (int s = 0; s <= nowstage; ++s) {
            const Real coef_a_ex_dt = dt_ * CoefAEx(nowstage + 1, s);
            const Real coef_a_im_dt = dt_ * CoefAIm(nowstage + 1, s);
            ForEachIndex(box, [&](auto... idx) {
                const std::size_t flat = FlatFieldIndex(idx...);
                q(idx...) += coef_a_ex_dt * TendBufExFlat(flat, varID, s) + coef_a_im_dt * TendBufImFlat(flat, varID, s);
            });
        }
    }
}

// ---------------------------------------------------------------------------
// Advance (variable list)
// ---------------------------------------------------------------------------

template <int Rank>
void TimeIntRK<Rank>::AdvanceVarlist(int nowstage, const std::vector<ArrayView<Real, Rank>>& vars,
                                    const std::vector<int>& varIDs, const BoxRange<Rank>& box)
{
    const int n = static_cast<int>(varIDs.size());
    if (low_storage_flag_) {
        int i = 0;
        while (i < n) {
            if (i + 1 < n) {
                AdvanceLowStoragePair(nowstage, vars[i], vars[i + 1], varIDs[i], varIDs[i + 1], box);
                i += 2;
            } else {
                AdvanceLowStorage(nowstage, vars[i], varIDs[i], box);
                i += 1;
            }
        }
    } else {
        for (int i = 0; i < n; ++i) {
            Advance(nowstage, vars[i], varIDs[i], box);
        }
    }
}

// ---------------------------------------------------------------------------
// Advance (density-weighted tracer variable)
// ---------------------------------------------------------------------------

template <int Rank>
void TimeIntRK<Rank>::AdvanceTracerVar(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box,
                                      ArrayView<const Real, Rank> ddens, ArrayView<const Real, Rank> ddens0,
                                      ArrayView<const Real, Rank> dens_hyd)
{
    if (imex_flag_) {
        throw std::logic_error("TimeIntRK::AdvanceTracerVar: IMEX schemes are not supported");
    }
    if (low_storage_flag_) {
        AdvanceTracerVarLowStorage(nowstage, q, varID, box, ddens, ddens0, dens_hyd);
    } else {
        AdvanceTracerVarGeneral(nowstage, q, varID, box, ddens, ddens0, dens_hyd);
    }
}

template <int Rank>
void TimeIntRK<Rank>::AdvanceTracerVarLowStorage(int nowstage, ArrayView<Real, Rank> q, int varID,
                                                const BoxRange<Rank>& box, ArrayView<const Real, Rank> ddens,
                                                ArrayView<const Real, Rank> ddens0,
                                                ArrayView<const Real, Rank> dens_hyd)
{
    const Real sig_ss = CoefSigEx(nowstage + 1, nowstage);
    const Real sig_Ns = CoefSigEx(nstage_, nowstage);
    const Real one_minus_sig_ss = 1.0 - sig_ss;
    const Real gam_ss = dt_ * CoefGamEx(nowstage + 1, nowstage);
    const Real gam_Ns = dt_ * CoefGamEx(nstage_, nowstage);
    const Real c_ssm1 = CoefCEx(nowstage);

    if (nowstage == nstage_ - 1) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            const Real dens_ssm1 = dens_hyd(idx...) + ddens0(idx...) + c_ssm1 * (ddens(idx...) - ddens0(idx...));
            q(idx...) = (VarTmpFlat(flat, varID) + sig_ss * q(idx...) * dens_ssm1 +
                       gam_ss * TendBufExFlat(flat, varID, 0)) /
                      (dens_hyd(idx...) + ddens(idx...));
        });
        return;
    }

    const Real c_ss = CoefCEx(nowstage + 1);

    if (nowstage == 0) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            Var0Flat(flat, varID) = q(idx...) * (dens_hyd(idx...) + ddens0(idx...));
            VarTmpFlat(flat, varID) = 0.0;
        });
    }

    if (std::abs(sig_Ns) > kProvisionalEps || std::abs(gam_Ns) > kProvisionalEps) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            const Real dens_ssm1 = dens_hyd(idx...) + ddens0(idx...) + c_ssm1 * (ddens(idx...) - ddens0(idx...));
            VarTmpFlat(flat, varID) += sig_Ns * q(idx...) * dens_ssm1 + gam_Ns * TendBufExFlat(flat, varID, 0);
        });
    }

    ForEachIndex(box, [&](auto... idx) {
        const std::size_t flat = FlatFieldIndex(idx...);
        const Real dens_ssm1 = dens_hyd(idx...) + ddens0(idx...) + c_ssm1 * (ddens(idx...) - ddens0(idx...));
        const Real dens_ss = dens_hyd(idx...) + ddens0(idx...) + c_ss * (ddens(idx...) - ddens0(idx...));
        q(idx...) = (one_minus_sig_ss * Var0Flat(flat, varID) + sig_ss * q(idx...) * dens_ssm1 +
                   gam_ss * TendBufExFlat(flat, varID, 0)) /
                  dens_ss;
    });
}

template <int Rank>
void TimeIntRK<Rank>::AdvanceTracerVarGeneral(int nowstage, ArrayView<Real, Rank> q, int varID,
                                             const BoxRange<Rank>& box, ArrayView<const Real, Rank> ddens,
                                             ArrayView<const Real, Rank> ddens0,
                                             ArrayView<const Real, Rank> dens_hyd)
{
    const int tintbuf_ind = tend_buf_indmap_[nowstage];

    if (nstage_ == 1) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            VarTmpFlat(flat, varID) = q(idx...) * (dens_hyd(idx...) + ddens0(idx...));
        });
    }

    if (nowstage == nstage_ - 1) {
        const Real coef_b_ex_dt = dt_ * CoefBEx(nowstage);
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            q(idx...) = (VarTmpFlat(flat, varID) + coef_b_ex_dt * TendBufExFlat(flat, varID, tintbuf_ind)) /
                      (dens_hyd(idx...) + ddens(idx...));
        });
        return;
    }

    const Real c_ss = CoefCEx(nowstage + 1);

    if (nowstage == 0 && !imex_flag_) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            Var0Flat(flat, varID) = q(idx...) * (dens_hyd(idx...) + ddens0(idx...));
            VarTmpFlat(flat, varID) = Var0Flat(flat, varID);
        });
    }

    {
        const Real coef_b_ex_dt = dt_ * CoefBEx(nowstage);
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            q(idx...) = Var0Flat(flat, varID);
            VarTmpFlat(flat, varID) += coef_b_ex_dt * TendBufExFlat(flat, varID, tintbuf_ind);
        });
    }

    if (tend_buf_size_ == 1 && !imex_flag_) {
        const Real coef_a_ex = CoefAEx(nowstage + 1, nowstage);
        const Real coef_b_ex_dt = dt_ * coef_a_ex;
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            const Real dens_ = (dens_hyd(idx...) + ddens0(idx...)) + coef_a_ex * (ddens(idx...) - ddens0(idx...));
            q(idx...) = (Var0Flat(flat, varID) + coef_b_ex_dt * TendBufExFlat(flat, varID, 0)) / dens_;
        });
    } else if (!imex_flag_) {
        for (int s = 0; s <= nowstage; ++s) {
            const Real coef_a_ex_dt = dt_ * CoefAEx(nowstage + 1, s);
            ForEachIndex(box, [&](auto... idx) {
                const std::size_t flat = FlatFieldIndex(idx...);
                q(idx...) += coef_a_ex_dt * TendBufExFlat(flat, varID, s);
            });
        }
        ForEachIndex(box, [&](auto... idx) {
            const Real dens_ = dens_hyd(idx...) + ddens0(idx...) + c_ss * (ddens(idx...) - ddens0(idx...));
            q(idx...) = q(idx...) / dens_;
        });
    }
    // imex_flag_ == true is unreachable here: AdvanceTracerVar() rejects it
    // up front, mirroring the Fortran source's Advance_trcvar.
}

// ---------------------------------------------------------------------------
// StoreVar0 / StoreImplicit
// ---------------------------------------------------------------------------

template <int Rank>
void TimeIntRK<Rank>::StoreVar0(ArrayView<const Real, Rank> q, int varID, const BoxRange<Rank>& box)
{
    ForEachIndex(box, [&](auto... idx) {
        Var0Flat(FlatFieldIndex(idx...), varID) = q(idx...);
    });
}

template <int Rank>
void TimeIntRK<Rank>::StoreImplicit(int nowstage, ArrayView<Real, Rank> q, int varID, const BoxRange<Rank>& box)
{
    if (!imex_flag_) {
        return;
    }

    const int tintbuf_ind = tend_buf_indmap_[nowstage];
    const Real coef_a_im_dt = dt_ * CoefAIm(nowstage, nowstage);

    if (nowstage == 0) {
        ForEachIndex(box, [&](auto... idx) {
            const std::size_t flat = FlatFieldIndex(idx...);
            Var0Flat(flat, varID) = q(idx...);
            VarTmpFlat(flat, varID) = q(idx...);
        });
    }

    ForEachIndex(box, [&](auto... idx) {
        const std::size_t flat = FlatFieldIndex(idx...);
        q(idx...) += coef_a_im_dt * TendBufImFlat(flat, varID, tintbuf_ind);
    });
}

// ---------------------------------------------------------------------------
// Explicit instantiation for the ranks used by SCALE-DG (1D/2D/3D fields)
// ---------------------------------------------------------------------------

template class TimeIntRK<1>;
template class TimeIntRK<2>;
template class TimeIntRK<3>;

}  // namespace common
}  // namespace FElib
