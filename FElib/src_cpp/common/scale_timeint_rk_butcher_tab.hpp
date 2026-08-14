/*
 * @par Description
 * Butcher tableau (and Shu-Osher form, for low-storage SSP schemes) for the
 * Runge-Kutta time integration schemes used by TimeIntRK.
 *
 * Notation rules (mirrors FElib/src/common/scale_timeint_rk_butcher_tab.F90):
 *   * ERK schemes
 *     ERK(_abbrev1_)MsPo(_abbrev2) : explicit RK scheme with M-stage and P-th order.
 *       ERK_1s1o (, ERK_Euler) : Euler scheme
 *       ERK_4s4o (, ERK_RK4)   : classical 4 stage and 4th-order RK scheme
 *       ERK_SSP_MsPo : strong stability preserving (SSP) (M,P) RK scheme
 *   * IMEX schemes (ERK for the non-stiff part, DIRK for the stiff part)
 *     IMEX_abbrev : abbreviation follows the table in Vogl et al. (2019).
 *
 * @par Reference
 *  - Vogl et al. 2019: Evaluation of implicit-explicit additive Runge-Kutta
 *    integrators for the HOMME-NH dynamical core. JAMES, 11, 4228-4244.
 *  - Higueras and Roldan, 2019: New third order low-storage SSP explicit
 *    Runge-Kutta methods. J. Sci. Comput., 79(3), 1882-1906.
 *
 * @author Yuta Kawai, Team SCALE
 */
#pragma once

#include <string_view>
#include <vector>

#include "utility/array_view.hpp"

namespace FElib {
namespace common {
namespace timeint_rk_butcher_tab {

using utility::ArrayView;
using Real = double;

//! Static shape/dispatch info for an RK scheme, looked up by name.
struct SchemeInfo {
    int nstage;
    int tend_buf_size;
    bool low_storage_flag;
    bool imex_flag;
};

//! @throws std::invalid_argument if scheme_name is not a supported scheme.
SchemeInfo GetInfo(std::string_view scheme_name);

//! Coefficient tables for an RK scheme: both the Butcher form (a, b, c) and,
//! for low-storage SSP schemes, the Shu-Osher form (sig, gam) it was derived
//! from. All indices are 0-based (row = stage index, col = stage index).
struct ButcherTable {
    int nstage = 0;

    std::vector<Real> coef_a_ex;    //< shape (nstage, nstage), explicit part
    std::vector<Real> coef_b_ex;    //< shape (nstage)
    std::vector<Real> coef_c_ex;    //< shape (nstage)

    std::vector<Real> coef_sig_ex;  //< shape (nstage+1, nstage), Shu-Osher form
    std::vector<Real> coef_gam_ex;  //< shape (nstage+1, nstage)

    std::vector<Real> coef_a_im;    //< shape (nstage, nstage), implicit part (only meaningful if imex)
    std::vector<Real> coef_b_im;    //< shape (nstage)
    std::vector<Real> coef_c_im;    //< shape (nstage)

    //! 0-based index into a stage's tendency-buffer slot, one per stage.
    std::vector<int> tend_buf_indmap;  //< shape (nstage)

    ArrayView<const Real, 2> CoefAEx() const { return ArrayView<const Real, 2>(coef_a_ex.data(), nstage, nstage); }
    ArrayView<const Real, 2> CoefSigEx() const { return ArrayView<const Real, 2>(coef_sig_ex.data(), nstage + 1, nstage); }
    ArrayView<const Real, 2> CoefGamEx() const { return ArrayView<const Real, 2>(coef_gam_ex.data(), nstage + 1, nstage); }
    ArrayView<const Real, 2> CoefAIm() const { return ArrayView<const Real, 2>(coef_a_im.data(), nstage, nstage); }
};

//! @throws std::invalid_argument if scheme_name is not a supported scheme.
ButcherTable Get(std::string_view scheme_name);

}  // namespace timeint_rk_butcher_tab
}  // namespace common
}  // namespace FElib
