#include "scale_timeint_rk_butcher_tab.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace FElib {
namespace common {
namespace timeint_rk_butcher_tab {

namespace {

// LAPACK's LU factorization/solve, as linked via SCALE_MATHLIB_LIBS.
// Used only by SolveLinEq() below to convert the Shu-Osher form into the
// Butcher form (ShuOsher2Butcher). scale_linalgebra.F90 (the Fortran module
// that provides this on the Fortran side) has not been ported to C++ yet, so
// this is a small private helper scoped to this file rather than a call into
// a shared FElib::common::linalgebra module -- see PORTING_STATUS.md.
extern "C" {
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a, const int* lda,
             const int* ipiv, double* b, const int* ldb, int* info);
}

//! Solve A X = B for X, where A is n x n and B/X are n x nrhs.
//! A, B, X all use the row-major ArrayView convention (A(i,j), not LAPACK's
//! native column-major layout), so this builds column-major scratch copies
//! for the LAPACK call rather than reinterpreting the row-major buffers
//! directly (which would silently solve a transposed system instead).
void SolveLinEq(ArrayView<const Real, 2> A, ArrayView<const Real, 2> B, ArrayView<Real, 2> X)
{
    const int n = static_cast<int>(A.extent(0));
    const int nrhs = static_cast<int>(B.extent(1));

    std::vector<Real> a_col(static_cast<std::size_t>(n) * n);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            a_col[static_cast<std::size_t>(j) * n + i] = A(i, j);
        }
    }

    std::vector<int> ipiv(n);
    int info = 0;
    dgetrf_(&n, &n, a_col.data(), &n, ipiv.data(), &info);
    if (info != 0) {
        throw std::runtime_error("timeint_rk_butcher_tab::SolveLinEq: matrix is singular (dgetrf info=" +
                                  std::to_string(info) + ")");
    }

    std::vector<Real> b_col(static_cast<std::size_t>(n) * nrhs);
    for (int j = 0; j < nrhs; ++j) {
        for (int i = 0; i < n; ++i) {
            b_col[static_cast<std::size_t>(j) * n + i] = B(i, j);
        }
    }

    dgetrs_("N", &n, &nrhs, a_col.data(), &n, ipiv.data(), b_col.data(), &n, &info);
    if (info != 0) {
        throw std::runtime_error("timeint_rk_butcher_tab::SolveLinEq: dgetrs failed (info=" +
                                  std::to_string(info) + ")");
    }

    for (int j = 0; j < nrhs; ++j) {
        for (int i = 0; i < n; ++i) {
            X(i, j) = b_col[static_cast<std::size_t>(j) * n + i];
        }
    }
}

//! Convert the Shu-Osher form (sig, gam) into the Butcher form (a, b).
//! For details of the strategy, see Higueras and Roldan (2019).
void ShuOsher2Butcher(ArrayView<const Real, 2> sig, ArrayView<const Real, 2> gam, int nstage,
                      ArrayView<Real, 2> a, std::vector<Real>& b)
{
    // L = I - SIG(0:nstage-1, :)
    std::vector<Real> L_buf(static_cast<std::size_t>(nstage) * nstage);
    ArrayView<Real, 2> L(L_buf.data(), nstage, nstage);
    for (int i = 0; i < nstage; ++i) {
        for (int j = 0; j < nstage; ++j) {
            L(i, j) = -sig(i, j);
        }
        L(i, i) += 1.0;
    }

    // GAM(0:nstage-1, :) -- the same row range used for L above.
    std::vector<Real> gam_top_buf(static_cast<std::size_t>(nstage) * nstage);
    ArrayView<Real, 2> gam_top(gam_top_buf.data(), nstage, nstage);
    for (int i = 0; i < nstage; ++i) {
        for (int j = 0; j < nstage; ++j) {
            gam_top(i, j) = gam(i, j);
        }
    }

    // A = L^-1 GAM
    SolveLinEq(ArrayView<const Real, 2>(L.data(), nstage, nstage),
               ArrayView<const Real, 2>(gam_top.data(), nstage, nstage), a);

    // b(n) = GAM(nstage, n) + sum_k SIG(nstage, k) * A(k, n)
    for (int n = 0; n < nstage; ++n) {
        Real s = gam(nstage, n);
        for (int k = 0; k < nstage; ++k) {
            s += sig(nstage, k) * a(k, n);
        }
        b[n] = s;
    }
}

}  // namespace

SchemeInfo GetInfo(std::string_view scheme_name)
{
    if (scheme_name == "ERK_1s1o" || scheme_name == "ERK_Euler") {
        return {1, 1, false, false};
    }
    if (scheme_name == "ERK_4s4o" || scheme_name == "ERK_RK4") {
        return {4, 1, false, false};
    }
    // Shu and Osher, 1998: Efficient implementation of essentially
    // nonoscillatory shock-capturing schemes, J. Comput. Phys.;
    // Gottlieb and Shu, 1998: Total variation diminishing Runge-Kutta
    // schemes, Math. Comput.
    if (scheme_name == "ERK_SSP_2s2o") {
        return {2, 1, true, false};
    }
    // Shu and Osher, 1998. SSP coefficient: 1, effective SSP coefficient Ceff: 1/3.
    if (scheme_name == "ERK_SSP_3s3o") {
        return {3, 1, true, false};
    }
    if (scheme_name == "ERK_SSP_4s3o") {
        return {4, 1, true, false};
    }
    // Higueras and Roldan, 2018: New third order low-storage SSP explicit
    // Runge-Kutta methods.
    if (scheme_name == "ERK_SSP_5s3o_2N2*") {
        return {5, 1, true, false};
    }
    // Ketcheson, 2008: Highly efficient strong stability-preserving
    // Runge-Kutta methods with low-storage implementations, SIAM J. Sci.
    // Comput. SSP coefficient: 6, effective SSP coefficient Ceff: 0.6.
    if (scheme_name == "ERK_SSP_10s4o_2N") {
        return {10, 1, true, false};
    }
    // Giraldo et al. (2013): Implicit-Explicit Formulations of a
    // Three-Dimensional Nonhydrostatic Unified Model of the Atmosphere
    // (NUMA), SIAM J. Sci. Comp.
    if (scheme_name == "IMEX_ARK232") {
        return {3, 3, false, true};
    }
    // Kennedy and Carpenter, 2003: Additive Runge-Kutta schemes for
    // convection-diffusion-reaction equations, Appl. Numer. Math.
    if (scheme_name == "IMEX_ARK324") {
        return {4, 4, false, true};
    }

    throw std::invalid_argument("timeint_rk_butcher_tab::GetInfo: unsupported scheme \"" +
                                 std::string(scheme_name) + "\"");
}

ButcherTable Get(std::string_view scheme_name)
{
    const SchemeInfo info = GetInfo(scheme_name);
    const int nstage = info.nstage;

    ButcherTable t;
    t.nstage = nstage;
    t.coef_a_ex.assign(static_cast<std::size_t>(nstage) * nstage, 0.0);
    t.coef_b_ex.assign(nstage, 0.0);
    t.coef_c_ex.assign(nstage, 0.0);
    t.coef_sig_ex.assign(static_cast<std::size_t>(nstage + 1) * nstage, 0.0);
    t.coef_gam_ex.assign(static_cast<std::size_t>(nstage + 1) * nstage, 0.0);
    t.coef_a_im.assign(static_cast<std::size_t>(nstage) * nstage, 0.0);
    t.coef_b_im.assign(nstage, 0.0);
    t.coef_c_im.assign(nstage, 0.0);
    t.tend_buf_indmap.assign(nstage, 0);

    ArrayView<Real, 2> a_ex(t.coef_a_ex.data(), nstage, nstage);
    ArrayView<Real, 2> sig_ex(t.coef_sig_ex.data(), nstage + 1, nstage);
    ArrayView<Real, 2> gam_ex(t.coef_gam_ex.data(), nstage + 1, nstage);
    ArrayView<Real, 2> a_im(t.coef_a_im.data(), nstage, nstage);

    bool call_shu_osher_to_butcher = false;

    if (scheme_name == "ERK_1s1o" || scheme_name == "ERK_Euler") {
        t.coef_b_ex[0] = 1.0;

        sig_ex(1, 0) = 1.0;
        gam_ex(1, 0) = 1.0;

    } else if (scheme_name == "ERK_4s4o" || scheme_name == "ERK_RK4") {
        a_ex(1, 0) = 0.5;
        a_ex(2, 1) = 0.5;
        a_ex(3, 2) = 1.0;
        t.coef_b_ex = {1.0 / 6.0, 2.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0};

    } else if (scheme_name == "ERK_SSP_2s2o") {
        sig_ex(1, 0) = 1.0;
        sig_ex(2, 0) = 0.5;
        sig_ex(2, 1) = 0.5;
        gam_ex(1, 0) = 1.0;
        gam_ex(2, 1) = 0.5;

        call_shu_osher_to_butcher = true;

    } else if (scheme_name == "ERK_SSP_3s3o") {
        sig_ex(1, 0) = 1.0;
        sig_ex(2, 0) = 3.0 / 4.0;
        sig_ex(2, 1) = 1.0 / 4.0;
        sig_ex(3, 0) = 1.0 / 3.0;
        sig_ex(3, 1) = 0.0;
        sig_ex(3, 2) = 2.0 / 3.0;
        gam_ex(1, 0) = 1.0;
        gam_ex(2, 1) = 1.0 / 4.0;
        gam_ex(3, 2) = 2.0 / 3.0;

        call_shu_osher_to_butcher = true;

    } else if (scheme_name == "ERK_SSP_4s3o") {
        sig_ex(1, 0) = 1.0;
        sig_ex(2, 0) = 0.0;
        sig_ex(2, 1) = 1.0;
        sig_ex(3, 0) = 2.0 / 3.0;
        sig_ex(3, 1) = 0.0;
        sig_ex(3, 2) = 1.0 / 3.0;
        sig_ex(4, 0) = 0.0;
        sig_ex(4, 1) = 0.0;
        sig_ex(4, 2) = 0.0;
        sig_ex(4, 3) = 1.0;
        gam_ex(1, 0) = 0.5;
        gam_ex(2, 1) = 0.5;
        gam_ex(3, 2) = 1.0 / 6.0;
        gam_ex(4, 3) = 0.5;

        call_shu_osher_to_butcher = true;

    } else if (scheme_name == "ERK_SSP_5s3o_2N2*") {
        sig_ex(1, 0) = 1.0;
        sig_ex(2, 0) = 0.0;
        sig_ex(2, 1) = 1.0;
        sig_ex(3, 0) = 0.682342861037239;
        sig_ex(3, 1) = 0.0;
        sig_ex(3, 2) = 0.317657138962761;
        sig_ex(4, 0) = 0.0;
        sig_ex(4, 1) = 0.0;
        sig_ex(4, 2) = 0.0;
        sig_ex(4, 3) = 1.0;
        sig_ex(5, 0) = 0.045230974482400;
        sig_ex(5, 1) = 0.0;
        sig_ex(5, 2) = 0.0;
        sig_ex(5, 3) = 0.0;
        sig_ex(5, 4) = 0.954769025517600;
        gam_ex(1, 0) = 0.465388589249323;
        gam_ex(2, 1) = 0.465388589249323;
        gam_ex(3, 2) = 0.124745797313998;
        gam_ex(4, 3) = 0.465388589249323;
        gam_ex(5, 4) = 0.154263303748666;

        call_shu_osher_to_butcher = true;

    } else if (scheme_name == "ERK_SSP_10s4o_2N") {
        for (int m : {0, 1, 2, 3}) {
            sig_ex(m + 1, m) = 1.0;
            gam_ex(m + 1, m) = 1.0 / 6.0;
        }
        sig_ex(5, 0) = 3.0 / 5.0;
        sig_ex(5, 1) = 0.0;
        sig_ex(5, 2) = 0.0;
        sig_ex(5, 3) = 0.0;
        sig_ex(5, 4) = 2.0 / 5.0;
        gam_ex(5, 4) = 1.0 / 15.0;
        for (int m : {5, 6, 7, 8}) {
            sig_ex(m + 1, m) = 1.0;
            gam_ex(m + 1, m) = 1.0 / 6.0;
        }
        sig_ex(10, 0) = 0.2 * 0.2;
        sig_ex(10, 4) = 1.8 * 0.2;
        sig_ex(10, 9) = 3.0 * 0.2;
        gam_ex(10, 4) = 1.8 / 30.0;
        gam_ex(10, 9) = 3.0 / 30.0;

        call_shu_osher_to_butcher = true;

    } else if (scheme_name == "IMEX_ARK232") {
        const Real alp = (3.0 + 2.0 * std::sqrt(2.0)) / 6.0;
        const Real gam = 1.0 - 1.0 / std::sqrt(2.0);
        const Real del = 1.0 / (2.0 * std::sqrt(2.0));

        a_ex(1, 0) = 2.0 * gam;
        a_ex(2, 0) = 1.0 - alp;
        a_ex(2, 1) = alp;
        a_ex(2, 2) = 0.0;
        t.coef_b_ex = {del, del, gam};

        a_im(1, 0) = gam;
        a_im(1, 1) = gam;
        a_im(1, 2) = 0.0;
        a_im(2, 0) = del;
        a_im(2, 1) = del;
        a_im(2, 2) = gam;
        t.coef_b_im = t.coef_b_ex;

        t.tend_buf_indmap = {0, 1, 2};

    } else if (scheme_name == "IMEX_ARK324") {
        a_ex(1, 0) = 1767732205903.0 / 2027836641118.0;
        a_ex(2, 0) = 5535828885825.0 / 10492691773637.0;
        a_ex(2, 1) = 788022342437.0 / 10882634858940.0;
        a_ex(3, 0) = 6485989280629.0 / 16251701735622.0;
        a_ex(3, 1) = -4246266847089.0 / 9704473918619.0;
        a_ex(3, 2) = 10755448449292.0 / 10357097424841.0;
        t.coef_b_ex[0] = 1471266399579.0 / 7840856788654.0;
        t.coef_b_ex[1] = -4482444167858.0 / 7529755066697.0;
        t.coef_b_ex[2] = 11266239266428.0 / 11593286722821.0;
        t.coef_b_ex[3] = 1767732205903.0 / 4055673282236.0;

        a_im(1, 0) = 1767732205903.0 / 4055673282236.0;
        a_im(1, 1) = 1767732205903.0 / 4055673282236.0;
        a_im(2, 0) = 2746238789719.0 / 10658868560708.0;
        a_im(2, 1) = -640167445237.0 / 6845629431997.0;
        a_im(2, 2) = 1767732205903.0 / 4055673282236.0;
        for (int j = 0; j < 4; ++j) {
            a_im(3, j) = t.coef_b_ex[j];
        }
        t.coef_b_im = t.coef_b_ex;

        t.tend_buf_indmap = {0, 1, 2, 3};

    } else {
        throw std::invalid_argument("timeint_rk_butcher_tab::Get: unsupported scheme \"" +
                                     std::string(scheme_name) + "\"");
    }

    if (call_shu_osher_to_butcher) {
        ShuOsher2Butcher(ArrayView<const Real, 2>(t.coef_sig_ex.data(), nstage + 1, nstage),
                          ArrayView<const Real, 2>(t.coef_gam_ex.data(), nstage + 1, nstage), nstage, a_ex,
                          t.coef_b_ex);
    }

    // The Fortran source recomputes coef_c_ex from coef_a_ex a second time
    // under `if (imex_flag)`, but that second computation is byte-for-byte
    // identical to this one (both read coef_a_ex, not coef_a_im), so it is
    // not repeated here. coef_c_im is consequently never actually populated
    // by this routine in either language -- it stays all-zero, which we
    // preserve rather than "fixing" (see
    // FElib/src/common/scale_timeint_rk_butcher_tab.F90).
    for (int n = 0; n < nstage; ++n) {
        Real s = 0.0;
        for (int j = 0; j <= n; ++j) {
            s += a_ex(n, j);
        }
        t.coef_c_ex[n] = s;
    }

    return t;
}

}  // namespace timeint_rk_butcher_tab
}  // namespace common
}  // namespace FElib
