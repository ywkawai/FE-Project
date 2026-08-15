// Unit test of scale_timeint_rk / scale_timeint_rk_butcher_tab
//
// Mirrors FElib/test/common/timeint_rk/test_timeint_rk.f90 (harmonic-oscillator
// convergence-rate checks for every ERK/IMEX scheme, plus the Shu-Osher ->
// Butcher conversion check), and adds coverage the Fortran test does not
// exercise: that the Rank-templated update (the C++ replacement for the
// erb-generated d=1,2,3 duplication) gives identical, per-point-independent
// results across Rank=1/2/3 and across the plain-Advance vs. AdvanceVarlist
// (paired low-storage) code paths, and an exact mass-conservation identity
// for AdvanceTracerVar.

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "cpp_test_harness.hpp"
#include "scale_timeint_rk.hpp"
#include "scale_timeint_rk_butcher_tab.hpp"

using FElib::common::TimeIntRK;
using FElib::utility::ArrayView;
using FElib::utility::BoxRange;
using FElib::utility::ForEachIndex;
using namespace FElib::test;

namespace {

template <int Rank, typename T>
ArrayView<T, Rank> MakeView(T* data, const std::array<std::size_t, Rank>& extents)
{
    return std::apply([&](auto... e) { return ArrayView<T, Rank>(data, e...); }, extents);
}

template <int Rank>
BoxRange<Rank> FullBox(const std::array<std::size_t, Rank>& extents)
{
    BoxRange<Rank> box;
    box.lo.fill(0);
    box.hi = extents;
    return box;
}

constexpr int kIdU = 0;
constexpr int kIdV = 1;

// ---------------------------------------------------------------------------
// A. Convergence-rate checks for the harmonic oscillator, ported from
//    test_timeint_rk.f90's calc_sol_ex / calc_sol_imex / check_terror.
// ---------------------------------------------------------------------------

constexpr double kEndTime = 2.0;
constexpr double kDt1 = 0.02;
constexpr int kSaveErrorDstep1 = 1;
constexpr double kDt2 = 0.0025;
constexpr int kSaveErrorDstep2 = 8;
constexpr int kErrorArraySize = 100;
const double kPi = std::acos(-1.0);

std::vector<double> CalcSolEx(std::string_view scheme, double dt, int save_error_dstep)
{
    TimeIntRK<1> tint;
    tint.Init(scheme, dt, 2, {1});

    const double omg = 2.0 * kPi;
    const int nstep = static_cast<int>(kEndTime / dt);

    double u = 1.0, v = 0.0;
    ArrayView<double, 1> uView(&u, 1);
    ArrayView<double, 1> vView(&v, 1);
    const BoxRange<1> box = FullBox<1>({1});

    std::vector<double> errors;
    for (int n = 1; n <= nstep; ++n) {
        for (int stage = 0; stage < tint.Nstage(); ++stage) {
            const int tb = tint.TendBufIndex(stage);
            tint.TendBufEx(kIdU, tb, {0}) = omg * v;
            tint.TendBufEx(kIdV, tb, {0}) = -omg * u;
            tint.Advance(stage, uView, kIdU, box);
            tint.Advance(stage, vView, kIdV, box);
        }
        if (n % save_error_dstep == 0) {
            const double ans = std::cos(omg * n * dt);
            errors.push_back(u - ans);
        }
    }
    return errors;
}

std::vector<double> CalcSolImex(std::string_view scheme, double dt, int save_error_dstep)
{
    TimeIntRK<1> tint;
    tint.Init(scheme, dt, 2, {1});

    const double omg = 2.0 * kPi;
    const double r = 0.1;
    const double omgg = std::sqrt(omg * omg - r * r / 4.0);
    const int nstep = static_cast<int>(kEndTime / dt);

    double u = 1.0, v = 0.0;
    ArrayView<double, 1> uView(&u, 1);
    ArrayView<double, 1> vView(&v, 1);
    const BoxRange<1> box = FullBox<1>({1});

    std::vector<double> errors;
    for (int n = 1; n <= nstep; ++n) {
        for (int stage = 0; stage < tint.Nstage(); ++stage) {
            const int tb = tint.TendBufIndex(stage);
            const double implicit_fac = tint.GetImplicitDiagFac(stage);
            const double coef = implicit_fac * omg;

            if (std::abs(implicit_fac) > 0.0) {
                const double ui = (u + coef * v) / (1.0 + coef * coef);
                const double vi = (v - coef * u) / (1.0 + coef * coef);
                tint.TendBufIm(kIdU, tb, {0}) = (ui - u) / implicit_fac;
                tint.TendBufIm(kIdV, tb, {0}) = (vi - v) / implicit_fac;
            } else {
                tint.TendBufIm(kIdU, tb, {0}) = omg * v;
                tint.TendBufIm(kIdV, tb, {0}) = -omg * u;
            }

            tint.StoreImplicit(stage, uView, kIdU, box);
            tint.StoreImplicit(stage, vView, kIdV, box);

            tint.TendBufEx(kIdU, tb, {0}) = -r * u;
            tint.TendBufEx(kIdV, tb, {0}) = 0.0;

            tint.Advance(stage, uView, kIdU, box);
            tint.Advance(stage, vView, kIdV, box);
        }
        if (n % save_error_dstep == 0) {
            const double ans =
                (std::cos(omgg * n * dt) - r / (2.0 * omgg) * std::sin(omgg * n * dt)) * std::exp(-0.5 * r * n * dt);
            errors.push_back(u - ans);
        }
    }
    return errors;
}

double Rms(const std::vector<double>& e)
{
    double s = 0.0;
    for (double v : e) s += v * v;
    return std::sqrt(s);
}

void CheckConvergenceOrder(const std::string& tname, const std::vector<double>& err1, const std::vector<double>& err2,
                          double dt_ratio, double check_ord)
{
    const double e1 = Rms(err1);
    const double e2 = Rms(err2);
    const double conv_rate = std::log(e1 / e2) / std::log(dt_ratio);
    std::cout << "  " << tname << ": conv rate=" << conv_rate << ", error1=" << e1 << ", error2=" << e2 << "\n";
    Check(conv_rate >= check_ord, tname + ": convergence rate >= " + std::to_string(check_ord));
}

void CheckTschemeEx(std::string_view scheme, double check_ord)
{
    const auto err1 = CalcSolEx(scheme, kDt1, kSaveErrorDstep1);
    const auto err2 = CalcSolEx(scheme, kDt2, kSaveErrorDstep2);
    CheckConvergenceOrder(std::string(scheme), err1, err2, kDt1 / kDt2, check_ord);
}

void CheckTschemeImex(std::string_view scheme, double check_ord)
{
    const auto err1 = CalcSolImex(scheme, kDt1, kSaveErrorDstep1);
    const auto err2 = CalcSolImex(scheme, kDt2, kSaveErrorDstep2);
    CheckConvergenceOrder(std::string(scheme), err1, err2, kDt1 / kDt2, check_ord);
}

void TestConvergenceRates()
{
    CheckTschemeEx("ERK_1s1o", 0.98);
    CheckTschemeEx("ERK_4s4o", 3.98);
    CheckTschemeEx("ERK_SSP_2s2o", 1.98);
    CheckTschemeEx("ERK_SSP_3s3o", 2.98);
    CheckTschemeEx("ERK_SSP_4s3o", 2.98);
    CheckTschemeEx("ERK_SSP_5s3o_2N2*", 2.98);
    CheckTschemeEx("ERK_SSP_10s4o_2N", 3.98);

    CheckTschemeImex("IMEX_ARK232", 1.98);
    CheckTschemeImex("IMEX_ARK324", 2.98);
}

// ---------------------------------------------------------------------------
// B. Shu-Osher -> Butcher conversion, ported from check_ShuOsher2Butcher /
//    check_ButcherMat.
// ---------------------------------------------------------------------------

void CheckButcherMat(const std::string& label, const FElib::common::timeint_rk_butcher_tab::ButcherTable& table,
                    const std::vector<double>& a_ans, const std::vector<double>& b_ans)
{
    const int nstage = table.nstage;
    double a_err = 0.0;
    for (int i = 0; i < nstage; ++i) {
        for (int j = 0; j < nstage; ++j) {
            const double d = table.CoefAEx()(i, j) - a_ans[static_cast<std::size_t>(i) * nstage + j];
            a_err += d * d;
        }
    }
    a_err = std::sqrt(a_err);
    Check(a_err <= 1.0e-14 * nstage * nstage, label + ": Butcher A matches Shu-Osher-derived A");

    double b_err = 0.0;
    for (int i = 0; i < nstage; ++i) {
        const double d = table.coef_b_ex[i] - b_ans[i];
        b_err += d * d;
    }
    b_err = std::sqrt(b_err);
    Check(b_err <= 1.0e-14 * nstage, label + ": Butcher b matches Shu-Osher-derived b");
}

void TestShuOsherToButcher()
{
    using FElib::common::timeint_rk_butcher_tab::Get;

    {
        const auto t = Get("ERK_SSP_3s3o");
        // clang-format off
        const std::vector<double> a_ans = {
            0.0,       0.0,       0.0,
            1.0,       0.0,       0.0,
            1.0/4.0,   1.0/4.0,   0.0,
        };
        // clang-format on
        const std::vector<double> b_ans = {1.0 / 6.0, 1.0 / 6.0, 4.0 / 6.0};
        CheckButcherMat("ERK_SSP_3s3o", t, a_ans, b_ans);
    }
    {
        const auto t = Get("ERK_SSP_4s3o");
        // clang-format off
        const std::vector<double> a_ans = {
            0.0,       0.0,       0.0,       0.0,
            1.0/2.0,   0.0,       0.0,       0.0,
            1.0/2.0,   1.0/2.0,   0.0,       0.0,
            1.0/6.0,   1.0/6.0,   1.0/6.0,   0.0,
        };
        // clang-format on
        const std::vector<double> b_ans = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 3.0 / 6.0};
        CheckButcherMat("ERK_SSP_4s3o", t, a_ans, b_ans);
    }
}

// ---------------------------------------------------------------------------
// C. Cross-Rank / AdvanceVarlist consistency.
//
// The harmonic oscillator has no spatial coupling, so if every grid point
// starts from the same (u, v) and is fed the same tendency, every point must
// stay numerically identical for the run's whole duration -- for any Rank,
// and whether Advance() is called once per variable or AdvanceVarlist() is
// used instead. This is exactly the property a Rank/indexing bug (e.g. a
// wrong stride in FlatFieldIndex, or points aliasing across the pair in
// AdvanceLowStoragePair) would break.
// ---------------------------------------------------------------------------

template <int Rank>
void SetUniform(TimeIntRK<Rank>& tint, int varID, int stage, double value, const BoxRange<Rank>& box)
{
    ForEachIndex(box, [&](auto... idx) {
        tint.TendBufEx(varID, stage, std::array<std::size_t, Rank>{static_cast<std::size_t>(idx)...}) = value;
    });
}

// Runs the (undamped) harmonic oscillator on a multi-point field where every
// point shares the same initial condition, using either Advance() per
// variable or AdvanceVarlist() for both variables together. Returns the
// final u at every point.
template <int Rank>
std::vector<double> RunOscillatorMultiPoint(std::string_view scheme, double dt, int nstep,
                                            const std::array<std::size_t, Rank>& extents, bool use_varlist)
{
    TimeIntRK<Rank> tint;
    tint.Init(scheme, dt, 2, extents);

    std::size_t n = 1;
    for (auto e : extents) n *= e;

    const double omg = 2.0 * kPi;
    std::vector<double> u(n, 1.0), v(n, 0.0);
    auto uView = MakeView<Rank>(u.data(), extents);
    auto vView = MakeView<Rank>(v.data(), extents);
    const BoxRange<Rank> box = FullBox<Rank>(extents);

    for (int step = 0; step < nstep; ++step) {
        for (int stage = 0; stage < tint.Nstage(); ++stage) {
            const int tb = tint.TendBufIndex(stage);
            // All points are identical by construction; point 0 is the representative scalar.
            SetUniform(tint, kIdU, tb, omg * v[0], box);
            SetUniform(tint, kIdV, tb, -omg * u[0], box);

            if (use_varlist) {
                tint.AdvanceVarlist(stage, {uView, vView}, {kIdU, kIdV}, box);
            } else {
                tint.Advance(stage, uView, kIdU, box);
                tint.Advance(stage, vView, kIdV, box);
            }
        }
    }
    return u;
}

void TestMultiPointConsistency()
{
    const std::vector<std::string> schemes = {"ERK_SSP_3s3o", "ERK_4s4o"};  // low-storage and general
    const double dt = 0.05;
    const int nstep = 20;

    for (const auto& scheme : schemes) {
        const auto ref = RunOscillatorMultiPoint<1>(scheme, dt, nstep, {1}, /*use_varlist=*/false);

        for (bool use_varlist : {false, true}) {
            const auto u2 = RunOscillatorMultiPoint<2>(scheme, dt, nstep, {2, 3}, use_varlist);
            const auto u3 = RunOscillatorMultiPoint<3>(scheme, dt, nstep, {2, 2, 2}, use_varlist);

            // Exact equality is not asserted here: the per-point arithmetic is
            // identical, but Rank=1 (a single-element loop) and Rank=2/3 (a
            // nested nested loop over several points) may be vectorized
            // differently by the compiler, which can shift the last few ULPs.
            const double tol = 1.0e-12 * std::max(1.0, std::abs(ref[0]));
            bool all_match = true;
            for (double x : u2) all_match = all_match && (std::abs(x - ref[0]) <= tol);
            for (double x : u3) all_match = all_match && (std::abs(x - ref[0]) <= tol);

            Check(all_match, scheme + ": Rank=2/3 multi-point results match the Rank=1 reference (varlist=" +
                              std::string(use_varlist ? "true" : "false") + ")");
        }
    }
}

// ForEachIndex only actually forks an OpenMP team when a box has at least
// FElib::utility::kOmpMinPoints points (see index_range.hpp) -- every box
// above is far smaller than that, so none of the checks so far have
// exercised that branch at all. This test uses boxes at/above the threshold
// so the parallel path itself (not just the serial fallback) is verified.
void TestLargeBoxParallelPath()
{
    const double dt = 0.05;
    const int nstep = 5;
    const std::string scheme = "ERK_SSP_3s3o";

    const auto ref = RunOscillatorMultiPoint<1>(scheme, dt, nstep, {1}, /*use_varlist=*/false);
    const double tol = 1.0e-9 * std::max(1.0, std::abs(ref[0]));

    const std::size_t n1 = FElib::utility::kOmpMinPoints * 2;
    const auto u1 = RunOscillatorMultiPoint<1>(scheme, dt, nstep, {n1}, /*use_varlist=*/false);
    bool all_match_1 = true;
    for (double x : u1) all_match_1 = all_match_1 && (std::abs(x - ref[0]) <= tol);
    Check(all_match_1, "Rank=1 box with " + std::to_string(n1) +
                       " points (>= kOmpMinPoints, exercises the OpenMP path) matches the serial reference");

    const std::array<std::size_t, 3> extents3 = {22, 22, 22};  // 10648 >= kOmpMinPoints
    const auto u3 = RunOscillatorMultiPoint<3>(scheme, dt, nstep, extents3, /*use_varlist=*/true);
    bool all_match_3 = true;
    for (double x : u3) all_match_3 = all_match_3 && (std::abs(x - ref[0]) <= tol);
    Check(all_match_3, "Rank=3 box with 22^3 points (>= kOmpMinPoints, exercises the collapse(2) OpenMP path) "
                       "matches the serial reference");
}

// ---------------------------------------------------------------------------
// D. AdvanceTracerVar mass conservation: with a zero tendency, the
//    density-weighted update degenerates to an exact identity
//    q_final * (dens_hyd + ddens) == q0 * (dens_hyd + ddens0)
//    for ANY scheme and ANY stage count (every coefficient-weighted zero
//    contribution vanishes), so this is an exact check, not asymptotic.
// ---------------------------------------------------------------------------

template <int Rank>
void TestTracerMassConservationForScheme(std::string_view scheme, const std::array<std::size_t, Rank>& extents)
{
    TimeIntRK<Rank> tint;
    const double dt = 0.1;
    tint.Init(scheme, dt, 1, extents);
    constexpr int ID_Q = 0;

    std::size_t n = 1;
    for (auto e : extents) n *= e;

    std::vector<double> q(n), dens_hyd(n), ddens0(n), ddens(n);
    for (std::size_t i = 0; i < n; ++i) {
        q[i] = 1.0 + 0.1 * static_cast<double>(i);
        dens_hyd[i] = 1.0;
        ddens0[i] = 0.05 * static_cast<double>(i);
        ddens[i] = 0.2 + 0.01 * static_cast<double>(i);
    }
    std::vector<double> mass0(n);
    for (std::size_t i = 0; i < n; ++i) mass0[i] = q[i] * (dens_hyd[i] + ddens0[i]);

    auto qView = MakeView<Rank>(q.data(), extents);
    const auto densHydView = MakeView<Rank, const double>(dens_hyd.data(), extents);
    const auto ddens0View = MakeView<Rank, const double>(ddens0.data(), extents);
    const auto ddensView = MakeView<Rank, const double>(ddens.data(), extents);
    const BoxRange<Rank> box = FullBox<Rank>(extents);

    for (int stage = 0; stage < tint.Nstage(); ++stage) {
        const int tb = tint.TendBufIndex(stage);
        SetUniform(tint, ID_Q, tb, 0.0, box);
        tint.AdvanceTracerVar(stage, qView, ID_Q, box, ddensView, ddens0View, densHydView);
    }

    double max_err = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double mass_final = q[i] * (dens_hyd[i] + ddens[i]);
        max_err = std::max(max_err, std::abs(mass_final - mass0[i]));
    }
    Check(max_err < 1.0e-12, std::string(scheme) + ": AdvanceTracerVar conserves mass under zero tendency");
}

}  // namespace

int main()
{
    return RunTestMain("test_timeint_rk_cpp", [] {
        TestConvergenceRates();
        TestShuOsherToButcher();
        TestMultiPointConsistency();
        TestLargeBoxParallelPath();
        TestTracerMassConservationForScheme<1>("ERK_SSP_3s3o", {5});
        TestTracerMassConservationForScheme<2>("ERK_4s4o", {2, 3});
    });
}
