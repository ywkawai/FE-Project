#include "scale_element_base.hpp"

#include <stdexcept>
#include <string>

#include "scale_linalgebra.hpp"

namespace FElib {
namespace element {

namespace {

void CheckSquare(ArrayView<const Real, 2> A, const char* name, const char* who)
{
    if (A.extent(0) != A.extent(1)) {
        throw std::invalid_argument(std::string(who) + ": " + name + " must be square");
    }
}

}  // namespace

void ElementBase::InitBase(bool lumped_mat_flag)
{
    lumped_mat_flag_ = lumped_mat_flag;

    const auto np = static_cast<std::size_t>(Np_);
    const auto nfp_tot = static_cast<std::size_t>(NfpTot_);

    M_.assign(np * np, 0.0);
    invM_.assign(np * np, 0.0);
    V_.assign(np * np, 0.0);
    invV_.assign(np * np, 0.0);
    Lift_.assign(np * nfp_tot, 0.0);
    IntWeightLgl_.assign(np, 0.0);
}

void ElementBase1D::InitBase1D(bool lumped_mat_flag)
{
    InitBase(lumped_mat_flag);

    const auto np = static_cast<std::size_t>(Np_);
    const auto nfp = static_cast<std::size_t>(Nfp_);
    const auto nfaces = static_cast<std::size_t>(Nfaces_);

    x1_.assign(np, 0.0);
    Fmask_.assign(nfp * nfaces, 0);
    Dx1_.assign(np * np, 0.0);
    Sx1_.assign(np * np, 0.0);
}

void ConstructMassMat(ArrayView<const Real, 2> V, ArrayView<Real, 2> mass_mat, ArrayView<Real, 2> inv_mass_mat)
{
    CheckSquare(V, "V", "ConstructMassMat");
    const std::size_t n = V.extent(0);
    if (mass_mat.extent(0) != n || mass_mat.extent(1) != n) {
        throw std::invalid_argument("ConstructMassMat: mass_mat's shape must match V's shape");
    }
    if (inv_mass_mat.extent(0) != n || inv_mass_mat.extent(1) != n) {
        throw std::invalid_argument("ConstructMassMat: inv_mass_mat's shape must match V's shape");
    }

    // invMassMat = V * V^T
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            Real s = 0.0;
            for (std::size_t k = 0; k < n; ++k) s += V(i, k) * V(j, k);
            inv_mass_mat(i, j) = s;
        }
    }

    common::linalgebra::Inv(ArrayView<const Real, 2>(inv_mass_mat.data(), n, n), mass_mat);
}

void ConstructStiffMat(ArrayView<const Real, 2> mass_mat, ArrayView<const Real, 2> inv_mass_mat,
                       ArrayView<const Real, 2> d_mat, ArrayView<Real, 2> stiff_mat)
{
    CheckSquare(mass_mat, "mass_mat", "ConstructStiffMat");
    const std::size_t n = mass_mat.extent(0);
    if (inv_mass_mat.extent(0) != n || inv_mass_mat.extent(1) != n) {
        throw std::invalid_argument("ConstructStiffMat: inv_mass_mat's shape must match mass_mat's shape");
    }
    if (d_mat.extent(0) != n || d_mat.extent(1) != n) {
        throw std::invalid_argument("ConstructStiffMat: d_mat's shape must match mass_mat's shape");
    }
    if (stiff_mat.extent(0) != n || stiff_mat.extent(1) != n) {
        throw std::invalid_argument("ConstructStiffMat: stiff_mat's shape must match mass_mat's shape");
    }

    // tmp1 = MassMat * DMat
    std::vector<Real> tmp1_buf(n * n);
    ArrayView<Real, 2> tmp1(tmp1_buf.data(), n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            Real s = 0.0;
            for (std::size_t k = 0; k < n; ++k) s += mass_mat(i, k) * d_mat(k, j);
            tmp1(i, j) = s;
        }
    }

    // StiffMat = InvMassMat * tmp1^T
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            Real s = 0.0;
            for (std::size_t k = 0; k < n; ++k) s += inv_mass_mat(i, k) * tmp1(j, k);
            stiff_mat(i, j) = s;
        }
    }
}

void ConstructLiftMat(ArrayView<const Real, 2> inv_mass_mat, ArrayView<const Real, 2> e_mat,
                      ArrayView<Real, 2> lift_mat)
{
    CheckSquare(inv_mass_mat, "inv_mass_mat", "ConstructLiftMat");
    const std::size_t n = inv_mass_mat.extent(0);
    const std::size_t nfp_tot = e_mat.extent(1);
    if (e_mat.extent(0) != n) {
        throw std::invalid_argument("ConstructLiftMat: e_mat's row count must match inv_mass_mat's size");
    }
    if (lift_mat.extent(0) != n || lift_mat.extent(1) != nfp_tot) {
        throw std::invalid_argument("ConstructLiftMat: lift_mat's shape must match (inv_mass_mat rows, e_mat cols)");
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < nfp_tot; ++j) {
            Real s = 0.0;
            for (std::size_t k = 0; k < n; ++k) s += inv_mass_mat(i, k) * e_mat(k, j);
            lift_mat(i, j) = s;
        }
    }
}

}  // namespace element
}  // namespace FElib
