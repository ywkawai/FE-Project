#include "scale_element_line.hpp"

#include <algorithm>
#include <cmath>

#include "scale_linalgebra.hpp"
#include "scale_polynomial.hpp"

namespace FElib {
namespace element {

using utility::ArrayView;

void LineElement::Init(int elem_order, bool lumped_mass_mat_flag)
{
    PolyOrder_ = elem_order;
    Nv_ = 2;
    Np_ = elem_order + 1;
    Nfp_ = 1;
    Nfaces_ = 2;
    NfpTot_ = Nfp_ * Nfaces_;

    InitBase1D(lumped_mass_mat_flag);
    ConstructElement();
}

void LineElement::ConstructElement()
{
    using common::polynomial::GenDLagrangePoly_lglpt;
    using common::polynomial::GenGaussLobattoPt;
    using common::polynomial::GenGaussLobattoPtIntWeight;
    using common::polynomial::GenLegendrePoly;

    const auto np = static_cast<std::size_t>(Np_);

    const std::vector<Real> lgl_pts = GenGaussLobattoPt(PolyOrder_);
    ArrayView<const Real, 1> lgl_pts_view(lgl_pts.data(), np);

    // P1D_ori(i, l) = value of the degree-l Legendre polynomial at lgl_pts[i].
    const std::vector<Real> P1D_ori = GenLegendrePoly(PolyOrder_, lgl_pts_view);
    ArrayView<const Real, 2> P1D_view(P1D_ori.data(), np, np);

    // DLagr1D(k, n) = derivative of the k-th Lagrange basis function at lgl_pts[n].
    // scale_element_line.F90 also computes polynomial_genDLegendrePoly(...) into a
    // local (DP1D_ori) that is never subsequently read -- that dead computation is
    // not reproduced here.
    const std::vector<Real> DLagr1D = GenDLagrangePoly_lglpt(PolyOrder_, lgl_pts_view);
    ArrayView<const Real, 2> DLagr1D_view(DLagr1D.data(), np, np);

    // Fmask(:, 0) = the single boundary node on the left face (node 0);
    // Fmask(:, 1) = the single boundary node on the right face (node Np-1).
    Fmask_[0] = 0;
    Fmask_[1] = static_cast<int>(np) - 1;

    ArrayView<Real, 2> V(V_.data(), np, np);
    ArrayView<Real, 2> Dx1(Dx1_.data(), np, np);
    for (std::size_t n = 0; n < np; ++n) {
        x1_[n] = lgl_pts[n];
        for (std::size_t l = 0; l < np; ++l) {
            V(n, l) = P1D_view(n, l) * std::sqrt(static_cast<Real>(l) + 0.5);
            Dx1(n, l) = DLagr1D_view(l, n);
        }
    }
    common::linalgebra::Inv(ArrayView<const Real, 2>(V_.data(), np, np), ArrayView<Real, 2>(invV_.data(), np, np));

    IntWeightLgl_ = GenGaussLobattoPtIntWeight(PolyOrder_);

    if (IsLumpedMatrix()) {
        std::fill(M_.begin(), M_.end(), 0.0);
        std::fill(invM_.begin(), invM_.end(), 0.0);
        ArrayView<Real, 2> M(M_.data(), np, np);
        ArrayView<Real, 2> invM(invM_.data(), np, np);
        for (std::size_t i = 0; i < np; ++i) {
            M(i, i) = IntWeightLgl_[i];
            invM(i, i) = 1.0 / IntWeightLgl_[i];
        }
    } else {
        ConstructMassMat(ArrayView<const Real, 2>(V_.data(), np, np), ArrayView<Real, 2>(M_.data(), np, np),
                         ArrayView<Real, 2>(invM_.data(), np, np));
    }

    ConstructStiffMat(ArrayView<const Real, 2>(M_.data(), np, np), ArrayView<const Real, 2>(invM_.data(), np, np),
                      ArrayView<const Real, 2>(Dx1_.data(), np, np), ArrayView<Real, 2>(Sx1_.data(), np, np));

    // Emat(Np, NfpTot): for each face f, the Nfp x Nfp identity block
    // (MassEdge) is scattered into rows Fmask(:,f), columns [f*Nfp, (f+1)*Nfp).
    // Since MassEdge is diagonal and Emat starts at zero, only the diagonal
    // entries need to be written (the Fortran source's full MassEdge loop --
    // which also assigns the always-zero off-diagonal entries -- collapses
    // to this without changing the result).
    const auto nfp = static_cast<std::size_t>(Nfp_);
    const auto nfaces = static_cast<std::size_t>(Nfaces_);
    const auto nfp_tot = static_cast<std::size_t>(NfpTot_);
    std::vector<Real> Emat_buf(np * nfp_tot, 0.0);
    ArrayView<Real, 2> Emat(Emat_buf.data(), np, nfp_tot);
    ArrayView<const int, 2> fmask(Fmask_.data(), nfp, nfaces);
    for (std::size_t f = 0; f < nfaces; ++f) {
        for (std::size_t k = 0; k < nfp; ++k) {
            Emat(static_cast<std::size_t>(fmask(k, f)), f * nfp + k) = 1.0;
        }
    }
    ConstructLiftMat(ArrayView<const Real, 2>(invM_.data(), np, np), ArrayView<const Real, 2>(Emat_buf.data(), np, nfp_tot),
                     ArrayView<Real, 2>(Lift_.data(), np, nfp_tot));
}

LineElement::IntrpResult LineElement::GenIntGaussLegendreIntrpMat(int interp_poly_order) const
{
    using common::polynomial::GenGaussLegendrePt;
    using common::polynomial::GenGaussLegendrePtIntWeight;
    using common::polynomial::GenLegendrePoly;

    IntrpResult result;
    result.points = GenGaussLegendrePt(interp_poly_order);
    result.weights = GenGaussLegendrePtIntWeight(interp_poly_order);

    const auto n_intrp = static_cast<std::size_t>(interp_poly_order);
    const auto np = static_cast<std::size_t>(Np_);

    ArrayView<const Real, 1> pts_view(result.points.data(), n_intrp);
    // P_int1D_ori(p1_, p1) = value of the degree-p1 Legendre polynomial at result.points[p1_].
    const std::vector<Real> P_int1D_ori = GenLegendrePoly(PolyOrder_, pts_view);
    ArrayView<const Real, 2> P_int1D_view(P_int1D_ori.data(), n_intrp, np);

    std::vector<Real> Vint_buf(n_intrp * np);
    ArrayView<Real, 2> Vint(Vint_buf.data(), n_intrp, np);
    for (std::size_t p1_ = 0; p1_ < n_intrp; ++p1_) {
        for (std::size_t p1 = 0; p1 < np; ++p1) {
            Vint(p1_, p1) = P_int1D_view(p1_, p1) * std::sqrt(static_cast<Real>(p1) + 0.5);
        }
    }

    // IntrpMat = Vint * invV
    result.matrix.assign(n_intrp * np, 0.0);
    ArrayView<Real, 2> intrp_mat(result.matrix.data(), n_intrp, np);
    ArrayView<const Real, 2> inv_v(invV_.data(), np, np);
    for (std::size_t i = 0; i < n_intrp; ++i) {
        for (std::size_t j = 0; j < np; ++j) {
            Real s = 0.0;
            for (std::size_t k = 0; k < np; ++k) s += Vint(i, k) * inv_v(k, j);
            intrp_mat(i, j) = s;
        }
    }

    return result;
}

}  // namespace element
}  // namespace FElib
