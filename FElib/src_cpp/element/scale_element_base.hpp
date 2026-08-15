/*
 * @par Description
 * Base classes for finite elements. FElib/src/element/scale_element_base.F90
 * uses Fortran type extension (ElementBase <- ElementBase1D/2D/3D); this is
 * the first FElib C++ module that needs a real class hierarchy rather than a
 * single class or a free-function namespace.
 *
 * Only ElementBase and ElementBase1D are ported so far -- the part needed by
 * FElib/src_cpp/element/scale_element_line.hpp's LineElement.
 * ElementBase2D/ElementBase3D (and their L2-projection / interpolation /
 * modal-truncation machinery) are not ported yet: with no concrete 2D/3D
 * element (Quadrilateral/Hexahedral) implemented, there is nothing to
 * exercise or test them against. See PORTING_STATUS.md.
 *
 * Design notes:
 *  - ElementBase/ElementBase1D model an inheritance-friendly base: their
 *    constructors are protected (so only a concrete leaf element such as
 *    LineElement can be instantiated), and the destructor is virtual (so a
 *    concrete element can safely be owned/destroyed through a
 *    `ElementBase&`/`ElementBase*`, once more than one concrete element type
 *    exists). Neither base class declares any *other* virtual function:
 *    nothing here is overridden polymorphically by a subclass (matching the
 *    Fortran source, where e.g. IsLumpedMatrix is inherited unchanged, never
 *    redefined), so ordinary (non-virtual) member functions are used
 *    everywhere else to keep accessors trivially inlinable.
 *  - The owned matrices/arrays are `protected` data (not private-plus-a
 *    separate mutable-accessor pair): only this class family's own
 *    construction code (e.g. LineElement's element-construction step) ever
 *    writes to them, exactly as in the Fortran source, so protected direct
 *    access costs nothing extra at runtime while keeping the data closed to
 *    unrelated external code. External callers only ever see the `public`
 *    read-only ArrayView accessors below.
 *
 * @author Yuta Kawai, Team SCALE
 */
#pragma once

#include <vector>

#include "utility/array_view.hpp"

namespace FElib {
namespace element {

using utility::ArrayView;
using Real = double;

//! Common base for every finite element: the Vandermonde matrix and the
//! matrices/weights derived from it (mass, lift, LGL integration weights).
class ElementBase {
public:
    virtual ~ElementBase() = default;

    int Np() const noexcept { return Np_; }             //< Number of nodes within the element.
    int Nfaces() const noexcept { return Nfaces_; }      //< Number of faces.
    int NfpTot() const noexcept { return NfpTot_; }       //< Total number of nodes on faces.
    int Nv() const noexcept { return Nv_; }              //< Number of vertices.
    bool IsLumpedMatrix() const noexcept { return lumped_mat_flag_; }

    ArrayView<const Real, 2> V() const noexcept { return {V_.data(), static_cast<std::size_t>(Np_), static_cast<std::size_t>(Np_)}; }
    ArrayView<const Real, 2> InvV() const noexcept { return {invV_.data(), static_cast<std::size_t>(Np_), static_cast<std::size_t>(Np_)}; }
    ArrayView<const Real, 2> M() const noexcept { return {M_.data(), static_cast<std::size_t>(Np_), static_cast<std::size_t>(Np_)}; }
    ArrayView<const Real, 2> InvM() const noexcept { return {invM_.data(), static_cast<std::size_t>(Np_), static_cast<std::size_t>(Np_)}; }
    ArrayView<const Real, 2> Lift() const noexcept { return {Lift_.data(), static_cast<std::size_t>(Np_), static_cast<std::size_t>(NfpTot_)}; }
    ArrayView<const Real, 1> IntWeightLgl() const noexcept { return {IntWeightLgl_.data(), static_cast<std::size_t>(Np_)}; }

protected:
    ElementBase() = default;

    //! Allocates M_/invM_/V_/invV_/Lift_/IntWeightLgl_. Np_/NfaceS_/NfpTot_/Nv_
    //! must already be set by the concrete element's own Init() before this
    //! is called (mirrors the Fortran source's ElementBase_Init, which reads
    //! elem%Np/elem%NfpTot to size the allocations).
    void InitBase(bool lumped_mat_flag);

    int Np_ = 0;
    int Nfaces_ = 0;
    int NfpTot_ = 0;
    int Nv_ = 0;
    bool lumped_mat_flag_ = false;

    std::vector<Real> V_;
    std::vector<Real> invV_;
    std::vector<Real> M_;
    std::vector<Real> invM_;
    std::vector<Real> Lift_;
    std::vector<Real> IntWeightLgl_;
};

//! Base for a 1D (line) reference element.
class ElementBase1D : public ElementBase {
public:
    int PolyOrder() const noexcept { return PolyOrder_; }
    int Nfp() const noexcept { return Nfp_; }

    //! Fmask(k, f) = the 0-based node index of the k-th node on face f (k in [0,Nfp), f in [0,Nfaces)).
    ArrayView<const int, 2> Fmask() const noexcept { return {Fmask_.data(), static_cast<std::size_t>(Nfp_), static_cast<std::size_t>(Nfaces_)}; }
    ArrayView<const Real, 1> x1() const noexcept { return {x1_.data(), static_cast<std::size_t>(Np_)}; }
    //! Elementwise differential matrix for the x1 direction (Dx1 = M^-1 Sx1).
    ArrayView<const Real, 2> Dx1() const noexcept { return {Dx1_.data(), static_cast<std::size_t>(Np_), static_cast<std::size_t>(Np_)}; }
    //! Elementwise stiffness matrix for the x1 direction.
    ArrayView<const Real, 2> Sx1() const noexcept { return {Sx1_.data(), static_cast<std::size_t>(Np_), static_cast<std::size_t>(Np_)}; }

protected:
    ElementBase1D() = default;

    //! Calls InitBase(), then allocates x1_/Fmask_/Dx1_/Sx1_. PolyOrder_/Nfp_
    //! (in addition to Np_/Nfaces_/NfpTot_/Nv_) must already be set.
    void InitBase1D(bool lumped_mat_flag);

    int PolyOrder_ = 0;
    int Nfp_ = 0;
    std::vector<int> Fmask_;

    std::vector<Real> x1_;
    std::vector<Real> Dx1_;
    std::vector<Real> Sx1_;
};

//! MassMat = inv(V V^T); InvMassMat = V V^T. Unlike the Fortran source
//! (where invMassMat is an optional output), this always computes and writes
//! both -- invMassMat is already computed as an internal step regardless of
//! whether the Fortran caller asks for it, so there is nothing to save by
//! making it optional (same simplification already applied to
//! LineElement::GenIntGaussLegendreIntrpMat's formerly-optional outputs).
//! @throws std::invalid_argument if V is not square, or mass_mat's/inv_mass_mat's shape doesn't match V's.
//! @throws std::runtime_error if V V^T is singular.
void ConstructMassMat(ArrayView<const Real, 2> V, ArrayView<Real, 2> mass_mat, ArrayView<Real, 2> inv_mass_mat);

//! StiffMat = InvMassMat * (MassMat * DMat)^T.
//! @throws std::invalid_argument on a shape mismatch between the arguments.
void ConstructStiffMat(ArrayView<const Real, 2> mass_mat, ArrayView<const Real, 2> inv_mass_mat,
                       ArrayView<const Real, 2> d_mat, ArrayView<Real, 2> stiff_mat);

//! LiftMat = InvMassMat * EMat.
//! @throws std::invalid_argument on a shape mismatch between the arguments.
void ConstructLiftMat(ArrayView<const Real, 2> inv_mass_mat, ArrayView<const Real, 2> e_mat,
                      ArrayView<Real, 2> lift_mat);

}  // namespace element
}  // namespace FElib
