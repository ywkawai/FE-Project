/*
 * @par Description
 * A 1D (line) finite element: LGL nodes, Vandermonde/mass/stiffness/lift
 * matrices, built from FElib::common::polynomial. Ported from
 * FElib/src/element/scale_element_line.F90.
 *
 * @author Yuta Kawai, Team SCALE
 */
#pragma once

#include <vector>

#include "scale_element_base.hpp"

namespace FElib {
namespace element {

class LineElement : public ElementBase1D {
public:
    //! @param elem_order polynomial order of the element (Np = elem_order + 1 LGL nodes).
    //! @param lumped_mass_mat_flag if true, use the diagonal (LGL-quadrature) mass matrix
    //!        instead of the exact one.
    void Init(int elem_order, bool lumped_mass_mat_flag);

    //! Result of GenIntGaussLegendreIntrpMat: a matrix interpolating nodal
    //! values (at this element's LGL nodes) onto `interp_poly_order`
    //! Gauss-Legendre points, plus those points and their integration
    //! weights. The Fortran source (LineElement_gen_IntGaussLegendreIntrpMat)
    //! makes `points`/`weights` optional outputs, but both are computed
    //! internally regardless of whether the caller wants them, so this
    //! always returns all three (same simplification as
    //! ElementBase's ConstructMassMat always returning inv_mass_mat).
    struct IntrpResult {
        std::vector<Real> matrix;   //< shape (interp_poly_order, Np)
        std::vector<Real> points;   //< length interp_poly_order
        std::vector<Real> weights;  //< length interp_poly_order
    };

    //! @throws std::invalid_argument if interp_poly_order < 1.
    IntrpResult GenIntGaussLegendreIntrpMat(int interp_poly_order) const;

private:
    void ConstructElement();
};

}  // namespace element
}  // namespace FElib
