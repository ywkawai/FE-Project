/*
 * @par Description
 * A collection of polynomial utilities used to build finite-element bases: Legendre polynomials, Gauss-Lobatto-Legendre / Gauss-Legendre quadrature points and weights, and Lagrange basis functions on GLL points.
 * 
 * @author Yuta Kawai, Team SCALE
 */
#pragma once

#include <cstddef>
#include <vector>

#include "utility/array_view.hpp"

namespace FElib {
namespace common {
namespace polynomial {

using utility::ArrayView;

// ---------------------------------------------------------------------------
// Legendre polynomials
// ---------------------------------------------------------------------------

//! P(i, n) = value of the degree-n Legendre polynomial at x(i), for n = 0..order.
//! @param P must have shape (x.extent(0), order+1).
//! @throws std::invalid_argument if order < 0 or P's shape does not match.
void GenLegendrePoly(int order, ArrayView<const double, 1> x, ArrayView<double, 2> P);

//! Convenience overload: allocates and returns a buffer with shape (x.extent(0), order+1).
std::vector<double> GenLegendrePoly(int order, ArrayView<const double, 1> x);

//! gradP(i, n) = d/dx of the degree-n Legendre polynomial at x(i).
//! @param P must already hold GenLegendrePoly(order, x) -- this mirrors the
//!        Fortran source, which reuses a precomputed P instead of recomputing it.
//! @param gradP must have the same shape as P.
//! @throws std::invalid_argument if order < 0 or P's/gradP's shape does not match.
void GenDLegendrePoly(int order, ArrayView<const double, 1> x, ArrayView<const double, 2> P,
                      ArrayView<double, 2> gradP);

std::vector<double> GenDLegendrePoly(int order, ArrayView<const double, 1> x,
                                      ArrayView<const double, 2> P);

// ---------------------------------------------------------------------------
// Gauss-Lobatto-Legendre (GLL) points and integration weights
// ---------------------------------------------------------------------------

//! order+1 Gauss-Lobatto-Legendre points on [-1,1], including both endpoints.
//! @param pts must have extent(0) == order+1.
//! @throws std::invalid_argument if order < 1 (order == 0 has no well-defined
//!         GLL rule; the Fortran version silently returns a degenerate point set).
void GenGaussLobattoPt(int order, ArrayView<double, 1> pts);
std::vector<double> GenGaussLobattoPt(int order);

//! Integration weights for the order+1 GLL points above.
//! @throws std::invalid_argument if order < 1 (the formula divides by order*(order+1)).
void GenGaussLobattoPtIntWeight(int order, ArrayView<double, 1> weights);
std::vector<double> GenGaussLobattoPtIntWeight(int order);

// ---------------------------------------------------------------------------
// Gauss-Legendre (GL) points and integration weights
// ---------------------------------------------------------------------------

//! order Gauss-Legendre points on [-1,1] (interior points only, no endpoints).
//! @param pts must have extent(0) == order.
//! @throws std::invalid_argument if order < 1.
void GenGaussLegendrePt(int order, ArrayView<double, 1> pts);
std::vector<double> GenGaussLegendrePt(int order);

//! Integration weights for the order GL points above.
//! @throws std::invalid_argument if order < 1.
void GenGaussLegendrePtIntWeight(int order, ArrayView<double, 1> weights);
std::vector<double> GenGaussLegendrePtIntWeight(int order);

// ---------------------------------------------------------------------------
// Lagrange basis functions on GLL points
// ---------------------------------------------------------------------------

//! L(i, n) = value of the n-th Lagrange basis function (the one equal to 1 at xLgl(n) and 0 at the other GLL points) evaluated at x(i).
//! @param xLgl must have extent(0) == order+1 (typically GenGaussLobattoPt(order)).
//! @param L must have shape (x.extent(0), order+1).
//! @throws std::invalid_argument if order < 1 (the formula divides by order*(order+1))
//!         or a shape mismatch is detected.
void GenLagrangePoly(int order, ArrayView<const double, 1> xLgl, ArrayView<const double, 1> x,
                     ArrayView<double, 2> L);
std::vector<double> GenLagrangePoly(int order, ArrayView<const double, 1> xLgl,
                                     ArrayView<const double, 1> x);

//! Lr(k, n) = derivative of the k-th Lagrange basis function evaluated at xLgl(n).
//! I.e. the row index is the basis-function index and the column index is the evaluation-point index: for nodal values f(xLgl(k)), the derivative of the interpolant at point xLgl(n) is sum_k Lr(k, n) * f(xLgl(k)).
//! @param Lr must have shape (order+1, order+1).
//! @throws std::invalid_argument if order < 0 or a shape mismatch is detected.
void GenDLagrangePoly_lglpt(int order, ArrayView<const double, 1> xLgl, ArrayView<double, 2> Lr);
std::vector<double> GenDLagrangePoly_lglpt(int order, ArrayView<const double, 1> xLgl);

}  // namespace polynomial
}  // namespace common
}  // namespace FElib
