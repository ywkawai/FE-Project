#include "scaleFElib.h"
#include "scale_c_binding.h"
module scale_polynomial_cbind
  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_cbinding_util, only: &
    logical_c2f => cbinding_util_logical_c2f, &
    arrayPtr_c2f => cbinding_util_ArrayPtr_c2f
  use iso_c_binding
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
contains
  subroutine Polynomial_GenLagrangePoly( Nord, x_lgl, x, Nx, &
    l ) bind(C, name="CPolynomial_GenLagrangePoly")
    use scale_polynominal, only: Polynominal_GenLagrangePoly
    implicit none
    integer(c_int), value :: Nord
    type(c_ptr), value :: x_lgl
    type(c_ptr), value :: x
    integer(c_int), value :: Nx
    type(c_ptr), value :: l

    real(c_double), pointer :: l_(:,:)
    !----------------------------
    call c_f_pointer(l, l_, [Nx,Nord+1])
    l_(:,:) = Polynominal_GenLagrangePoly( &
        Nord, arrayptr_c2f(x_lgl, Nord+1), arrayptr_c2f(x, Nx) )
    return
  end subroutine Polynomial_GenLagrangePoly

  subroutine Polynomial_GenDLagrangePoly_lglpt( Nord, x_lgl, x, Nx, &
    l ) bind(C, name="CPolynomial_GenDLagrangePoly_lglpt")
    use scale_polynominal, only: Polynominal_GenLagrangePoly
    implicit none
    integer(c_int), value :: Nord
    type(c_ptr), value :: x_lgl
    type(c_ptr), value :: x
    integer(c_int), value :: Nx
    type(c_ptr), value :: l

    real(c_double), pointer :: l_(:,:)
    !----------------------------
    call c_f_pointer(l, l_, [Nx,Nord+1])
    l_(:,:) = Polynominal_GenLagrangePoly( &
        Nord, arrayptr_c2f(x_lgl, Nord+1), arrayptr_c2f(x, Nx) )
    return
  end subroutine Polynomial_GenDLagrangePoly_lglpt

  subroutine Polynomial_GenLegendrePoly( Nord, x, Nx, &
    P ) bind(C, name="CPolynomial_GenLegendrePoly")
    use scale_polynominal, only: Polynominal_GenLegendrePoly
    implicit none
    integer(c_int), value :: Nord
    type(c_ptr), value :: x
    integer(c_int), value :: Nx
    type(c_ptr), value :: P

    real(c_double), pointer :: P_(:,:)
    !----------------------------
    call c_f_pointer(P, P_, [Nx,Nord+1])
    P_(:,:) = Polynominal_GenLegendrePoly( &
        Nord,  arrayptr_c2f(x, Nx) )
    return
  end subroutine Polynomial_GenLegendrePoly

!******
!***** Getter
!******

end module scale_polynomial_cbind
