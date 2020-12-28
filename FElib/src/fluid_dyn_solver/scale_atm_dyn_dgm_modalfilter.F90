!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!      Modal filter for Atmospheric dynamical process. 
!!      The modal filter surpresses the numerical instability due to the aliasing errors. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_modalfilter
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: ElementBase
  use scale_localmesh_base, only: LocalMeshBase
  use scale_element_modalfilter, only: ModalFilter

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_modalfilter_apply

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

contains

!OCL SERIAL
  subroutine atm_dyn_dgm_modalfilter_apply(  &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,     & ! (inout)
    lmesh, elem, filter )                      ! (in)

    implicit none

    class(LocalMeshBase), intent(in) :: lmesh
    class(ElementBase), intent(in) :: elem    
    real(RP), intent(inout)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DRHOT_(elem%Np,lmesh%NeA)
    class(ModalFilter), intent(in) :: filter
    
    integer :: ke
    !------------------------------------

    !$omp parallel do
    do ke=1, lmesh%Ne
      DDENS_(:,ke) = matmul(filter%FilterMat,DDENS_(:,ke))
      MOMX_ (:,ke) = matmul(filter%FilterMat,MOMX_ (:,ke))
      MOMY_ (:,ke) = matmul(filter%FilterMat,MOMY_ (:,ke))
      MOMZ_ (:,ke) = matmul(filter%FilterMat,MOMZ_ (:,ke))
      DRHOT_(:,ke) = matmul(filter%FilterMat,DRHOT_(:,ke))
    end do

    return
  end subroutine atm_dyn_dgm_modalfilter_apply

end module scale_atm_dyn_dgm_modalfilter
