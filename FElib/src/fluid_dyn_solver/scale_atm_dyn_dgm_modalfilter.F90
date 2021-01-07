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
    real(RP) :: tmp(elem%Np,5)
    integer :: ii, kk
    real(RP) :: Mik
    !------------------------------------

    !$omp parallel do private( tmp, ii, kk, Mik )
    do ke=lmesh%NeS, lmesh%NeE

      tmp(:,:) = 0.0_RP
      do ii=1, elem%Np
      do kk=1, elem%Np
        Mik = filter%FilterMat(ii,kk)
        tmp(ii,1) = tmp(ii,1) + Mik * DDENS_(kk,ke)
        tmp(ii,2) = tmp(ii,2) + Mik * MOMX_ (kk,ke)
        tmp(ii,3) = tmp(ii,3) + Mik * MOMY_ (kk,ke)
        tmp(ii,4) = tmp(ii,4) + Mik * MOMZ_ (kk,ke)
        tmp(ii,5) = tmp(ii,5) + Mik * DRHOT_(kk,ke)
      end do
      end do      
      DDENS_(:,ke) = tmp(:,1)
      MOMX_ (:,ke) = tmp(:,2)
      MOMY_ (:,ke) = tmp(:,3)
      MOMZ_ (:,ke) = tmp(:,4)
      DRHOT_(:,ke) = tmp(:,5)
    end do

    return
  end subroutine atm_dyn_dgm_modalfilter_apply

end module scale_atm_dyn_dgm_modalfilter
