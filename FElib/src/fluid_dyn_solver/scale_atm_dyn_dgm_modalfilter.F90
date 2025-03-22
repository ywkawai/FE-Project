!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Common / Modal filter
!!
!! @par Description
!!      Modal filter for Atmospheric dynamical process. 
!!      The modal filter suppresses the numerical instability due to the aliasing errors. 
!!
!! @author Yuta Kawai, Team SCALE
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
  use scale_element_operation_base, only: ElementOperationBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_modalfilter_apply
  public :: atm_dyn_dgm_tracer_modalfilter_apply

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
    lmesh, elem, elem_operation,             & ! (in) 
    do_weight_Gsqrt )                          ! (in)
    implicit none

    class(LocalMeshBase), intent(in) :: lmesh
    class(ElementBase), intent(in) :: elem    
    real(RP), intent(inout)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DRHOT_(elem%Np,lmesh%NeA)
    class(ElementOperationBase3D), intent(in) :: elem_operation
    logical, intent(in), optional :: do_weight_Gsqrt
    
    integer :: ke
    real(RP) :: tmp(elem%Np,5)
    real(RP) :: work(elem%Np)
    real(RP) :: tmp_out(elem%Np,5)
    integer :: kk
    logical :: do_weight_Gsqrt_
    real(RP) :: RGsqrt(elem%Np)
    !------------------------------------

    if ( present( do_weight_Gsqrt ) ) then
      do_weight_Gsqrt_ = do_weight_Gsqrt
    else
      do_weight_Gsqrt_ = .false.
    end if

    if ( do_weight_Gsqrt_ ) then
      !$omp parallel do private( tmp, tmp_out, work, kk, RGsqrt )
      do ke=lmesh%NeS, lmesh%NeE

        do kk=1, elem%Np
          tmp(kk,1) = lmesh%Gsqrt(kk,ke) * DDENS_(kk,ke)
          tmp(kk,2) = lmesh%Gsqrt(kk,ke) * MOMX_(kk,ke)
          tmp(kk,3) = lmesh%Gsqrt(kk,ke) * MOMY_(kk,ke)
          tmp(kk,4) = lmesh%Gsqrt(kk,ke) * MOMZ_(kk,ke)
          tmp(kk,5) = lmesh%Gsqrt(kk,ke) * DRHOT_(kk,ke)
        end do

        call elem_operation%ModalFilter_var5( tmp, work, &
          tmp_out )

        RGsqrt(:) = 1.0_RP / lmesh%Gsqrt(:,ke)
        DDENS_(:,ke) = tmp_out(:,1) * RGsqrt(:)
        MOMX_ (:,ke) = tmp_out(:,2) * RGsqrt(:)
        MOMY_ (:,ke) = tmp_out(:,3) * RGsqrt(:)
        MOMZ_ (:,ke) = tmp_out(:,4) * RGsqrt(:)
        DRHOT_(:,ke) = tmp_out(:,5) * RGsqrt(:)
      end do    

    else

      !$omp parallel do private( tmp, work, kk )
      do ke=lmesh%NeS, lmesh%NeE

        do kk=1, elem%Np
          tmp(kk,1) = DDENS_(kk,ke)
          tmp(kk,2) = MOMX_(kk,ke)
          tmp(kk,3) = MOMY_(kk,ke)
          tmp(kk,4) = MOMZ_(kk,ke)
          tmp(kk,5) = DRHOT_(kk,ke)
        end do
        
        call elem_operation%ModalFilter_var5( tmp, work, &
          tmp_out )

        DDENS_(:,ke) = tmp(:,1)
        MOMX_ (:,ke) = tmp(:,2)
        MOMY_ (:,ke) = tmp(:,3)
        MOMZ_ (:,ke) = tmp(:,4)
        DRHOT_(:,ke) = tmp(:,5)
      end do

    end if

    return
  end subroutine atm_dyn_dgm_modalfilter_apply

!OCL SERIAL
  subroutine atm_dyn_dgm_tracer_modalfilter_apply(  &
    QTRC_,                                   & ! (inout)
    DENS_hyd_, DDENS_,                       & ! (inout)
    lmesh, elem, elem_operation              ) ! (in)

    implicit none

    class(LocalMeshBase), intent(in) :: lmesh
    class(ElementBase), intent(in) :: elem    
    real(RP), intent(inout)  :: QTRC_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DENS_hyd_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    class(ElementOperationBase3D), intent(in) :: elem_operation

    integer :: ke
    real(RP) :: tmp(elem%Np)
    real(RP) :: work(elem%Np)
    real(RP) :: tmp_out(elem%Np)
    real(RP) :: weight(elem%Np)
    integer :: kk
    !------------------------------------

    !$omp parallel do private( tmp, tmp_out, work, kk, weight )
    do ke=lmesh%NeS, lmesh%NeE

      do kk=1, elem%Np
        weight(kk) = lmesh%Gsqrt(kk,ke) &
                   * ( DENS_hyd_(kk,ke) + DDENS_(kk,ke) ) 
        tmp(kk) = weight(kk) * QTRC_(kk,ke)
      end do

      call elem_operation%ModalFilter_tracer( tmp, work, &
        tmp_out )

      QTRC_(:,ke) = tmp_out(:) / weight(:)
    end do    

    return
  end subroutine atm_dyn_dgm_tracer_modalfilter_apply

end module scale_atm_dyn_dgm_modalfilter
