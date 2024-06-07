!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Common / Modal filter
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
    lmesh, elem, filter, do_weight_Gsqrt     ) ! (in)

    implicit none

    class(LocalMeshBase), intent(in) :: lmesh
    class(ElementBase), intent(in) :: elem    
    real(RP), intent(inout)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DRHOT_(elem%Np,lmesh%NeA)
    class(ModalFilter), intent(in) :: filter
    logical, intent(in), optional :: do_weight_Gsqrt
    
    integer :: ke
    real(RP) :: tmp(elem%Np,5)
    integer :: ii, kk
    real(RP) :: Mik
    logical :: do_weight_Gsqrt_
    real(RP) :: RGsqrt(elem%Np)
    !------------------------------------

    if ( present( do_weight_Gsqrt ) ) then
      do_weight_Gsqrt_ = do_weight_Gsqrt
    else
      do_weight_Gsqrt_ = .false.
    end if

    if ( do_weight_Gsqrt_ ) then
      !$omp parallel do private( tmp, ii, kk, Mik, RGsqrt )
      do ke=lmesh%NeS, lmesh%NeE

        tmp(:,:) = 0.0_RP
        do ii=1, elem%Np
        do kk=1, elem%Np
          Mik = filter%FilterMat(ii,kk) * lmesh%Gsqrt(kk,ke)

          tmp(ii,1) = tmp(ii,1) + Mik * DDENS_(kk,ke)
          tmp(ii,2) = tmp(ii,2) + Mik * MOMX_ (kk,ke)
          tmp(ii,3) = tmp(ii,3) + Mik * MOMY_ (kk,ke)
          tmp(ii,4) = tmp(ii,4) + Mik * MOMZ_ (kk,ke)
          tmp(ii,5) = tmp(ii,5) + Mik * DRHOT_(kk,ke)
        end do
        end do

        RGsqrt(:) = 1.0_RP / lmesh%Gsqrt(:,ke)
        DDENS_(:,ke) = tmp(:,1) * RGsqrt(:)
        MOMX_ (:,ke) = tmp(:,2) * RGsqrt(:)
        MOMY_ (:,ke) = tmp(:,3) * RGsqrt(:)
        MOMZ_ (:,ke) = tmp(:,4) * RGsqrt(:)
        DRHOT_(:,ke) = tmp(:,5) * RGsqrt(:)
      end do    

    else

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

    end if

    return
  end subroutine atm_dyn_dgm_modalfilter_apply

!OCL SERIAL
  subroutine atm_dyn_dgm_tracer_modalfilter_apply(  &
    QTRC_,                                   & ! (inout)
    DENS_hyd_, DDENS_,                       & ! (inout)
    lmesh, elem, filter                      ) ! (in)

    implicit none

    class(LocalMeshBase), intent(in) :: lmesh
    class(ElementBase), intent(in) :: elem    
    real(RP), intent(inout)  :: QTRC_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DENS_hyd_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    class(ModalFilter), intent(in) :: filter

    integer :: ke
    real(RP) :: tmp(elem%Np)
    integer :: ii, kk
    real(RP) :: Mik

    real(RP) :: weight(elem%Np)
    !------------------------------------

    !$omp parallel do private( tmp, ii, kk, Mik, weight )
    do ke=lmesh%NeS, lmesh%NeE

      tmp(:) = 0.0_RP
      weight(:) = lmesh%Gsqrt(:,ke) &
                * ( DENS_hyd_(:,ke) + DDENS_(:,ke) )

      do ii=1, elem%Np
      do kk=1, elem%Np
        Mik = filter%FilterMat(ii,kk) * weight(kk)

        tmp(ii) = tmp(ii) + Mik * QTRC_(kk,ke)
      end do
      end do

      QTRC_(:,ke) = tmp(:) / weight(:)
    end do    

    return
  end subroutine atm_dyn_dgm_tracer_modalfilter_apply

end module scale_atm_dyn_dgm_modalfilter
