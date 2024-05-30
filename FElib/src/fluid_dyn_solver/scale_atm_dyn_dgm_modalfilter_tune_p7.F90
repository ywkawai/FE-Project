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
  private :: apply_filter_xyz_direction
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
    real(RP) :: FilterMat_tr(8,8)
    integer :: ii, kk
    real(RP) :: Mik
    logical :: do_weight_Gsqrt_
    real(RP) :: RGsqrt(elem%Np)
    !------------------------------------
    FilterMat_tr(:,:) = transpose(filter%FilterMat(:,:))


    if ( present( do_weight_Gsqrt ) ) then
      do_weight_Gsqrt_ = do_weight_Gsqrt
    else
      do_weight_Gsqrt_ = .false.
    end if

    if ( do_weight_Gsqrt_ ) then
      !$omp parallel do private( tmp, ii, kk, Mik, RGsqrt )
      do ke=lmesh%NeS, lmesh%NeE

        tmp(:,:) = 0.0_RP
        !do ii=1, elem%Np
        !do kk=1, elem%Np
        !  Mik = filter%FilterMat(ii,kk) * lmesh%Gsqrt(kk,ke)

        !  tmp(ii,1) = tmp(ii,1) + Mik * DDENS_(kk,ke)
        !  tmp(ii,2) = tmp(ii,2) + Mik * MOMX_ (kk,ke)
        !  tmp(ii,3) = tmp(ii,3) + Mik * MOMY_ (kk,ke)
        !  tmp(ii,4) = tmp(ii,4) + Mik * MOMZ_ (kk,ke)
        !  tmp(ii,5) = tmp(ii,5) + Mik * DRHOT_(kk,ke)
        !end do
        !end do

        do ii=1, elem%Np
          DDENS_(ii,ke) = lmesh%Gsqrt(ii,ke) * DDENS_(ii,ke)
          MOMX_(ii,ke) = lmesh%Gsqrt(ii,ke) * MOMX_(ii,ke)
          MOMY_(ii,ke) = lmesh%Gsqrt(ii,ke) * MOMY_(ii,ke)
          MOMZ_(ii,ke) = lmesh%Gsqrt(ii,ke) * MOMZ_(ii,ke)
          DRHOT_(ii,ke) = lmesh%Gsqrt(ii,ke) * DRHOT_(ii,ke)
        end do

        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, DDENS_(:,ke), MOMX_(:,ke), MOMY_(:,ke), MOMZ_(:,ke), DRHOT_(:,ke), tmp)

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
        !do ii=1, elem%Np
        !do kk=1, elem%Np
        !  Mik = filter%FilterMat(ii,kk)
        !  tmp(ii,1) = tmp(ii,1) + Mik * DDENS_(kk,ke)
        !  tmp(ii,2) = tmp(ii,2) + Mik * MOMX_ (kk,ke)
        !  tmp(ii,3) = tmp(ii,3) + Mik * MOMY_ (kk,ke)
        !  tmp(ii,4) = tmp(ii,4) + Mik * MOMZ_ (kk,ke)
        !  tmp(ii,5) = tmp(ii,5) + Mik * DRHOT_(kk,ke)
        !end do
        !end do      

        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, DDENS_(:,ke), MOMX_(:,ke), MOMY_(:,ke), MOMZ_(:,ke), DRHOT_(:,ke), tmp)

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

!OCL SERIAL
  subroutine apply_filter_xyz_direction(filterMat, filterMat_tr, DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, tmp )
    implicit none

    real(RP), intent(in) :: filterMat(8, 8)
    real(RP), intent(in) :: filterMat_tr(8, 8)
    real(RP), intent(inout) :: DDENS_(8, 8, 8)
    real(RP), intent(inout) :: MOMX_(8, 8, 8)
    real(RP), intent(inout) :: MOMY_(8, 8, 8)
    real(RP), intent(inout) :: MOMZ_(8, 8, 8)
    real(RP), intent(inout) :: DRHOT_(8, 8, 8)
    real(RP), intent(inout) :: tmp(8, 8, 8, 5)
    
    integer :: i, j, k

    !-- x direction
    do k=1, 8
    do j=1, 8
    do i=1, 8

      tmp(i,j,k,1) = filterMat(i,1) * DDENS_(1,j,k) + &
                     filterMat(i,2) * DDENS_(2,j,k) + &
                     filterMat(i,3) * DDENS_(3,j,k) + & 
                     filterMat(i,4) * DDENS_(4,j,k) + & 
                     filterMat(i,5) * DDENS_(5,j,k) + & 
                     filterMat(i,6) * DDENS_(6,j,k) + & 
                     filterMat(i,7) * DDENS_(7,j,k) + & 
                     filterMat(i,8) * DDENS_(8,j,k)  

      tmp(i,j,k,2) = filterMat(i,1) * MOMX_(1,j,k) + &
                     filterMat(i,2) * MOMX_(2,j,k) + &
                     filterMat(i,3) * MOMX_(3,j,k) + & 
                     filterMat(i,4) * MOMX_(4,j,k) + & 
                     filterMat(i,5) * MOMX_(5,j,k) + & 
                     filterMat(i,6) * MOMX_(6,j,k) + & 
                     filterMat(i,7) * MOMX_(7,j,k) + & 
                     filterMat(i,8) * MOMX_(8,j,k)  

      tmp(i,j,k,3) = filterMat(i,1) * MOMY_(1,j,k) + &
                     filterMat(i,2) * MOMY_(2,j,k) + &
                     filterMat(i,3) * MOMY_(3,j,k) + & 
                     filterMat(i,4) * MOMY_(4,j,k) + & 
                     filterMat(i,5) * MOMY_(5,j,k) + & 
                     filterMat(i,6) * MOMY_(6,j,k) + & 
                     filterMat(i,7) * MOMY_(7,j,k) + & 
                     filterMat(i,8) * MOMY_(8,j,k)  

      tmp(i,j,k,4) = filterMat(i,1) * MOMZ_(1,j,k) + &
                     filterMat(i,2) * MOMZ_(2,j,k) + &
                     filterMat(i,3) * MOMZ_(3,j,k) + & 
                     filterMat(i,4) * MOMZ_(4,j,k) + & 
                     filterMat(i,5) * MOMZ_(5,j,k) + & 
                     filterMat(i,6) * MOMZ_(6,j,k) + & 
                     filterMat(i,7) * MOMZ_(7,j,k) + & 
                     filterMat(i,8) * MOMZ_(8,j,k)  

      tmp(i,j,k,5) = filterMat(i,1) * DRHOT_(1,j,k) + &
                     filterMat(i,2) * DRHOT_(2,j,k) + &
                     filterMat(i,3) * DRHOT_(3,j,k) + & 
                     filterMat(i,4) * DRHOT_(4,j,k) + & 
                     filterMat(i,5) * DRHOT_(5,j,k) + & 
                     filterMat(i,6) * DRHOT_(6,j,k) + & 
                     filterMat(i,7) * DRHOT_(7,j,k) + & 
                     filterMat(i,8) * DRHOT_(8,j,k)  

    end do
    end do
    end do

    !-- y direction
    do k=1, 8
    do j=1, 8
    do i=1, 8
      DDENS_(i,j,k) = tmp(i,1,k,1) * filterMat_tr(1,j) + &
                      tmp(i,2,k,1) * filterMat_tr(2,j) + &
                      tmp(i,3,k,1) * filterMat_tr(3,j) + &
                      tmp(i,4,k,1) * filterMat_tr(4,j) + &
                      tmp(i,5,k,1) * filterMat_tr(5,j) + &
                      tmp(i,6,k,1) * filterMat_tr(6,j) + &
                      tmp(i,7,k,1) * filterMat_tr(7,j) + &
                      tmp(i,8,k,1) * filterMat_tr(8,j) 

      MOMX_(i,j,k) = tmp(i,1,k,2) * filterMat_tr(1,j) + &
                     tmp(i,2,k,2) * filterMat_tr(2,j) + &
                     tmp(i,3,k,2) * filterMat_tr(3,j) + &
                     tmp(i,4,k,2) * filterMat_tr(4,j) + &
                     tmp(i,5,k,2) * filterMat_tr(5,j) + &
                     tmp(i,6,k,2) * filterMat_tr(6,j) + &
                     tmp(i,7,k,2) * filterMat_tr(7,j) + &
                     tmp(i,8,k,2) * filterMat_tr(8,j) 

      MOMY_(i,j,k) = tmp(i,1,k,3) * filterMat_tr(1,j) + &
                     tmp(i,2,k,3) * filterMat_tr(2,j) + &
                     tmp(i,3,k,3) * filterMat_tr(3,j) + &
                     tmp(i,4,k,3) * filterMat_tr(4,j) + &
                     tmp(i,5,k,3) * filterMat_tr(5,j) + &
                     tmp(i,6,k,3) * filterMat_tr(6,j) + &
                     tmp(i,7,k,3) * filterMat_tr(7,j) + &
                     tmp(i,8,k,3) * filterMat_tr(8,j) 

      MOMZ_(i,j,k) = tmp(i,1,k,4) * filterMat_tr(1,j) + &
                     tmp(i,2,k,4) * filterMat_tr(2,j) + &
                     tmp(i,3,k,4) * filterMat_tr(3,j) + &
                     tmp(i,4,k,4) * filterMat_tr(4,j) + &
                     tmp(i,5,k,4) * filterMat_tr(5,j) + &
                     tmp(i,6,k,4) * filterMat_tr(6,j) + &
                     tmp(i,7,k,4) * filterMat_tr(7,j) + &
                     tmp(i,8,k,4) * filterMat_tr(8,j) 

      DRHOT_(i,j,k) = tmp(i,1,k,5) * filterMat_tr(1,j) + &
                      tmp(i,2,k,5) * filterMat_tr(2,j) + &
                      tmp(i,3,k,5) * filterMat_tr(3,j) + &
                      tmp(i,4,k,5) * filterMat_tr(4,j) + &
                      tmp(i,5,k,5) * filterMat_tr(5,j) + &
                      tmp(i,6,k,5) * filterMat_tr(6,j) + &
                      tmp(i,7,k,5) * filterMat_tr(7,j) + &
                      tmp(i,8,k,5) * filterMat_tr(8,j) 

    end do
    end do
    end do

    !-- z direction
    do k=1, 8
    do j=1, 8
    do i=1, 8
      tmp(i,j,k,1) = DDENS_(i,j,1) * filterMat_tr(1,k) + &
                     DDENS_(i,j,2) * filterMat_tr(2,k) + &
                     DDENS_(i,j,3) * filterMat_tr(3,k) + &
                     DDENS_(i,j,4) * filterMat_tr(4,k) + &
                     DDENS_(i,j,5) * filterMat_tr(5,k) + &
                     DDENS_(i,j,6) * filterMat_tr(6,k) + &
                     DDENS_(i,j,7) * filterMat_tr(7,k) + &
                     DDENS_(i,j,8) * filterMat_tr(8,k) 

      tmp(i,j,k,2) = MOMX_(i,j,1) * filterMat_tr(1,k) + &
                     MOMX_(i,j,2) * filterMat_tr(2,k) + &
                     MOMX_(i,j,3) * filterMat_tr(3,k) + &
                     MOMX_(i,j,4) * filterMat_tr(4,k) + &
                     MOMX_(i,j,5) * filterMat_tr(5,k) + &
                     MOMX_(i,j,6) * filterMat_tr(6,k) + &
                     MOMX_(i,j,7) * filterMat_tr(7,k) + &
                     MOMX_(i,j,8) * filterMat_tr(8,k) 

      tmp(i,j,k,3) = MOMY_(i,j,1) * filterMat_tr(1,k) + &
                     MOMY_(i,j,2) * filterMat_tr(2,k) + &
                     MOMY_(i,j,3) * filterMat_tr(3,k) + &
                     MOMY_(i,j,4) * filterMat_tr(4,k) + &
                     MOMY_(i,j,5) * filterMat_tr(5,k) + &
                     MOMY_(i,j,6) * filterMat_tr(6,k) + &
                     MOMY_(i,j,7) * filterMat_tr(7,k) + &
                     MOMY_(i,j,8) * filterMat_tr(8,k) 

      tmp(i,j,k,4) = MOMZ_(i,j,1) * filterMat_tr(1,k) + &
                     MOMZ_(i,j,2) * filterMat_tr(2,k) + &
                     MOMZ_(i,j,3) * filterMat_tr(3,k) + &
                     MOMZ_(i,j,4) * filterMat_tr(4,k) + &
                     MOMZ_(i,j,5) * filterMat_tr(5,k) + &
                     MOMZ_(i,j,6) * filterMat_tr(6,k) + &
                     MOMZ_(i,j,7) * filterMat_tr(7,k) + &
                     MOMZ_(i,j,8) * filterMat_tr(8,k) 

      tmp(i,j,k,5) = DRHOT_(i,j,1) * filterMat_tr(1,k) + &
                     DRHOT_(i,j,2) * filterMat_tr(2,k) + &
                     DRHOT_(i,j,3) * filterMat_tr(3,k) + &
                     DRHOT_(i,j,4) * filterMat_tr(4,k) + &
                     DRHOT_(i,j,5) * filterMat_tr(5,k) + &
                     DRHOT_(i,j,6) * filterMat_tr(6,k) + &
                     DRHOT_(i,j,7) * filterMat_tr(7,k) + &
                     DRHOT_(i,j,8) * filterMat_tr(8,k) 

    end do
    end do
    end do
  end subroutine apply_filter_xyz_direction

end module scale_atm_dyn_dgm_modalfilter
