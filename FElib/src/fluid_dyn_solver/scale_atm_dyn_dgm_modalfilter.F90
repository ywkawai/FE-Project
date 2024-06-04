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
    real(RP) :: FilterMat_tr(12,12)
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

        do ii=1, elem%Np
          DDENS_(ii,ke) = lmesh%Gsqrt(ii,ke) * DDENS_(ii,ke)
          MOMX_(ii,ke) = lmesh%Gsqrt(ii,ke) * MOMX_(ii,ke)
          MOMY_(ii,ke) = lmesh%Gsqrt(ii,ke) * MOMY_(ii,ke)
          MOMZ_(ii,ke) = lmesh%Gsqrt(ii,ke) * MOMZ_(ii,ke)
          DRHOT_(ii,ke) = lmesh%Gsqrt(ii,ke) * DRHOT_(ii,ke)
        end do

        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, DDENS_(:,ke), tmp(:,1))
        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, MOMX_(:,ke),  tmp(:,2))
        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, MOMY_(:,ke),  tmp(:,3))
        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, MOMZ_(:,ke),  tmp(:,4))
        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, DRHOT_(:,ke), tmp(:,5))

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

        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, DDENS_(:,ke), tmp(:,1))
        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, MOMX_(:,ke),  tmp(:,2))
        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, MOMY_(:,ke),  tmp(:,3))
        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, MOMZ_(:,ke),  tmp(:,4))
        call apply_filter_xyz_direction(filter%FilterMat, FilterMat_tr, DRHOT_(:,ke), tmp(:,5))

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
  subroutine apply_filter_xyz_direction(filterMat, filterMat_tr, q_in, q_tmp )
    implicit none

    real(RP), intent(in) :: filterMat(12, 12)
    real(RP), intent(in) :: filterMat_tr(12, 12)
    real(RP), intent(inout) :: q_in(12, 12, 12)
    real(RP), intent(inout) :: q_tmp(12, 12, 12)
    
    real(RP) :: tmp1, tmp2, tmp3
    integer :: i, j, k

    !-- x direction
    do k=1, 12
    do j=1, 12
    do i=1, 12

      tmp1 = filterMat(i,1)  * q_in(1,j,k) + &
             filterMat(i,2)  * q_in(2,j,k) + &
             filterMat(i,3)  * q_in(3,j,k) + & 
             filterMat(i,4)  * q_in(4,j,k)
              
      tmp2 = filterMat(i,5)  * q_in(5,j,k) + & 
             filterMat(i,6)  * q_in(6,j,k) + & 
             filterMat(i,7)  * q_in(7,j,k) + & 
             filterMat(i,8)  * q_in(8,j,k) 

      tmp3 = filterMat(i,9)  * q_in(9,j,k) + & 
             filterMat(i,10) * q_in(10,j,k) + & 
             filterMat(i,11) * q_in(11,j,k) + & 
             filterMat(i,12) * q_in(12,j,k)  

      q_tmp(i,j,k) = tmp1 + tmp2 + tmp3

    end do
    end do
    end do

    !-- y direction
    do k=1, 12
    do j=1, 12
    do i=1, 12

      tmp1 = q_tmp(i,1,k)  * filterMat_tr(1,j) + &
             q_tmp(i,2,k)  * filterMat_tr(2,j) + &
             q_tmp(i,3,k)  * filterMat_tr(3,j) + &
             q_tmp(i,4,k)  * filterMat_tr(4,j)

      tmp2 = q_tmp(i,5,k)  * filterMat_tr(5,j) + &
             q_tmp(i,6,k)  * filterMat_tr(6,j) + &
             q_tmp(i,7,k)  * filterMat_tr(7,j) + &
             q_tmp(i,8,k)  * filterMat_tr(8,j)

      tmp3 = q_tmp(i,9,k)  * filterMat_tr(9,j) + &
             q_tmp(i,10,k) * filterMat_tr(10,j) + &
             q_tmp(i,11,k) * filterMat_tr(11,j) + &
             q_tmp(i,12,k) * filterMat_tr(12,j) 

      q_in(i,j,k) = tmp1 + tmp2 + tmp3

    end do
    end do
    end do

    !-- z direction
    do k=1, 12
    do j=1, 12
    do i=1, 12

      tmp1 = q_in(i,j,1)  * filterMat_tr(1,k) + &
             q_in(i,j,2)  * filterMat_tr(2,k) + &
             q_in(i,j,3)  * filterMat_tr(3,k) + &
             q_in(i,j,4)  * filterMat_tr(4,k)

      tmp2 = q_in(i,j,5)  * filterMat_tr(5,k) + &
             q_in(i,j,6)  * filterMat_tr(6,k) + &
             q_in(i,j,7)  * filterMat_tr(7,k) + &
             q_in(i,j,8)  * filterMat_tr(8,k)

      tmp3 = q_in(i,j,9)  * filterMat_tr(9,k) + &
             q_in(i,j,10) * filterMat_tr(10,k) + &
             q_in(i,j,11) * filterMat_tr(11,k) + &
             q_in(i,j,12) * filterMat_tr(12,k) 

      q_tmp(i,j,k) = tmp1 + tmp2 + tmp3

    end do
    end do
    end do
  end subroutine apply_filter_xyz_direction

end module scale_atm_dyn_dgm_modalfilter
