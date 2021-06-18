!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!      Sponge layer for Atmospheric dynamical process. 
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_spongelayer
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: &
    ElementBase, ElementBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_spongelayer_add_tend

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  private :: calc_wdampcoef

contains

  subroutine atm_dyn_dgm_spongelayer_add_tend( MOMX_dt, MOMY_dt, MOMZ_dt, &
    MOMX_, MOMY_, MOMZ_, wdamp_tau, wdamp_height, hveldamp_flag,          &
    lmesh, elem   )

    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(inout) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: wdamp_tau
    real(RP), intent(in) :: wdamp_height
    logical, intent(in) :: hveldamp_flag

    integer :: ke
    integer :: ke_x, ke_y, ke_z
    integer :: keZtop
    real(RP) :: wdamp_coef(elem%Np)
    real(RP) :: zTop(elem%Nnode_h1D**2)
    real(RP) :: s
    !-----------------------------------------------------------------

    if ( hveldamp_flag ) then
      s = 1.0_RP
    else
      s = 0.0_RP
    end if

    !$omp parallel do collapse(3) private(ke,keZtop,zTop,wdamp_coef)
    do ke_z = 1, lmesh%NeZ
    do ke_y = 1, lmesh%NeY
    do ke_x = 1, lmesh%NeX
      ke = ke_x + (ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      keZtop =  ke_x + (ke_y-1)*lmesh%NeX + (lmesh%NeZ-1)*lmesh%NeX*lmesh%NeY
      zTop(:) = lmesh%pos_en(elem%Hslice(:,elem%Nnode_v),keZtop,3)

      call calc_wdampcoef( &
        wdamp_tau, wdamp_height, lmesh%pos_en(:,ke,3), zTop(:), &
        elem%Nnode_h1D, elem%Nnode_v,                           &
        wdamp_coef(:) )

      MOMX_dt(:,ke) = MOMX_dt(:,ke) - s * wdamp_coef(:) * MOMX_(:,ke)
      MOMY_dt(:,ke) = MOMY_dt(:,ke) - s * wdamp_coef(:) * MOMY_(:,ke)
      MOMZ_dt(:,ke) = MOMZ_dt(:,ke) - wdamp_coef(:) * MOMZ_(:,ke)
    end do
    end do
    end do

    return
  end subroutine atm_dyn_dgm_spongelayer_add_tend

!-- private ------------------------------

!OCL SERIAL
  subroutine calc_wdampcoef( &
    wdamp_tau, wdamp_height, z, zTop, Nnode_h1D, Nnode_v, & ! (in)
    wdamp_coef                                            ) ! (out)

    use scale_const, only: &
      PI => CONST_PI
    implicit none

    integer, intent(in) :: Nnode_h1D
    integer, intent(in) :: Nnode_v
    real(RP), intent(out) :: wdamp_coef(Nnode_h1D**2,Nnode_v)
    real(RP), intent(in) :: wdamp_tau
    real(RP), intent(in) :: wdamp_height
    real(RP), intent(in) :: z(Nnode_h1D**2,Nnode_v)
    real(RP), intent(in) :: zTop(Nnode_h1D**2)

    integer :: p_z
    real(RP) :: sw(Nnode_h1D**2)
    real(RP) :: r_wdamp_tau
    !-----------------------------------------------------------------

    r_wdamp_tau = 1.0_RP / wdamp_tau
    do p_z=1, Nnode_v
      wdamp_coef(:,p_z) = 0.25_RP * r_wdamp_tau                                        &
        * ( 1.0_RP + sign( 1.0_RP, z(:,p_z) - wdamp_height )                         ) &
        * ( 1.0_RP - cos( PI * (z(:,p_z) - wdamp_height)/(zTop(:) - wdamp_height) )  )
    end do

    return
  end subroutine calc_wdampcoef

end module scale_atm_dyn_dgm_spongelayer