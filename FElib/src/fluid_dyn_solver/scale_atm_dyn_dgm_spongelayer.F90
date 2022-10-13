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
  use scale_mesh_base3d, only: MeshBase3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  type, public :: AtmDynSpongeLayer
    real(RP) :: wdamp_tau
    real(RP) :: wdamp_height
    logical  :: hveldamp_flag
  contains
    procedure :: Init => atm_dyn_dgm_spongelayer_Init
    procedure :: Final => atm_dyn_dgm_spongelayer_Final
    procedure :: AddTend => atm_dyn_dgm_spongelayer_add_tend
  end type AtmDynSpongeLayer

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

!OCL SERIAL
  subroutine atm_dyn_dgm_spongelayer_Init( this, mesh3D, dtsec )
    use scale_prc, only: PRC_abort
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
    implicit none

    class(AtmDynSpongeLayer), target, intent(inout) :: this
    class(MeshBase3D), target, intent(in) :: mesh3D
    real(RP), intent(in) :: dtsec

    real(RP) :: SL_WDAMP_TAU        = -1.0_RP ! the maximum tau for Rayleigh damping of w [s]
    real(RP) :: SL_WDAMP_HEIGHT     = -1.0_RP ! the height to start apply Rayleigh damping [m]
    integer  :: SL_WDAMP_LAYER      = -1      ! the vertical number of finite element to start apply Rayleigh damping [num]
    logical  :: SL_HORIVELDAMP_FLAG = .false. ! Is the horizontal velocity damped? 
    
    namelist /PARAM_ATMOS_DYN_SPONGELAYER/ &
      SL_WDAMP_TAU,                        &                
      SL_WDAMP_HEIGHT,                     &
      SL_WDAMP_LAYER,                      &
      SL_HORIVELDAMP_FLAG
    
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D
  
    integer :: NeGZ
    integer :: ierr
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_SPONGELAYER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup_spongelayer",*) 'Not found namelist. Default used.'
    else if( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup_spongelayer",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_SPONGELAYER. Check!'
      call PRC_abort
    end if
    LOG_NML(PARAM_ATMOS_DYN_SPONGELAYER)

    this%wdamp_tau    = SL_WDAMP_TAU 
    this%wdamp_height = SL_WDAMP_HEIGHT

    lcmesh3D => mesh3D%lcmesh_list(1)
    elem3D => lcmesh3D%refElem3D
    select type(mesh3D)
    type is (MeshCubeDom3D)
      NeGZ = mesh3D%NeGZ
    type is (MeshCubedSphereDom3D)
      NeGZ = mesh3D%NeGZ
    end select

    if ( SL_WDAMP_LAYER > NeGZ ) then
      LOG_ERROR("ATMOS_DYN_setup_spongelayer",*) 'SL_wdamp_layer should be less than total of vertical elements (NeGZ). Check!'
      call PRC_abort
    else if( SL_WDAMP_LAYER > 0 ) then
      this%wdamp_height = lcmesh3D%pos_en(1,1+(SL_WDAMP_LAYER-1)*lcmesh3D%NeX*lcmesh3D%NeY,3)
    end if
    if ( this%wdamp_tau < 0.0_RP ) then
      this%wdamp_tau = dtsec * 10.0_RP
    else if ( this%wdamp_tau < dtsec ) then
      LOG_ERROR("ATMOS_DYN_setup_spongelayer",*) 'SL_wdamp_tau should be larger than TIME_DT (ATMOS_DYN). Check!'
      call PRC_abort
    end if
    
    this%hveldamp_flag = SL_HORIVELDAMP_FLAG

    return
  end subroutine atm_dyn_dgm_spongelayer_Init

!OCL SERIAL
  subroutine atm_dyn_dgm_spongelayer_Final( this )
    implicit none
    class(AtmDynSpongeLayer), target, intent(inout) :: this
    !---------------------------------------------------------------

    return
  end subroutine atm_dyn_dgm_spongelayer_Final

  subroutine atm_dyn_dgm_spongelayer_add_tend( this, MOMX_dt, MOMY_dt, MOMZ_dt, &
    MOMX_, MOMY_, MOMZ_,                                                        &
    lmesh, elem   )

    implicit none

    class(AtmDynSpongeLayer), intent(in) :: this
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lmesh%NeA)

    integer :: ke
    integer :: ke_x, ke_y, ke_z
    integer :: keZtop
    real(RP) :: wdamp_coef(elem%Np)
    real(RP) :: zTop(elem%Nnode_h1D**2)
    real(RP) :: s
    !-----------------------------------------------------------------

    if ( this%hveldamp_flag ) then
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
        this%wdamp_tau, this%wdamp_height, lmesh%pos_en(:,ke,3), zTop(:), &
        elem%Nnode_h1D, elem%Nnode_v,                                     &
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