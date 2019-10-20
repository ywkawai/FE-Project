!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_vars
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_meshfield_base, only: MeshField2D
  use scale_element_base, only: ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: ATMOS_VARS_setup
  public :: ATMOS_VARS_finalize
  public :: ATMOS_VARS_output

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type(MeshField2D), public, save, target :: DDENS
  type(MeshField2D), public, save, target :: MOMX
  type(MeshField2D), public, save,target :: MOMZ
  type(MeshField2D), public, save,target :: DRHOT
  
  type(MeshField2D), public, save, target :: U
  type(MeshField2D), public, save, target :: W
  type(MeshField2D), public, save, target :: DPRES
  type(MeshField2D), public, save, target :: TEMP
  type(MeshField2D), public, save, target :: DTHETA
  type(MeshField2D), public, save, target :: GxU, GzU
  type(MeshField2D), public, save, target :: GxW, GzW
  type(MeshField2D), public, save, target :: GxPT, GzPT

  type(MeshField2D), public, save, target :: PRES_hydro
  type(MeshField2D), public, save, target :: DENS_hydro


  integer, public, parameter :: VARS_DDENS_ID = 1
  integer, public, parameter :: VARS_MOMX_ID  = 2
  integer, public, parameter :: VARS_MOMZ_ID  = 3
  integer, public, parameter :: VARS_DRHOT_ID = 4
  integer, public, parameter :: PROG_VARS_NUM = 4

  integer, public, parameter :: VARS_U_ID       = 5
  integer, public, parameter :: VARS_W_ID       = 6
  integer, public, parameter :: VARS_DPRES_ID   = 7
  integer, public, parameter :: VARS_TEMP_ID    = 8
  integer, public, parameter :: VARS_DTHETA_ID  = 9
  integer, public, parameter :: DIAG_VARS_NUM   = 5

  integer, public, parameter :: VARS_PRES_HYDRO_ID = 10
  integer, public, parameter :: VARS_DENS_HYDRO_ID = 11
  integer, public, parameter :: AUX_VARS_NUM       = 2

  integer, public, parameter :: VARS_GxU_ID       = 1
  integer, public, parameter :: VARS_GzU_ID       = 2
  integer, public, parameter :: VARS_GxW_ID       = 3
  integer, public, parameter :: VARS_GzW_ID       = 4
  integer, public, parameter :: VARS_GxPT_ID      = 5
  integer, public, parameter :: VARS_GzPT_ID      = 6
  integer, public, parameter :: AUX_DIFFVARS_NUM  = 6


  integer, public, parameter :: VARS_TOT_NUM  = PROG_VARS_NUM + DIAG_VARS_NUM + AUX_VARS_NUM + 3

  type(MeshFieldCommRectDom2D), public, save :: PROG_VARS_comm
  type(MeshFieldContainer), public, save :: PROG_VARS_list(PROG_VARS_NUM)

  type(MeshFieldCommRectDom2D), public, save :: AUX_DIFFVARS_comm
  type(MeshFieldContainer), public, save :: AUX_DIFFVARS_list(AUX_DIFFVARS_NUM)

  type(MeshField2D), public, save :: DxMOMX
  type(MeshField2D), public, save  :: DzMOMZ
  type(MeshField2D), public, save :: LiftDDENS
  integer, public, parameter :: VARS_DxMOMX_ID       = 12
  integer, public, parameter :: VARS_DzMOMZ_ID       = 13
  integer, public, parameter :: VARS_LiftDDENS_ID       = 14
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  integer :: HST_ID(VARS_TOT_NUM)


contains  
  subroutine ATMOS_VARS_setup()
    use scale_file_history, only: FILE_HISTORY_reg
    use mod_atmos_mesh, only: mesh
    implicit none

    integer :: n, k
    type(LocalMesh2D), pointer :: lcmesh
    !-------------------------------------------------------------------------
  
    call DxMOMX%Init( "DxMOMX", "", mesh)
    call DzMOMZ%Init( "DzMOMZ", "", mesh)
    call LiftDDENS%Init( "LiftDDENS", "", mesh)
    
    call DDENS%Init( "DDENS", "kg/m3", mesh )
    call MOMX%Init( "MOMX", "kg/m2/s", mesh )
    call MOMZ%Init( "MOMZ", "kg/m2/s", mesh )
    call DRHOT%Init( "DRHOT", "kg/m3.K", mesh )

    call U%Init( "U", "m/s", mesh )
    call W%Init( "W", "m/s", mesh )
    call DPRES%Init( "DPRES", "kg.s-2.m-1", mesh )
    call TEMP%Init( "TEMP", "K", mesh )    
    call DTHETA%Init( "DTHETA", "K", mesh )    

    call GxU%Init( "GxU", "s-1", mesh )
    call GzU%Init( "GzU", "s-1", mesh )
    call GxW%Init( "GxW", "s-1", mesh )
    call GzW%Init( "GzW", "s-1", mesh )
    call GxPT%Init( "GxPT", "K/m", mesh )
    call GzPT%Init( "GzPT", "K/m", mesh )
    do n=1, mesh%LOCAL_MESH_NUM
      GxU%local(n)%val(:,:) = 0.0_RP
      GzU%local(n)%val(:,:) = 0.0_RP
      GxW%local(n)%val(:,:) = 0.0_RP
      GzW%local(n)%val(:,:) = 0.0_RP
      GxPT%local(n)%val(:,:) = 0.0_RP
      GzPT%local(n)%val(:,:) = 0.0_RP
    end do

    call DENS_hydro%Init( "DENS_hydro", "kg/m3", mesh )
    call PRES_hydro%Init( "PRES_hydro", "kg.s-2.m-1", mesh )

    call PROG_VARS_comm%Init(PROG_VARS_NUM, 0, mesh)
    PROG_VARS_list(VARS_DDENS_ID)%field2d => DDENS
    PROG_VARS_list(VARS_MOMX_ID)%field2d => MOMX
    PROG_VARS_list(VARS_MOMZ_ID)%field2d => MOMZ
    PROG_VARS_list(VARS_DRHOT_ID)%field2d => DRHOT

    call AUX_DIFFVARS_comm%Init(AUX_DIFFVARS_NUM, 0, mesh)
    AUX_DIFFVARS_list(VARS_GxU_ID)%field2d => GxU
    AUX_DIFFVARS_list(VARS_GzU_ID)%field2d => GzU
    AUX_DIFFVARS_list(VARS_GxW_ID)%field2d => GxW
    AUX_DIFFVARS_list(VARS_GzW_ID)%field2d => GzW
    AUX_DIFFVARS_list(VARS_GxPT_ID)%field2d => GxPT
    AUX_DIFFVARS_list(VARS_GzPT_ID)%field2d => GzPT

    call FILE_HISTORY_reg( DDENS%varname, "deviation of density", DDENS%unit, HST_ID(VARS_DDENS_ID), dim_type='XY')
    call FILE_HISTORY_reg( MOMX%varname, "momentum (X)", MOMX%unit, HST_ID(VARS_MOMX_ID), dim_type='XY')
    call FILE_HISTORY_reg( MOMZ%varname, "momentum (Z)", MOMZ%unit, HST_ID(VARS_MOMZ_ID), dim_type='XY')
    call FILE_HISTORY_reg( DRHOT%varname, "deviation of rho * theta", DRHOT%unit, HST_ID(VARS_DRHOT_ID), dim_type='XY')

    call FILE_HISTORY_reg( U%varname, "velocity (X)", U%unit, HST_ID(VARS_U_ID), dim_type='XY')
    call FILE_HISTORY_reg( W%varname, "velocity (Z)", W%unit, HST_ID(VARS_W_ID), dim_type='XY')
    call FILE_HISTORY_reg( DPRES%varname, "deviation of pressure", DPRES%unit, HST_ID(VARS_DPRES_ID), dim_type='XY')
    call FILE_HISTORY_reg( TEMP%varname, "temperature", TEMP%unit, HST_ID(VARS_TEMP_ID), dim_type='XY')
    call FILE_HISTORY_reg( DTHETA%varname, "deviation of potential temperature", DTHETA%unit, HST_ID(VARS_DTHETA_ID), dim_type='XY')

    call FILE_HISTORY_reg( PRES_hydro%varname, "hydrostatic pressure", DPRES%unit, HST_ID(VARS_PRES_HYDRO_ID), dim_type='XY')
    call FILE_HISTORY_reg( DENS_hydro%varname, "hydrostatic density", TEMP%unit, HST_ID(VARS_DENS_HYDRO_ID), dim_type='XY')

    call FILE_HISTORY_reg( DxMOMX%varname, "DxMOMX", DxMOMX%unit, HST_ID(VARS_DxMOMX_ID), dim_type='XY' )
    call FILE_HISTORY_reg( DzMOMZ%varname, "DzMOMZ", DzMOMZ%unit, HST_ID(VARS_DzMOMZ_ID), dim_type='XY' )
    call FILE_HISTORY_reg( LiftDDENS%varname, "LiftDDENS", LiftDDENS%unit, HST_ID(VARS_LiftDDENS_ID), dim_type='XY' )
    return
  end subroutine ATMOS_VARS_setup

  subroutine ATMOS_VARS_finalize()
    implicit none
    !-------------------------------------------------------------------------

    call PROG_VARS_comm%Final()
    call AUX_DIFFVARS_comm%Final()

    call DDENS%Final()
    call MOMX%Final()
    call MOMZ%Final()
    call DRHOT%Final()
    
    call U%Final()
    call W%Final()
    call DPRES%Final()
    call TEMP%Final()
    call DTHETA%Final()

    call GxU%Final()
    call GzU%Final()
    call GxW%Final()
    call GzW%Final()
    call GxPT%Final()
    call GzPT%Final()

    call DENS_hydro%Final()
    call PRES_hydro%Final()

    return
  end subroutine ATMOS_VARS_finalize


  subroutine ATMOS_VARS_output( tsec_ )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_put,   &
      FILE_HISTORY_meshfield_write
    use mod_atmos_mesh, only: mesh
    implicit none

    real(RP), intent(in) :: tsec_
    
    integer :: n
    type(LocalMesh2D), pointer :: lcmesh
    !-------------------------------------------------------------------------

    call FILE_HISTORY_meshfield_put(HST_ID(VARS_DDENS_ID), DDENS)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_MOMX_ID), MOMX)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_MOMZ_ID), MOMZ)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_DRHOT_ID), DRHOT)

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)

      call vars_calc_diagnoseVar( &
        U%local(n)%val, W%local(n)%val, DPRES%local(n)%val, TEMP%local(n)%val, DTHETA%local(n)%val,   &
        DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,                 &
        PRES_hydro%local(n)%val, DENS_hydro%local(n)%val,                                             &
        lcmesh, lcmesh%refElem2D )
    end do
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_U_ID), U)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_W_ID), W)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_DTHETA_ID), DTHETA)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_DPRES_ID), DPRES)

    call FILE_HISTORY_meshfield_put(HST_ID(VARS_DENS_HYDRO_ID), DENS_hydro)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_PRES_HYDRO_ID), PRES_hydro)

    call FILE_HISTORY_meshfield_put(HST_ID(VARS_DxMOMX_ID), DxMOMX)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_DzMOMZ_ID), DzMOMZ)
    call FILE_HISTORY_meshfield_put(HST_ID(VARS_LiftDDENS_ID), LiftDDENS)

    call FILE_HISTORY_meshfield_write()   

    return
  end subroutine ATMOS_VARS_output

  subroutine vars_calc_diagnoseVar( &
    U_, W_, DPRES_, TEMP_, DTHETA_,                      &
    DDENS_, MOMX_, MOMZ_, DRHOT_, PRES_hyd, DENS_hyd,    &
    lcmesh, elem )

    use scale_const, only: &
      GRAV => CONST_GRAV,  &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    implicit none

    type(LocalMesh2D), intent(in) :: lcmesh
    type(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: U_(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: W_(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DPRES_(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: TEMP_(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DTHETA_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT_(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)

    integer :: k
    real(RP) :: RHOT(elem%Np), DENS(elem%Np), PRES(elem%Np), THETA(elem%Np)

    !-------------------------------------------------------------------------

    do k=1, lcmesh%Ne
      DENS(:) = DDENS_(:,k) + DENS_hyd(:,k)

      U_(:,k) = MOMX_(:,k)/DENS(:)
      W_(:,k) = MOMZ_(:,k)/DENS(:)

      RHOT(:) = PRES00/Rdry * (PRES_hyd(:,k)/PRES00)**(CVdry/CPdry) + DRHOT_(:,k)
      
      PRES(:) = PRES00 * (Rdry*RHOT(:)/PRES00)**(CPdry/Cvdry)
      THETA(:) = RHOT(:)/DENS(:)

      DPRES_(:,k) = PRES(:) - PRES_hyd(:,k)
      TEMP_(:,k) = PRES(:)/(Rdry*DENS(:))

      DTHETA_(:,k) = THETA(:) - PRES_hyd(:,k)/(Rdry*DENS_hyd(:,k)) * (PRES00/PRES_hyd(:,k))**(Rdry/CPdry)
    end do

    return
  end subroutine vars_calc_diagnoseVar

end module mod_atmos_vars