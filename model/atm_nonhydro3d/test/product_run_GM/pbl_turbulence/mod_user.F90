!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of planetary boudary layer turbulence
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort  

  use scale_const, only: &
    PI => CONST_PI,        &
    GRAV => CONST_GRAV,    &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    Pstd   => CONST_Pstd  
  
  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: &
    ElementBase3D, ElementBase2D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_meshfield_base, only: MeshField3D

  use scale_sparsemat, only: &
    SparseMat, SparseMat_matmul
  
  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public, extends(UserBase) :: User
  contains
    procedure :: mkinit_ => USER_mkinit
    generic :: mkinit => mkinit_
    procedure :: setup_ => USER_setup
    generic :: setup => setup_
    procedure :: calc_tendency => USER_calc_tendency
  end type User

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical :: is_PREShyd_ref_set

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('idealized_pbl_turbulence')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_pbl_turblence )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL
  subroutine USER_setup( this, atm )
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do                   = .false. !< do user step?
    namelist / PARAM_USER / &
       USER_do

    integer :: ierr    
    !------------------------------------------


    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    call this%UserBase%Setup( atm, USER_do )
    !-

    is_PREShyd_ref_set = .false. 

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRESHYD_VID => AUXVAR_PRESHYDRO_ID,        &
      PRESHYD_REF_VID => AUXVAR_PRESHYDRO_REF_ID    
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(MeshField3D), pointer :: PRES_hyd
    class(MeshField3D), pointer :: PRES_hyd_ref
    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: domid, ke
    !------------------------------------------


    ! Set reference hydrostatic pressure
    if ( .not. is_PREShyd_ref_set ) then
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_VID, PRES_hyd )
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_REF_VID, PRES_hyd_ref )

      do domid=1, PRES_hyd_ref%mesh%LOCAL_MESH_NUM
        lcmesh3D => PRES_hyd_ref%mesh%lcmesh_list(domid)
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          PRES_hyd_ref%local(domid)%val(:,ke) = PRES_hyd%local(domid)%val(:,ke)
        end do
      end do
      is_PREShyd_ref_set = .true.
    end if

    return
  end subroutine USER_calc_tendency

  !------

!OCL SERIAL
  subroutine exp_SetInitCond_pbl_turblence( this,                          &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_random, only: &
      RANDOM_uniform
    use scale_atm_dyn_dgm_hydrostatic, only:  &
      hydrostatic_calc_basicstate_constPTLAPS, &
      hydrostaic_build_rho_XYZ 
    
    use mod_experiment, only: &
      TracerLocalMeshField_ptr 
         
    implicit none

    class(Experiment), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
    type(TracerLocalMeshField_ptr), intent(inout) :: tracer_field_list(:)    
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    real(RP), intent(in) :: dom_zmin, dom_zmax
    
    real(RP) :: ENV_PRES_SFC    
    real(RP) :: ENV_U          = 5.0_RP
    real(RP) :: ENV_V          = 0.0_RP
    real(RP) :: ENV_THETA_SFC  = 298.0_RP 
    real(RP) :: ENV_THETA_LAPS = 4.0E-3_RP
    real(RP) :: RANDOM_THETA   = 1.0_RP
    real(RP) :: RANDOM_U       = 0.0_RP
    real(RP) :: RANDOM_V       = 0.0_RP        

    namelist /PARAM_EXP/ &
      ENV_U,            &
      ENV_THETA_SFC,    &
      ENV_THETA_LAPS,   &
      ENV_PRES_SFC,     &
      RANDOM_THETA

    real(RP) :: rndm(elem%Np) 
    real(RP) :: POT (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: DENS(elem%Np)

    integer :: ke, ke2D
    integer :: ke_x, ke_y, ke_z
    integer :: ierr

    real(RP), allocatable :: bnd_SFC_PRES(:,:)
    real(RP) :: sfc_rhot(elem%Nfp_v)

    !-----------------------------------------------------------------------------

    ENV_PRES_SFC = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("PBL_TURBULENCE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("PBL_TURBULENCE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    
    call hydrostatic_calc_basicstate_constPTLAPS( DENS_hyd, PRES_hyd,     &
      ENV_THETA_LAPS, ENV_THETA_SFC, ENV_PRES_SFC,                        &
      lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3),   &
      lcmesh, elem )
    
    !---

    allocate( bnd_SFC_PRES(elem%Nfp_v,lcmesh%Ne2DA) )

    !$omp parallel do private( ke, ke2D, rndm, sfc_rhot )
    do ke_z=1, lcmesh%NeZ
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX
      ke2D = ke_x + (ke_y-1)*lcmesh%NeX
      ke = ke2D + (ke_z-1)*lcmesh%NeX*lcmesh%NeY

      call RANDOM_uniform( rndm )
      POT(:,ke_z,ke_x,ke_y) = ENV_THETA_SFC + ENV_THETA_LAPS * z(:,ke)    &
                            + ( rndm(:) * 2.0_RP - 1.0_RP ) * RANDOM_THETA
      
      if ( ke_z == 1 ) then
        sfc_rhot(:) = DENS_hyd(elem%Hslice(:,1),ke2D) * POT(elem%Hslice(:,1),ke_z,ke_x,ke_y)        
        bnd_SFC_PRES(:,ke2D) = PRES00 * ( Rdry * sfc_rhot(:) / PRES00 )**( CPdry/CVdry )
      end if
    end do
    end do
    end do
    
    call hydrostaic_build_rho_XYZ( DDENS, &
      DENS_hyd, PRES_hyd, POT,                                            &
      lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3),   &
      lcmesh, elem, bnd_SFC_PRES=bnd_SFC_PRES                             )

    !$parallel do private( ke, DENS, rndm, ke_x, ke_y )
    do ke_z=1, lcmesh%NeZ
    do ke_y=1, lcmesh%NeY
    do ke_x=1, lcmesh%NeX  
      ke = ke_x + (ke_y-1)*lcmesh%NeX + (ke_z-1)*lcmesh%NeX*lcmesh%NeY

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      DRHOT(:,ke) = DENS(:) * POT(:,ke_z,ke_x,ke_y)                             &
                  - DENS_hyd(:,ke) * ( ENV_THETA_SFC + ENV_THETA_LAPS * z(:,ke) )

      call RANDOM_uniform( rndm )
      MOMX(:,ke) = DENS(:) * (ENV_U + (rndm(:) * 2.0_RP - 1.0_RP ) * RANDOM_U )
      MOMY(:,ke) = DENS(:) * (ENV_V + (rndm(:) * 2.0_RP - 1.0_RP ) * RANDOM_V )    
    end do
    end do
    end do

    return
  end subroutine exp_SetInitCond_pbl_turblence

end module mod_user
