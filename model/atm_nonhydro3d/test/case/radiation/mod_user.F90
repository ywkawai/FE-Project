!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of radiation scheme.
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_user

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort  

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D   
  use scale_meshfield_base, only: MeshField3D

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

  type(MeshField3D), private :: PRES_diff

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('radiation')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_radiation )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL
  subroutine USER_setup( this, atm )
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do        = .false. !< do user
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
    if ( USER_do ) call PRES_diff%Init( 'PRES_diff', 'Pa', atm%mesh%ptr_mesh )
  
    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    implicit none

    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    if ( this%USER_do ) then
      call atm%vars%Calc_diagVar( 'PRES_diff', PRES_diff )
      call FILE_HISTORY_meshfield_in( PRES_diff, "perturbation of PRES" )
    end if

    return
  end subroutine USER_calc_tendency

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_radiation( this, &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    use scale_tracer, only: &
      TRACER_inq_id
    use scale_const, only: &
      PI => CONST_PI,       &
      GRAV => CONST_GRAV,   &
      Rdry => CONST_Rdry,   &
      Rvap => CONST_Rvap,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    use scale_atmos_saturation, only: &
      ATMOS_SATURATION_psat_all,      &
      ATMOS_SATURATION_pres2qsat_all
    use scale_atmos_hydrometeor, only: &
      CV_VAPOR, &
      CV_WATER, &
      CP_VAPOR, &
      CP_WATER    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT, &
      hydrostatic_build_rho_XYZ
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
    
    real(RP) :: TEMP0    = 250.0_RP !< Initial temperature [K]
    real(RP) :: SFC_PRES = 1.0E5_RP !< Surface pressure [Pa]
    real(RP) :: ENV_RH   =  0.0_RP  !< Relative Humidity of environment [%]
    integer :: NITER_RH  = 3
    namelist /PARAM_EXP/ &
      TEMP0,  &
      ENV_RH, &
      NITER_RH
    integer :: ierr

    integer :: iq_QV    
    
    integer :: ke, ke2D
    integer :: ke_x, ke_y, ke_z
    integer :: p, p3, p2D

    real(RP) :: PT_tmp(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: Rtot  (elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: CPtot(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: CPtot_ov_CVtot(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)
    real(RP) :: bnd_SFC_PRES(elem%Nnode_h1D**2,lcmesh%lcmesh2D%NeA)
    real(RP) :: QV(elem%Np)
    real(RP) :: PRES(elem%Np)
    real(RP) :: psat0

    integer :: itr
    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("RADIATION_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("RADIATION_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    call hydrostatic_calc_basicstate_constT( &
      DENS_hyd, PRES_hyd,                    & ! (out)
      TEMP0, SFC_PRES, x, y, z, lcmesh, elem ) ! (in)

    if ( ENV_RH > 0.0_RP ) then
      LOG_INFO("RADIATION_setup",*) 'Calculate QV from RH'
      LOG_INFO("RADIATION_setup",*) 'ENV_RH = ', ENV_RH, ' [%]'
      LOG_INFO("RADIATION_setup",*) 'NITER_RH = ', NITER_RH

      call TRACER_inq_id( "QV", iq_QV )
      call ATMOS_SATURATION_psat_all( TEMP0, psat0 )

      do itr=1, NITER_RH
        LOG_INFO("RADIATION_setup",*) 'RH iteration: ', itr

        !$omp parallel do collapse(3) private(ke_z,ke_x,ke_y,ke,ke2D,p3,p2D,p, QV,PRES)
        do ke_y=1, lcmesh%NeY
        do ke_x=1, lcmesh%NeX
        do ke_z=1, lcmesh%NeZ
          ke2D = ke_x + (ke_y-1)*lcmesh%NeX
          ke = ke2D + (ke_z-1)*lcmesh%NeX*lcmesh%NeY

          do p3=1, elem%Nnode_v
          do p2D=1, elem%Nnode_h1D**2
            p = p2D + (p3 - 1)*elem%Nnode_h1D**2
            QV(p) = ENV_RH * 1.0E-2_RP * psat0 / ( ( DENS_hyd(p,ke) + DDENS(p,ke) ) * Rvap * TEMP0 )
          end do
          end do
          tracer_field_list(iq_QV)%ptr%val(:,ke) = QV(:)

          Rtot          (:,ke_z,ke_x,ke_y) = Rdry  * ( 1.0_RP - QV(:) ) + Rvap     * QV(:)
          CPtot         (:,ke_z,ke_x,ke_y) = CPdry * ( 1.0_RP - QV(:) ) + CP_VAPOR * QV(:)
          CPtot_ov_CVtot(:,ke_z,ke_x,ke_y) = CPtot(:,ke_z,ke_x,ke_y)                         &
                                          / ( CVdry * ( 1.0_RP - QV(:) ) + CV_VAPOR * QV(:) )

          PRES(:) = ( DENS_hyd(:,ke) + DDENS(:,ke) ) * Rtot(:,ke_z,ke_x,ke_y) * TEMP0
          PT_tmp(:,ke_z,ke_x,ke_y) = TEMP0 * ( PRES00 / PRES(:) )**( Rtot(:,ke_z,ke_x,ke_y) / CPtot(:,ke_z,ke_x,ke_y) )
          if ( ke_z == 1 ) then
            bnd_SFC_PRES(:,ke2D) = PRES(elem%Hslice(:,1))
          end if
        end do
        end do
        end do

        call hydrostatic_build_rho_XYZ( DDENS, & ! (out)
          DENS_hyd, PRES_hyd, PT_tmp,         & ! (in)
          Rtot, CPtot_ov_CVtot,               & ! (in)
          x, y, z, lcmesh, elem,              & ! (in)
          bnd_SFC_PRES                        ) ! (in)
      end do ! End of RH iteration

      !$omp parallel do collapse(3) private(ke_z,ke_x,ke_y,ke,ke2D)
      do ke_y=1, lcmesh%NeY
      do ke_x=1, lcmesh%NeX
      do ke_z=1, lcmesh%NeZ
        ke2D = ke_x + (ke_y-1)*lcmesh%NeX
        ke = ke2D + (ke_z-1)*lcmesh%NeX*lcmesh%NeY
        DRHOT(:,ke) = ( DENS_hyd(:,ke) + DDENS(:,ke) ) * PT_tmp(:,ke_z,ke_x,ke_y) &
          - DENS_hyd(:,ke) * TEMP0 * ( PRES00 / PRES_hyd(:,ke) )**( Rdry / CPdry )
      end do
      end do
      end do

    end if


    return
  end subroutine exp_SetInitCond_radiation

end module mod_user
