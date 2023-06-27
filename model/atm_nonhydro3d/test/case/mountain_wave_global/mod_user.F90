!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of mountain wave in global model
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

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D

  use mod_user_base, only: UserBase
  use mod_experiment, only: Experiment

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
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

  type(MeshField3D), private :: PT_diff

  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('mountain_wave_global')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_mountain_wave )
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
    if ( this%USER_do ) call PT_diff%Init( 'PT_diff', 'K', atm%mesh%ptr_mesh )

    return
  end subroutine USER_setup

!OCL SERIAL  
  subroutine USER_calc_tendency( this, atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in

    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CpDry => CONST_CPdry,    &
      RPlanet => CONST_RADIUS, &
      PI => CONST_PI

    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use scale_cubedsphere_coord_cnv, only: CubedSphereCoordCnv_LonLat2CSVec
  
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      MOMX_p  => PHYTEND_MOMX_ID, &
      MOMY_p  => PHYTEND_MOMY_ID, &
      MOMZ_p  => PHYTEND_MOMZ_ID, &
      RHOH_p  => PHYTEND_RHOH_ID

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars

    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(ElementBase3D), pointer :: elem
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot

    real(RP), parameter :: rtau = 1.0_RP / 1800.0_RP

    integer :: n
    integer :: ke

    real(RP), allocatable :: DENS(:)
    real(RP), allocatable :: sfac(:)
    real(RP), allocatable :: UmetOvCosLat(:,:), Vmet(:,:)
    real(RP), allocatable :: U0(:,:), V0(:,:)
    !------------------------------------------

    if ( this%USER_do ) then
      call atm%vars%Calc_diagVar( 'PT_diff', PT_diff )
      call FILE_HISTORY_meshfield_in( PT_diff, "perturbation of potential temperature" )
    end if

    return
    
    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )     
      
      call AtmosVars_GetLocalMeshPhyAuxVars( n,  atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                          )
      
      elem => lcmesh%refElem3D
      allocate( DENS(elem%Np), sfac(elem%Np) )
      allocate( U0(elem%Np,lcmesh%Ne), V0(elem%Np,lcmesh%Ne) )
      allocate( UmetOvCosLat(elem%Np,lcmesh%Ne), Vmet(elem%Np,lcmesh%Ne) )

      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeE
        UmetOvCosLat(:,ke) = 40.0_RP
        Vmet(:,ke)         = 0.0_RP
      end do
      
      call CubedSphereCoordCnv_LonLat2CSVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), elem%Np * lcmesh%Ne, RPlanet, &
        UmetOvCosLat(:,:), Vmet(:,:),                                                             &
        U0(:,:), V0(:,:)                                                                          )

      !$omp parallel do private(DENS, sfac)
      do ke=lcmesh%NeS, lcmesh%NeE
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        sfac(:) = rtau * 0.5_RP * ( 1.0_RP - cos( PI * ( lcmesh%pos_en(:,ke,3) - 25.E3_RP ) / (40.E3_RP - 25.E3_RP) ) ) 

        where ( lcmesh%pos_en(:,ke,3) > 25.E3_RP ) 
          atm%vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMX_p)%local(n)%val(:,ke)   &
            - sfac(:) * ( MOMX%val(:,ke) - DENS(:) * U0(:,ke) )
          atm%vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMY_p)%local(n)%val(:,ke)   &
            - sfac(:) * ( MOMY%val(:,ke) - DENS(:) * V0(:,ke) )
          atm%vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(MOMZ_p)%local(n)%val(:,ke)   &
            - sfac(:) * MOMZ%val(:,ke)

          atm%vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke) = atm%vars%PHY_TEND(RHOH_p)%local(n)%val(:,ke)   &
            - DENS(:) * sfac(:) * CpDry * ( PRES%val(:,ke) / DENS(:) - PRES_hyd%val(:,ke) / DENS_hyd%val(:,ke) ) / Rdry
        end where
      end do
      deallocate( DENS, sfac )
      deallocate( UmetOvCosLat, Vmet, U0, V0 )
    end do

    return
  end subroutine USER_calc_tendency

  !------

!OCL SERIAL  
  subroutine exp_SetInitCond_mountain_wave( this,                 &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
     
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      OHM => CONST_OHM,      &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_RADIUS
    
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constBVFreq
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
  
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
    
    real(RP) :: THETA0 = 300.0_RP
    real(RP) :: BruntVaisalaFreq = 0.0187_RP  
    real(RP) :: U0     = 40.0_RP
    namelist /PARAM_EXP/ &
      U0,                &
      BruntVaisalaFreq
    integer, parameter :: IntrpPolyOrder_h = 8
    integer, parameter :: IntrpPolyOrder_v = 8

    real(RP) :: p_lat(elem%Nfp_v,lcmesh%lcmesh2D%Ne)
    real(RP) :: G_(elem%Nfp_v), G0(elem%Nfp_v), lat0(elem%Nfp_v)

    real(RP) :: cos_lat(elem%Np)
    real(RP) :: T(elem%Np)
    real(RP) :: MOMX_met_ov_coslat(elem%Np,lcmesh%Ne)
    real(RP) :: MOMY_met          (elem%Np,lcmesh%Ne)

    integer :: ke, ke2d
    integer :: ierr

    type(LocalMesh2D), pointer :: lmesh2D
    real(RP) :: H0_pres
    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MOUNTAIN_WAVE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MOUNTAIN_WAVE_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)
    !---

    call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
      BruntVaisalaFreq, THETA0, PRES00,                               & ! (in)
      x, y, lcmesh%zlev, lcmesh, elem                                 ) ! (in)
    
    !---

    H0_pres = Grav * Rdry / CPdry / BruntVaisalaFreq**2

    lmesh2D => lcmesh%lcmesh2D
    !$omp parallel do private(G0, G_, lat0)
    do ke2D=lmesh2D%NeS, lmesh2D%NeE
      lat0(:) = 0.0_RP
      call calc_G( G0(:), &
        U0, OHM, RPlanet, Grav, H0_pres, lat0, elem%Nfp_v )
      call calc_G( G_(:), &
        U0, OHM, RPlanet, Grav, H0_pres, lcmesh%lat2D(:,ke2D), elem%Nfp_v )

      p_lat(:,ke2d) = PRES00 * ( G0(:) / G_(:) )**( RPlanet * 0.25_RP / H0_pres )
    end do

    !$omp parallel do private( ke2d, cos_lat, T )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2d = lcmesh%EMap3Dto2D(ke)
      cos_lat(:) = cos(lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D))

      PRES_hyd(:,ke) = p_lat(elem%IndexH2Dto3D,ke2d) * exp( - lcmesh%zlev(:,ke) / H0_pres )

      T(:) = H0_pres / Rdry * ( Grav - U0 * cos_lat * ( U0 * cos_lat / RPlanet + 2.0_RP * OHM * cos_lat ) )
      DENS_hyd(:,ke) = p_lat(elem%IndexH2Dto3D,ke2d) / ( Rdry * T(:) ) * exp( - lcmesh%zlev(:,ke) / H0_pres )

      MOMX_met_ov_coslat(:,ke) = DENS_hyd(:,ke) * U0
      MOMY_met(:,ke)           = 0.0_RP
    end do

    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), elem%Np * lcmesh%Ne, RPlanet, &
      MOMX_met_ov_coslat(:,:), MOMY_met(:,:),                                                   &
      MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE)                              )

    return
  contains
    subroutine calc_G( G, U, OMG, a, Grav, H0p, lat, np )
      implicit none
      integer, intent(in) :: np
      real(RP), intent(out) :: G(np)
      real(RP), intent(in) :: U
      real(RP), intent(in) :: OMG
      real(RP), intent(in) :: a
      real(RP), intent(in) :: Grav
      real(RP), intent(in) :: H0p
      real(RP), intent(in) :: lat(np)

      real(RP) :: cos_2lat(np)
      real(RP) :: cos_4lat(np)
      real(RP) :: g1(np), g2
      !----------------

      cos_2lat(:) = cos(2.0_RP * lat(:))
      cos_4lat(:) = cos(4.0_RP * lat(:))

      g1(:) = ( 3.0_RP + 4.0_RP * cos_2lat(:) + cos_4lat(:) ) * 2.0_RP * U**2 * ( U**2 + 4.0_RP * a * OMG * ( U + a * OMG ) ) &
            - ( 1.0_RP + cos_2lat(:) ) * 16.0_RP * U * Grav * a * ( U + 2.0_RP * a * OMG )                                             &
            + 16.0_RP * ( a * Grav )**2 
      g2    = U**2 * ( U**2 + 4.0_RP * a * OMG * ( U + a * OMG ) )

      G(:) = g1(:) / g2
  
      return
    end subroutine 
  end subroutine exp_SetInitCond_mountain_wave

end module mod_user
