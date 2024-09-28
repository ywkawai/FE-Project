!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of two-dimensional linear advection in global domain
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
    PI  => CONST_PI,        &
    P00 => CONST_PRE00,     &
    Grav => CONST_GRAV,     &
    Rdry => CONST_Rdry,     &
    RPlanet => CONST_RADIUS, &
    EPS => CONST_EPS
  use scale_tracer, only: &
    TRACER_inq_id

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D   
  use scale_mesh_base, only: MeshBase
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_meshfield_base, only: MeshField3D
  use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_in

  use scale_cubedsphere_coord_cnv, only: &
    CubedSphereCoordCnv_LonLat2CSVec
  
  use mod_atmos_component, only: &
    AtmosComponent

  use scale_meshfield_analysis_numerror, only: MeshFieldAnalysisNumerror3D

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
    procedure :: final => User_final
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

  ! Surface state
  real(RP)               :: SFC_PRES           ! surface pressure [Pa]
  ! Environment state
  real(RP)               :: ENV_TEMP = 300.0_RP ! potential temperature of environment [K]


  character(len=H_MID) :: FLOW_TYPE        = 'SOLID_BODY_ROTATION_FLOW' ! DEFORMATION_FLOW
  real(RP) :: SOLID_BODY_ROT_TAU  = 86400.0_RP * 12.0_RP
  real(RP) :: SOLID_BODY_ROT_ALPH = 0.0_RP

  character(len=H_MID) :: INIT_TRACER_PROF = 'GAUSSIAN'            ! DCMIP2012_q1
  real(RP) :: GAUSSIAN_rh
  real(RP) :: lonc   = 0.0_RP
  real(RP) :: latc   = 0.0_RP
  integer :: InitGP_PolyOrder_h = 11
  integer :: InitGP_PolyOrder_v = 11


  type(MeshFieldAnalysisNumerror3D) :: numerror_analysis
  integer :: numerror_vid(1)

  type(MeshField3D) :: PTracer_exact
  type(MeshField3D) :: lon, lat

  real(RP) :: lonc_now, latc_now

  !-----------------------------------------------------------------------------
contains
  subroutine USER_mkinit ( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('tracer_advection_global')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_tracer_advection )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()
    
    return
  end subroutine USER_mkinit

  subroutine USER_setup( this, atm )
    use scale_tracer, only: &
       TRACER_regist    
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    integer :: PolyOrderErrorCheck
    character(len=H_MID) :: LOG_OUT_BASENAME = 'LOG_NUMERROR'
    integer :: LOG_STEP_INTERVAL

    logical :: USER_do        = .false. !< do user step?
    namelist / PARAM_USER / &
       USER_do,              &
       polyOrderErrorCheck,  &
       LOG_OUT_BASENAME,     &
       LOG_STEP_INTERVAL    


    class(MeshBase), pointer :: ptr_mesh
    class(MeshCubedSphereDom3D), pointer :: mesh

    class(ElementBase3D), pointer :: refElem
    type(HexahedralElement) :: refElem3D      

    integer :: ierr    
    integer :: iq        
    !------------------------------------------

    PolyOrderErrorCheck = 6
    LOG_STEP_INTERVAL   = 5

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

    !--- read namelist
    call read_exp_params()

    !-
    call this%UserBase%Setup( atm, USER_do )

    !-
    call atm%mesh_gm%GetModelMesh( ptr_mesh )
    select type( ptr_mesh )
    type is (MeshCubedSphereDom3D)
      mesh => ptr_mesh
    end select

    !-
    call refElem3D%Init( mesh%refElem3D%PolyOrder_h, mesh%refElem3D%PolyOrder_v, .false. )
    call numerror_analysis%Init( &
      PolyOrderErrorCheck, LOG_OUT_BASENAME, LOG_STEP_INTERVAL, refElem3D )
    call refElem3D%Final()

    call numerror_analysis%Regist( "PTracer", "1", numerror_vid(1) )

    !-
    call TRACER_REGIST( iq,                   & ! [OUT]
                        1,                    & ! [IN]
                        (/'PTracer'/),        & ! [IN]
                        (/'Passive tracer'/), & ! [IN]
                        (/'1'/)               ) ! [IN]
        
    !--
    call PTracer_exact%Init( "PTracer_exact", "kg/m3", mesh )
    call lon%Init( "lon", "rad", mesh )
    call lat%Init( "lat", "rad", mesh )

    return
  end subroutine USER_setup

!OCL SERIAL
  subroutine USER_final( this )
    implicit none
    class(User), intent(inout) :: this
    !------------------------------------------
    
    call numerror_analysis%Final()
    return
  end subroutine USER_final

!OCL SERIAL  
  subroutine USER_calc_tendency( this, atm )
    use scale_time_manager, only:  TIME_NOWSTEP
    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use scale_geographic_coord_cnv, only: &
      GeographicCoordCnv_orth_to_geo_pos, &
      GeographicCoordCnv_geo_to_orth_pos, &
      GeographicCoordCnv_rotateY,         &
      GeographicCoordCnv_rotateZ    
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars   
    use scale_prc 
    implicit none

    class(User), intent(inout) :: this 
    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    type(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D
    integer :: iq
    real(RP) :: time

    real(RP), allocatable :: svec(:,:,:)
    real(RP), allocatable :: W(:)
    real(RP), allocatable :: lon3D(:)
    real(RP), allocatable :: lat3D(:)  

    class(MeshBase), pointer :: ptr_mesh
    class(MeshCubedSphereDom3D), pointer :: mesh
    real(RP) :: tsec


    real(RP) :: center_pos_geo(3,1), center_pos_orth(3,1)
    real(RP) :: center_pos_tmp(3,1), center_pos_onXY(3)    
    !------------------------------------------

    call this%UserBase%calc_tendency( atm )

    call atm%mesh_gm%GetModelMesh( ptr_mesh )
    select type( ptr_mesh )
    type is (MeshCubedSphereDom3D)
      mesh => ptr_mesh
    end select
    tsec = atm%time_manager%dtsec * real( TIME_NOWSTEP - 1, kind=RP )
    !if (PRC_myrank==0) write(*,*) "time=", time
    
    !--
    call TRACER_inq_id( "PTracer", iq )

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      
      elem3D => lcmesh%refElem3D

      !$omp parallel do private(ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)
        lon%local(n)%val(:,ke) = lcmesh%lon2D(elem3D%IndexH2Dto3D(:),ke2D)
        lat%local(n)%val(:,ke) = lcmesh%lat2D(elem3D%IndexH2Dto3D(:),ke2D)
      end do

      select case( FLOW_TYPE )
      case( 'SOLID_BODY_ROTATION_FLOW' )
        cycle
      end select  

      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      
      elem3D => lcmesh%refElem3D

      allocate( svec(elem3D%Np,lcmesh%Ne,2), W(elem3D%Np) )
      allocate( lon3D(elem3D%Np), lat3D(elem3D%Np) )

      !$omp parallel do private( ke2D, W, lon3D, lat3D )
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)
        lon3D(:) = lcmesh%lon2D(elem3D%IndexH2Dto3D(:),ke2D)
        lat3D(:) = lcmesh%lat2D(elem3D%IndexH2Dto3D(:),ke2D)
  
        select case( FLOW_TYPE )
        case( 'DEFORMATION_FLOW' )
          call calc_DCMIP2012_deformation_flow( &
            svec(:,ke,1), svec(:,ke,2), W(:),        &
            lon3D(:), lat3D(:), lcmesh%zlev(:,ke),   &
            time, elem3D%Np                          )  
        end select          
        svec    (:,ke,1) = DENS_hyd%val(:,ke) * svec(:,ke,1)
        svec    (:,ke,2) = DENS_hyd%val(:,ke) * svec(:,ke,2)
        MOMZ%val(:,ke  ) = DENS_hyd%val(:,ke) * W(:)
      end do

      call CubedSphereCoordCnv_LonLat2CSVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),    & ! (in)
        lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), elem3D%Np * lcmesh%Ne,    & ! (in)
        svec(:,:,1), svec(:,:,2),                                      & ! (in)
        MOMX%val(:,lcmesh%NeS:lcmesh%NeE), MOMY%val(:,lcmesh%NeS:lcmesh%NeE) ) ! (out)
        
      deallocate( svec )
      deallocate( lon3D, lat3D )
    end do

    select case( FLOW_TYPE )
    case( 'SOLID_BODY_ROTATION_FLOW' )
      center_pos_geo(:,1) = (/ lonc, latc, RPlanet /)
      call GeographicCoordCnv_geo_to_orth_pos( center_pos_geo, 1, & ! (in) 
        center_pos_orth ) ! (out)

      call GeographicCoordCnv_rotateY( center_pos_orth(:,1), - SOLID_BODY_ROT_ALPH, & ! (in) 
        center_pos_tmp(:,1)   ) ! (out)
      call GeographicCoordCnv_rotateZ( center_pos_tmp(:,1), 2.0_RP * PI / SOLID_BODY_ROT_TAU * tsec, & ! (in) 
        center_pos_onXY       ) ! (out)

      call GeographicCoordCnv_rotateY( center_pos_onXY, + SOLID_BODY_ROT_ALPH, & ! (in)
        center_pos_tmp(:,1) ) ! (out)
      call GeographicCoordCnv_orth_to_geo_pos( center_pos_tmp, 1, & ! (in) 
        center_pos_geo ) ! (out)

      lonc_now = center_pos_geo(1,1); latc_now = center_pos_geo(2,1)
      LOG_INFO("USER_calc_tendency",*) "time=", tsec, "pos=", lonc_now, latc_now
      
    end select

    !--- Evaulate  numerical error

    call PROF_rapstart( 'USER_calc_tendency_numanal_eval', 2)
    call numerror_analysis%Evaluate( &
      TIME_NOWSTEP, tsec, mesh, evaluate_error_lc )
    call PROF_rapend( 'USER_calc_tendency_numanal_eval', 2)
    
    call FILE_HISTORY_meshfield_in( PTracer_exact, "PTracer_exact" )
    call FILE_HISTORY_meshfield_in( lon, "lon" )
    call FILE_HISTORY_meshfield_in( lat, "lat" )

    return

  contains 
    subroutine evaluate_error_lc( this, q, qexact, qexact_intrp, lcmesh, elem3D, intrp_epos )
      use scale_localmeshfield_base, only: LocalMeshFieldBase
      use mod_atmos_vars, only: &
        AtmosVars_GetLocalMeshPrgVars,    &
        AtmosVars_GetLocalMeshPhyAuxVars

      use mod_mkinit_util, only: &
        mkinitutil_calc_cosinebell_global, &
        mkinitutil_GalerkinProjection_global

      implicit none
      class(MeshFieldAnalysisNumerror3D), intent(in) :: this
      class(LocalMesh3D), intent(in) :: lcmesh
      class(ElementBase3D) :: elem3D
      real(RP), intent(out) :: q(elem3D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact(elem3D%Np,lcmesh%Ne,this%var_num)
      real(RP), intent(out) :: qexact_intrp(this%intrp_np,lcmesh%Ne,this%var_num)
      real(RP), intent(in) :: intrp_epos(this%intrp_np,3)

      class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
      class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
      class(LocalMeshFieldBase), pointer :: PRES, PT
      class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
  
      integer :: n
      integer :: ke
      real(RP) :: DENS(elem3D%Np)
      !---------------------------------------------

      call PROF_rapstart( 'USER_calc_tendency_get_exactsol', 2)

      n = lcmesh%lcdomID

      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot                   )      
      
      call AtmosVars_GetLocalMeshPhyAuxVars( n,  atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                          )
      
      call evaluate_error_lc_core( q(:,:,:), qexact(:,:,:), qexact_intrp(:,:,:), &
        PTracer_exact%local(n)%val,                                                 &
        tsec, atm%vars%QTRC_VARS(iq)%local(n)%val,                                  &
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                         &
        DENS_hyd%val, PRES_hyd%val,                                                 &
        mesh%zmin_gl, mesh%zmax_gl,                                                 &
        lcmesh, elem3D, intrp_epos, numerror_analysis%intrp_np                      )
            
      call PROF_rapend( 'USER_calc_tendency_get_exactsol', 2)

      return
    end subroutine evaluate_error_lc

  end subroutine USER_calc_tendency

!OCL SERIAL
  subroutine evaluate_error_lc_core( q, qexact, qexact_intrp, PTracer_exact, &
    tsec, PTracer, DDENS, MOMX, MOMY, MOMZ, DRHOT, DENS_hyd, PRES_hyd, &
    dom_zmin, dom_zmax,                                                &
    lcmesh, elem3D, epos_intrp, intrp_np )

    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    integer, intent(in) :: intrp_np
    real(RP), intent(out) :: q(elem3D%Np,lcmesh%Ne,numerror_analysis%var_num)
    real(RP), intent(out) :: qexact(elem3D%Np,lcmesh%Ne,numerror_analysis%var_num)
    real(RP), intent(out) :: qexact_intrp(numerror_analysis%intrp_np,lcmesh%Ne,numerror_analysis%var_num)
    real(RP), intent(out) :: PTracer_exact(elem3D%Np,lcmesh%NeA)
    real(DP), intent(in) :: tsec
    real(RP), intent(in) :: PTracer(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: dom_zmin, dom_zmax
    real(RP), intent(in) :: epos_intrp(intrp_np,3)

    integer :: ke
    !------------------------------------------------------------------------

    call PROF_rapstart( 'USER_calc_exactsol_core_1', 2)

    call calc_tracer_profile_now( qexact_intrp(:,:,1),    & ! (out)
      INIT_TRACER_PROF, tsec,                                   & ! (in)
      lcmesh, epos_intrp(:,1),epos_intrp(:,2), epos_intrp(:,3), & ! (in)
      intrp_np, elem3D%Nv                                       ) ! (in)
  
    call calc_tracer_profile_now( qexact(:,:,1), & ! (out)
      INIT_TRACER_PROF, tsec,                    & ! (in)
      lcmesh, elem3D%x1, elem3D%x2, elem3D%x3,   & ! (in)
      elem3D%Np, elem3D%Nv                       ) ! (in)

    !$omp parallel do
    do ke=lcmesh%NeS, lcmesh%NeE
      q(:,ke,1) = PTracer(:,ke)
      PTracer_exact(:,ke) = qexact(:,ke,1)
    end do

    call PROF_rapend( 'USER_calc_exactsol_core_1', 2)

  end subroutine evaluate_error_lc_core

  !------

  subroutine exp_SetInitCond_tracer_advection( this,                       &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )

    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT
    use mod_experiment, only: &
      TracerLocalMeshField_ptr
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2CartPos, &
      CubedSphereCoordCnv_Cart2CSVec
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

    integer :: ke
    integer :: ke2D
    integer :: iq

    real(RP) :: svec(elem%Np,lcmesh%Ne,2)
    real(RP) :: W(elem%Np)
    real(RP) :: lon3D(elem%Np), lat3D(elem%Np)

    integer :: p

    real(RP) :: Cart_pos(elem%Np,lcmesh%Ne,3)
    real(RP) :: Cart_vec(elem%Np,lcmesh%Ne,3)
    !-----------------------------------------------------------------------------

    call read_exp_params()

    !----
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd,                            &
      ENV_TEMP, SFC_PRES, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3), &
      lcmesh, elem )

    !---

    call TRACER_inq_id( "PTracer", iq )
    call calc_tracer_profile_GP( tracer_field_list(iq)%ptr%val, & ! (out)
      INIT_TRACER_PROF, lcmesh, elem,                           & ! (in) 
      InitGP_PolyOrder_h, InitGP_PolyOrder_v                    )
  
    call CubedSphereCoordCnv_CS2CartPos( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), & ! (in)
      elem%Np * lcmesh%Ne,                              & ! (in)
      Cart_pos(:,:,1), Cart_pos(:,:,2), Cart_pos(:,:,3) ) ! (out)
    
    !$omp parallel do private(ke2D, lon3D, lat3D, W)
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)
      lon3D(:) = lcmesh%lon2D(elem%IndexH2Dto3D(:),ke2D)
      lat3D(:) = lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D)

      select case (FLOW_TYPE) 
      case( 'SOLID_BODY_ROTATION_FLOW' )
        ! call calc_solid_body_rot_flow( &
        !   svec(:,ke,1), svec(:,ke,2), W(:),        &
        !   lon3D(:), lat3D(:), lcmesh%zlev(:,ke),   &
        !   0.0_RP, elem%Np                          )

        call calc_solid_body_rot_flow_cartvec( &
          Cart_vec(:,ke,1), Cart_vec(:,ke,2), Cart_vec(:,ke,3),                 & ! (out)
          Cart_pos(:,ke,1), Cart_pos(:,ke,2), Cart_pos(:,ke,3), 0.0_RP, elem%Np ) ! (in)
        
        Cart_vec(:,ke,1) = DENS_hyd(:,ke) * Cart_vec(:,ke,1)
        Cart_vec(:,ke,2) = DENS_hyd(:,ke) * Cart_vec(:,ke,2)
        Cart_vec(:,ke,3) = DENS_hyd(:,ke) * Cart_vec(:,ke,3)
        W(:) = 0.0_RP 
      case( 'DEFORMATION_FLOW' )
        call calc_DCMIP2012_deformation_flow( &
          svec(:,ke,1), svec(:,ke,2), W(:),        &
          lon3D(:), lat3D(:), lcmesh%zlev(:,ke),   &
          0.0_RP, elem%Np                          )

        svec(:,ke,1) = DENS_hyd(:,ke) * svec(:,ke,1)
        svec(:,ke,2) = DENS_hyd(:,ke) * svec(:,ke,2)
      end select

      MOMZ(:,ke) = DENS_hyd(:,ke) * W(:)
    end do
        
    select case (FLOW_TYPE) 
    case( 'SOLID_BODY_ROTATION_FLOW' )
      
      call CubedSphereCoordCnv_Cart2CSVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),    & ! (in)
        lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), elem%Np * lcmesh%Ne,      & ! (in)
        Cart_vec(:,:,1), Cart_vec(:,:,2), Cart_vec(:,:,3),             & ! (in)
        MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE) )   ! (out)

    case( 'DEFORMATION_FLOW' )

      call CubedSphereCoordCnv_LonLat2CSVec( &
        lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),    & ! (in)
        lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), elem%Np * lcmesh%Ne,      & ! (in)
        svec(:,:,1), svec(:,:,2),                                      & ! (in)
        MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE) )   ! (out)
      
    end select


    return
  end subroutine exp_SetInitCond_tracer_advection

!- private --

!OCL SERIAL
  subroutine read_exp_params()
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      CVdry => CONST_CVdry,  &
      PRES00 => CONST_PRE00, &
      Pstd   => CONST_Pstd
    implicit none
   
    namelist /PARAM_EXP/ &
      SFC_PRES,            &
      ENV_TEMP,            &
      FLOW_TYPE,           &
      SOLID_BODY_ROT_TAU,  &
      SOLID_BODY_ROT_ALPH, &
      INIT_TRACER_PROF,    &
      GAUSSIAN_rh,         &
      lonc, latc,          &
      InitGP_PolyOrder_h,  &
      InitGP_PolyOrder_v

    integer :: ierr
    !-----------------------------------------------------------------------------

    SFC_PRES = Pstd
    ENV_TEMP = 300.0_RP      

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("TRACER_ADVECTION_global_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("TRACER_ADVECTION_global_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)
  
    lonc_now = lonc
    latc_now = latc

    return
  end subroutine read_exp_params

!OCL SERIAL
  subroutine calc_tracer_profile_GP( qtrc, & ! (inout) 
    trc_prof_typename, lcmesh, elem,       &
    GP_PolyOrder_h, GP_PolyOrder_v         ) ! (in)

    use scale_const, only: &
      RPlanet => CONST_RADIUS    
    use mod_mkinit_util, only: &
      mkinitutil_calc_cosinebell_global, &
      mkinitutil_GalerkinProjection_global
        
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: qtrc(elem%Np,lcmesh%NeA)
    character(len=*), intent(in) :: trc_prof_typename
    integer, intent(in), optional :: GP_PolyOrder_h
    integer, intent(in), optional :: GP_PolyOrder_v
    !---------------------------------------------

    select case( trim(trc_prof_typename) )
    case( 'GAUSSIAN' )
      call mkinitutil_GalerkinProjection_global( qtrc,             &  ! (out)
        calc_gaussian, GP_PolyOrder_h, GP_PolyOrder_v,             &  ! (in)
        lcmesh, elem, RPlanet                                      )  ! (in)
    case( 'COSBELL' )
      call mkinitutil_GalerkinProjection_global( qtrc,             &  ! (out)
        calc_cosbell, GP_PolyOrder_h, GP_PolyOrder_v,             &  ! (in)
        lcmesh, elem, RPlanet                                      )  ! (in)        
    case( 'DCMIP2012_q1' )
      call mkinitutil_GalerkinProjection_global( qtrc,              &  ! (out)
        calc_DCMIP2012_tracer_q1, GP_PolyOrder_h, GP_PolyOrder_v,   &  ! (in)
        lcmesh, elem, RPlanet                                       )  ! (in)
    case default
      LOG_ERROR("calc_tracer_profile",*) trim(trc_prof_typename)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine calc_tracer_profile_GP

!OCL SERIAL
  subroutine calc_tracer_profile_now( qtrc,   & ! (inout) 
    trc_prof_typename, tsec,                  & ! (in)
    lcmesh3D, elem_x1, elem_x2, elem_x3, Np, Nv ) ! (in)

    use scale_const, only: &
      RPlanet => CONST_RADIUS
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2LonLatPos
    use mod_mkinit_util, only: &
      mkinitutil_calc_cosinebell_global, &
      mkinitutil_GalerkinProjection_global
        
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    integer, intent(in) :: Np, Nv
    real(RP), intent(out) :: qtrc(Np,lcmesh3D%NeA)
    character(len=*), intent(in) :: trc_prof_typename
    real(RP), intent(in) :: tsec
    real(RP), intent(in) :: elem_x1(Np)
    real(RP), intent(in) :: elem_x2(Np)
    real(RP), intent(in) :: elem_x3(Np)

    real(RP) :: vx(Nv), vy(Nv), vz(Nv)
    real(RP) :: alpha(Np,lcmesh3D%Ne), beta(Np,lcmesh3D%Ne), z(Np,lcmesh3D%Ne)
    real(RP) :: gam(Np,lcmesh3D%Ne)
    real(RP) :: lon(Np,lcmesh3D%Ne), lat(Np,lcmesh3D%Ne)

    integer :: ke
    !---------------------------------------------

    !$omp parallel do private(vx, vy, vz)
    do ke=lcmesh3D%NeS, lcmesh3D%NeE
      vx(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),1)
      vy(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),2)
      vz(:) = lcmesh3D%pos_ev(lcmesh3D%EToV(ke,:),3)
      alpha(:,ke) = vx(1) + 0.5_RP * ( elem_x1(:) + 1.0_RP ) * ( vx(2) - vx(1) ) 
      beta(:,ke) = vy(1) + 0.5_RP * ( elem_x2(:) + 1.0_RP ) * ( vy(4) - vy(1) )
      z(:,ke) = vz(1) + 0.5_RP * ( elem_x3(:) + 1.0_RP ) * ( vz(5) - vz(1) )
      gam(:,ke) = 1.0_RP ! + z(:,ke) /  RPlanet
    end do

    call CubedSphereCoordCnv_CS2LonLatPos( &
      lcmesh3D%panelID, alpha, beta, gam, Np * lcmesh3D%Ne, & ! (in)
      lon(:,:), lat(:,:) ) ! (out)
    
    select case( trim(trc_prof_typename) )
    case( 'GAUSSIAN' )
      !$omp parallel do
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        call calc_gaussian_core( qtrc(:,ke), &
          lon(:,ke), lat(:,ke), z(:,ke), Np, RPlanet )
      end do
    case( 'COSBELL' )
      !$omp parallel do
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        call calc_cosbell_core( qtrc(:,ke), &
          lon(:,ke), lat(:,ke), z(:,ke), Np, RPlanet )
      end do      
    case( 'DCMIP2012_q1' )
      !$omp parallel do
      do ke=lcmesh3D%NeS, lcmesh3D%NeE
        call calc_DCMIP2012_tracer_q1_core( qtrc(:,ke), &
          lon(:,ke), lat(:,ke), z(:,ke), Np, RPlanet )
      end do      
    case default
      LOG_ERROR("calc_tracer_profile",*) trim(trc_prof_typename)//' is not supported. Check!'
      call PRC_abort
    end select

    return
  end subroutine calc_tracer_profile_now

!OCL SERIAL
  subroutine calc_gaussian(  qtrc_intrp,   &
    lon_intrp, lat_intrp, z_intrp, elem_intrp, RPlanet )
    implicit none
    class(ElementBase3D), intent(in) :: elem_intrp
    real(RP), intent(out) :: qtrc_intrp(elem_intrp%Np)
    real(RP), intent(in) :: lon_intrp(elem_intrp%Np)
    real(RP), intent(in) :: lat_intrp(elem_intrp%Np)
    real(RP), intent(in) :: z_intrp(elem_intrp%Np)
    real(RP), intent(in) :: RPlanet

    !---------------------------------------------
    call calc_gaussian_core(  qtrc_intrp,   &
      lon_intrp, lat_intrp, z_intrp, elem_intrp%Np, RPlanet )

    return
  end subroutine calc_gaussian

!OCL SERIAL
  subroutine calc_gaussian_core(  qtrc_intrp,   &
    lon_intrp, lat_intrp, z_intrp, Np, RPlanet )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: qtrc_intrp(Np)
    real(RP), intent(in) :: lon_intrp(Np)
    real(RP), intent(in) :: lat_intrp(Np)
    real(RP), intent(in) :: z_intrp(Np)
    real(RP), intent(in) :: RPlanet

    real(RP) :: r_intrp(Np)
    !---------------------------------------------

    r_intrp(:) = RPlanet / GAUSSIAN_rh * acos( sin(latc_now) * sin(lat_intrp(:)) + cos(latc_now) * cos(lat_intrp(:)) * cos(lon_intrp(:) - lonc_now) )
    qtrc_intrp(:) = exp( - r_intrp(:)**2 )
    return
  end subroutine calc_gaussian_core

!OCL SERIAL
  subroutine calc_cosbell(  qtrc_intrp,   &
    lon_intrp, lat_intrp, z_intrp, elem_intrp, RPlanet )
    implicit none
    class(ElementBase3D), intent(in) :: elem_intrp
    real(RP), intent(out) :: qtrc_intrp(elem_intrp%Np)
    real(RP), intent(in) :: lon_intrp(elem_intrp%Np)
    real(RP), intent(in) :: lat_intrp(elem_intrp%Np)
    real(RP), intent(in) :: z_intrp(elem_intrp%Np)
    real(RP), intent(in) :: RPlanet

    !---------------------------------------------
    call calc_cosbell_core(  qtrc_intrp,   &
      lon_intrp, lat_intrp, z_intrp, elem_intrp%Np, RPlanet )

    return
  end subroutine calc_cosbell

!OCL SERIAL
  subroutine calc_cosbell_core(  qtrc_intrp,   &
    lon_intrp, lat_intrp, z_intrp, Np, RPlanet )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: qtrc_intrp(Np)
    real(RP), intent(in) :: lon_intrp(Np)
    real(RP), intent(in) :: lat_intrp(Np)
    real(RP), intent(in) :: z_intrp(Np)
    real(RP), intent(in) :: RPlanet

    real(RP) :: r_intrp(Np)
    !---------------------------------------------

    r_intrp(:) = RPlanet / GAUSSIAN_rh * acos( sin(latc_now) * sin(lat_intrp(:)) + cos(latc) * cos(lat_intrp(:)) * cos(lon_intrp(:) - lonc_now) )
    where( r_intrp(:) <= 1.0_RP ) 
      qtrc_intrp(:) = ( 0.5_RP * (1.0_RP + cos( PI * r_intrp(:) ) ) )
    elsewhere
      qtrc_intrp(:) = 0.0_RP
    end where

    return
  end subroutine calc_cosbell_core

!OCL SERIAL
  subroutine calc_DCMIP2012_tracer_q1( q,  &
    lon, lat, z, elem_intrp, RPlanet       )

    implicit none
    class(ElementBase3D), intent(in) :: elem_intrp
    real(RP), intent(out) :: q(elem_intrp%Np)
    real(RP), intent(in) :: lon(elem_intrp%Np)
    real(RP), intent(in) :: lat(elem_intrp%Np)
    real(RP), intent(in) :: z(elem_intrp%Np)
    real(RP), intent(in) :: RPlanet
    !----------------------------------------------

    call calc_DCMIP2012_tracer_q1_core( q,  &
      lon, lat, z, elem_intrp%Np, RPlanet   )

    return
  end subroutine calc_DCMIP2012_tracer_q1

!OCL SERIAL
  subroutine calc_DCMIP2012_tracer_q1_core( q,  &
    lon, lat, z, Np, RPlanet       )

    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: q(Np)
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: z(Np)
    real(RP), intent(in) :: RPlanet
    
    integer :: i
    real(RP) :: d(Np,2)
    real(RP) :: r(Np)

    real(RP) :: Rt                           !< Horizontal half-width of tracers
    real(RP), parameter :: Zt = 1000.0_RP    !< Vertical half-width of tracers
    real(RP) :: lon_c(2)                     !< Initial longitude of first and second tracers
    real(RP), parameter :: lat_c = 0.0_RP    !< Initial latitude of tracers
    real(RP), parameter :: z_c   = 5000.0_RP !< Initial altitude of tracers
    !----------------------------------------------

    Rt = 0.5_RP * RPlanet 
    lon_c(:) = (/ 5.0_RP, 7.0_RP /) * PI / 6.0_RP

    do i=1, 2
      r(:) = RPlanet * acos(sin(lat_c)*sin(lat(:)) + cos(lat_c)*cos(lat(:))*cos(lon(:)-lon_c(i)))
      d(:,i) = min( 1.0_RP, (r(:) / Rt)**2 + ((z(:) - z_c) / Zt)**2 )
    end do
    q(:) = 1.0_RP + 0.5_RP * ( cos(PI * d(:,1)) + cos(PI * d(:,2)) )

    return
  end subroutine calc_DCMIP2012_tracer_q1_core

!OCL SERIAL
  subroutine calc_solid_body_rot_flow( U, V, W, &
    lon, lat, z, time, Np )

    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: U(Np)
    real(RP), intent(out) :: V(Np)
    real(RP), intent(out) :: W(Np)
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: z(Np)
    real(RP), intent(in) :: time

    real(RP) :: U0
    !----------------------------------------------

    U0 = 2.0_RP * PI * RPlanet / SOLID_BODY_ROT_TAU

    U(:) =   U0 * ( cos(SOLID_BODY_ROT_ALPH) * cos(lat(:)) + sin(SOLID_BODY_ROT_ALPH) * cos(lon(:)) * sin(lat(:)) )
    V(:) = - U0 * sin(SOLID_BODY_ROT_ALPH) * sin(lon(:)) 
    W(:) = 0.0_RP

    return
  end subroutine calc_solid_body_rot_flow

!OCL SERIAL
  subroutine calc_solid_body_rot_flow_cartvec( Vec_x, Vec_y, Vec_z, &
    x, y, z, time, Np )

    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: Vec_x(Np)
    real(RP), intent(out) :: Vec_y(Np)
    real(RP), intent(out) :: Vec_z(Np)
    real(RP), intent(in) :: x(Np)
    real(RP), intent(in) :: y(Np)
    real(RP), intent(in) :: z(Np)
    real(RP), intent(in) :: time

    real(RP) :: U0
    real(RP) :: OMG_vec(3)
    !----------------------------------------------

    U0 = 2.0_RP * PI * RPlanet / SOLID_BODY_ROT_TAU
    OMG_vec(:) = 2.0_RP * PI / SOLID_BODY_ROT_TAU &
               * (/ sin(SOLID_BODY_ROT_ALPH), 0.0_RP, cos(SOLID_BODY_ROT_ALPH) /)

    Vec_x(:) = OMG_vec(2) * z(:) - y(:) * OMG_vec(3)
    Vec_y(:) = OMG_vec(3) * x(:) - z(:) * OMG_vec(1)
    Vec_z(:) = OMG_vec(1) * y(:) - x(:) * OMG_vec(2)

    return
  end subroutine calc_solid_body_rot_flow_cartvec

!OCL SERIAL
  subroutine calc_DCMIP2012_deformation_flow( U, V, W, &
    lon, lat, z, time, Np )

    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: U(Np)
    real(RP), intent(out) :: V(Np)
    real(RP), intent(out) :: W(Np)
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: z(Np)
    real(RP), intent(in) :: time

    real(RP) :: lon2(Np)
    real(RP) :: omg(Np)
    real(RP) :: pres(Np)
    real(RP) :: ud(Np), ua(Np)
    real(RP) :: vd(Np), va(Np)

    real(RP), parameter :: tau = 1036800.0_RP  !> Period of motion [sec]
    real(RP) :: OMG0                           !> Maximum of the vertical pressure velocity in units Pa/s
    real(RP), parameter :: b    = 0.2_RP       !> Normalized pressure depth of the divergent layer
    real(RP), parameter :: ptop = 254.944E2_RP
    real(RP), parameter :: T0   = 300.0_RP     !> Isothermal atmospheric temperature [K]
    !----------------------------------------------

    lon2(:) = lon(:) - 2.0_RP * PI * time / tau
    OMG0 = 23000.0_RP * PI / tau

    pres(:) = P00 * exp(- Grav * z(:) / ( Rdry * T0 ) )
    omg(:) = OMG0 * sin(lon2(:)) * cos(lat(:)) * cos(2.0_RP * PI * time / tau)             &
      * ( 1.0_RP + exp( (ptop - P00) / (b * ptop) ) - exp( (pres(:) - P00) / (b * ptop) )  &
                                                    - exp( (ptop - pres(:)) / (b * ptop) ) )

    ua(:) = 10.0_RP * RPlanet / tau * sin(lon2(:))**2 * sin(2.0_RP*lat(:)) * cos(PI * time / tau) &
          + 2.0_RP * PI * RPlanet / tau *  cos(lat(:))
    va(:) = 10.0_RP * RPlanet / tau * sin(2.0_RP * lon2(:)) * cos(lat(:)) * cos(PI * time / tau)

    ud(:) = OMG0 * RPlanet / (b * ptop) * cos(lon2(:)) * cos(lat(:))**2 * cos(2.0_RP * PI * time / tau) &
            * ( - exp( (pres(:) - P00) / (b * ptop) ) + exp( (ptop - pres(:)) / (b * ptop) )  )
    vd(:) = 0.0_RP

    U(:) = ua(:) + ud(:)
    V(:) = va(:) + vd(:)
    W(:) = - omg(:) / (Grav * pres(:) / ( Rdry * T0 ) )

    return
  end subroutine calc_DCMIP2012_deformation_flow

end module mod_user
