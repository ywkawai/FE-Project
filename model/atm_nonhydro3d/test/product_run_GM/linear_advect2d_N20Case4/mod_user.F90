!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of tracer advection in global domain
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
    RPlanet => CONST_RADIUS

  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D   
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_cubedsphere_coord_cnv, only: &
    CubedSphereCoordCnv_LonLat2CSVec
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D

  use scale_sparsemat, only: &
    SparseMat, sparsemat_matmul
  use mod_atmos_component, only: &
    AtmosComponent


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
  type(MeshField2D) :: lon2D, lat2D
  type(sparsemat), pointer :: Dx, Dy
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
    
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

  subroutine USER_setup( this, atm )
    use scale_tracer, only: &
       TRACER_regist    
    implicit none
    class(User), intent(inout) :: this    
    class(AtmosComponent), intent(inout) :: atm

    logical :: USER_do        = .false. !< do user step?
    namelist / PARAM_USER / &
       USER_do

    integer :: ierr    
    integer :: iq        

    integer :: n, ke, ke2D
    class(LocalMesh3D), pointer :: lcmesh3D
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase3D), pointer :: mesh3D
    class(ElementBase3D), pointer :: elem3D
    real(RP), allocatable :: psi(:), W(:)
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
    call TRACER_REGIST( iq,                   & ! [OUT]
                        1,                    & ! [IN]
                        (/'PTracer'/),        & ! [IN]
                        (/'Passive tracer'/), & ! [IN]
                        (/'1'/)               ) ! [IN]

    !-
    call atm%dyn_proc%trcadv_driver%Set_massflux_func( calc_deformation_flow )

    mesh3D => atm%mesh%ptr_mesh
    call atm%mesh%ptr_mesh%GetMesh2D(mesh2D)

    call lon2D%Init( "lon2D", "rad", mesh2D )
    call lat2D%Init( "lat2D", "rad", mesh2D )

    Dx => atm%mesh%DOptrMat(1)
    Dy => atm%mesh%DOptrMat(2)
    return
  end subroutine USER_setup

  subroutine USER_calc_tendency( this, atm )
    use scale_time_manager, only:  TIME_NOWSTEP
    use scale_localmeshfield_base, only: LocalMeshFieldBase
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars   
    use scale_prc 
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
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
    real(RP) :: time
    !------------------------------------------

    time = atm%time_manager%dtsec * real( TIME_NOWSTEP - 1, kind=RP )
    LOG_INFO("USER_calc_tendency",*) TIME_NOWSTEP, "time=", time
    
    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      
      elem3D => lcmesh%refElem3D

      call calc_deformation_flow( time, &
        MOMX%val, MOMY%val, MOMZ%val,   &
        DENS_hyd%val, lcmesh, elem3D )

      lon2D%local(n)%val(:,1:lcmesh%Ne2D) = lcmesh%lon2D(:,:)
      lat2D%local(n)%val(:,1:lcmesh%Ne2D) = lcmesh%lat2D(:,:)
    end do

    call FILE_HISTORY_meshfield_in( lon2D, "longitude" )
    call FILE_HISTORY_meshfield_in( lat2D, "latitude" )

    return
  end subroutine USER_calc_tendency

  subroutine calc_deformation_flow( time, &
    MFLX_x, MFLX_y, MFLX_z, DENS_hyd, lcmesh, elem )
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSVec
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(in) :: time
    real(RP), intent(out) :: MFLX_x(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MFLX_y(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MFLX_z(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)

    real(RP) :: svec(elem%Np,lcmesh%Ne,2), psi(elem%Np)
    real(RP) :: lon3D(elem%Np), lat3D(elem%Np)
    integer :: ke, ke2D
    real(RP) :: Fx(elem%Np), Fy(elem%Np)
    !-----------------------------------------------------------
  
    !$omp parallel do private( ke2D, lon3D, lat3D, Fx, Fy, psi )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)
      lon3D(:) = lcmesh%lon2D(elem%IndexH2Dto3D(:),ke2D)
      lat3D(:) = lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D)

      call DCMIP2012_deformation_flow( &
        svec(:,ke,1), svec(:,ke,2), psi(:),       &
        lon3D(:), lat3D(:),                &
        time, elem%Np                    )  
        
      svec    (:,ke,1) = DENS_hyd(:,ke) * svec(:,ke,1)
      svec    (:,ke,2) = DENS_hyd(:,ke) * svec(:,ke,2)

      ! call sparsemat_matmul( Dx, psi(:), Fx )
      ! call sparsemat_matmul( Dy, psi(:), Fy )
      ! MFLX_x(:,ke) = - DENS_hyd(:,ke) * lcmesh%Escale(:,ke,2,2) * Fy(:) / lcmesh%Gsqrt(:,ke)
      ! MFLX_y(:,ke) = + DENS_hyd(:,ke) * lcmesh%Escale(:,ke,1,1) * Fx(:) / lcmesh%Gsqrt(:,ke)

      MFLX_z(:,ke) = 0.0_RP
    end do

    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),        & ! (in)
      lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), lcmesh%Ne * elem%Np,          & ! (in)
      svec(:,:,1), svec(:,:,2),                                          & ! (in)
      MFLX_x(:,lcmesh%NeS:lcmesh%NeE), MFLX_y(:,lcmesh%NeS:lcmesh%NeE) ) ! (out)

    ! if ( lcmesh%panelID > 4 ) then
    !   !$omp parallel do
    !   do ke=lcmesh%NeS, lcmesh%NeE
    !     where ( sqrt(tan(lcmesh%pos_en(:,ke,1))**2 + tan(lcmesh%pos_en(:,ke,2))**2) < 1E-16 )
    !       MFLX_x(:,ke) = 0.0_RP; MFLX_y(:,ke) = 0.0_RP
    !     end where
    !   end do
    ! end if

    return
  end subroutine calc_deformation_flow

  !------

  subroutine exp_SetInitCond_tracer_advection( this,                       &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,        &
      GRAV => CONST_GRAV,    &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      CVdry => CONST_CVdry,  &
      PRES00 => CONST_PRE00, &
      Pstd   => CONST_Pstd
    use scale_tracer, only: &
      TRACER_inq_id

    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT
    use mod_experiment, only: &
      TracerLocalMeshField_ptr
  
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2CartPos

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
    
    ! Surface state
     real(RP)               :: SFC_PRES           ! surface pressure [Pa]
    ! Environment state
    real(RP)               :: ENV_TEMP = 300.0_RP ! potential temperature of environment [K]
    character(len=H_SHORT) :: qtrc_name = 'cosine_bells'
   
    namelist /PARAM_EXP/ &
      SFC_PRES, &
      ENV_TEMP, &
      qtrc_name 

    integer :: ke
    integer :: ke2D
    integer :: ierr
    integer :: iq

    real(RP) :: svec(elem%Np,lcmesh%Ne,2), psi(elem%Np)
    real(RP) :: lon3D(elem%Np), lat3D(elem%Np)
    real(RP) :: X_cartC(elem%Np,lcmesh%Ne), Y_cartC(elem%Np,lcmesh%Ne), Z_cartC(elem%Np,lcmesh%Ne)
    !-----------------------------------------------------------------------------

    SFC_PRES = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("TRACER_ADVECTION_global_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("TRACER_ADVECTION_global_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    
    ENV_TEMP = 300.0_RP      
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd,                            &
      ENV_TEMP, SFC_PRES, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3), &
      lcmesh, elem )
    PRES_hyd(:,:) = SFC_PRES
    DENS_hyd(:,:) = SFC_PRES / ( Rdry * ENV_TEMP )

    !---
    call CubedSphereCoordCnv_CS2CartPos( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), &
      lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), elem%Np * lcmesh%Ne,   &
      X_cartC, Y_cartC, Z_cartC )

    call TRACER_inq_id( "PTracer", iq )
    !$omp parallel do private(ke2D, lon3D, lat3D, psi)
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)
      lon3D(:) = lcmesh%lon2D(elem%IndexH2Dto3D(:),ke2D)
      lat3D(:) = lcmesh%lat2D(elem%IndexH2Dto3D(:),ke2D)

      call DCMIP2012_tracer( &
        tracer_field_list(iq)%ptr%val(:,ke),         &
        qtrc_name, lon3D(:), lat3D(:),               &
        X_cartC(:,ke), Y_cartC(:,ke), Z_cartC(:,ke), &
        elem%Np   )
      
      call DCMIP2012_deformation_flow( &
        svec(:,ke,1), svec(:,ke,2), psi(:), &
        lon3D(:), lat3D(:),                 &
        0.0_RP, elem%Np                )
      
      svec(:,ke,1) = DENS_hyd(:,ke) * svec(:,ke,1)
      svec(:,ke,2) = DENS_hyd(:,ke) * svec(:,ke,2)
      MOMZ(:,ke) = 0.0_RP
    end do
        
    call CubedSphereCoordCnv_LonLat2CSVec( &
      lcmesh%panelID, lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),  & ! (in)
      lcmesh%gam(:,lcmesh%NeS:lcmesh%NeE), lcmesh%Ne * elem%Np,    & ! (in)
      svec(:,:,1), svec(:,:,2),                                    & ! (in)
      MOMX(:,lcmesh%NeS:lcmesh%NeE), MOMY(:,lcmesh%NeS:lcmesh%NeE) ) ! (out)

    call calc_normalized_factor( tracer_field_list(iq)%ptr%val(:,ke), &
      lcmesh, elem )
    return
  end subroutine exp_SetInitCond_tracer_advection

  subroutine calc_normalized_factor( ptracer, lcmesh, elem )
    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD 
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(in) :: ptracer(elem%Np,lcmesh%NeA)

    integer :: ke
    real(RP) :: sum_tmp(3), sum_out(3)
    real(RP) :: int_weight(elem%Np)
    integer :: ierr
    !-------------------------------

    sum_tmp(:) = 0.0_RP
    !$omp parallel do private(ke, int_weight) reduction(+:sum_tmp)
    do ke=lcmesh%NeS, lcmesh%NeE
      int_weight(:) = lcmesh%Gsqrt(:,ke) * lcmesh%J(:,ke) * elem%IntWeight_lgl(:)
      sum_tmp(1) = sum_tmp(1) + sum( int_weight(:) )
      sum_tmp(2) = sum_tmp(2) + sum( int_weight(:) * abs(ptracer(:,ke)) )
      sum_tmp(3) = sum_tmp(3) + sum( int_weight(:) * ptracer(:,ke)**2 )
    end do

    call MPI_AllReduce( sum_tmp, sum_out, size(sum_tmp), &
      MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr       )

    LOG_INFO("Normalized_Factor",*) sum_out(1), sum_out(2)/sum_out(1), sqrt(sum_out(3)/sum_out(1))

    return
  end subroutine calc_normalized_factor

!- private --

!OCL SERIAL
  subroutine DCMIP2012_tracer( q, &
    qtrcname, lon, lat, X, Y, Z, Np )
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_CS2CartPos
    implicit none

    integer, intent(in) :: Np
    real(RP), intent(out) :: q(Np)
    character(*), intent(in) :: qtrcname
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: X(Np), Y(Np), Z(Np)

    integer :: i
    real(RP) :: d(Np,2)
    real(RP) :: r(Np)

    real(RP) :: Rt                           !< Horizontal half-width of tracers
    real(RP) :: lon_c(2)                     !< Initial longitude of first and second tracers
    real(RP), parameter :: lat_c = 0.0_RP    !< Initial latitude of tracers
!    real(RP), parameter :: lat_c = 0.5_RP * PI    !< Initial latitude of tracers

    real(RP) :: X_c(2), Y_c(2), Z_c(2)
    !----------------------------------------------

    lon_c(:) = (/ 5.0_RP, 7.0_RP /) * PI / 6.0_RP

    X_c(:) = RPlanet * cos(lat_c) * cos(lon_c(:)) 
    Y_c(:) = RPlanet * cos(lat_c) * sin(lon_c(:)) 
    Z_c(:) = RPlanet * sin(lat_c)

    do i=1, 2
      select case(qtrcname)
      case('cosine_bells')
        Rt = 0.5_RP 
        r(:) = acos(sin(lat_c)*sin(lat(:)) + cos(lat_c)*cos(lat(:))*cos(lon(:)-lon_c(i)))
        d(:,i) = min( 1.0_RP, r(:) / Rt )
      case('gaussian_hills')
        d(:,i) = ( (X(:) - X_c(i))**2 + (Y(:) - Y_c(i))**2 + (Z(:) - Z_c(i))**2 ) / RPlanet**2
      case('slotted_cylinders')
        Rt = 0.5_RP 
        r(:) = acos(sin(lat_c)*sin(lat(:)) + cos(lat_c)*cos(lat(:))*cos(lon(:)-lon_c(i)))
        d(:,i) = r(:) / Rt
      end select
    end do

    select case(qtrcname)
    case('cosine_bells')
      q(:) = 0.9_RP * ( 1.0_RP + 0.5_RP * ( cos(PI * d(:,1)) + cos(PI * d(:,2)) ) ) + 0.1_RP
    case('gaussian_hills')
      q(:) = 0.95_RP * ( exp( - 5.0_RP * d(:,1) ) + exp( - 5.0_RP * d(:,2) ) )  
    case('slotted_cylinders')
      q(:) = 0.0_RP
      where ( ( d(:,1) <= 1.0_RP ) .and. ( abs(lon(:)-lon_c(1)) >= Rt / 6.0_RP ) )
        q(:) = 1.0_RP
      elsewhere ( ( d(:,2) <= 1.0_RP ) .and. ( abs(lon(:)-lon_c(2)) >= Rt/6.0_RP ) )
        q(:) = 1.0_RP
      elsewhere ( ( d(:,1) <= 1.0_RP ) .and. ( abs(lon(:)-lon_c(1)) < Rt/6.0_RP ) .and. ( lat(:) - lat_c < - 5.0_RP/12.0_RP*Rt ) )
        q(:) = 1.0_RP
      elsewhere ( ( d(:,2) <= 1.0_RP ) .and. ( abs(lon(:)-lon_c(2)) < Rt/6.0_RP ) .and. ( lat(:) - lat_c > 5.0_RP/12.0_RP*Rt ) )
        q(:) = 1.0_RP
      end where
    case default
      LOG_ERROR('DCMIP2012_tracer',*) trim(qtrcname)//' is not supported. Check!'
      call PRC_abort
    end select
    return
  end subroutine DCMIP2012_tracer

!OCL SERIAL
  subroutine DCMIP2012_deformation_flow( U, V, psi, &
    lon, lat, time, Np )

    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: U(Np)
    real(RP), intent(out) :: V(Np)
    real(RP), intent(out) :: psi(Np)
    real(RP), intent(in) :: lon(Np)
    real(RP), intent(in) :: lat(Np)
    real(RP), intent(in) :: time

    real(RP) :: lon2(Np)
    real(RP), parameter :: tau = 1036800.0_RP  !> Period of motion [sec]
    real(RP), parameter :: sw = 1.0_RP
    real(RP) :: k
    !----------------------------------------------

    k = 10.0_RP * RPlanet / tau
    lon2(:) = lon(:) - sw * 2.0_RP * PI * time / tau !+ 0.5_RP*PI
    U(:) = k * sin(lon2(:))**2 * sin(2.0_RP*lat(:)) * cos(PI * time / tau) &
          + sw * 2.0_RP * PI * RPlanet / tau *  cos(lat(:))
    V(:) = k * sin(2.0_RP * lon2(:)) * cos(lat(:)) * cos(PI * time / tau)
    Psi(:) = RPlanet * ( &
      k * sin(lon2(:))**2 * cos(lat(:))**2 * cos(PI * time / tau) &
      - sw * 2.0_RP * PI * RPlanet / tau *  sin(lat(:)) )

    return
  end subroutine DCMIP2012_deformation_flow

end module mod_user
