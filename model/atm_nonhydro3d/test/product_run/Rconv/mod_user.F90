!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a test case of Rayleigh–Bénard convection
!!
!! @author Yuta Kawai, Team SCALE
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

  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D   
  use scale_localmesh_3d, only: LocalMesh3D   
  use scale_localmeshfield_base, only: LocalMeshFieldBase
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

  integer :: IntrpPolyOrder_h = 8
  integer :: IntrpPolyOrder_v = 8

  real(RP), private :: BTM_FIXED_HEAT_FLUX  = 15.88_RP
  real(RP), private :: TOP_FIXED_HEAT_FLUX  = 15.88_RP  ! ztop=1.6 km

  integer, parameter :: BC_HEAT_FIXED_TEMP    = 1
  integer, parameter :: BC_HEAT_FIXED_FLUX    = 2
  integer, parameter :: BC_HEAT_USER_MOD_NOSPEC = -1

  integer, parameter :: BC_MOM_NoSlip         = 1
  integer, parameter :: BC_MOM_CONST_BULKCOEF = 2
  integer, parameter :: BC_MOM_USER_MOD_NOSPEC = -1
  real(RP) :: MOM_CONST_BULKCOEF

  integer :: BTM_BC_MOM_TYPEID  = BC_MOM_USER_MOD_NOSPEC
  integer :: TOP_BC_MOM_TYPEID  = BC_MOM_USER_MOD_NOSPEC
  integer :: BTM_BC_HEAT_TYPEID = BC_HEAT_USER_MOD_NOSPEC ! 2: FIXED FLUX
  integer :: TOP_BC_HEAT_TYPEID = BC_HEAT_USER_MOD_NOSPEC ! 2: FIXED FLUX

  logical :: is_PREShyd_ref_set

  real(RP) :: StabCoef_bnd = 0.0_RP
  real(RP) :: StabCoef_bnd2 = 0.0_RP

  real(RP) :: U0
  !-----------------------------------------------------------------------------
contains
!OCL SERIAL
  subroutine USER_mkinit( this, atm )
    implicit none
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    type(Experiment) :: exp_manager
    !------------------------------------------

    call exp_manager%Init('Rayleigh--Bénard_convection')
    call exp_manager%Regist_SetInitCond( exp_SetInitCond_RB_convection )
    call this%UserBase%mkinit( atm, exp_manager )
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL
  subroutine USER_setup( this, atm )
    implicit none
    
    class(User), intent(inout) :: this
    class(AtmosComponent), intent(inout) :: atm

    character(H_SHORT) :: BTM_BC_TYPE_HEAT
    character(H_SHORT) :: TOP_BC_TYPE_HEAT
    character(H_SHORT) :: BTM_BC_TYPE_MOM
    character(H_SHORT) :: TOP_BC_TYPE_MOM

    logical :: USER_do = .false. !< do user step?
    namelist / PARAM_USER / &
       USER_do, &
       BTM_BC_TYPE_MOM,     &
       TOP_BC_TYPE_MOM,     &
       BTM_BC_TYPE_HEAT,    &
       TOP_BC_TYPE_HEAT,    &
       BTM_FIXED_HEAT_FLUX, &
       TOP_FIXED_HEAT_FLUX, &
       MOM_CONST_BULKCOEF,  &
       StabCoef_bnd,        &
       StabCoef_bnd2,       &
       U0

    integer :: ierr    
    !------------------------------------------


    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    BTM_BC_TYPE_MOM = "NoSpec"
    TOP_BC_TYPE_MOM = "NoSpec"

    BTM_BC_TYPE_HEAT = "NoSpec"
    TOP_BC_TYPE_HEAT = "NoSpec"

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
    BTM_BC_MOM_TYPEID = get_bctype_heat_id(BTM_BC_TYPE_HEAT)
    TOP_BC_MOM_TYPEID = get_bctype_heat_id(TOP_BC_TYPE_HEAT)

    BTM_BC_HEAT_TYPEID = get_bctype_heat_id(BTM_BC_TYPE_HEAT)
    TOP_BC_HEAT_TYPEID = get_bctype_heat_id(TOP_BC_TYPE_HEAT)

    !-
    is_PREShyd_ref_set = .false. 

    !-
    return
  end subroutine USER_setup

  function get_bctype_mom_id(bctype) result(type_id)
    character(*), intent(in) :: bctype
    integer :: type_id
    !------------------
    select case( trim(bctype) )
    case('CONST_BULKCOEF')
      type_id = BC_MOM_CONST_BULKCOEF
    case('NoSpec')
      type_id = BC_MOM_USER_MOD_NOSPEC
    case default
      LOG_ERROR("USER_setup",*) 'Not appropriate boundary condition. Check!', trim(bctype)
      call PRC_abort
    end select
    return
  end function get_bctype_mom_id  

  function get_bctype_heat_id(bctype) result(type_id)
    character(*), intent(in) :: bctype
    integer :: type_id
    !------------------
    select case( trim(bctype) )
    case('FixedFlux')
      type_id = BC_HEAT_FIXED_FLUX
    case('NoSpec')
      type_id = BC_HEAT_USER_MOD_NOSPEC
    case default
      LOG_ERROR("USER_setup",*) 'Not appropriate boundary condition. Check!', trim(bctype)
      call PRC_abort
    end select
    return
  end function get_bctype_heat_id  

  !------

!OCL SERIAL
  subroutine USER_calc_tendency( this, atm )
    use scale_time_manager, only:  TIME_NOWSTEP
    use scale_prc 
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in    
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,   &
      AtmosVars_GetLocalMeshPhyTends,  &
      AtmosVars_GetLocalMeshPhyAuxVars      
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRESHYD_VID => AUXVAR_PRESHYDRO_ID,         &
      THERMHYD_VID => AUXVAR_THERMHYDRO_ID,       &
      PRESHYD_REF_VID => AUXVAR_PRESHYDRO_REF_ID
      
    implicit none

    class(User), intent(inout) :: this 
    class(AtmosComponent), intent(inout) :: atm

    integer :: n
    class(LocalMesh3D), pointer :: lcmesh
    type(ElementBase3D), pointer :: elem3D

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(MeshField3D), pointer :: THERM_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p

    class(MeshField3D), pointer :: PRES_hyd_ref
    class(MeshField3D), pointer :: PRES_hyd_

    real(RP) :: dt
    integer :: ke
    !----------------------------

    dt = atm%time_manager%dtsec

    if ( .not. is_PREShyd_ref_set ) then
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_VID, PRES_hyd_ )
      call atm%vars%AUXVARS_manager%Get3D( PRESHYD_REF_VID, PRES_hyd_ref )

      do n=1, PRES_hyd_ref%mesh%LOCAL_MESH_NUM
        call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )        

        !$omp parallel do
        do ke=lcmesh%NeS, lcmesh%NeE
          PRES_hyd_ref%local(n)%val(:,ke) = PRES_hyd_%local(n)%val(:,ke)
        end do
      end do
      call atm%vars%AUXVARS_manager%MeshFieldComm_Exchange()
      call atm%dyn_proc%dyncore_driver%Update_phyd_hgrad( PRES_hyd_, PRES_hyd_ref, &
        PRES_hyd_%mesh, atm%mesh%element3D_operation )

      is_PREShyd_ref_set = .true.
    end if

    call atm%vars%AUXVARS_manager%Get3D( THERMHYD_VID, THERM_hyd )

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )     
       
      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )

      call AtmosVars_GetLocalMeshPhyTends( n,  atm%mesh%ptr_mesh, & 
        atm%vars%PHYTENDS_manager,                     &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp,   &
        RHOH_p  )

      !--
      call cal_tend_from_sfcflx( &
        DENS_tp%val, MOMX_tp%val, MOMY_tp%val, MOMZ_tp%val, RHOT_tp%val, &
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                  &
        DENS_hyd%val, PRES_hyd%val, THERM_hyd%local(n)%val,                  &
        PRES%val, CPtot%val,                                                 &
        atm%mesh%DOptrMat(3), atm%mesh%SOptrMat(3), atm%mesh%LiftOptrMat,    & 
        lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D )
    end do

    return
  end subroutine USER_calc_tendency


!OCL SERIAL  
  subroutine cal_tend_from_sfcflx( &
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp,  &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,               &
    DENS_hyd, PRES_hyd, THERM_hyd,                &
    PRES, CPtot,                                  &
    Dz, Sz, Lift, lcmesh, elem, lcmesh2D, elem2D  )

    use scale_sparsemat, only: &
      SparseMat, &
      sparsemat_matmul        
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(inout) :: DENS_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMZ_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOT_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: THERM_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem%Np,lcmesh%NeA)
    type(SparseMat), intent(in) :: Dz, Sz, Lift

    integer :: ke, ke2D, ke_top, ke_btm
    integer :: ij
    real(RP) :: SFLX_MOMX_top(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: SFLX_MOMX_btm(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: SFLX_MOMY_top(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: SFLX_MOMY_btm(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: SFLX_SH_top(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: SFLX_SH_btm(elem2D%Np,lcmesh2D%Ne)
    real(RP) :: LiftDelFlx(elem%Np)
    real(RP) :: del_flux_grad(elem%NfpTot,lcmesh%Ne,3)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%Ne,3)

    logical :: cal_tend_flag

    real(RP) :: Kh_top(elem%Np,lcmesh%Ne2D)
    real(RP) :: Kh_btm(elem%Np,lcmesh%Ne2D)

    real(RP) :: DzPT(elem%Np,lcmesh%NeA)
    real(RP) :: DzU(elem%Np,lcmesh%NeA)
    real(RP) :: DzV(elem%Np,lcmesh%NeA)

    real(RP) :: DENS(elem%Np)
    real(RP) :: PT(elem%Np)
    real(RP) :: Fz(elem%Np)
    real(RP) :: lz

    real(RP), parameter :: Pr = 0.7_RP
    !---------------------------------------------

    cal_tend_flag = .false.

    if ( BTM_BC_HEAT_TYPEID == BC_HEAT_FIXED_FLUX ) then
      SFLX_SH_btm(:,:) = BTM_FIXED_HEAT_FLUX
      cal_tend_flag = .true.
    else
      SFLX_SH_btm(:,:) = 0.0_RP
    end if

    if ( TOP_BC_HEAT_TYPEID == BC_HEAT_FIXED_FLUX ) then
      SFLX_SH_top(:,:) = TOP_FIXED_HEAT_FLUX
      cal_tend_flag = .true.
    else
      SFLX_SH_top(:,:) = 0.0_RP
    end if

    if ( cal_tend_flag ) then
      lz = 1.5_RP * ( lcmesh%zmax - lcmesh%zmin ) / real( lcmesh%NeZ * elem%Nnode_v, kind=RP )

      call cal_del_flux_grad( del_flux_grad, &
        lcmesh%normal_fn(:,:,3),                                   &
        DENS_hyd, THERM_hyd, DDENS, MOMX, MOMY, DRHOT,             &
        lcmesh%VMapM, lcmesh%VMapP, lcmesh, elem, lcmesh2D, elem2D )

      !$omp parallel do private(ke2D, ke_top, ke_btm, Fz, LiftDelFlx, DENS, PT)
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
        ke_btm = ke2D
        ke_top = ke2D + (lcmesh%NeZ-1)*lcmesh2D%Ne
        Kh_btm(:,ke2D) = StabCoef_bnd * exp( - (lcmesh%zlev(:,ke_btm) / lz)**2 )
        Kh_top(:,ke2D) = StabCoef_bnd * exp( - ( ( lcmesh%zmax - lcmesh%zlev(:,ke_top) ) / lz)**2 )

        !-
        DENS(:) = DENS_hyd(:,ke_btm) + DDENS(:,ke_btm)
        PT(:) = ( THERM_hyd(:,ke_btm) + DRHOT(:,ke_btm) ) / DENS(:)

        call sparsemat_matmul( Dz, MOMX(:,ke_btm) / DENS(:), Fz)
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_btm)*del_flux_grad(:,ke_btm,1), LiftDelFlx )
        DzU(:,ke_btm) = lcmesh%Escale(:,ke_btm,3,3) * Fz(:) + LiftDelFlx(:)

        call sparsemat_matmul( Dz, MOMY(:,ke_btm) / DENS(:), Fz)
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_btm)*del_flux_grad(:,ke_btm,2), LiftDelFlx )
        DzV(:,ke_btm) = lcmesh%Escale(:,ke_btm,3,3) * Fz(:) + LiftDelFlx(:)

        call sparsemat_matmul( Dz, PT(:), Fz)
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_btm)*del_flux_grad(:,ke_btm,3), LiftDelFlx )
        DzPT(:,ke_btm) = lcmesh%Escale(:,ke_btm,3,3) * Fz(:) + LiftDelFlx(:)

        !-
        DENS(:) = DENS_hyd(:,ke_top) + DDENS(:,ke_top)
        PT(:) = ( THERM_hyd(:,ke_top) + DRHOT(:,ke_top) ) / DENS(:)

        call sparsemat_matmul( Dz, MOMX(:,ke_top) / DENS(:), Fz)
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_btm)*del_flux_grad(:,ke_top,1), LiftDelFlx )
        DzU(:,ke_top) = lcmesh%Escale(:,ke_top,3,3) * Fz(:) + LiftDelFlx(:)

        call sparsemat_matmul( Dz, MOMY(:,ke_top) / DENS(:), Fz)
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_top)*del_flux_grad(:,ke_top,2), LiftDelFlx )
        DzV(:,ke_top) = lcmesh%Escale(:,ke_top,3,3) * Fz(:) + LiftDelFlx(:)

        call sparsemat_matmul( Dz, PT(:), Fz)
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_top)*del_flux_grad(:,ke_top,3), LiftDelFlx )
        DzPT(:,ke_top) = lcmesh%Escale(:,ke_top,3,3) * Fz(:) + LiftDelFlx(:)
      end do
      
      call cal_del_flux( del_flux, &
        SFLX_MOMX_btm, SFLX_MOMX_top, SFLX_MOMY_btm, SFLX_MOMY_top, &
        SFLX_SH_btm, SFLX_SH_top, DzPT, lcmesh%normal_fn(:,:,3),   &
        DENS_hyd, THERM_hyd, DDENS, MOMX, MOMY, DRHOT,             &
        lcmesh%VMapM, lcmesh%VMapP, lcmesh, elem, lcmesh2D, elem2D )

      !$omp parallel do private(ke2D, ke_top, ke_btm, Fz, LiftDelFlx, DENS)
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
        ke_btm = ke2D
        ke_top = ke2D + (lcmesh%NeZ-1)*lcmesh2D%Ne
        
        !-
        DENS(:) = DENS_hyd(:,ke_btm) + DDENS(:,ke_btm)

        ! call sparsemat_matmul( Sz, DENS(:) * Pr * Kh_btm(:,ke2D) * DzU(:,ke_btm), Fz)
        ! call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_btm)*del_flux(:,ke_btm,1), LiftDelFlx )
        ! MOMX_tp (:,ke_btm) = MOMX_tp (:,ke_btm) &
        !   - lcmesh%Escale(:,ke_btm,3,3) * Fz(:) !&
        !   ! - LiftDelFlx(:)

        ! call sparsemat_matmul( Sz, DENS(:) * Pr * Kh_btm(:,ke2D) * DzV(:,ke_btm), Fz)
        ! call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_btm)*del_flux(:,ke_btm,2), LiftDelFlx )
        ! MOMY_tp (:,ke_btm) = MOMY_tp (:,ke_btm) &
        !   - lcmesh%Escale(:,ke_btm,3,3) * Fz(:) !&
        !   ! - LiftDelFlx(:)

        call sparsemat_matmul( Sz, DENS(:) * Kh_btm(:,ke2D) * DzPT(:,ke_btm), Fz)
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_btm)*del_flux(:,ke_btm,3), LiftDelFlx )
        RHOT_tp (:,ke_btm) = RHOT_tp (:,ke_btm) &
          - lcmesh%Escale(:,ke_btm,3,3) * Fz(:) &
          - LiftDelFlx(:) / CPtot(:,ke_btm)

        !-
        DENS(:) = DENS_hyd(:,ke_top) + DDENS(:,ke_top)

        ! call sparsemat_matmul( Sz, DENS(:) * Pr * Kh_top(:,ke2D) * DzU(:,ke_top), Fz)
        ! call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_top)*del_flux(:,ke_top,1), LiftDelFlx )
        ! MOMX_tp (:,ke_top) = MOMX_tp (:,ke_top) &
        !   - lcmesh%Escale(:,ke_top,3,3) * Fz(:) !&
        !   ! - LiftDelFlx(:)

        ! call sparsemat_matmul( Sz, DENS(:) * Pr * Kh_top(:,ke2D) * DzV(:,ke_top), Fz)
        ! call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_top)*del_flux(:,ke_top,2), LiftDelFlx )
        ! MOMY_tp (:,ke_top) = MOMY_tp (:,ke_top) &
        !   - lcmesh%Escale(:,ke_top,3,3) * Fz(:) !&
        !   ! - LiftDelFlx(:)

        call sparsemat_matmul( Sz, DENS(:) * Kh_top(:,ke2D) * DzPT(:,ke_top), Fz)
        call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke_top)*del_flux(:,ke_top,3), LiftDelFlx )
        RHOT_tp (:,ke_top) = RHOT_tp (:,ke_top) &
          - lcmesh%Escale(:,ke_top,3,3) * Fz(:) &
          - LiftDelFlx(:) / CPtot(:,ke_top)
      end do
    end if

    return
  end subroutine cal_tend_from_sfcflx

!OCL SERIAL  
  subroutine cal_del_flux( del_flux,      &
    sflx_momx_btm, sflx_momx_top,         &
    sflx_momy_btm, sflx_momy_top,         &
    sflx_sh_btm,  sflx_sh_top,            &
    dz_pt,                                &
    nz,                                   &
    DENS_hyd, THERM_hyd, DDENS, MOMX, MOMY, DRHOT,      &
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D          )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D   
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(in) :: sflx_momx_btm(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: sflx_momx_top(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: sflx_momy_btm(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: sflx_momy_top(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: sflx_sh_btm(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: sflx_sh_top(elem2D%Np,lmesh2D%Ne)
    real(RP), intent(in) :: dz_pt(elem%Np*lmesh2D%NeA)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: THERM_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DDENS(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: MOMX(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: MOMY(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem%Np*lmesh%NeA)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)

    integer :: ke, ke2D, p
    integer :: i, i2
    integer :: iP, iM

    real(RP) :: DENS_M, DENS_P
    real(RP) :: MOMX_P, MOMY_P
    !------------------------------------------------------------------------

    !$omp parallel private(i, i2, iP, iM, ke, ke2D, p, DENS_M, DENS_P, MOMX_P, MOMY_P)
    !$omp workshare
    del_flux(:,:) = 0.0_RP
    !$omp end workshare

    !$omp do collapse(2)
    do ke2D=1, lmesh2D%Ne
    do p=1, elem2D%Np
      i  = elem%Nfaces_h*elem%Nfp_h + p              + (ke2D-1)*elem%NfpTot
      i2 = elem%Nfaces_h*elem%Nfp_h + elem%Nfp_v + p + (ke2D-1)*elem%NfpTot

      iM = vmapM(i); iP = vmapP(i)
      DENS_M = DDENS(iM) + DENS_hyd(iM)
      DENS_P = DDENS(iP) + DENS_hyd(iP)
      MOMX_P = 2.0_RP * DENS_M * U0 - MOMX(iM)
      MOMY_P =                      - MOMY(iM)

      del_flux(i,1) = &
!         sflx_momx_btm(p,ke2D) * nz(i) &
        - StabCoef_bnd2 * ( MOMX_P - MOMX(iM) )
      del_flux(i,2) = &
!          sflx_momy_btm(p,ke2D) * nz(i) &
        - StabCoef_bnd2 * ( MOMY_P - MOMY(iM) )
      del_flux(i,3) = sflx_sh_btm(p,ke2D) * nz(i)
    end do
    end do
    !$omp do collapse(2)
    do ke2D=1, lmesh2D%Ne
    do p=1, elem2D%Np
      ke = ke2D + (lmesh%NeZ-1)*lmesh%Ne2D
      i  = (ke-1)*elem%NfpTot + elem%Nfaces_h*elem%Nfp_h + elem%Nfp_v + p
      i2 = (ke-1)*elem%NfpTot + elem%Nfaces_h*elem%Nfp_h              + p
      
      !--
      iM = vmapM(i); iP = vmapP(i)
      DENS_M = DDENS(iM) + DENS_hyd(iM)
      DENS_P = DDENS(iP) + DENS_hyd(iP)
      MOMX_P = 2.0_RP * DENS_M * U0 - MOMX(iM)
      MOMY_P =                      - MOMY(iM)

      del_flux(i,1) = &
!          sflx_momx_top(p,ke2D) * nz(i) &
        - StabCoef_bnd2 * ( MOMX_P - MOMX(iM) )
      del_flux(i,2) = &
!          sflx_momy_top(p,ke2D) * nz(i)  &
        - StabCoef_bnd2 * ( MOMY_P - MOMY(iM) )
      del_flux(i,3) = sflx_sh_top(p,ke2D) * nz(i)
    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine cal_del_flux

!OCL SERIAL  
  subroutine cal_del_flux_grad( del_flux_grad,     &
    nz,                                            &
    DENS_hyd, THERM_hyd, DDENS, MOMX, MOMY, DRHOT, &
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D     )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D   
    real(RP), intent(out) ::  del_flux_grad(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: THERM_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DDENS(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: MOMX(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: MOMY(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem%Np*lmesh%NeA)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)

    integer :: ke, ke2D, p
    integer :: i, i2
    integer :: iP, iM

    real(RP) :: DENS_M, DENS_P
    real(RP) :: MOMX_P, MOMY_P

    real(RP) :: PT_M, PT_P
    !------------------------------------------------------------------------

    !$omp parallel private(i, i2, iP, iM, ke, ke2D, p, DENS_M, DENS_P, MOMX_P, MOMY_P, PT_M, PT_P)
    !$omp workshare
    del_flux_grad(:,:) = 0.0_RP
    !$omp end workshare

    !$omp do collapse(2)
    do ke2D=1, lmesh2D%Ne
    do p=1, elem2D%Np
      i  = elem%Nfaces_h*elem%Nfp_h + p              + (ke2D-1)*elem%NfpTot
      i2 = elem%Nfaces_h*elem%Nfp_h + elem%Nfp_v + p + (ke2D-1)*elem%NfpTot

      iM = vmapM(i); iP = vmapP(i)
      DENS_M = DDENS(iM) + DENS_hyd(iM)
      DENS_P = DDENS(iP) + DENS_hyd(iP)
      PT_M = ( THERM_hyd(iM) + DRHOT(iM) ) / DENS_M
      PT_P = ( THERM_hyd(iP) + DRHOT(iP) ) / DENS_P
      MOMX_P = 2.0_RP * DENS_M * U0 - MOMX(iM)
      MOMY_P =                      - MOMY(iM)

      del_flux_grad(i,1) = ( U0 - MOMX(iM) / DENS_M ) * nz(i)
      del_flux_grad(i,2) =      - MOMY(iM) / DENS_M   * nz(i)
      del_flux_grad(i,3) = 0.5_RP * ( PT_P - PT_M ) * nz(i)

      !---
      iM = vmapM(i2); iP = vmapP(i2)
      DENS_M = DDENS(iM) + DENS_hyd(iM)
      DENS_P = DDENS(iP) + DENS_hyd(iP)
      PT_M = ( THERM_hyd(iM) + DRHOT(iM) ) / DENS_M
      PT_P = ( THERM_hyd(iP) + DRHOT(iP) ) / DENS_P
      del_flux_grad(i2,1) = ( MOMX(iP) / DENS_P - MOMX(iM) / DENS_M ) * nz(i2)
      del_flux_grad(i2,2) = ( MOMY(iP) / DENS_P - MOMY(iM) / DENS_M ) * nz(i2)
      del_flux_grad(i2,3) = 0.5_RP * ( PT_P - PT_M ) * nz(i2)
    end do
    end do
    !$omp do collapse(2)
    do ke2D=1, lmesh2D%Ne
    do p=1, elem2D%Np
      ke = ke2D + (lmesh%NeZ-1)*lmesh%Ne2D
      i  = (ke-1)*elem%NfpTot + elem%Nfaces_h*elem%Nfp_h + elem%Nfp_v + p
      i2 = (ke-1)*elem%NfpTot + elem%Nfaces_h*elem%Nfp_h              + p
      
      !--
      iM = vmapM(i); iP = vmapP(i)
      DENS_M = DDENS(iM) + DENS_hyd(iM)
      DENS_P = DDENS(iP) + DENS_hyd(iP)
      PT_M = ( THERM_hyd(iM) + DRHOT(iM) ) / DENS_M
      PT_P = ( THERM_hyd(iP) + DRHOT(iP) ) / DENS_P
      MOMX_P = 2.0_RP * DENS_M * U0 - MOMX(iM)
      MOMY_P =                      - MOMY(iM)

      del_flux_grad(i,1) = ( U0 - MOMX(iM) / DENS_M ) * nz(i)
      del_flux_grad(i,2) =      - MOMY(iM) / DENS_M   * nz(i)
      del_flux_grad(i,3) = 0.5_RP * ( PT_P - PT_M ) * nz(i)

      !---
      iM = vmapM(i2); iP = vmapP(i2)
      DENS_M = DDENS(iM) + DENS_hyd(iM)
      DENS_P = DDENS(iP) + DENS_hyd(iP)
      del_flux_grad(i2,1) = ( MOMX(iP) / DENS_P - MOMX(iM) / DENS_M ) * nz(i2)
      del_flux_grad(i2,2) = ( MOMY(iP) / DENS_P - MOMY(iM) / DENS_M ) * nz(i2)
      del_flux_grad(i2,3) = 0.5_RP * ( PT_P - PT_M ) * nz(i2)
    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine cal_del_flux_grad

!OCL SERIAL  
  subroutine exp_SetInitCond_RB_convection( this,                   &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use scale_const, only: &
      PI => CONST_PI,       &
      GRAV => CONST_GRAV,   &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    
    use scale_random, only: &
      RANDOM_uniform
    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constBVFreq, &
      hydrostaic_build_rho_XYZ 
    
    use mod_mkinit_util, only: &
      mkinitutil_gen_GPMat
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
    
    real(RP) :: THETA0
    real(RP) :: DTHETA
    real(RP) :: BruntVaisalaFreq = 0.0E-2_RP ! [s-1]
    real(RP) :: U0

    real(RP) :: x_c, y_c, z_c
    real(RP) :: r_x, r_y, r_z
    character(H_SHORT) :: PERTURB_TYPE ! BUBBLE, RANDOM

    namelist /PARAM_EXP/ &
      BruntVaisalaFreq,         &
      U0,                       &
      THETA0, DTHETA,           &
      PERTURB_TYPE,             &
      x_c, y_c, z_c,            &
      r_x, r_y, r_z,            &
      IntrpPolyOrder_h,         &
      IntrpPolyOrder_v

    type(HexahedralElement) :: elem_intrp
    real(RP), allocatable :: IntrpMat(:,:)
    real(RP), allocatable :: x_intrp(:), y_intrp(:), z_intrp(:)
    real(RP) :: vx(elem%Nv), vy(elem%Nv), vz(elem%Nv)    
    
    real(RP) :: RovCp
    real(RP) :: PT  (elem%Np)
    real(RP) :: DENS(elem%Np)
    real(RP) :: exner_sfc
    real(RP) :: EXNER(elem%Np)

    real(RP) :: PT_hyd(elem%Np,lcmesh%NeA), PT_pertub(elem%Np,lcmesh%NeA)
    real(RP) :: PT_zxy(elem%Np,lcmesh%NeZ,lcmesh%NeX,lcmesh%NeY)

    real(RP) :: sfc_rhot(elem%Nfp_v)
    real(RP) :: bnd_SFC_PRES(elem%Nnode_h1D**2,lcmesh%Ne2DA)
    real(RP) :: rndm(elem%Np)

    integer :: ke, ke2D, ke_x, Ke_y, ke_z
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.0_RP; y_c = 0.0_RP; z_c = 5.E3_RP
    r_x = 5.0E3_RP; r_y = 5.0E3_RP; r_z = 1E2_RP

    U0      = 0.0_RP
    THETA0  = 300.0_RP
    DTHETA  = 0.01_RP
    PERTURB_TYPE = "BUBBLE"

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("RB_convection_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("RB_convection_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    call elem_intrp%Init( IntrpPolyOrder_h, IntrpPolyOrder_v, .false. )
    
    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    call mkinitutil_gen_GPMat( IntrpMat, elem_intrp, elem )

    allocate( x_intrp(elem_intrp%Np), y_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )

    RovCp = Rdry / CPdry
    exner_sfc = (PRES00 / PRES00)**RovCP
    
    if ( abs(BruntVaisalaFreq) < 1E-16_RP ) then
      !$omp parallel do private(PT, exner)
      do ke=lcmesh%NeS, lcmesh%NeE
        ! d exner / dz = - g / ( Cp * PT0 ) * exp (- N2/g * z)
        ! exner = exner(zs) - g^2 / (Cp * N^2) [ 1/PT (z) - 1/PT(zs) ] 
        PT(:) = THETA0
        exner(:) = exner_sfc - Grav / ( CpDry * THETA0 ) * lcmesh%zlev(:,ke)
  
        PRES_hyd(:,ke) = PRES00 * exner(:)**(CPdry/Rdry)
        DENS_hyd(:,ke) =  PRES_hyd(:,ke) / ( Rdry * exner(:) * PT(:) )
      end do      
    else
      call hydrostatic_calc_basicstate_constBVFreq( DENS_hyd, PRES_hyd, & ! (out)
        BruntVaisalaFreq, THETA0, PRES00, x, y, z, lcmesh, elem         ) ! (in)
    end if
    
    !--- Generate perturbation

    select case( PERTURB_TYPE )
    case ( "BUBBLE" )
      !$omp parallel do private( vx, vy, vz, x_intrp, y_intrp, z_intrp, PT, DENS )
      do ke=lcmesh%NeS, lcmesh%NeE
        vx(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),1)
        vy(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),2)
        vz(:) = lcmesh%pos_ev(lcmesh%EToV(ke,:),3)
        x_intrp(:) = vx(1) + 0.5_RP*(elem_intrp%x1(:) + 1.0_RP)*(vx(2) - vx(1))
        y_intrp(:) = vy(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vy(4) - vy(1))
        z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x3(:) + 1.0_RP)*(vz(5) - vz(1))
  
        PT(:) = matmul( IntrpMat, & 
            THETA0 * exp( BruntVaisalaFreq**2 / Grav * z_intrp(:) ) &
          + DTHETA / ( 1.0_RP + ((x_intrp(:) - x_c)/r_x)**2 + ((y_intrp(:) - y_c)/r_y)**2   + ((z_intrp(:) - z_c)/r_z)**2 ) )

        DENS(:) = PRES_hyd(:,ke) / ( Rdry * PT(:) * (PRES_hyd(:,ke)/PRES00)**(RovCp) )
        DDENS(:,ke) = DENS(:) - DENS_hyd(:,ke)
      end do
    case ( "RANDOM" )
      !$omp parallel do private( ke, ke2D, rndm, EXNER, sfc_rhot ) collapse(3)
      do ke_z=1, lcmesh%NeZ
      do ke_y=1, lcmesh%NeY
      do ke_x=1, lcmesh%NeX
        call RANDOM_uniform( rndm )
        ke2D = ke_x + (ke_y-1)*lcmesh%NeX
        ke = ke2D + (ke_z-1)*lcmesh%NeX*lcmesh%NeY

        EXNER(:) = ( PRES_hyd(:,ke) / PRES00 )**RovCp 
        PT_hyd(:,ke) = PRES00 / ( DENS_hyd(:,ke) * Rdry ) * EXNER(:)**(CVdry/Rdry)
        PT_zxy(:,ke_z,ke_x,ke_y) = PT_hyd(:,ke) &
            + ( rndm(:) * 2.0_RP - 1.0_RP ) * DTHETA

        DDENS(:,ke) = 0.0_RP

        if ( ke_z==1 ) then
          sfc_rhot(:) = DENS_hyd(elem%Hslice(:,1),ke2D) * PT_zxy(elem%Hslice(:,1),ke_z,ke_x,ke_y)        
          bnd_SFC_PRES(:,ke2D) = PRES00 * ( Rdry * sfc_rhot(:) / PRES00 )**( CPdry/CVdry )
        end if
      end do
      end do
      end do

      call hydrostaic_build_rho_XYZ( DDENS, &
        DENS_hyd, PRES_hyd, PT_zxy,                                         &
        lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2), lcmesh%pos_en(:,:,3),   &
        lcmesh, elem, bnd_SFC_PRES=bnd_SFC_PRES                             )

      !$omp parallel do private(ke2D, ke) collapse(3)
      do ke_z=1, lcmesh%NeZ
      do ke_y=1, lcmesh%NeY
      do ke_x=1, lcmesh%NeX
        ke2D = ke_x + (ke_y-1)*lcmesh%NeX
        ke = ke2D + (ke_z-1)*lcmesh%NeX*lcmesh%NeY
        DRHOT(:,ke) = &
            ( DENS_hyd(:,ke) + DDENS(:,ke) ) * PT_zxy(:,ke_z,ke_x,Ke_y) &
          - DENS_hyd(:,ke) * PT_hyd(:,ke)

        MOMX(:,ke) = ( DENS_hyd(:,ke) + DDENS(:,ke) ) * U0
      end do
      end do
      end do
    end select

    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_RB_convection

end module mod_user
