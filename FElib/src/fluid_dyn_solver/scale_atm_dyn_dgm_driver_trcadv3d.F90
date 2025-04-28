!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / DGM driver (tracer advection)
!!
!! @par Description
!!      Driver module for tracer advection based on DGM 
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_driver_trcadv3d
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,  &
    Rdry => CONST_Rdry,  &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00  
  use scale_tracer, only: QA

  use scale_timeint_rk, only: TimeInt_RK
  use scale_sparsemat, only: SparseMat

  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_element_operation_base, only: ElementOperationBase3D

  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField2D, MeshField3D
  use scale_localmeshfield_base, only: &
      LocalMeshFieldBaseList

  use scale_element_modalfilter, only: &
    ModalFilter

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo
  use scale_model_mesh_manager, only: ModelMesh3D

  use scale_atm_dyn_dgm_driver_base, only: &
    AtmDynDGMDriver_base3D, AtmDynDGMDriver_base3D_Init, AtmDynDGMDriver_base3D_Final

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    PRGVAR_NUM, &
    DENS_VID => PRGVAR_DDENS_ID, THERM_VID => PRGVAR_THERM_ID,&
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    AUXVAR_NUM, &
    PRESHYD_VID => AUXVAR_PRESHYDRO_ID, DENSHYD_VID => AUXVAR_DENSHYDRO_ID,                 &
    CPTOT_VID => AUXVAR_CPtot_ID, CVTOT_VID => AUXVAR_CVtot_ID, RTOT_VID => AUXVAR_Rtot_ID, &
    PRES_VID => AUXVAR_PRES_ID, &
    PHYTEND_NUM1 => PHYTEND_NUM
  
  use scale_atm_dyn_dgm_bnd, only: AtmDynBnd

  use scale_atm_dyn_dgm_trcadvect3d_heve, only: &
    atm_dyn_dgm_trcadvect3d_heve_Init,                 &
    atm_dyn_dgm_trcadvect3d_heve_Final,                &    
    atm_dyn_dgm_trcadvect3d_heve_calc_fct_coef,        &
    atm_dyn_dgm_trcadvect3d_heve_cal_tend,             &
    atm_dyn_dgm_trcadvect3d_TMAR,                      &
    atm_dyn_dgm_trcadvect3d_save_massflux,             &
    atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_advtest

  use scale_atm_dyn_dgm_driver_nonhydro3d, only: &
    AtmDynDGMDriver_nonhydro3d
    
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  type, extends(AtmDynDGMDriver_base3D), public :: AtmDynDGMDriver_trcadv3d
    integer :: EQS_TYPEID
    logical :: ONLY_TRACERADV_FLAG
    type(SparseMat) :: FaceIntMat    

    ! Limiter
    logical :: disable_limiter
    type(MeshField3D) :: alphaDensM, alphaDensP ! coeffcient with stabilizion terms in numerical flux (element boundary data)

    ! 
    type(MeshField3D), allocatable :: TRCVARS3D(:)
    type(ModelVarManager) :: TRCVAR3D_manager
    integer :: TRCVAR3D_commid

    type(MeshField3D), allocatable :: AUX_TRCVARS3D(:)
    type(ModelVarManager) :: AUXTRCVAR3D_manager
    integer :: AUXTRCVAR3D_commid

    type(MeshField3D), allocatable :: AUXTRC_FLUX_VARS3D(:)
    type(ModelVarManager) :: AUXTRC_FLUX_VAR3D_manager
    integer :: AUXTRC_FLUX_VAR3D_commid    

    ! element-wise modal filter
    logical :: MODALFILTER_FLAG
    type(ModalFilter) :: modal_filter_3d

    ! boundary_condition
    type(AtmDynBnd), pointer :: boundary_cond

  contains
    procedure :: Init => AtmDynDGMDriver_trcadv3d_Init
    procedure :: Final => AtmDynDGMDriver_trcadv3d_Final
    procedure :: Update => AtmDynDGMDriver_trcadv3d_update
  end type AtmDynDGMDriver_trcadv3d

  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: EQS_TYPEID_TRCADV3D_HEVE             = 1
  
  !-
  integer, public, parameter :: AUXTRCVARS3D_NUM         = 1
  integer, public, parameter :: AUXTRCVARS3D_FCTCOEF_ID  = 1

  type(VariableInfo), public :: ATMOS_DYN_AUXTRCVARS3D_VINFO(AUXTRCVARS3D_NUM)
  DATA ATMOS_DYN_AUXTRCVARS3D_VINFO / &
    VariableInfo( AUXTRCVARS3D_FCTCOEF_ID, 'TRCADV_FCTCOEF', '',  &
                  '1',  3, 'XYZ',  ''                          )  / 

  !-
  integer, public, parameter :: TRCVARS3D_NUM         = 3
  integer, public, parameter :: TRCVARS3D_TRCADV_ID   = 1
  integer, public, parameter :: TRCVARS3D_DENS_ID     = 2
  integer, public, parameter :: TRCVARS3D_DENS0_ID    = 3

  type(VariableInfo), public :: ATMOS_DYN_TRCVARS3D_VINFO(TRCVARS3D_NUM)
  DATA ATMOS_DYN_TRCVARS3D_VINFO / &
    VariableInfo( TRCVARS3D_TRCADV_ID, 'TRCADV', '',     '1',  3, 'XYZ',  '' ), &
    VariableInfo( TRCVARS3D_DENS_ID  ,   'DENS', '', 'kg/m3',  3, 'XYZ',  '' ), & 
    VariableInfo( TRCVARS3D_DENS0_ID  ,  'DENS0', '', 'kg/m3',  3, 'XYZ',  '' ) / 

  !-
  integer, public, parameter :: MASS_FLUX_NUM   = 3
  integer, public, parameter :: MASSFLX_Z_ID    = 1  
  integer, public, parameter :: MASSFLX_X_ID    = 2
  integer, public, parameter :: MASSFLX_Y_ID    = 3

  type(VariableInfo), public :: ATMOS_DYN_MASS_FLUX_VINFO(MASS_FLUX_NUM)
  DATA ATMOS_DYN_MASS_FLUX_VINFO / &
    VariableInfo( MASSFLX_Z_ID, 'MASSFLX_Z', 'flux in z-direction',  &
                  'kg/s/m2',  3, 'XYZ',  ''                     ),   &
    VariableInfo( MASSFLX_X_ID, 'MASSFLX_X', 'flux in x-direction',  &
                  'kg/s/m2',  3, 'XYZ',  ''                     ),   & 
    VariableInfo( MASSFLX_Y_ID, 'MASSFLX_Y', 'flux in y-direction',  &
                  'kg/s/m2',  3, 'XYZ',  ''                     )    / 

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: setup_modalfilter

contains
!OCL SERIAL  
  subroutine AtmDynDGMDriver_trcadv3d_Init( this, &
    eqs_type_name, tint_type_name, dtsec,               &
    modal_filter_flag, disable_limiter,                 &
    model_mesh3D, boundary_cond, ONLY_TRACERADV_FLAG    )

    use scale_localmeshfield_base, only: LOCAL_MESHFIELD_TYPE_NODES_FACEVAL
    implicit none

    class(AtmDynDGMDriver_trcadv3d), intent(inout) :: this
    character(len=*), intent(in) :: eqs_type_name
    character(len=*), intent(in) :: tint_type_name
    real(DP), intent(in) :: dtsec
    logical, intent(in) :: modal_filter_flag
    logical, intent(in) :: disable_limiter
    class(ModelMesh3D), intent(inout), target :: model_mesh3D
    class(AtmDynBnd), intent(in), target :: boundary_cond
    logical, intent(in) :: ONLY_TRACERADV_FLAG

    class(MeshBase3D), pointer :: mesh3D
    class(HexahedralElement), pointer :: refElem3D
    class(ElementBase), pointer :: refElem

    integer :: iv
    !-----------------------------------------------------------------------------

    mesh3D => model_mesh3D%ptr_mesh

    call AtmDynDGMDriver_base3D_Init( this, &
      1,                                    &
      tint_type_name, dtsec,                &
      mesh3D                                )

    select case(eqs_type_name)
    !-- HEVE ------------------
    case("TRCADV3D_HEVE")
    case default
        LOG_ERROR("AtmDynDGMDriver_trcadv3d_Init",*) 'Invalid EQS_TYPE in namelist PARAM_ATMOS_DYN. Check!'
        call PRC_abort  
    end select
    
    this%ONLY_TRACERADV_FLAG = ONLY_TRACERADV_FLAG
    this%disable_limiter = disable_limiter

    !- Initialize a module for tracer equations
    call atm_dyn_dgm_trcadvect3d_heve_Init( mesh3D, this%FaceIntMat ) 


    !- Initialize 3D auxiliary variables for tracer advection

    call this%TRCVAR3D_manager%Init()
    allocate( this%TRCVARS3D(TRCVARS3D_NUM) )

    do iv = 1, TRCVARS3D_NUM
      call this%TRCVAR3D_manager%Regist( &
        ATMOS_DYN_TRCVARS3D_VINFO(iv), mesh3D,     & ! (in) 
        this%TRCVARS3D(iv),                        & ! (inout)
        .false., fill_zero=.true.                  ) ! (in)
    end do

    iv = TRCVARS3D_TRCADV_ID
    call model_mesh3D%Create_communicator( &
      1, 0, 0,                         & ! (in) 
      this%TRCVAR3D_manager,           & ! (in)
      this%TRCVARS3D(iv:iv),           & ! (in)
      this%TRCVAR3D_commid             ) ! (out)
      
    !- Initialize variables to store time-averaged 3D mass flux

    call this%AUXTRC_FLUX_VAR3D_manager%Init()
    allocate( this%AUXTRC_FLUX_VARS3D(MASS_FLUX_NUM) )
    
    do iv = 1, MASS_FLUX_NUM
      call this%AUXTRC_FLUX_VAR3D_manager%Regist(  &
        ATMOS_DYN_MASS_FLUX_VINFO(iv), mesh3D,     & ! (in) 
        this%AUXTRC_FLUX_VARS3D(iv),               & ! (inout)
        .false., fill_zero=.true.                  ) ! (in)
    end do

    call model_mesh3D%Create_communicator( &
      1, 1, 0,                           & ! (in) 
      this%AUXTRC_FLUX_VAR3D_manager,    & ! (in)
      this%AUXTRC_FLUX_VARS3D(:),        & ! (in)
      this%AUXTRC_FLUX_VAR3D_commid      ) ! (out)

    !- Initialize 3D auxiliary variables for preserving nonnegativity in tracer advection

    call this%AUXTRCVAR3D_manager%Init()
    allocate( this%AUX_TRCVARS3D(AUXTRCVARS3D_NUM) )

    do iv = 1, AUXTRCVARS3D_NUM
      call this%AUXTRCVAR3D_manager%Regist( &
        ATMOS_DYN_AUXTRCVARS3D_VINFO(iv), mesh3D,  & ! (in) 
        this%AUX_TRCVARS3D(iv),                    & ! (inout)
        .false., fill_zero=.true.                  ) ! (in)
    end do

    call model_mesh3D%Create_communicator( &
      1, 0, 0,                         & ! (in) 
      this%AUXTRCVAR3D_manager,        & ! (in)
      this%AUX_TRCVARS3D(:),           & ! (in)
      this%AUXTRCVAR3D_commid          ) ! (out)

    !-    
    call this%alphaDensM%Init( "alphaDensM", "kg/m3.m/s", mesh3D, LOCAL_MESHFIELD_TYPE_NODES_FACEVAL )
    call this%alphaDensP%Init( "alphaDensP", "kg/m3.m/s", mesh3D, LOCAL_MESHFIELD_TYPE_NODES_FACEVAL )
  
    !- initialize an object to manage boundary conditions
    this%boundary_cond => boundary_cond


    !- initialize an object to manage modal filter
    this%MODALFILTER_FLAG = modal_filter_flag
    if (this%MODALFILTER_FLAG) then
      refElem => mesh3D%refElem
      select type(refElem)
      class is (HexahedralElement)
        refElem3D => refElem
      end select

      call setup_modalfilter( this, refElem3D, model_mesh3D%element3D_operation )
    end if

    return
  end subroutine AtmDynDGMDriver_trcadv3d_Init

!OCL SERIAL
  subroutine AtmDynDGMDriver_trcadv3d_update( this, &
    TRC_VARS, PROG_VARS, AUX_VARS, PHYTENDS,             &
    element_operation,                                   &
    Dx, Dy, Dz, Sx, Sy, Sz, Lift, mesh3D,                &
    dyn_driver                                           )

    use scale_tracer, only: &
      TRACER_ADVC, TRACER_NAME
    use scale_atm_dyn_dgm_modalfilter, only: &
      atm_dyn_dgm_tracer_modalfilter_apply
    implicit none

    class(AtmDynDGMDriver_trcadv3d), intent(inout) :: this
    class(MeshBase3D), intent(in), target :: mesh3D
    class(ModelVarManager), intent(inout) :: TRC_VARS
    class(ModelVarManager), intent(inout) :: PROG_VARS
    class(ModelVarManager), intent(inout) :: AUX_VARS
    class(ModelVarManager), intent(inout) :: PHYTENDS
    class(ElementOperationBase3D), intent(in) :: element_operation
    type(SparseMat), intent(in) :: Dx, Dy, Dz
    type(SparseMat), intent(in) :: Sx, Sy, Sz
    type(SparseMat), intent(in) :: Lift
    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: dyn_driver

    integer :: rkstage
    integer :: tintbuf_ind
    real(RP) :: dt
    real(RP) :: dttmp_trc

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    integer :: ke

    class(MeshField3D), pointer :: DDENS, MOMX, MOMY, MOMZ, THERM
    class(MeshField3D), pointer :: DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, PRES
    class(MeshField3D), pointer :: QTRC, RHOQ_tp

    class(MeshField3D), pointer :: QTRC_tmp, DDENS_TRC, DDENS0_TRC
    class(MeshField3D), pointer :: MFLX_x_tavg, MFLX_y_tavg, MFLX_z_tavg

    integer :: iq    
    !-----------------------------------------------------------------------------
    
    !-- prepairation

    call PROF_rapstart( 'ATM_DYN_trc_update_pre', 2)
    call PROG_VARS%Get3D(DENS_VID , DDENS)
    call PROG_VARS%Get3D(THERM_VID, THERM)

    call AUX_VARS%Get3D( PRESHYD_VID, PRES_hyd )
    call AUX_VARS%Get3D( DENSHYD_VID, DENS_hyd )
    call AUX_VARS%Get3D( PRES_VID, PRES )
    call AUX_VARS%Get3D( RTOT_VID, Rtot )
    call AUX_VARS%Get3D( CVTOT_VID, CVtot )
    call AUX_VARS%Get3D( CPTOT_VID, CPtot )

    call this%TRCVAR3D_manager%Get3D( TRCVARS3D_TRCADV_ID, QTRC_tmp )
    call this%TRCVAR3D_manager%Get3D(  TRCVARS3D_DENS_ID,  DDENS_TRC )
    call this%TRCVAR3D_manager%Get3D( TRCVARS3D_DENS0_ID, DDENS0_TRC )

    call this%AUXTRC_FLUX_VAR3D_manager%Get3D( MASSFLX_Z_ID, MFLX_z_tavg )
    call this%AUXTRC_FLUX_VAR3D_manager%Get3D( MASSFLX_X_ID, MFLX_x_tavg )
    call this%AUXTRC_FLUX_VAR3D_manager%Get3D( MASSFLX_Y_ID, MFLX_y_tavg )

    if ( this%ONLY_TRACERADV_FLAG ) then
      call PROG_VARS%Get3D(MOMZ_VID , MOMZ)
      call PROG_VARS%Get3D(MOMX_VID , MOMX)
      call PROG_VARS%Get3D(MOMY_VID , MOMY)

      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(n)

        !$omp parallel do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          DDENS0_TRC %local(n)%val(:,ke) = DDENS%local(n)%val(:,ke) 
          DDENS_TRC  %local(n)%val(:,ke) = DDENS%local(n)%val(:,ke)
          MFLX_z_tavg%local(n)%val(:,ke) = MOMZ%local(n)%val(:,ke)
          MFLX_x_tavg%local(n)%val(:,ke) = MOMX%local(n)%val(:,ke)
          MFLX_y_tavg%local(n)%val(:,ke) = MOMY%local(n)%val(:,ke)
        end do
        call atm_dyn_dgm_trcadvect3d_heve_cal_alphdens_advtest( &
          this%alphaDensM%local(n)%face_val, this%alphaDensP%local(n)%face_val,            & ! (inout)
          DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val,     & ! (in)
          DENS_hyd%local(n)%val,                                                           & ! (in)
          lcmesh3D%Gsqrt,                                                                  & ! (in)
          lcmesh3D%normal_fn(:,:,1), lcmesh3D%normal_fn(:,:,2), lcmesh3D%normal_fn(:,:,3), & ! (in)
          lcmesh3D%VMapM, lcmesh3D%VMapP, lcmesh3D, lcmesh3D%refElem3D                     ) ! (in)
      end do
    end if

    !* Exchange halo data of mass flux

    call PROF_rapstart( 'ATM_DYN_exchange_mflx', 3)
    call this%AUXTRC_FLUX_VAR3D_manager%MeshFieldComm_Exchange()
    call PROF_rapend( 'ATM_DYN_exchange_mflx', 3)

    call PROF_rapstart( 'ATM_DYN_applyBC_mflux', 3)
    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)

      call this%boundary_cond%ApplyBC_PROGVARS_lc( n, & ! (in)
        DDENS_TRC%local(n)%val(:,:), MFLX_x_tavg%local(n)%val, MFLX_y_tavg%local(n)%val, MFLX_z_tavg%local(n)%val, THERM%local(n)%val, & ! (inout)
        DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                                                                                  & ! (in)
        lcmesh3D%Gsqrt(:,:), lcmesh3D%GsqrtH(:,:), lcmesh3D%GIJ(:,:,1,1), lcmesh3D%GIJ(:,:,1,2), lcmesh3D%GIJ(:,:,2,2),                & ! (in)
        lcmesh3D%GI3(:,:,1), lcmesh3D%GI3(:,:,2),                                                                                      & ! (in)
        lcmesh3D%normal_fn(:,:,1), lcmesh3D%normal_fn(:,:,2), lcmesh3D%normal_fn(:,:,3),                                               & ! (in)
        lcmesh3D%vmapM, lcmesh3D%vmapP, lcmesh3D%vmapB, lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D   ) ! (in)
    end do
    call PROF_rapend( 'ATM_DYN_applyBC_mflux', 3)

    call PROF_rapend( 'ATM_DYN_trc_update_pre', 2)

    !-- Tracer advection ------------------------------------------------

    call PROF_rapstart( 'ATM_DYN_trc_update', 2)

    do iq = 1, QA
    
      call TRC_VARS%Get3D( iq, QTRC )
      call PHYTENDS%Get3D( PHYTEND_NUM1 + iq, RHOQ_tp )

      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(n)

        !$omp parallel do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          QTRC_tmp%local(n)%val(:,ke) = QTRC%local(n)%val(:,ke)
        end do            
      end do

      do rkstage=1, this%tint(1)%nstage
        if ( TRACER_ADVC(iq) ) then
          call PROF_rapstart( 'ATM_DYN_exchange_trc', 3)
          call this%TRCVAR3D_manager%MeshFieldComm_Exchange()
          call PROF_rapend( 'ATM_DYN_exchange_trc', 3)


          do n=1, mesh3D%LOCAL_MESH_NUM
            lcmesh3D => mesh3D%lcmesh_list(n)

            dt = this%tint(n)%Get_deltime()
            dttmp_trc = dt * this%tint(n)%coef_gam_ex(rkstage+1,rkstage) &
                         / this%tint(n)%coef_sig_ex(rkstage+1,rkstage)

            call atm_dyn_dgm_trcadvect3d_heve_calc_fct_coef( &
              this%AUX_TRCVARS3D(AUXTRCVARS3D_FCTCOEF_ID)%local(n)%val,                     & ! (out)
              QTRC_tmp%local(n)%val,                                                        & ! (in)
              MFLX_x_tavg%local(n)%val, MFLX_y_tavg%local(n)%val, MFLX_z_tavg%local(n)%val, & ! (in) 
              RHOQ_tp%local(n)%val,                                                         & ! (in)
              this%alphaDensM%local(n)%face_val, this%alphaDensP%local(n)%face_val,         & ! (in)
              DENS_hyd%local(n)%val, DDENS_TRC%local(n)%val, DDENS0_TRC%local(n)%val,       & ! (in)
              this%tint(n)%coef_c_ex(rkstage), dttmp_trc,                                   & ! (in) 
              Dx, Dy, Dz, Sx, Sy, Sz, Lift, this%FaceIntMat,                                & ! (in)
              lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D, & ! (in)
              this%disable_limiter                                                          ) ! (in)
          end do

          call PROF_rapstart( 'ATM_DYN_exchange_trc', 3)
          call this%AUXTRCVAR3D_manager%MeshFieldComm_Exchange()
          call PROF_rapend( 'ATM_DYN_exchange_trc', 3)
        end if

        do n=1, mesh3D%LOCAL_MESH_NUM
          lcmesh3D => mesh3D%lcmesh_list(n)

          tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)

          call PROF_rapstart( 'ATM_DYN_update_caltend_ex_trc', 3)  
          if ( TRACER_ADVC(iq) ) then 
            call atm_dyn_dgm_trcadvect3d_heve_cal_tend( &        
              this%tint(n)%tend_buf2D_ex(:,:,1,tintbuf_ind),                                & ! (out)
              QTRC_tmp%local(n)%val,                                                        & ! (in)
              MFLX_x_tavg%local(n)%val, MFLX_y_tavg%local(n)%val, MFLX_z_tavg%local(n)%val, & ! (in) 
              this%alphaDensM%local(n)%face_val, this%alphaDensP%local(n)%face_val,         & ! (in)
              this%AUX_TRCVARS3D(AUXTRCVARS3D_FCTCOEF_ID)%local(n)%val,                     & ! (out)
              RHOQ_tp%local(n)%val,                                                         & ! (in)
              Dx, Dy, Dz, Sx, Sy, Sz, Lift, this%FaceIntMat,                                & ! (in)
              lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D  ) ! (in)
          else
            !$omp parallel do
            do ke=lcmesh3D%NeS, lcmesh3D%NeE
              this%tint(n)%tend_buf2D_ex(:,ke,1,tintbuf_ind) = RHOQ_tp%local(n)%val(:,ke)
            end do
          end if
          call PROF_rapend( 'ATM_DYN_update_caltend_ex_trc', 3)

          call PROF_rapstart( 'ATM_DYN_update_advance_trc', 3)                
          call this%tint(n)%Advance_trcvar( &
            rkstage, QTRC_tmp%local(n)%val, 1,                                     &
            1, lcmesh3D%refElem%Np, lcmesh3D%NeS, lcmesh3D%NeE,                    &
            DDENS_TRC%local(n)%val, DDENS0_TRC%local(n)%val, DENS_hyd%local(n)%val ) 
          call PROF_rapend( 'ATM_DYN_update_advance_trc', 3)

          if ( rkstage == this%tint(1)%nstage .and. this%MODALFILTER_FLAG ) then
            call PROF_rapstart( 'ATM_DYN_update_qtrc_modalfilter', 3)
            call atm_dyn_dgm_tracer_modalfilter_apply( &
              QTRC_tmp%local(n)%val,                              & ! (inout)
              DENS_hyd%local(n)%val, DDENS_TRC%local(n)%val,      & ! (in)
              lcmesh3D, lcmesh3D%refElem3D, element_operation     ) ! (in)
            call PROF_rapend( 'ATM_DYN_update_qtrc_modalfilter', 3)
          end if

          if ( TRACER_ADVC(iq)                   &
            .and. rkstage == this%tint(1)%nstage &
            .and. this%ONLY_TRACERADV_FLAG       &            
            .and. ( .not. this%disable_limiter ) ) then

            call PROF_rapstart( 'ATM_DYN_update_trc_TMAR', 3)             
            call atm_dyn_dgm_trcadvect3d_TMAR( &
              QTRC_tmp%local(n)%val,                                                       & ! (inout)
              DENS_hyd%local(n)%val, DDENS_TRC%local(n)%val,                               & ! (in)
              lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D ) ! (in)
            call PROF_rapend( 'ATM_DYN_update_trc_TMAR', 3)  
          end if

        end do
      end do ! end for RK loop

      do n=1, mesh3D%LOCAL_MESH_NUM
        !$omp parallel do
        do ke=lcmesh3D%NeS, lcmesh3D%NeE
          QTRC%local(n)%val(:,ke) = ( DENS_hyd%local(n)%val(:,ke) + DDENS_TRC%local(n)%val(:,ke) ) &
                                  / ( DENS_hyd%local(n)%val(:,ke) + DDENS%local(n)%val(:,ke) )     &
                                  * QTRC_tmp%local(n)%val(:,ke)
        end do            
      end do

    end do ! end do for iq

    ! Update pressure
    call update_pressure_specific_heat( PRES, Rtot, CVtot, CPtot, & ! (inout)
      TRC_VARS, PROG_VARS, AUX_VARS, dyn_driver, mesh3D           ) ! (in)

    ! Negative fixer
    if (         .not. this%disable_limiter       & 
         .and. ( .not. this%ONLY_TRACERADV_FLAG ) ) then

      call PROF_rapstart( 'ATM_DYN_trc_negative_fixer', 2)          
      call fix_negative_val( &
        TRC_VARS, DDENS, THERM, PRES, CVtot, CPtot, Rtot,               & ! (inout)
        DENS_hyd, PRES_hyd, dt, mesh3D, dyn_driver%Is_THERMVAR_RHOT()   ) ! (in)
      call PROF_rapend( 'ATM_DYN_trc_negative_fixer', 2)
    end if

    call PROF_rapend( 'ATM_DYN_trc_update', 2)

    return
  end subroutine AtmDynDGMDriver_trcadv3d_update

!OCL SERIAL
  subroutine AtmDynDGMDriver_trcadv3d_Final( this )
    implicit none

    class(AtmDynDGMDriver_trcadv3d), intent(inout) :: this
    !-----------------------------------------------------------------------------

    call atm_dyn_dgm_trcadvect3d_heve_Final()

    nullify( this%boundary_cond )
    
    call this%TRCVAR3D_manager%Final()
    deallocate( this%TRCVARS3D )

    call this%AUXTRCVAR3D_manager%Final()
    deallocate( this%AUX_TRCVARS3D )

    call this%AUXTRC_FLUX_VAR3D_manager%Final()
    deallocate( this%AUXTRC_FLUX_VARS3D )

    call this%alphaDensM%Final()
    call this%alphaDensP%Final()

    if ( this%MODALFILTER_FLAG ) then
      call this%modal_filter_3d%Final()
    end if
    
    call this%FaceIntMat%Final()
    call AtmDynDGMDriver_base3D_Final( this )   

    return
  end subroutine AtmDynDGMDriver_trcadv3d_Final

!-- private
!OCL SERIAL
  subroutine fix_negative_val( &
    TRC_VARS, DDENS, THERM, PRES, CVtot, CPtot, Rtot, & ! (inout)
    DENS_hyd, PRES_hyd, dt, mesh3D, IS_THERMVAR_RHOT  ) ! (in)
    use scale_atmos_hydrometeor, only: &
      QLA, QIA
    use scale_atm_phy_mp_dgm_common, only: &
      atm_phy_mp_dgm_common_negative_fixer
    implicit none
    class(ModelVarManager), intent(inout) :: TRC_VARS
    class(MeshField3D), intent(inout) :: DDENS
    class(MeshField3D), intent(inout) :: THERM
    class(MeshField3D), intent(inout) :: PRES
    class(MeshField3D), intent(inout) :: CVtot
    class(MeshField3D), intent(inout) :: CPtot
    class(MeshField3D), intent(inout) :: Rtot
    class(MeshField3D), intent(inout) :: DENS_hyd
    class(MeshField3D), intent(inout) :: PRES_hyd
    real(RP), intent(in) :: dt
    class(MeshBase3D), intent(in), target :: mesh3D
    logical, intent(in) :: IS_THERMVAR_RHOT

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n

    integer :: iq

    integer :: trcid_list(QA)
    type(LocalMeshFieldBaseList) :: lc_qtrc(QA)
    !---------------------------------------------------------------

    do iq=1, QA 
      trcid_list(iq) = iq
    end do

    do n=1, mesh3D%LOCAL_MESH_NUM
      call TRC_VARS%GetLocalMeshFieldList( trcid_list, n, lc_qtrc )
      lcmesh3D => mesh3D%lcmesh_list(n)

      if ( IS_THERMVAR_RHOT ) then
        call atm_phy_mp_dgm_common_negative_fixer( &
          lc_qtrc, DDENS%local(n)%val, PRES%local(n)%val,            & ! (inout)
          CVtot%local(n)%val, CPtot%local(n)%val, Rtot%local(n)%val, & ! (inout)
          DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,              & ! (in)
          dt, lcmesh3D, lcmesh3D%refElem3D, QA, QLA, QIA,            & ! (in)
          DRHOT=THERM%local(n)%val                                   ) ! (inout)
      else
        call atm_phy_mp_dgm_common_negative_fixer( &
          lc_qtrc, DDENS%local(n)%val, PRES%local(n)%val,            & ! (inout)
          CVtot%local(n)%val, CPtot%local(n)%val, Rtot%local(n)%val, & ! (inout)
          DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,              & ! (in)
          dt, lcmesh3D, lcmesh3D%refElem3D, QA, QLA, QIA             ) ! (in)
      end if     
    end do

    return
  end subroutine fix_negative_val

  subroutine update_pressure_specific_heat( &
    PRES, Rtot, CVtot, CPtot,                         &
    TRC_VARS, PROG_VARS, AUX_VARS, dyn_driver, mesh3D )

    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      atm_dyn_dgm_nonhydro3d_common_calc_pressure      
    implicit none

    class(MeshField3D), intent(inout) :: PRES
    class(MeshField3D), intent(inout) :: Rtot
    class(MeshField3D), intent(inout) :: CVtot
    class(MeshField3D), intent(inout) :: CPtot
    class(ModelVarManager), intent(inout) :: TRC_VARS
    class(ModelVarManager), intent(inout) :: PROG_VARS
    class(ModelVarManager), intent(inout) :: AUX_VARS
    class(AtmDynDGMDriver_nonhydro3d), intent(inout) :: dyn_driver
    class(MeshBase3D), intent(in), target :: mesh3D

    integer :: n
    integer :: ke
    integer :: iq

    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    real(RP), allocatable :: q_tmp(:,:)
    real(RP), allocatable :: qdry(:)

    integer :: trcid_list(QA)
    type(LocalMeshFieldBaseList) :: lc_qtrc(QA)
    !---------------------------------------------------------------

    do iq=1, QA 
      trcid_list(iq) = iq
    end do

    do n=1, mesh3D%LOCAL_MESH_NUM
      call TRC_VARS%GetLocalMeshFieldList( trcid_list, n, lc_qtrc )
      
      lcmesh3D => mesh3D%lcmesh_list(n)
      elem3D => lcmesh3D%refElem3D
      allocate( q_tmp(elem3D%Np,QA), qdry(elem3D%Np) )

      !$omp parallel do private(ke, iq, q_tmp, qdry)
      do ke = lcmesh3D%NeS, lcmesh3D%NeE
        do iq = 1, QA
          q_tmp(:,iq) = lc_qtrc(iq)%ptr%val(:,ke)
        end do
        call ATMOS_THERMODYN_specific_heat( &
          elem3D%Np, 1, elem3D%Np, QA,                                                      & ! (in)
          q_tmp(:,:), TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:),              & ! (in)
          qdry, Rtot%local(n)%val(:,ke), CVtot%local(n)%val(:,ke), CPtot%local(n)%val(:,ke) ) ! (out)
      end do
      deallocate(q_tmp, qdry)      
    end do

    call dyn_driver%calc_pressure( PRES, PROG_VARS, AUX_VARS )

    return
  end subroutine update_pressure_specific_heat

  !-- Setup modal filter
!OCL SERIAL
  subroutine setup_modalfilter( this, refElem3D, element_operation )
    use scale_element_line, only: LineElement
    implicit none

    class(AtmDynDGMDriver_trcadv3d), target, intent(inout) :: this
    class(HexahedralElement), target, intent(in) :: refElem3D
    class(ElementOperationBase3D), intent(inout) :: element_operation

    real(RP) :: MF_ETAC_h  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_h = 36.0_RP
    integer  :: MF_ORDER_h = 16
    real(RP) :: MF_ETAC_v  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_v = 36.0_RP
    integer  :: MF_ORDER_v = 16

    namelist /PARAM_ATMOS_DYN_TRACER_MODALFILTER/ &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    

    integer :: ierr

    type(LineElement) :: elemV1D
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_TRACER_MODALFILTER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_TRCADV3D_setup_modalfilter",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_TRCADV3D_setup_modalfilter",*) 'Invalid names in namelist PARAM_ATMOS_DYN_TRACER_MODALFILTER. Check!'
      call PRC_abort
    endif

    LOG_NML(PARAM_ATMOS_DYN_TRACER_MODALFILTER)

    call this%modal_filter_3d%Init( &
        refElem3D,                           & ! (in)
        MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
        MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)

    call element_operation%Setup_ModalFilter_tracer( &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)

    return
  end subroutine setup_modalfilter

end module scale_atm_dyn_dgm_driver_trcadv3d