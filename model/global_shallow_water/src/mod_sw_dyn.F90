!-------------------------------------------------------------------------------
!> module  dynamical process
!!
!! @par Description
!!          Module for shallow water equation
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_sw_dyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof
  use scale_const, only: &
    UNDEF8 => CONST_UNDEF8

  use scale_sparsemat, only: SparseMat
  use scale_timeint_rk, only: TimeInt_RK
  
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_element_base, only: &
    ElementBase, ElementBase2D
    
  use scale_meshfield_base, only: MeshFieldBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc

  use scale_atm_dyn_dgm_globalsw, only: &
    atm_dyn_dgm_globalsw_Init,          &
    atm_dyn_dgm_globalsw_Final,         &
    atm_dyn_dgm_globalsw_cal_tend

  use scale_element_modalfilter, only: ModalFilter

  use mod_sw_mesh, only: SWMesh
  use mod_sw_vars, only: &
    SWVars_GetLocalMeshPrgVar,  &
    SWVars_GetLocalMeshPrgVars, &
    SW_PROGVARS_NUM,            &
    h_ID  => SW_PROGVARS_h_ID,  &
    U_ID  => SW_PROGVARS_U_ID,  &
    V_ID  => SW_PROGVARS_V_ID,  &
    hs_ID => SW_AUXVARS_hs_ID,  &
    u1_ID => SW_AUXVARS_u1_ID,  &
    u2_ID => SW_AUXVARS_u2_ID

  use mod_sw_bnd, only:SWBnd

  use mod_sw_dyn_vars, only: &
    SWDynVars, &
    SWDynAuxVars_GetLocalMeshFields


  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  abstract interface    
    subroutine atm_dyn_globalsw_cal_tend_ex( &
      h_dt, U_dt, V_dt,                      & ! (out)
      h_, U_, V_, hs_, u1_, u2_, CORIOLIS,   & ! (in)
      Dx, Dy, Sx, Sy, Lift, lmesh, elem      )

      import RP
      import LocalMesh2D
      import elementbase2D
      import SparseMat
      implicit none

      class(LocalMesh2D), intent(in) :: lmesh
      class(elementbase2D), intent(in) :: elem
      type(SparseMat), intent(in) :: Dx, Dy, Sx, Sy, Lift
      real(RP), intent(out) :: h_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: U_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: V_dt(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: h_ (elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: U_ (elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: V_ (elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: hs_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: u1_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: u2_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CORIOLIS(elem%Np,lmesh%NeA)  
    end subroutine atm_dyn_globalsw_cal_tend_ex
  end interface

  type, extends(ModelComponentProc), public :: SWDyn
    integer :: EQS_TYPEID
    type(TimeInt_RK), allocatable :: tint(:)
    type(SWBnd) :: boundary_cond
    type(SWDynVars) :: dyn_vars

    procedure (atm_dyn_globalsw_cal_tend_ex), pointer, nopass :: cal_tend_ex => null()

    ! explicit numerical diffusion
    logical :: CALC_NUMDIFF_FLAG
    integer  :: ND_LAPLACIAN_NUM
    real(RP) :: ND_COEF_H

    ! element-wise modal filter
    logical :: MODALFILTER_FLAG
    type(ModalFilter) :: modal_filter

  contains
    procedure, public :: setup => SWDyn_setup 
    procedure, public :: calc_tendency => SWDyn_calc_tendency
    procedure, public :: update => SWDyn_update
    procedure, public :: finalize => SWDyn_finalize
  end type SWDyn

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------
  
  integer, public, parameter :: EQS_TYPEID_GLOBAL_SHALLOW_WATER    = 1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  ! private :: cal_numfilter_tend
  ! private :: add_phy_tend

  private :: setup_modalfilter
  private :: setup_numdiff
  private :: setup_coriolis_parameter

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

contains

  subroutine SWDyn_setup( this, model_mesh, tm_parent_comp )
    use mod_sw_mesh, only: SWMesh
    use mod_sw_vars, only: SW_PROGVARS_NUM
    use scale_time_manager, only: TIME_manager_component

    implicit none

    class(SWDyn), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    character(len=H_MID) :: EQS_TYPE      = "GLOBAL_SHALLOW_WATER"
    character(len=H_SHORT) :: TINTEG_TYPE = 'RK_TVD_3'
    real(DP) :: TIME_DT                             = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  
    
    logical  :: MODALFILTER_FLAG = .false.
    logical :: NUMDIFF_FLAG      = .false.

    namelist / PARAM_SW_DYN /       &
      EQS_TYPE,                               &
      TINTEG_TYPE,                            &
      TIME_DT,                                &
      TIME_DT_UNIT,                           &
      MODALFILTER_FLAG,                       &
      NUMDIFF_FLAG                           
    
    class(SWMesh), pointer     :: sw_mesh
    class(MeshBase), pointer      :: ptr_mesh
    class(LocalMeshBase), pointer :: ptr_lcmesh
    class(ElementBase2D), pointer :: elem2D
    integer :: n
    real(DP) :: dtsec

    integer :: ierr
    !--------------------------------------------------

    if (.not. this%IsActivated()) return
    LOG_INFO('SWDyn_setup',*)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SW_DYN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("SW_DYN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("SW_DYN_setup",*) 'Not appropriate names in namelist PARAM_SW_DYN. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_SW_DYN)
    
    !- get mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    type is (SWMesh)
      sw_mesh => model_mesh
    end select

    !- Setup the temporal integrator

    call tm_parent_comp%Regist_process( 'SW_DYN', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                 ! (out)

    dtsec = tm_parent_comp%process_list(this%tm_process_id)%dtsec
    
    allocate( this%tint(ptr_mesh%LOCAL_MESH_NUM) )
    do n = 1, ptr_mesh%LOCAL_MESH_NUM
      call ptr_mesh%GetLocalMesh( n, ptr_lcmesh )
      call this%tint(n)%Init( TINTEG_TYPE, dtsec, SW_PROGVARS_NUM, 2, &
        (/ ptr_mesh%refElem%Np, ptr_lcmesh%NeA /) )
    end do

    !- initialize an object to manage boundary conditions
    call this%boundary_cond%Init()

    !- initialize the variables 
    call this%dyn_vars%Init( model_mesh )
    call setup_coriolis_parameter( this%dyn_vars, sw_mesh )

    !- Initialize a module for 3D dynamical core 

    select case(EQS_TYPE)
    case("GLOBAL_SHALLOW_WATER")
      this%EQS_TYPEID = EQS_TYPEID_GLOBAL_SHALLOW_WATER
      call atm_dyn_dgm_globalsw_Init( sw_mesh%mesh )
      this%cal_tend_ex => atm_dyn_dgm_globalsw_cal_tend
    case default
      LOG_ERROR("SW_DYN_setup",*) 'Not appropriate names in namelist PARAM_SW_DYN. Check!'
      call PRC_abort
    end select    

    !- Setup the numerical diffusion
    this%CALC_NUMDIFF_FLAG = NUMDIFF_FLAG
    if( NUMDIFF_FLAG ) call setup_numdiff( this, sw_mesh )

    !- Setup the modal filter
    this%MODALFILTER_FLAG = MODALFILTER_FLAG
    if ( MODALFILTER_FLAG ) call setup_modalfilter( this, sw_mesh )

    return
  end subroutine SWDyn_setup

  subroutine SWDyn_calc_tendency( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
    implicit none
    
    class(SWDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    !LOG_INFO('AtmosDyn_tendency',*)

    return  
  end subroutine SWDyn_calc_tendency

!OCL SERIAL
  subroutine SWDyn_update( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )

    use scale_atm_dyn_dgm_modalfilter, only: &
      atm_dyn_dgm_modalfilter_apply

    implicit none

    class(SWDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    integer :: rkstage
    integer :: tintbuf_ind

    class(MeshBase), pointer :: mesh
    class(MeshBase2D), pointer :: mesh2D    
    class(LocalMesh2D), pointer :: lcmesh
    integer :: n
    integer :: ke

    class(LocalMeshFieldBase), pointer :: h, U, V
    class(LocalMeshFieldBase), pointer :: hs, u1, u2
    class(LocalMeshFieldBase), pointer :: Coriolis

    integer :: vid
    real(RP) :: implicit_fac
    real(RP) :: dt
    character(len=H_SHORT) :: labl
    !--------------------------------------------------
    
    call PROF_rapstart( 'SW_DYN_update', 1)   

    call model_mesh%GetModelMesh( mesh )

    !-
    do rkstage=1, this%tint(1)%nstage

      !* Exchange halo data
      call PROF_rapstart( 'SW_DYN_exchange_prgv', 2)
      call prgvars_list%MeshFieldComm_Exchange()
      call PROF_rapend( 'SW_DYN_exchange_prgv', 2)
  
      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart( 'SW_DYN_get_localmesh_ptr', 2)         
        call SWVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list, &
          h, U, V, hs, u1, u2, lcmesh       )
        call PROF_rapend( 'SW_DYN_get_localmesh_ptr', 2)         
        
        !* Apply boundary conditions
        call PROF_rapstart( 'SW_DYN_applyBC_prgv', 2)
        call PROF_rapend( 'SW_DYN_applyBC_prgv', 2)
      end do

      do n=1, mesh%LOCAL_MESH_NUM
        tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'SW_DYN_get_localmesh_ptr', 2)         
        call SWVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list, &
          h, U, V, hs, u1, u2, lcmesh       )
      
        call SWDynAuxVars_GetLocalMeshFields( n, &
          mesh, this%dyn_vars%AUXVARS2D_manager, &
          Coriolis )
        call PROF_rapend( 'SW_DYN_get_localmesh_ptr', 2)

        call PROF_rapstart( 'SW_DYN_update_caltend_ex', 2)
        call this%cal_tend_ex( &
          this%tint(n)%tend_buf2D_ex(:,:,h_ID,tintbuf_ind),   &
          this%tint(n)%tend_buf2D_ex(:,:,U_ID ,tintbuf_ind),  &
          this%tint(n)%tend_buf2D_ex(:,:,V_ID ,tintbuf_ind),  &
          h%val, U%val, V%val,                                &
          hs%val, u1%val, u2%val,                             &
          Coriolis%val,                                       &
          model_mesh%DOptrMat(1), model_mesh%DOptrMat(2),     &
          model_mesh%SOptrMat(1), model_mesh%SOptrMat(2),     &
          model_mesh%LiftOptrMat,                             &
          lcmesh, lcmesh%refElem2D                            ) 
        call PROF_rapend( 'SW_DYN_update_caltend_ex', 2)
        
        call PROF_rapstart( 'SW_DYN_update_add_tp', 2)
        ! call add_phy_tend( &
        !   this, this%tint(n)%tend_buf2D_ex(:,:,:,tintbuf_ind), & ! (inout)
        !   forcing_list,                                        & ! (in)
        !   mesh, n, lcmesh, lcmesh%refElem2D                    ) ! (in)
        call PROF_rapend( 'SW_DYN_update_add_tp', 2)

        call PROF_rapstart( 'SW_DYN_update_advance', 2)      
        call this%tint(n)%Advance( rkstage, h%val, h_ID, &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call this%tint(n)%Advance( rkstage, U%val, U_ID,         &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call this%tint(n)%Advance( rkstage, V%val, V_ID,         &
                    1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        !$omp parallel do
        do ke=lcmesh%NeS, lcmesh%NeE
          u1%val(:,ke) = lcmesh%G_ij(:,ke,1,1) * U%val(:,ke) + lcmesh%G_ij(:,ke,1,2) * V%val(:,ke)
          u2%val(:,ke) = lcmesh%G_ij(:,ke,2,1) * U%val(:,ke) + lcmesh%G_ij(:,ke,2,2) * V%val(:,ke)
        end do
        call PROF_rapend( 'SW_DYN_update_advance', 2)

        !------------------------------------------------------------------------------
      end do
    end do

    !-- modal filter
    if ( this%MODALFILTER_FLAG ) then
      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart( 'SW_DYN_get_localmesh_ptr', 2)         
        call SWVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list, &
          h, U, V, hs, u1, u2, lcmesh       )
        call PROF_rapend( 'SW_DYN_get_localmesh_ptr', 2)

        call PROF_rapstart( 'SW_DYN_update_expfilter', 2)
        call atm_dyn_dgm_modalfilter_apply(              & ! (inout)
          h%val, U%val, V%val, u1%val, u2%val,           & ! (in)
          lcmesh, lcmesh%refElem2D, this%modal_filter    ) ! (in)
        call PROF_rapend( 'SW_DYN_update_expfilter', 2)
      end do
    end if

    !-- numerical diffusion
    if ( this%CALC_NUMDIFF_FLAG ) then
      
      call PROF_rapstart( 'ATM_DYN_numfilter', 1)
      call prgvars_list%MeshFieldComm_Exchange()

      do vid = 1, SW_PROGVARS_NUM
        ! call cal_numfilter_tend( this, model_mesh, prgvars_list, auxvars_list, vid )
      end do

      do n=1, mesh%LOCAL_MESH_NUM
        dt = this%tint(n)%Get_deltime()

        call SWVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list, &
          h, U, V, hs, u1, u2, lcmesh       )
        
        !$omp parallel do
        do ke=1, lcmesh%Ne
         h%val(:,ke) = h%val(:,ke) + dt * this%tint(n)%tend_buf2D_ex(:,ke,h_ID,1)
         U%val(:,ke) = U%val(:,ke) + dt * this%tint(n)%tend_buf2D_ex(:,ke,U_ID,1)
         V%val(:,ke) = V%val(:,ke) + dt * this%tint(n)%tend_buf2D_ex(:,ke,V_ID,1)
        end do
      end do

      call PROF_rapend( 'SW_DYN_numfilter', 2)
    end if

    call PROF_rapend( 'SW_DYN_update', 1)

    return  
  end subroutine SWDyn_update

  subroutine SWDyn_finalize( this )
    implicit none
    class(SWDyn), intent(inout) :: this

    integer :: n
    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    LOG_INFO('SWDyn_finalize',*)

    select case(this%EQS_TYPEID)
    case(EQS_TYPEID_GLOBAL_SHALLOW_WATER)
      call atm_dyn_dgm_globalsw_Final()
    end select 

    ! if (this%CALC_NUMDIFF_FLAG) then
    !   call atm_dyn_dgm_nonhydro3d_numdiff_Final()
    ! end if

    if (this%MODALFILTER_FLAG) then
      call this%modal_filter%Final()
    end if
    
    do n = 1, size(this%tint)
      call this%tint(n)%Final()
    end do
    deallocate( this%tint )

    call this%boundary_cond%Final()
    call this%dyn_vars%Final()

    return  
  end subroutine SWDyn_finalize  

  !--- private ---------------

!OCL SERIAL
  ! subroutine add_phy_tend( this,      & ! (in)
  !   dyn_tends,                        & ! (inout)
  !   phytends_list,                    & ! (in)
  !   mesh, domID, lcmesh, elem2D       ) ! (in)

  !   use scale_const, only: &
  !     Rdry => CONST_Rdry,   &
  !     CPdry => CONST_CPdry, &
  !     CVdry => CONST_CVdry, &
  !     PRES00 => CONST_PRE00
    
  !   use mod_sw_vars, only: &
  !     SWVars_GetLocalMeshPhyTends

  !   implicit none

  !   class(SWDyn), intent(inout) :: this
  !   class(LocalMesh2D), intent(in) :: lcmesh
  !   class(elementbase2D), intent(in) :: elem2D
  !   real(RP), intent(inout) :: dyn_tends(elem2D%Np,lcmesh%NeA,SW_PROGVARS_NUM)
  !   class(ModelVarManager), intent(inout) :: phytends_list
  !   class(MeshBase), intent(in) :: mesh
  !   integer, intent(in) :: domID

  !   class(LocalMeshFieldBase), pointer :: h_tp, U_tp, V_tp
  !   integer :: ke
  !   !---------------------------------------------------------------------------------

  !   call SWVars_GetLocalMeshPhyTends( domID, mesh, phytends_list, & ! (in)
  !     DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p            ) ! (out)

  !   !$omp parallel do          &
  !   !$Omp private( RHOT, EXNER )
  !   do ke=lcmesh%NeS, lcmesh%NeE
  !     RHOT(:) = P0ovR * (PRES_hyd(:,ke) * rP0)**rgamm + DRHOT(:,ke)
  !     EXNER(:) = (RovP0 * RHOT(:))**RovCv

  !     dyn_tends(:,ke,DDENS_ID) = dyn_tends(:,ke,DDENS_ID) + DENS_tp%val(:,ke)
  !     dyn_tends(:,ke,MOMX_ID ) = dyn_tends(:,ke,MOMX_ID ) + MOMX_tp%val(:,ke)
  !     dyn_tends(:,ke,MOMY_ID ) = dyn_tends(:,ke,MOMY_ID ) + MOMY_tp%val(:,ke)
  !     dyn_tends(:,ke,MOMZ_ID ) = dyn_tends(:,ke,MOMZ_ID ) + MOMZ_tp%val(:,ke)
  !     dyn_tends(:,ke,DRHOT_ID) = dyn_tends(:,ke,DRHOT_ID) + RHOT_tp%val(:,ke) &
  !                              + RHOH_p %val(:,ke) / ( CpDry * EXNER(:) )
  !   end do

  !   return
  ! end subroutine add_phy_tend

!OCL SERIAL
  ! subroutine cal_numfilter_tend( this, model_mesh, prgvars_list, auxvars_list, varid )

  !   use mod_atmos_dyn_vars, only: &
  !     AtmosDynAuxVars_GetLocalMeshFields,     &
  !     AtmosDynNumDiffFlux_GetLocalMeshFields, &
  !     AtmosDynNumDiffTend_GetLocalMeshFields
    
  !   use scale_atm_dyn_dgm_nonhydro3d_numdiff, only: &
  !     atm_dyn_dgm_nonhydro3d_numdiff_tend,          &
  !     atm_dyn_dgm_nonhydro3d_numdiff_cal_laplacian, &
  !     atm_dyn_dgm_nonhydro3d_numdiff_cal_flx

  !   implicit none
          
  !   class(AtmosDyn), intent(inout) :: this
  !   class(ModelMeshBase), intent(in) :: model_mesh
  !   class(ModelVarManager), intent(inout) :: prgvars_list
  !   class(ModelVarManager), intent(inout) :: auxvars_list
  !   integer, intent(in) :: varid

  !   class(LocalMeshFieldBase), pointer :: var
  !   class(LocalMeshFieldBase), pointer :: ND_flx_x, ND_flx_y, ND_flx_z
  !   class(LocalMeshFieldBase), pointer :: ND_lapla_h, ND_lapla_v
  !   class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
  !   class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd

  !   class(MeshBase), pointer :: mesh
  !   class(LocalMesh3D), pointer :: lcmesh
  !   integer :: n
  !   integer :: ke
  !   integer :: nd_itr
  !   real(RP) :: nd_sign
  !   logical :: dens_weight_flag
  !   logical, allocatable :: is_bound(:,:)

  !   !-----------------------------------------

  !   nd_sign = (-1)**(mod(this%ND_LAPLACIAN_NUM+1,2))
  !   dens_weight_flag = (varid /= DDENS_ID)

  !   call model_mesh%GetModelMesh( mesh )

  !   do n=1, mesh%LOCAL_MESH_NUM
  !     call AtmosVars_GetLocalMeshPrgVar( n, mesh, prgvars_list, auxvars_list, &
  !       varid, var,                                                           &
  !       DENS_hyd, PRES_hyd, lcmesh                                            )
  !     call AtmosVars_GetLocalMeshPrgVars( n, mesh, prgvars_list, auxvars_list, &
  !       DDENS, MOMX, MOMY, MOMZ, DRHOT,                                        &
  !       DENS_hyd, PRES_hyd                                                     )
  !     call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
  !       ND_flx_x, ND_flx_y, ND_flx_z )
      
  !     allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
  !     call this%boundary_cond%ApplyBC_numdiff_even_lc( var%val, is_bound, varid, n, &
  !       MOMX%val, MOMY%val, MOMZ%val, DENS_hyd%val, PRES_hyd%val,                    &
  !       lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),   &
  !       lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )
      
  !     call atm_dyn_dgm_nonhydro3d_numdiff_cal_flx( ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, &
  !       var%val, var%val, DDENS%val, DENS_hyd%val,                                       &
  !       model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),          &
  !       model_mesh%LiftOptrMat,                                                          &
  !       lcmesh, lcmesh%refElem3D, is_bound, dens_weight_flag ) 

  !     deallocate( is_bound )
  !   end do

  !   !* Exchange halo data
  !   call this%dyn_vars%NUMDIFF_FLUX_manager%MeshFieldComm_Exchange()

  !   do nd_itr=1, this%ND_LAPLACIAN_NUM-1
  !     do n = 1, mesh%LOCAL_MESH_NUM
  !       call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
  !         ND_flx_x, ND_flx_y, ND_flx_z, lcmesh)
  !       call AtmosDynNumDiffTend_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_TEND_manager, &
  !         ND_lapla_h, ND_lapla_v )
          
  !       allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
  !       call this%boundary_cond%ApplyBC_numdiff_odd_lc( &
  !         ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, is_bound, varid, n,              &
  !         lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
  !         lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )

  !       call atm_dyn_dgm_nonhydro3d_numdiff_cal_laplacian( ND_lapla_h%val, ND_lapla_v%val, &
  !         ND_flx_x%val, ND_flx_y%val, ND_flx_z%val,                                        &
  !         model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),          &
  !         model_mesh%LiftOptrMat,                                                          &
  !         lcmesh, lcmesh%refElem3D, is_bound )
        
  !       deallocate( is_bound )
  !     end do
  !     !* Exchange halo data
  !     call this%dyn_vars%NUMDIFF_TEND_manager%MeshFieldComm_Exchange()

  !     do n = 1, mesh%LOCAL_MESH_NUM
  !       call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
  !         ND_flx_x, ND_flx_y, ND_flx_z, lcmesh)
  !       call AtmosDynNumDiffTend_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_TEND_manager, &
  !         ND_lapla_h, ND_lapla_v )    
  !       call AtmosVars_GetLocalMeshPrgVars( n, mesh, prgvars_list, auxvars_list, &
  !         DDENS, MOMX, MOMY, MOMZ, DRHOT,                                        &
  !         DENS_hyd, PRES_hyd                                                     )
          
  !       allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
  !       call this%boundary_cond%ApplyBC_numdiff_even_lc( &
  !         ND_lapla_h%val, is_bound, varid, n,                                        &
  !         MOMX%val, MOMY%val, MOMZ%val, DENS_hyd%val, PRES_hyd%val,                  &
  !         lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
  !         lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )
          
  !       call atm_dyn_dgm_nonhydro3d_numdiff_cal_flx( ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, &
  !         ND_lapla_h%val, ND_lapla_v%val, DDENS%val, DENS_hyd%val,                             &
  !         model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),              &
  !         model_mesh%LiftOptrMat,                                                              &
  !         lcmesh, lcmesh%refElem3D, is_bound, .false. ) 

  !       deallocate( is_bound )
  !     end do
  !     !* Exchange halo data
  !     call this%dyn_vars%NUMDIFF_FLUX_manager%MeshFieldComm_Exchange()
  !   end do
    
  !   do n = 1, mesh%LOCAL_MESH_NUM

  !     call AtmosDynNumDiffFlux_GetLocalMeshFields( n, mesh, this%dyn_vars%NUMDIFF_FLUX_manager, &
  !       ND_flx_x, ND_flx_y, ND_flx_z, lcmesh)  
  !     call AtmosVars_GetLocalMeshPrgVar( n, mesh, prgvars_list, auxvars_list,  &
  !       varid, var,                                                            &
  !       DENS_hyd, PRES_hyd, lcmesh                                             )
  !     call AtmosVars_GetLocalMeshPrgVar( n, mesh, prgvars_list, auxvars_list,  &
  !       DDENS_ID, DDENS                                                        )

  !     allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
  !     call this%boundary_cond%ApplyBC_numdiff_odd_lc(                              &
  !       ND_flx_x%val, ND_flx_y%val, ND_flx_z%val, is_bound, varid, n,              &
  !       lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
  !       lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )

  !     call atm_dyn_dgm_nonhydro3d_numdiff_tend( this%tint(n)%tend_buf2D_ex(:,:,varid,1),  &
  !       ND_flx_x%val, ND_flx_y%val, ND_flx_z%val,                                         &
  !       DDENS%val, DENS_hyd%val, nd_sign * this%ND_COEF_H, nd_sign * this%ND_COEF_V,      &
  !       model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),           &
  !       model_mesh%LiftOptrMat,                                                           &
  !       lcmesh, lcmesh%refElem3D, is_bound, dens_weight_flag ) 

  !     deallocate( is_bound )
  !   end do

  !   return
  ! end subroutine cal_numfilter_tend

  !-- Setup modal filter
!OCL SERIAL
  subroutine setup_modalfilter( this, sw_mesh )
    implicit none

    class(SWDyn), target, intent(inout) :: this
    class(SWMesh), target, intent(in) :: sw_mesh

    real(RP) :: MF_ETAC  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA = 36.0_RP
    integer  :: MF_ORDER = 16

    namelist /PARAM_SW_DYN_MODALFILTER/ &
      MF_ETAC, MF_ALPHA, MF_ORDER

    integer :: ierr
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SW_DYN_MODALFILTER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("SW_DYN_setup_modalfilter",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("SW_DYN_setup_modalfilter",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_MODALFILTER. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_SW_DYN_MODALFILTER)      

    call sw_mesh%Construct_ModalFilter2D( &
      this%modal_filter,                  & ! (inout)
      MF_ETAC, MF_ALPHA, MF_ORDER         ) ! (in)

    return
  end subroutine setup_modalfilter

  !-- Setup explicit numerical diffusion
!OCL SERIAL
  subroutine setup_numdiff( this, sw_mesh )
    implicit none

    class(SWDyn), target, intent(inout) :: this
    class(SWMesh), target, intent(in) :: sw_mesh

    integer ::  ND_LAPLACIAN_NUM = 1
    real(RP) :: ND_COEF_h        = 0.0_RP

    namelist /PARAM_SW_DYN_NUMDIFF/ &
      ND_LAPLACIAN_NUM,                &
      ND_COEF_h

    integer :: ierr
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SW_DYN_NUMDIFF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_SW_setup_numdiff",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_SW_setup_numdiff",*) 'Not appropriate names in namelist PARAM_SW_DYN_NUMDIFF. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_SW_DYN_NUMDIFF)      

    this%ND_LAPLACIAN_NUM = ND_LAPLACIAN_NUM
    this%ND_COEF_H = ND_COEF_h
!    call atm_dyn_dgm_nonhydro3d_numdiff_Init( sw_mesh%mesh )

    return
  end subroutine setup_numdiff

  !-- Setup Coriolis parameter

!OCL SERIAL
  subroutine setup_coriolis_parameter( this, sw_mesh )

    use scale_coriolis_param, only: &
      get_coriolis_parameter
    implicit none

    class(SWDynVars), target, intent(inout) :: this
    class(SWMesh), target, intent(in) :: sw_mesh

    class(LocalMeshFieldBase), pointer :: coriolis
    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: n
    !---------------------------------------------------------------

    do n = 1, sw_mesh%mesh%LOCAL_MESH_NUM
      call SWDynAuxVars_GetLocalMeshFields( n, sw_mesh%mesh, this%AUXVARS2D_manager, &
        coriolis, lcmesh2D )

      call get_coriolis_parameter( &
        coriolis%val(:,lcmesh2D%NeS:lcmesh2D%NeE),                       & ! (out)
        "SPHERE", lcmesh2D%refElem2D%Np * lcmesh2D%Ne,                   & ! (in)
        lat=lcmesh2D%lat(:,:)                                            ) ! (in)
    end do

    return
  end subroutine setup_coriolis_parameter

!--------

!   subroutine cal_MOMZ_tend( &
!     MOMZ_t, MOMZ_t_advx, MOMZ_t_advY, MOMZ_t_advZ, MOMZ_t_lift, MOMZ_t_buoy,     & ! (out)
!      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                    & ! (in)
!      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D )
 
!      use scale_element_base
!      use scale_sparsemat
!      use scale_const, only: &
!       GRAV => CONST_GRAV,  &
!       Rdry => CONST_Rdry,  &
!       CPdry => CONST_CPdry, &
!       CVdry => CONST_CVdry, &
!       PRES00 => CONST_PRE00
!      use scale_atm_dyn_nonhydro3d, only: IntrpMat_VPOrdM1
!      implicit none
 
!      class(LocalMesh3D), intent(in) :: lmesh
!      class(elementbase3D), intent(in) :: elem
!      class(LocalMesh2D), intent(in) :: lmesh2D
!      class(elementbase2D), intent(in) :: elem2D
!      type(SparseMat), intent(in) :: Dx, Dy, Dz, Sx, Sy, Sz, Lift
!      real(RP), intent(out) :: MOMZ_t(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_advx(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_advy(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_advz(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_lift(elem%Np,lmesh%NeA)
!      real(RP), intent(out) :: MOMZ_t_buoy(elem%Np,lmesh%NeA)

!      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
!      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
 
!      real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
!      real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
!      real(RP) :: dens_(elem%Np), RHOT_(elem%Np), dpres_(elem%Np)
!      real(RP) :: pres_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np)
 
!      integer :: ke
!      !------------------------------------------------------------------------
 
!      call cal_del_flux_dyn( del_flux,                                          & ! (out)
!        DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,                & ! (in)
!        lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
!        lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
!        lmesh, elem )                                                             ! (in)
  
!      !-----
!      !$omp parallel do private(RHOT_,pres_,dpres_,dens_,u_,v_,w_,Fx,Fy,Fz,LiftDelFlx)
!      do ke = lmesh%NeS, lmesh%NeE
!        !--
 
!        RHOT_(:) = PRES00/Rdry * (PRES_hyd(:,ke)/PRES00)**(CVdry/CPdry) + DRHOT_(:,ke)
!        pres_(:) = PRES00 * (Rdry*RHOT_(:)/PRES00)**(CPdry/Cvdry)
!        dpres_(:) = pres_(:) - PRES_hyd(:,ke)
!        dens_(:) = DDENS_(:,ke) + DENS_hyd(:,ke)
 
!        u_(:) = MOMX_(:,ke)/dens_(:)
!        v_(:) = MOMY_(:,ke)/dens_(:)
!        w_(:) = MOMZ_(:,ke)/dens_(:)
 
!        !-- MOMZ
!        call sparsemat_matmul(Dx, u_(:)*MOMZ_(:,ke), Fx)
!        call sparsemat_matmul(Dy, v_(:)*MOMZ_(:,ke), Fy)
!        call sparsemat_matmul(Dz, w_(:)*MOMZ_(:,ke), Fz)
!        MOMZ_t_advx(:,ke) = - lmesh%Escale(:,ke,1,1) * Fx(:)
!        MOMZ_t_advy(:,ke) = - lmesh%Escale(:,ke,2,2) * Fy(:)
!        MOMZ_t_advz(:,ke) = - lmesh%Escale(:,ke,3,3) * Fz(:)

!        call sparsemat_matmul(Dz, dpres_(:), Fz)
!        MOMZ_t_buoy(:,ke) = - lmesh%Escale(:,ke,3,3) * Fz(:) &
!                            - matmul(IntrpMat_VPOrdM1, DDENS_(:,ke)) * Grav

!        call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux(:,ke), LiftDelFlx)
!        MOMZ_t_lift(:,ke) = - LiftDelFlx(:)

!        MOMZ_t(:,ke) = MOMZ_t_advx(:,ke) +  MOMZ_t_advy(:,ke) +  MOMZ_t_advz(:,ke) &
!                     + MOMZ_t_lift(:,ke) + MOMZ_t_buoy(:,ke)
!      end do
 
!      return
!  end subroutine cal_MOMZ_tend


!  subroutine cal_del_flux_dyn( del_flux, &
!    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,   &
!    nx, ny, nz, vmapM, vmapP, lmesh, elem )

!    use scale_const, only: &
!     GRAV => CONST_GRAV,  &
!     Rdry => CONST_Rdry,  &
!     CPdry => CONST_CPdry, &
!     CVdry => CONST_CVdry, &
!     PRES00 => CONST_PRE00
  
!    implicit none

!    class(LocalMesh3D), intent(in) :: lmesh
!    class(elementbase3D), intent(in) :: elem  
!    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
!    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
!    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
!    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
!    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
!    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
!    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
!    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
!    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
!    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
!    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
!    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
!    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
   
!    integer :: i, iP, iM
!    real(RP) :: VelP, VelM, alpha
!    real(RP) :: uM, uP, vM, vP, wM, wP, presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P
!    real(RP) :: gamm, rgamm
!    !------------------------------------------------------------------------

!    gamm = CpDry/CvDry
!    rgamm = CvDry/CpDry

!    !$omp parallel do private( &
!    !$omp iM, iP, uM, VelP, VelM, alpha, &
!    !$omp uP, vM, vP, wM, wP, presM, presP, dpresM, dpresP, densM, densP, rhotM, rhotP, rhot_hyd_M, rhot_hyd_P)
!    do i=1, elem%NfpTot*lmesh%Ne
!      iM = vmapM(i); iP = vmapP(i)

!      rhot_hyd_M = PRES00/Rdry * (PRES_hyd(iM)/PRES00)**rgamm
!      rhot_hyd_P = PRES00/Rdry * (PRES_hyd(iP)/PRES00)**rgamm
     
!      rhotM = rhot_hyd_M + DRHOT_(iM)
!      presM = PRES00 * (Rdry*rhotM/PRES00)**gamm
!      dpresM = presM - PRES_hyd(iM)*abs(nz(i))

!      rhotP = rhot_hyd_P + DRHOT_(iP) 
!      presP = PRES00 * (Rdry*rhotP/PRES00)**gamm
!      dpresP = presP - PRES_hyd(iP)*abs(nz(i))

!      densM = DDENS_(iM) + DENS_hyd(iM)
!      densP = DDENS_(iP) + DENS_hyd(iP)

!      VelM = (MOMX_(iM)*nx(i) + MOMY_(iM)*ny(i) + MOMZ_(iM)*nz(i))/densM
!      VelP = (MOMX_(iP)*nx(i) + MOMY_(iP)*ny(i) + MOMZ_(iP)*nz(i))/densP

!      alpha = max( sqrt(gamm*presM/densM) + abs(VelM), sqrt(gamm*presP/densP) + abs(VelP)  )


!      del_flux(i) = 0.5_RP*(                &
!                    ( MOMZ_(iP)*VelP - MOMZ_(iM)*VelM)   &
!                    + ( dpresP - dpresM )*nz(i)          &                    
!                    - alpha*(MOMZ_(iP) - MOMZ_(iM))      )
!    end do

!    return
!  end subroutine cal_del_flux_dyn
end module mod_sw_dyn