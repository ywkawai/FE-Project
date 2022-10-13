!> module ATMOSPHERE phyics / turbulence process
!!
!! @par Description
!!          Module for turbulence process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_tb
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
  use scale_tracer, only: QA

  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
    
  use scale_meshfield_base, only: MeshFieldBase
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc

  use scale_atm_dyn_dgm_bnd, only: AtmDynBnd
  
  use scale_atm_phy_tb_dgm_smg, only: &
    atm_phy_tb_dgm_smg_Init,          &
    atm_phy_tb_dgm_smg_Final,         &
    atm_phy_tb_dgm_smg_cal_grad,      &
    atm_phy_tb_dgm_smg_cal_tend,      &
    atm_phy_tb_dgm_smg_cal_grad_qtrc, &
    atm_phy_tb_dgm_smg_cal_tend_qtrc

  use mod_atmos_phy_tb_vars, only: AtmosPhyTbVars

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, extends(ModelComponentProc), public :: AtmosPhyTb
    integer :: TB_TYPEID
    type(AtmosPhyTbVars) :: vars
    type(AtmDynBnd), pointer :: dyn_bnd
  contains
    procedure, public :: setup => AtmosPhyTb_setup 
    procedure, public :: calc_tendency => AtmosPhyTb_calc_tendency
    procedure, public :: update => AtmosPhyTb_update
    procedure, public :: finalize => AtmosPhyTb_finalize
    procedure, public :: SetDynBC => AtmosPhyTB_SetDynBC
  end type AtmosPhyTb

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  integer, parameter :: TB_TYPEID_SMAGORINSKY  = 1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
    
contains
  subroutine AtmosPhyTb_setup( this, model_mesh, tm_parent_comp )
    use mod_atmos_mesh, only: AtmosMesh
    use scale_time_manager, only: TIME_manager_component

    implicit none
    class(AtmosPhyTb), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  

    character(len=H_MID) :: TB_TYPE = 'SMAGORINSKY'
    namelist /PARAM_ATMOS_PHY_TB/ &
      TIME_DT,          &
      TIME_DT_UNIT,     &
      TB_TYPE
    
    class(AtmosMesh), pointer     :: atm_mesh
    class(MeshBase), pointer      :: ptr_mesh

    integer :: ierr
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_TB_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_TB_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB)

    !- get mesh --------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )
    select type(model_mesh)
    class is (AtmosMesh)
      atm_mesh => model_mesh
    end select

    !--- Regist this compoent in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_TB', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                       ! (out)

    !- initialize the variables 
    call this%vars%Init( model_mesh )      

    !--- Set the type of turbulence scheme

    select case( TB_TYPE )
    case ('SMAGORINSKY')
      this%TB_TYPEID = TB_TYPEID_SMAGORINSKY
      call atm_phy_tb_dgm_smg_Init( atm_mesh%ptr_mesh )
    case default
      LOG_ERROR("ATMOS_PHY_TB_setup",*) 'Not appropriate names of TB_TYPE in namelist PARAM_ATMOS_PHY_TB. Check!'
      call PRC_abort
    end select

    !--
    this%dyn_bnd => null()

    return
  end subroutine AtmosPhyTb_setup


  subroutine AtmosPhyTb_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )


    use scale_tracer, only: &
      TRACER_ADVC   

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshQTRCVar,    &
      AtmosVars_GetLocalMeshPhyAuxVars, &
      AtmosVars_GetLocalMeshPhyTends
    use mod_atmos_phy_tb_vars, only:          &
      AtmosPhyTbVars_GetLocalMeshFields_tend,    &
      AtmosPhyTbVars_GetLocalMeshFields_aux,     &
      AtmosPhyTbVars_GetLocalMeshFields_aux_qtrc
    
    implicit none
    
    class(AtmosPhyTb), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, THERM, QTRC
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_P
    type(LocalMeshFieldBaseList) :: RHOQ_tp(QA)
    class(LocalMeshFieldBase), pointer :: tb_MOMX_t, tb_MOMY_t, tb_MOMZ_t, tb_RHOT_t
    class(LocalMeshFieldBase), pointer :: tb_RHOQ_t
    type(LocalMeshFieldBaseList) :: tb_RHOQ_t_list(QA)
    class(LocalMeshFieldBase), pointer :: S11, S12, S22, S23, S31, S33, TKE
    class(LocalMeshFieldBase), pointer :: dPTdx, dPTdy, dPTdz
    class(LocalMeshFieldBase), pointer :: dQTdx, dQTdy, dQTdz
    class(LocalMeshFieldBase), pointer :: Nu, Kh

    integer :: n
    integer :: ke
    integer :: iq

    type DYN_BNDInfo
      logical, allocatable :: is_bound(:,:)
    end type
    type(DYN_BNDInfo), allocatable :: bnd_info(:)

    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    !LOG_INFO('AtmosDyn_tendency',*)

    LOG_PROGRESS(*) 'atmosphere / physics / turbulence'

    call model_mesh%GetModelMesh( mesh )

    if (is_update) then
      allocate( bnd_info(mesh%LOCAL_MESH_NUM) )

      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart('ATM_PHY_TB_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,       &
          DDENS, MOMX, MOMY, MOMZ, THERM,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )      
        call AtmosVars_GetLocalMeshPhyAuxVars( n,      &
          mesh, auxvars_list,                          &
          PRES, PT )
        call AtmosPhyTbVars_GetLocalMeshFields_aux( n, &
          mesh, this%vars%auxvars_manager,             &
          S11, S12, S22, S23, S31, S33, TKE,           &
          dPTdx, dPTdy, dPTdz, Nu, Kh                  )
        call PROF_rapend('ATM_PHY_TB_get_localmesh_ptr', 2)   

        call PROF_rapstart('ATM_PHY_TB_inquire_bnd', 2)
        allocate( bnd_info(n)%is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
        call this%dyn_bnd%Inquire_bound_flag(  bnd_info(n)%is_bound, &
          n, lcmesh%VMapM, lcmesh%VMapP, lcmesh%VMapB,               &
          lcmesh, lcmesh%refElem3D                                   )
        call PROF_rapend('ATM_PHY_TB_inquire_bnd', 2)

        call PROF_rapstart('ATM_PHY_TB_cal_grad', 2)
        select case( this%TB_TYPEID )
        case (TB_TYPEID_SMAGORINSKY)
          call atm_phy_tb_dgm_smg_cal_grad( &
              S11%val, S12%val, S22%val, S23%val, S31%val, S33%val, TKE%val,          &
              dPTdx%val, dPTdy%val, dPTdz%val, Nu%val, Kh%val,                        &
              DDENS%val, MOMX%val, MOMY%val, MOMZ%val, THERM%val,                     &
              DENS_hyd%val, PRES_hyd%val, PRES%val, PT%val,                           &
              model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3), &
              model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3), &
              model_mesh%LiftOptrMat,                                                 &
              lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D,   &
              bnd_info(n)%is_bound                                                    )
        end select
        call PROF_rapend('ATM_PHY_TB_cal_grad', 2)
      end do

      !* Exchange halo data
      call PROF_rapstart('ATM_PHY_TB_exchange_prgv', 2)
      call this%vars%auxvars_manager%MeshFieldComm_Exchange()
      call PROF_rapend('ATM_PHY_TB_exchange_prgv', 2)

      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart('ATM_PHY_TB_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshPrgVars( n, &
          mesh, prgvars_list, auxvars_list,       &
          DDENS, MOMX, MOMY, MOMZ, THERM,         &
          DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
          lcmesh                                  )
        call AtmosVars_GetLocalMeshPhyAuxVars( n,      &
          mesh, auxvars_list,                          &
          PRES, PT )

        call AtmosPhyTbVars_GetLocalMeshFields_tend( n,          &
          mesh, this%vars%tends_manager,                         &
          tb_MOMX_t, tb_MOMY_t, tb_MOMZ_t, tb_RHOT_t             )

        call AtmosPhyTbVars_GetLocalMeshFields_aux( n, &
          mesh, this%vars%auxvars_manager,             &
          S11, S12, S22, S23, S31, S33, TKE,           &
          dPTdx, dPTdy, dPTdz, Nu, Kh                  )
        call PROF_rapend('ATM_PHY_TB_get_localmesh_ptr', 2)   
  
        call PROF_rapstart('ATM_PHY_TB_cal_tend', 2)
        select case( this%TB_TYPEID )
        case (TB_TYPEID_SMAGORINSKY)
          call atm_phy_tb_dgm_smg_cal_tend( &
              tb_MOMX_t%val, tb_MOMY_t%val, tb_MOMZ_t%val, tb_RHOT_t%val,                &
              S11%val, S12%val, S22%val, S23%val, S31%val, S33%val, TKE%val,             &
              dPTdx%val, dPTdy%val, dPTdz%val, Nu%val, Kh%val,                           &
              DDENS%val, MOMX%val, MOMY%val, MOMZ%val, THERM%val,                        &
              DENS_hyd%val, PRES_hyd%val, PRES%val, PT%val,                              &
              model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),    &
              model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3),    &
              model_mesh%LiftOptrMat,                                                    &
              lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D,      &
              bnd_info(n)%is_bound                                                       )
        end select
        call PROF_rapend('ATM_PHY_TB_cal_tend', 2)
      end do

      do iq = 1, QA
        if ( .not. TRACER_ADVC(iq) ) cycle
        
        do n=1, mesh%LOCAL_MESH_NUM

          call PROF_rapstart('ATM_PHY_TB_get_localmesh_ptr', 2)                 
          call AtmosVars_GetLocalMeshQTRCVar( n,       &
            mesh, trcvars_list, iq,                    &
            QTRC, lcmesh                               )
          call AtmosPhyTbVars_GetLocalMeshFields_aux_qtrc( n, &
            mesh, this%vars%auxvars_manager, this%vars%auxtrcvars_manager, &
            dQTdx, dQTdy, dQTdz, Kh                                        )
          call PROF_rapend('ATM_PHY_TB_get_localmesh_ptr', 2)         

          call PROF_rapstart('ATM_PHY_TB_cal_grad_qtrc', 2)
          call atm_phy_tb_dgm_smg_cal_grad_qtrc( &
            dQTdx%val, dQTdy%val, dQTdz%val,                                           & ! (out)
            QTRC%val,                                                                  & ! (in) 
            model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),    & ! (in)
            model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3),    & ! (in)
            model_mesh%LiftOptrMat,                                                    & ! (in)
            lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D,      & ! (in)
            bnd_info(n)%is_bound )                                                       ! (in)
          call PROF_rapend('ATM_PHY_TB_cal_grad_qtrc', 2)
        end do

        call this%vars%auxtrcvars_manager%MeshFieldComm_Exchange()

        do n=1, mesh%LOCAL_MESH_NUM
          call PROF_rapstart('ATM_PHY_TB_get_localmesh_ptr', 2)         
          call AtmosVars_GetLocalMeshPrgVars( n, &
            mesh, prgvars_list, auxvars_list,       &
            DDENS, MOMX, MOMY, MOMZ, THERM,         &
            DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, &
            lcmesh                                  )
          call AtmosPhyTbVars_GetLocalMeshFields_aux_qtrc( n, &
            mesh, this%vars%auxvars_manager, this%vars%auxtrcvars_manager, &
            dQTdx, dQTdy, dQTdz, Kh,                                       &
            this%vars%tends_manager, iq, tb_RHOQ_t,                        &
            lcmesh                                                         )
          call PROF_rapend('ATM_PHY_TB_get_localmesh_ptr', 2)

          call PROF_rapstart('ATM_PHY_TB_cal_tend_qtrc', 2)
          call atm_phy_tb_dgm_smg_cal_tend_qtrc( tb_RHOQ_t%val,                        & ! (out)
            dQTdx%val, dQTdy%val, dQTdz%val,                                           & ! (in)
            Kh%val, DDENS%val, DENS_hyd%val,                                           & ! (in) 
            model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),    & ! (in)
            model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3),    & ! (in)
            model_mesh%LiftOptrMat,                                                    & ! (in)
            lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D,      & ! (in)
            bnd_info(n)%is_bound )                                                       ! (in)
          call PROF_rapend('ATM_PHY_TB_cal_tend_qtrc', 2)
        end do        
      end do

      do n = 1, mesh%LOCAL_MESH_NUM
        deallocate( bnd_info(n)%is_bound )
      end do
    end if
      
    call PROF_rapstart('ATM_PHY_TB_add_tend', 2)
    do n=1, mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p, RHOQ_tp  )

      call AtmosPhyTbVars_GetLocalMeshFields_tend( n,               &
        mesh, this%vars%tends_manager,                              &
        tb_MOMX_t, tb_MOMY_t, tb_MOMZ_t, tb_RHOT_t, tb_RHOQ_t_list, &
        lcmesh                                                      )


      !$omp parallel private(ke, iq)
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        MOMX_tp%val(:,ke) = MOMX_tp%val(:,ke) + tb_MOMX_t%val(:,ke)
        MOMY_tp%val(:,ke) = MOMY_tp%val(:,ke) + tb_MOMY_t%val(:,ke)
        MOMZ_tp%val(:,ke) = MOMZ_tp%val(:,ke) + tb_MOMZ_t%val(:,ke)
        RHOT_tp%val(:,ke) = RHOT_tp%val(:,ke) + tb_RHOT_t%val(:,ke)
      end do
      !$omp end do    
      do iq = 1, QA
        if ( .not. TRACER_ADVC(iq) ) cycle
        !$omp do
        do ke = lcmesh%NeS, lcmesh%NeE
          RHOQ_tp(iq)%ptr%val(:,ke) = RHOQ_tp(iq)%ptr%val(:,ke)        &
                                    + tb_RHOQ_t_list(iq)%ptr%val(:,ke)
        end do
        !$omp end do
      end do
      !$omp end parallel
    end do
    call PROF_rapend('ATM_PHY_TB_add_tend', 2)

    return  
  end subroutine AtmosPhyTb_calc_tendency

  subroutine AtmosPhyTb_update( this, model_mesh, prgvars_list, trcvars_list, auxvars_list, forcing_list, is_update )
    implicit none

    class(AtmosPhyTb), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------

    return
  end subroutine AtmosPhyTb_update

  subroutine AtmosPhyTb_finalize( this )
    implicit none
    class(AtmosPhyTb), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    select case ( this%TB_TYPEID )
    case( TB_TYPEID_SMAGORINSKY )
      call atm_phy_tb_dgm_smg_Final()
    end select

    call this%vars%Final()

    return
  end subroutine AtmosPhyTb_finalize

  subroutine AtmosPhyTb_setDynBC( this, dyn_bnd )
    implicit none
    class(AtmosPhyTb), intent(inout) :: this
    type(AtmDynBnd), intent(in), target :: dyn_bnd
    !--------------------------------------------------

    this%dyn_bnd => dyn_bnd

    return
  end subroutine AtmosPhyTb_setDynBC
  
!- private ------------------------------------------------

end module mod_atmos_phy_tb