!> module ATMOSPHERE phyics / surface process
!!
!! @par Description
!!          Module for surface process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_sfc
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

  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
    
  use scale_meshfield_base, only: MeshFieldBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc
  
  use mod_atmos_phy_sfc_vars, only: AtmosPhySfcVars

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, extends(ModelComponentProc), public :: AtmosPhySfc
    integer :: SFCFLX_TYPEID
    type(AtmosPhySfcVars) :: vars
  contains
    procedure, public :: setup => AtmosPhySfc_setup 
    procedure, public :: calc_tendency => AtmosPhySfc_calc_tendency
    procedure, public :: update => AtmosPhySfc_update
    procedure, public :: finalize => AtmosPhySfc_finalize
  end type AtmosPhySfc

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !-----------------------------------------------------------------------------

  integer, parameter :: SFCFLX_TYPEID_CONST           = 1

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
  subroutine AtmosPhySfc_setup( this, model_mesh, tm_parent_comp )
    use mod_atmos_mesh, only: AtmosMesh
    use mod_atmos_vars, only: ATMOS_PROGVARS_NUM
    use scale_time_manager, only: TIME_manager_component

    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_setup

    implicit none
    class(AtmosPhySfc), intent(inout) :: this
    class(ModelMeshBase), target, intent(in) :: model_mesh
    class(TIME_manager_component), intent(inout) :: tm_parent_comp

    real(DP) :: TIME_DT                             = UNDEF8
    character(len=H_SHORT) :: TIME_DT_UNIT          = 'SEC'  

    character(len=H_MID) :: SFCFLX_TYPE = "CONST"
    namelist /PARAM_ATMOS_PHY_SFC/ &
      TIME_DT,          &
      TIME_DT_UNIT,     &
      SFCFLX_TYPE
    
    integer :: ierr
    !--------------------------------------------------

    if (.not. this%IsActivated()) return

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_SFC_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SFC,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_SFC_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_PHY_SFC_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_SFC. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_SFC)

    !--- Regist this compoent in the time manager
    
    call tm_parent_comp%Regist_process( 'ATMOS_PHY_SFC', TIME_DT, TIME_DT_UNIT, & ! (in)
      this%tm_process_id )                                                        ! (out)

    !- initialize the variables 
    call this%vars%Init( model_mesh )      

    !--- Set the type of surface flux scheme

    select case( SFCFLX_TYPE )
    case ('CONST')
      call ATMOS_PHY_SF_const_setup()
      this%SFCFLX_TYPEID = SFCFLX_TYPEID_CONST
    case default
      LOG_ERROR("ATMOS_PHY_SFC_setup",*) 'Not appropriate names of SFCFLX_TYPE in namelist PARAM_PHY_SFC. Check!'
      call PRC_abort
    end select

    return
  end subroutine AtmosPhySfc_setup


  subroutine AtmosPhySfc_calc_tendency( &
    this, model_mesh, prgvars_list, trcvars_list, &
    auxvars_list, forcing_list, is_update         )

    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,   &
      AtmosVars_GetLocalMeshPhyTends,  &
      AtmosVars_GetLocalMeshPhyAuxVars
    use mod_atmos_phy_sfc_vars, only: &
      AtmosPhySfcVars_GetLocalMeshFields
    
    implicit none
    
    class(AtmosPhySfc), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update

    class(MeshBase), pointer :: mesh
    class(LocalMesh3D), pointer :: lcmesh

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOH_p
    class(LocalMeshFieldBase), pointer :: SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH, SFLX_LH

    integer :: n
    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    !LOG_INFO('AtmosDyn_tendency',*)

    call model_mesh%GetModelMesh( mesh )
    do n=1, mesh%LOCAL_MESH_NUM
      call PROF_rapstart( 'ATM_PHY_SFC_get_localmesh_ptr', 2)         
      call AtmosVars_GetLocalMeshPrgVars( n,  &
        mesh, prgvars_list, auxvars_list,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,       &
        DENS_hyd, PRES_hyd, lcmesh            )
      
      call AtmosVars_GetLocalMeshPhyTends( n,        &
        mesh, forcing_list,                          &
        DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, &
        RHOH_p  )
      call AtmosVars_GetLocalMeshPhyAuxVars( n,      &
        mesh, auxvars_list,                          &
        PRES, PT )

      call AtmosPhySfcVars_GetLocalMeshFields( n,    &
        mesh, this%vars%SFCFLX_manager,              &
        SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH, SFLX_LH  )
      call PROF_rapend( 'ATM_PHY_SFC_get_localmesh_ptr', 2)   

      call PROF_rapstart( 'ATM_PHY_SFC_cal_tend', 2)
      call cal_tend_from_sfcflx( this, is_update,                               &
        DENS_tp%val, MOMX_tp%val, MOMY_tp%val, MOMZ_tp%val, RHOH_p%val,         &
        SFLX_MU%val, SFLX_MV%val, SFLX_MW%val, SFLX_SH%val, SFLX_LH%val,        &
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                     &
        DENS_hyd%val, PRES_hyd%val,                                             &
        PRES%val, PT%val,                                                       &
        model_mesh%DOptrMat(3), model_mesh%SOptrMat(3), model_mesh%LiftOptrMat, &
        lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D    )
      call PROF_rapend( 'ATM_PHY_SFC_cal_tend', 2)
    end do

    return  
  end subroutine AtmosPhySfc_calc_tendency

  subroutine AtmosPhySfc_update( this, model_mesh,           &
    prgvars_list, trcvars_list, auxvars_list,  forcing_list, &
    is_update )
    implicit none

    class(AtmosPhySfc), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: trcvars_list    
    class(ModelVarManager), intent(inout) :: auxvars_list
    class(ModelVarManager), intent(inout) :: forcing_list
    logical, intent(in) :: is_update
    !--------------------------------------------------

    return
  end subroutine AtmosPhySfc_update

  subroutine AtmosPhySfc_finalize( this )
    implicit none
    class(AtmosPhySfc), intent(inout) :: this

    !--------------------------------------------------
    if (.not. this%IsActivated()) return

    call this%vars%Final()

    return
  end subroutine AtmosPhySfc_finalize

!- private ------------------------------------------------

!OCL SERIAL  
  subroutine cal_tend_from_sfcflx( this, is_update_sflx, &
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOH_p,          &
    SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH, SFLX_LH,         &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,                      &
    DENS_hyd, PRES_hyd,                                  &
    PRES, PT,                                            &
    Dz, Sz, Lift,                                        &
    lcmesh, elem, lcmesh2D, elem2D )

    use scale_const, only: &
      Rdry => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
    
    use scale_atmos_phy_sf_const, only: &
       ATMOS_PHY_SF_const_flux
 
    use scale_sparsemat  

    implicit none

    class(AtmosPhySfc), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(elementbase2D), intent(in) :: elem2D
    logical, intent(in) :: is_update_sflx
    real(RP), intent(inout) :: DENS_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMZ_tp(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOH_p (elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: SFLX_MU(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_MV(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_MW(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_SH(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(inout) :: SFLX_LH(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PT(elem%Np,lcmesh%NeA)
    type(SparseMat), intent(in) :: Dz, Sz, Lift

    real(RP) :: ATM_W   (elem2D%Np,lcmesh2D%NeA)
    real(RP) :: ATM_U   (elem2D%Np,lcmesh2D%NeA)
    real(RP) :: ATM_V   (elem2D%Np,lcmesh2D%NeA)
    real(RP) :: Z1      (elem2D%Np,lcmesh2D%NeA)
    real(RP) :: SFC_TEMP(elem2D%Np,lcmesh2D%NeA)
    real(RP) :: SFC_DENS(elem2D%Np,lcmesh2D%NeA)

    ! dummy
    real(RP) :: SFLX_QV(elem2D%Np,lcmesh2D%NeA)
    real(RP) :: U10(elem2D%Np,lcmesh2D%NeA)
    real(RP) :: V10(elem2D%Np,lcmesh2D%NeA)

    integer :: ke
    integer :: ke2D
    integer :: hsliceZ0, hsliceZ1
    integer :: ij
    real(RP) :: dens
    real(RP) :: LiftDelFlx(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lcmesh%Ne,4)

    !-------------------------------------------------

    if (is_update_sflx) then
      !$omp parallel do collapse(2) private( &
      !$omp ke, hSliceZ0, hsliceZ1,          &
      !$omp dens                             )
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
      do ij=1, elem2D%Np

        ke = ke2D
        hsliceZ0 = elem%Hslice(ij,1)
        hsliceZ1 = elem%Hslice(ij,2)

        dens = DENS_hyd(hsliceZ1,ke) + DDENS(hsliceZ1,ke)
        ATM_U(ij,ke2D) = MOMX(hsliceZ1,ke) / dens
        ATM_V(ij,ke2D) = MOMY(hsliceZ1,ke) / dens
        ATM_W(ij,ke2D) = 0.0_RP

        SFC_DENS(ij,ke2D) = DENS_hyd(hsliceZ0,ke) + DDENS(hsliceZ0,ke)
        SFC_TEMP(ij,ke2D) = PRES(hsliceZ0,ke) / ( Rdry * SFC_DENS(ij,ke2D) )

        Z1(ij,ke2D) = lcmesh%pos_en(hsliceZ1,ke,3)
      end do
      end do

      select case ( this%SFCFLX_TYPEID )
      case ( SFCFLX_TYPEID_CONST )
        call ATMOS_PHY_SF_const_flux( &
          elem2D%Np, 1, elem2D%Np, lcmesh2D%NeA, 1, lcmesh2D%Ne, & ! [IN]
          ATM_W(:,:), ATM_U(:,:), ATM_V(:,:), SFC_TEMP(:,:),     & ! [IN]
          Z1(:,:), SFC_DENS(:,:),                                & ! [IN]
          SFLX_MW(:,:), SFLX_MU(:,:), SFLX_MV(:,:),              & ! [OUT]
          SFLX_SH(:,:), SFLX_LH(:,:), SFLX_QV(:,:),              & ! [OUT]
          U10(:,:), V10(:,:)                                     ) ! [OUT]
      end select

    end if

    !- Add the tendency due to the surface flux

    call cal_del_flux( del_flux,            &
      SFLX_MU, SFLX_MV, SFLX_MW, SFLX_SH,   &
      lcmesh%normal_fn(:,:,3),              &
      lcmesh, elem, lcmesh2D, elem2D        )
    
    !$omp parallel do private(ke, LiftDelFlx)
    do ke2D=1, lcmesh2D%Ne
      ke = ke2D

      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,1), LiftDelFlx )
      MOMX_tp(:,ke) = MOMX_tp(:,ke) - LiftDelFlx(:) 

      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,2), LiftDelFlx )
      MOMY_tp(:,ke) = MOMY_tp(:,ke) - LiftDelFlx(:) 

      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,3), LiftDelFlx )
      MOMZ_tp(:,ke) = MOMZ_tp(:,ke) - LiftDelFlx(:) 
      
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke)*del_flux(:,ke,4), LiftDelFlx )
      RHOH_p (:,ke) = RHOH_p (:,ke) - LiftDelFlx(:) 
    end do

  end subroutine cal_tend_from_sfcflx

!OCL SERIAL  
  subroutine cal_del_flux( del_flux,      &
    sflx_mu, sflx_mv, sflx_mw, sflx_sh,   &
    nz,                                   &
    lmesh, elem, lmesh2D, elem2D          )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(elementbase2D), intent(in) :: elem2D   
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne,4)
    real(RP), intent(in) :: sflx_mu(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: sflx_mv(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: sflx_mw(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: sflx_sh(elem2D%Np,lmesh2D%NeA)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)

    integer :: ke2D, p
    integer :: i
    !------------------------------------------------------------------------

    !$omp parallel private(i)
    !$omp workshare
    del_flux(:,:) = 0.0_RP
    !$omp end workshare

    !$omp do collapse(2)
    do ke2D=1, lmesh2D%Ne
    do p=1, elem2D%Np
      i = elem%Nfaces_h*elem%Nfp_h + p + (ke2D-1)*elem%NfpTot
      del_flux(i,1) = sflx_mu(p,ke2D) * nz(i)
      del_flux(i,2) = sflx_mv(p,ke2D) * nz(i)
      del_flux(i,3) = sflx_mw(p,ke2D) * nz(i)
      del_flux(i,4) = sflx_sh(p,ke2D) * nz(i)
    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine cal_del_flux

end module mod_atmos_phy_sfc