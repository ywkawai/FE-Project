!-------------------------------------------------------------------------------
!> module FElib / Data / Statistics / numerical error
!!
!! @par Description
!!           This module provides classes to evaluate numerical error for 1D, 2D, and 3D fields
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_meshfield_analysis_numerror
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: &
    PRC_ismaster, PRC_abort
  use scale_time_manager, only: &
    TIME_NSTEP

  use scale_element_base, only: &
    ElementBase,                &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_element_line, only: LineElement
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_meshfield_base, only: &
    MeshFieldBase, MeshField1D, MeshField2D, MeshField3D
    
  use scale_meshfield_analysis_numerror_base, only: &
    MeshFieldAnalysisNumerrorBase,   &
    MeshFieldAnalysisNumerrorInfoBase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !- 1D

  !> Derived type for numerical error analysis of 1D field
  type, public, extends(MeshFieldAnalysisNumerrorBase) :: MeshFieldAnalysisNumerror1D
    class(MeshBase1D), pointer :: mesh1D
    procedure(set_data_lc_1d), pointer :: set_data_lc
  contains
    procedure :: Init_1D => meshfield_analysis_numerror_1D_Init
    generic :: Init => Init_1D
    procedure :: Evaluate => meshfield_analysis_numerror_1D_evaluate
  end type MeshFieldAnalysisNumerror1D
  interface 
    subroutine set_data_lc_1d( this_, q, qexact, qexact_intrp, lcmesh, elem1D, intrp_epos, tsec )
      import RP
      import ElementBase1D
      import LocalMesh1D
      import MeshFieldAnalysisNumerror1D
      class(MeshFieldAnalysisNumerror1D), intent(in) :: this_
      class(LocalMesh1D), intent(in) :: lcmesh
      class(ElementBase1D), intent(in) :: elem1D
      real(RP), intent(out) :: q(elem1D%Np,lcmesh%Ne,this_%var_num)
      real(RP), intent(out) :: qexact(elem1D%Np,lcmesh%Ne,this_%var_num)
      real(RP), intent(out) :: qexact_intrp(this_%intrp_np,lcmesh%Ne,this_%var_num)
      real(RP), intent(in) :: intrp_epos(this_%intrp_np,this_%ndim)
      real(RP), intent(in) :: tsec
    end subroutine set_data_lc_1d
  end interface

  !- 2D
  !> Derived type for numerical error analysis of 2D field
  type, public, extends(MeshFieldAnalysisNumerrorBase) :: MeshFieldAnalysisNumerror2D
    class(MeshBase2D), pointer :: mesh2D
    procedure(set_data_lc_2d), pointer :: set_data_lc
  contains
    procedure :: Init_2D => meshfield_analysis_numerror_2D_Init
    generic :: Init => Init_2D
    procedure :: Evaluate => meshfield_analysis_numerror_2D_evaluate
  end type MeshFieldAnalysisNumerror2D
  interface 
    subroutine set_data_lc_2d( this_, q, qexact, qexact_intrp, lcmesh, elem2D, intrp_epos, tsec )
      import RP
      import ElementBase2D
      import LocalMesh2D
      import MeshFieldAnalysisNumerror2D
      class(MeshFieldAnalysisNumerror2D), intent(in) :: this_
      class(LocalMesh2D), intent(in) :: lcmesh
      class(ElementBase2D), intent(in) :: elem2D
      real(RP), intent(out) :: q(elem2D%Np,lcmesh%Ne,this_%var_num)
      real(RP), intent(out) :: qexact(elem2D%Np,lcmesh%Ne,this_%var_num)
      real(RP), intent(out) :: qexact_intrp(this_%intrp_np,lcmesh%Ne,this_%var_num)
      real(RP), intent(in) :: intrp_epos(this_%intrp_np,this_%ndim)
      real(RP), intent(in) :: tsec
    end subroutine set_data_lc_2d
  end interface

  !- 3D
  !> Derived type for numerical error analysis of 3D field
  type, public, extends(MeshFieldAnalysisNumerrorBase) :: MeshFieldAnalysisNumerror3D
    class(MeshBase3D), pointer :: mesh3D
    procedure(set_data_lc_3d), pointer :: set_data_lc
  contains
    procedure :: Init_3D => meshfield_analysis_numerror_3D_Init
    generic :: Init => Init_3D
    procedure :: Evaluate => meshfield_analysis_numerror_3D_evaluate
  end type MeshFieldAnalysisNumerror3D
  interface 
    subroutine set_data_lc_3d( this_, q, qexact, qexact_intrp, lcmesh, elem3D, intrp_epos, tsec )
      import RP
      import ElementBase3D
      import LocalMesh3D
      import MeshFieldAnalysisNumerror3D
      class(MeshFieldAnalysisNumerror3D), intent(in) :: this_
      class(LocalMesh3D), intent(in) :: lcmesh
      class(ElementBase3D), intent(in) :: elem3D
      real(RP), intent(out) :: q(elem3D%Np,lcmesh%Ne,this_%var_num)
      real(RP), intent(out) :: qexact(elem3D%Np,lcmesh%Ne,this_%var_num)
      real(RP), intent(out) :: qexact_intrp(this_%intrp_np,lcmesh%Ne,this_%var_num)
      real(RP), intent(in) :: intrp_epos(this_%intrp_np,this_%ndim)
      real(RP), intent(in) :: tsec
    end subroutine set_data_lc_3d
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains

!-- 1D --

  !-----------------------------------------------------------------------------
  !> Initialize a object to evaluate numerical error for 1D field
!OCL SERIAL
  subroutine meshfield_analysis_numerror_1D_Init( this,    &
    porder_error_check, log_fname_base, log_step_interval, &
    mesh, refElem1D,                                       &
    set_data_lc, numerror_analysis_info )
    implicit none
    class(MeshFieldAnalysisNumerror1D), intent(inout) :: this
    integer, intent(in) :: porder_error_check
    character(len=*), intent(in) :: log_fname_base
    integer, intent(in) :: log_step_interval
    class(MeshBase1D), intent(in), target :: mesh
    type(LineElement), intent(in) :: refElem1D
    procedure(set_data_lc_1d) :: set_data_lc      
    class(MeshFieldAnalysisNumerrorInfoBase), intent(in), target :: numerror_analysis_info
    !---------------------------------------------------------------------------

    call this%MeshFieldAnalysisNumerrorBase%Init( &
      porder_error_check, 1, refElem1D%Np, porder_error_check**1,     &
      log_fname_base, log_step_interval, mesh, numerror_analysis_info )

    this%IntrpMat(:,:) = refElem1D%GenIntGaussLegendreIntrpMat( &
      this%PolyOrderErrorCheck,                                 &  ! (in)
      this%intw_intrp, this%epos_intrp(:,1)                     )  ! (out)

    this%mesh1D => mesh
    this%set_data_lc => set_data_lc
    return
  end subroutine meshfield_analysis_numerror_1D_Init
  
!OCL SERIAL
  subroutine meshfield_analysis_numerror_1D_evaluate( this, &
    tstep, tsec )
    implicit none
    class(MeshFieldAnalysisNumerror1D), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP), intent(in) :: tsec
    !---------------------------------------------------------------------------

    call this%Evaluate_base( tstep, tsec, this%mesh1D%dom_vol, evaluate_error_core, calc_covariance_core )
    return
  end subroutine meshfield_analysis_numerror_1D_evaluate

!-- 2D --

  !-----------------------------------------------------------------------------
  !> Initialize a object to evaluate numerical error for 2D field
!OCL SERIAL
  subroutine meshfield_analysis_numerror_2D_Init( this, &
    porder_error_check, log_fname_base, log_step_interval, &
    mesh, refElem2D,                                       &
    set_data_lc, numerror_analysis_info )
    implicit none
    class(MeshFieldAnalysisNumerror2D), intent(inout) :: this
    integer, intent(in) :: porder_error_check
    character(len=*), intent(in) :: log_fname_base
    integer, intent(in) :: log_step_interval
    class(MeshBase2D), intent(in), target :: mesh
    type(QuadrilateralElement), intent(in) :: refElem2D
    class(MeshFieldAnalysisNumerrorInfoBase), intent(in), target :: numerror_analysis_info
    procedure(set_data_lc_2d) :: set_data_lc
    !---------------------------------------------------------------------------

    call this%MeshFieldAnalysisNumerrorBase%Init( &
      porder_error_check, 2, refElem2D%Np, porder_error_check**2,     &
      log_fname_base, log_step_interval, mesh, numerror_analysis_info )

    this%IntrpMat(:,:) = refElem2D%GenIntGaussLegendreIntrpMat( &
      this%PolyOrderErrorCheck,                                   &  ! (in)
      this%intw_intrp, this%epos_intrp(:,1), this%epos_intrp(:,2) )  ! (out)

    this%mesh2D => mesh
    this%set_data_lc => set_data_lc
    return
  end subroutine meshfield_analysis_numerror_2D_Init

!OCL SERIAL
  subroutine meshfield_analysis_numerror_2D_evaluate( this, &
    tstep, tsec )
    implicit none
    class(MeshFieldAnalysisNumerror2D), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP), intent(in) :: tsec
    !---------------------------------------------------------------------------

    call this%Evaluate_base( tstep, tsec, this%mesh2D%dom_vol, evaluate_error_core, calc_covariance_core )
    return
  end subroutine meshfield_analysis_numerror_2D_evaluate

!-- 3D --

  !-----------------------------------------------------------------------------
  !> Initialize a object to evaluate numerical error for 3D field
!OCL SERIAL
  subroutine meshfield_analysis_numerror_3D_Init( this,    &
    porder_error_check, log_fname_base, log_step_interval, &
    mesh, refElem3D,                                       &
    set_data_lc, numerror_analysis_info )
    implicit none
    class(MeshFieldAnalysisNumerror3D), intent(inout) :: this
    integer, intent(in) :: porder_error_check
    character(len=*), intent(in) :: log_fname_base
    integer, intent(in) :: log_step_interval
    class(MeshBase3D), intent(in), target :: mesh
    type(HexahedralElement), intent(in) :: refElem3D
    class(MeshFieldAnalysisNumerrorInfoBase), intent(in), target :: numerror_analysis_info
    procedure(set_data_lc_3d), pointer :: set_data_lc
    !---------------------------------------------------------------------------

    call this%MeshFieldAnalysisNumerrorBase%Init( &
      porder_error_check, 3, refElem3D%Np, porder_error_check**3,     &
      log_fname_base, log_step_interval, mesh, numerror_analysis_info )


    this%IntrpMat(:,:) = refElem3D%GenIntGaussLegendreIntrpMat( &
      this%PolyOrderErrorCheck,                                                         &  ! (in)
      this%intw_intrp, this%epos_intrp(:,1), this%epos_intrp(:,2), this%epos_intrp(:,3) )  ! (out)

    this%mesh3D => mesh
    this%set_data_lc => set_data_lc
    return
  end subroutine meshfield_analysis_numerror_3D_Init

!OCL SERIAL
  subroutine meshfield_analysis_numerror_3D_evaluate( &
    this, tstep, tsec, mesh  )
    implicit none
    class(MeshFieldAnalysisNumerror3D), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP), intent(in) :: tsec
    class(MeshBase3D), target, intent(in) :: mesh
    procedure(set_data_lc_3d) :: set_data_lc
    !---------------------------------------------------------------------------

    call this%Evaluate_base( tstep, tsec, this%mesh3D%dom_vol, evaluate_error_core, calc_covariance_core )
    return
  end subroutine meshfield_analysis_numerror_3D_evaluate

!--- private ---------------------

!OCL SERIAL
  subroutine evaluate_error_core( base, tsec, &
    num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
    numsol_mean_lc, exactsol_mean_lc                     )
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
    real(RP), intent(in) :: tsec
    real(RP), intent(inout) :: num_error_l1_lc(base%var_num)
    real(RP), intent(inout) :: num_error_l2_lc(base%var_num)
    real(RP), intent(inout) :: num_error_linf_lc(base%var_num)
    real(RP), intent(inout) :: numsol_mean_lc(base%var_num) 
    real(RP), intent(inout) :: exactsol_mean_lc(base%var_num) 

    real(RP), allocatable :: q(:,:,:)
    real(RP), allocatable :: qexact(:,:,:)
    real(RP), allocatable :: qexact_intrp(:,:,:)

    integer :: lcdomid      
    integer :: Np, Ne
    class(ElementBase), pointer :: elem
    class(LocalMeshBase), pointer :: lcmesh
    !---------------------------------------------------------------------------

    do lcdomid=1, base%mesh%LOCAL_MESH_NUM
      call prepare_lc_data( base, lcdomid, tsec,      & ! (in)
        q, qexact, qexact_intrp, lcmesh, elem, Np, Ne ) ! (out)

      call base%Evaluate_error_lc( &
        num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
        numsol_mean_lc, exactsol_mean_lc,                    &
        q, qexact, qexact_intrp, lcmesh, elem                )

      deallocate( q, qexact, qexact_intrp )
    end do
    return
  end subroutine evaluate_error_core

!OCL SERIAL
  subroutine calc_covariance_core( base, tsec,&
    cov_numsol_numsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc, &
    numsol_mean, exactsol_mean                      )
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
    real(RP), intent(in) :: tsec
    real(RP), intent(inout) :: cov_numsol_numsol_lc(base%var_num)
    real(RP), intent(inout) :: cov_numsol_exactsol_lc(base%var_num)
    real(RP), intent(inout) :: cov_exactsol_exactsol_lc(base%var_num)
    real(RP), intent(in) :: numsol_mean(base%var_num) 
    real(RP), intent(in) :: exactsol_mean(base%var_num) 

    real(RP), allocatable :: q(:,:,:)
    real(RP), allocatable :: qexact(:,:,:)
    real(RP), allocatable :: qexact_intrp(:,:,:)

    integer :: lcdomid      
    integer :: Np, Ne
    class(LocalMeshBase), pointer :: lcmesh
    class(ElementBase), pointer :: elem
    !---------------------------------------------------------------------------

    do lcdomid=1, base%mesh%LOCAL_MESH_NUM
      call prepare_lc_data( base, lcdomid, tsec,      & ! (in)
        q, qexact, qexact_intrp, lcmesh, elem, Np, Ne ) ! (out)

      call base%Evaluate_covariance_lc( &
        cov_numsol_numsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc, &
        q, qexact, qexact_intrp, numsol_mean, exactsol_mean, lcmesh, elem       )
      
      deallocate( q, qexact, qexact_intrp )
    end do
    return
  end subroutine calc_covariance_core    

!OCL SERIAL
  subroutine prepare_lc_data( base, lcdomid, tsec, &
    q, qexact, qexact_intrp, lcmesh, elem, Np, Ne  )
    implicit none
    class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
    integer, intent(in) :: lcdomid
    real(RP), intent(in) :: tsec
    real(RP), allocatable, intent(out) :: q(:,:,:)
    real(RP), allocatable, intent(out) :: qexact(:,:,:)
    real(RP), allocatable, intent(out) :: qexact_intrp(:,:,:)
    class(LocalMeshBase), pointer, intent(out) :: lcmesh
    class(ElementBase),   pointer, intent(out) :: elem
    integer, intent(out) :: Np, Ne

    class(LocalMesh1D), pointer :: lcmesh1D
    class(LocalMesh2D), pointer :: lcmesh2D
    class(LocalMesh3D), pointer :: lcmesh3D
    !-------------------------------------------------------------------

    !- Decide Np/Ne and pointers
    select type(base)
    class is (MeshFieldAnalysisNumerror1D)
      lcmesh1D => base%mesh1D%lcmesh_list(lcdomid)
      Np = lcmesh1D%refElem1D%Np;  Ne = lcmesh1D%Ne
      lcmesh => lcmesh1D;         elem => lcmesh1D%refElem1D
    class is (MeshFieldAnalysisNumerror2D)
      lcmesh2D => base%mesh2D%lcmesh_list(lcdomid)
      Np = lcmesh2D%refElem2D%Np;  Ne = lcmesh2D%Ne
      lcmesh => lcmesh2D;         elem => lcmesh2D%refElem2D
    class is (MeshFieldAnalysisNumerror3D)
      lcmesh3D => base%mesh3D%lcmesh_list(lcdomid)
      Np = lcmesh3D%refElem3D%Np;  Ne = lcmesh3D%Ne
      lcmesh => lcmesh3D;         elem => lcmesh3D%refElem3D
    end select

    allocate( q(Np,Ne,base%var_num) )
    allocate( qexact(Np,Ne,base%var_num) )
    allocate( qexact_intrp(base%intrp_np,Ne,base%var_num) )

    !- Fill arrays via set_data_lc
    select type(base)
    class is (MeshFieldAnalysisNumerror1D)
      call base%set_data_lc( q, qexact, qexact_intrp,       &
        lcmesh1D, lcmesh1D%refElem1D, base%epos_intrp, tsec )
    class is (MeshFieldAnalysisNumerror2D)
      call base%set_data_lc( q, qexact, qexact_intrp,       &
        lcmesh2D, lcmesh2D%refElem2D, base%epos_intrp, tsec )
    class is (MeshFieldAnalysisNumerror3D)
      call base%set_data_lc( q, qexact, qexact_intrp,       &
        lcmesh3D, lcmesh3D%refElem3D, base%epos_intrp, tsec )
    end select
    return
  end subroutine prepare_lc_data

end module scale_meshfield_analysis_numerror