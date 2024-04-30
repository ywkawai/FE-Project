!-------------------------------------------------------------------------------
!> module meshfield_analysis_numerror_base
!!
!! @par Description
!!           This module provides classes to evaluate numerical error for 1D, 2D, and 3D fields
!!
!! @author Team SCALE
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
    MeshFieldAnalysisNumerrorBase

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  type, public, extends(MeshFieldAnalysisNumerrorBase) :: MeshFieldAnalysisNumerror1D
  contains
    procedure :: Init_1D => meshfield_analysis_numerror_1D_Init
    generic :: Init => Init_1D
    procedure :: Evaluate => meshfield_analysis_numerror_1D_evaluate
  end type MeshFieldAnalysisNumerror1D

  type, public, extends(MeshFieldAnalysisNumerrorBase) :: MeshFieldAnalysisNumerror2D
  contains
    procedure :: Init_2D => meshfield_analysis_numerror_2D_Init
    generic :: Init => Init_2D
    procedure :: Evaluate => meshfield_analysis_numerror_2D_evaluate
  end type MeshFieldAnalysisNumerror2D

  type, public, extends(MeshFieldAnalysisNumerrorBase) :: MeshFieldAnalysisNumerror3D
  contains
    procedure :: Init_3D => meshfield_analysis_numerror_3D_Init
    generic :: Init => Init_3D
    procedure :: Evaluate => meshfield_analysis_numerror_3D_evaluate
  end type MeshFieldAnalysisNumerror3D

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
    refElem1D )
    implicit none
    class(MeshFieldAnalysisNumerror1D), intent(inout) :: this
    integer, intent(in) :: porder_error_check
    character(len=*), intent(in) :: log_fname_base
    integer, intent(in) :: log_step_interval
    type(LineElement), intent(in) :: refElem1D
    !---------------------------------------------------------------------------

    call this%MeshFieldAnalysisNumerrorBase%Init( &
      porder_error_check, 1, refElem1D%Np, porder_error_check**1, &
      log_fname_base, log_step_interval                           )

    this%IntrpMat(:,:) = refElem1D%GenIntGaussLegendreIntrpMat( &
      this%PolyOrderErrorCheck,                                 &  ! (in)
      this%intw_intrp, this%epos_intrp(:,1)                     )  ! (out)

    return
  end subroutine meshfield_analysis_numerror_1D_Init

!OCL SERIAL
  subroutine meshfield_analysis_numerror_1D_evaluate( &
    this, tstep, tsec, mesh, set_data_lc              )
    implicit none
    class(MeshFieldAnalysisNumerror1D), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP) :: tsec
    class(MeshBase1D), target, intent(in) :: mesh
    interface 
      subroutine set_data_lc( this_, q, qexact, qexact_intrp, lcmesh, elem1D, intrp_epos )
        import RP
        import ElementBase1D
        import LocalMesh1D
        import MeshFieldAnalysisNumerror1D
        class(MeshFieldAnalysisNumerror1D), intent(in) :: this_
        class(LocalMesh1D), intent(in) :: lcmesh
        class(ElementBase1D) :: elem1D
        real(RP), intent(out) :: q(elem1D%Np,lcmesh%Ne,this_%var_num)
        real(RP), intent(out) :: qexact(elem1D%Np,lcmesh%Ne,this_%var_num)
        real(RP), intent(out) :: qexact_intrp(this_%intrp_np,lcmesh%Ne,this_%var_num)
        real(RP), intent(in) :: intrp_epos(this_%intrp_np,this_%ndim)
      end subroutine set_data_lc
    end interface
    !---------------------------------------------------------------------------

    call this%Evaluate_base( tstep, tsec, mesh%dom_vol, evaluate_error_core, calc_covariance_core )

    return
  contains 
!OCL SERIAL
    subroutine evaluate_error_core( base, &
      num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
      numsol_mean_lc, exactsol_mean_lc                     )
      implicit none
      class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
      real(RP), intent(inout) :: num_error_l1_lc(base%var_num)
      real(RP), intent(inout) :: num_error_l2_lc(base%var_num)
      real(RP), intent(inout) :: num_error_linf_lc(base%var_num)
      real(RP), intent(inout) :: numsol_mean_lc(base%var_num) 
      real(RP), intent(inout) :: exactsol_mean_lc(base%var_num) 

      integer :: lcdomid      
      class(LocalMesh1D), pointer :: lcmesh
      real(RP), allocatable :: q(:,:,:)
      real(RP), allocatable :: qexact(:,:,:)
      real(RP), allocatable :: qexact_intrp(:,:,:)
     !---------------------------------------------------------------------------

      do lcdomid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(lcdomid)
        allocate( q(lcmesh%refElem1D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact(lcmesh%refElem1D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num) )
  
        call set_data_lc( this, q, qexact, qexact_intrp,       & ! (in)
          lcmesh, lcmesh%refElem1D, base%epos_intrp            ) ! (out)
  
        call this%Evaluate_error_lc( &
          num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
          numsol_mean_lc, exactsol_mean_lc,                    &
          q, qexact, qexact_intrp, lcmesh, lcmesh%refElem1D    )
  
        deallocate( q, qexact, qexact_intrp )
      end do

      return
    end subroutine evaluate_error_core

!OCL SERIAL
    subroutine calc_covariance_core( base, &
      cov_numsol_numsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc, &
      numsol_mean, exactsol_mean                      )
      implicit none
      class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
      real(RP), intent(inout) :: cov_numsol_numsol_lc(base%var_num)
      real(RP), intent(inout) :: cov_numsol_exactsol_lc(base%var_num)
      real(RP), intent(inout) :: cov_exactsol_exactsol_lc(base%var_num)
      real(RP), intent(in) :: numsol_mean(base%var_num) 
      real(RP), intent(in) :: exactsol_mean(base%var_num) 

      integer :: lcdomid      
      class(LocalMesh1D), pointer :: lcmesh
      real(RP), allocatable :: q(:,:,:)
      real(RP), allocatable :: qexact(:,:,:)
      real(RP), allocatable :: qexact_intrp(:,:,:)
     !---------------------------------------------------------------------------

      do lcdomid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(lcdomid)
        allocate( q(lcmesh%refElem1D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact(lcmesh%refElem1D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num) )
  
        call set_data_lc( this, q, qexact, qexact_intrp, & ! (in)
          lcmesh, lcmesh%refElem1D, base%epos_intrp      ) ! (out)
  
        call this%Evaluate_covariance_lc( &
          cov_numsol_exactsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc,     &
          q, qexact, qexact_intrp, numsol_mean, exactsol_mean, lcmesh, lcmesh%refElem1D )
        
        deallocate( q, qexact, qexact_intrp )
      end do

      return
    end subroutine calc_covariance_core    
  end subroutine meshfield_analysis_numerror_1D_evaluate

!-- 2D --

  !-----------------------------------------------------------------------------
  !> Initialize a object to evaluate numerical error for 2D field
!OCL SERIAL
  subroutine meshfield_analysis_numerror_2D_Init( this, &
    porder_error_check, log_fname_base, log_step_interval, &
    refElem2D )
    implicit none
    class(MeshFieldAnalysisNumerror2D), intent(inout) :: this
    integer, intent(in) :: porder_error_check
    character(len=*), intent(in) :: log_fname_base
    integer, intent(in) :: log_step_interval
    type(QuadrilateralElement), intent(in) :: refElem2D
    !---------------------------------------------------------------------------

    call this%MeshFieldAnalysisNumerrorBase%Init( &
      porder_error_check, 2, refElem2D%Np, porder_error_check**2, &
      log_fname_base, log_step_interval                           )

    this%IntrpMat(:,:) = refElem2D%GenIntGaussLegendreIntrpMat( &
      this%PolyOrderErrorCheck,                                   &  ! (in)
      this%intw_intrp, this%epos_intrp(:,1), this%epos_intrp(:,2) )  ! (out)

    return
  end subroutine meshfield_analysis_numerror_2D_Init

!OCL SERIAL
  subroutine meshfield_analysis_numerror_2D_evaluate( &
    this, tstep, tsec, mesh, set_data_lc              )
    implicit none
    class(MeshFieldAnalysisNumerror2D), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP) :: tsec
    class(MeshBase2D), target, intent(in) :: mesh
    interface 
      subroutine set_data_lc( this_, q, qexact, qexact_intrp, lcmesh, elem2D, intrp_epos )
        import RP
        import ElementBase2D
        import LocalMesh2D
        import MeshFieldAnalysisNumerror2D
        class(MeshFieldAnalysisNumerror2D), intent(in) :: this_
        class(LocalMesh2D), intent(in) :: lcmesh
        class(ElementBase2D) :: elem2D
        real(RP), intent(out) :: q(elem2D%Np,lcmesh%Ne,this_%var_num)
        real(RP), intent(out) :: qexact(elem2D%Np,lcmesh%Ne,this_%var_num)
        real(RP), intent(out) :: qexact_intrp(this_%intrp_np,lcmesh%Ne,this_%var_num)
        real(RP), intent(in) :: intrp_epos(this_%intrp_np,this_%ndim)
      end subroutine set_data_lc
    end interface
    !---------------------------------------------------------------------------

    call this%Evaluate_base( tstep, tsec, mesh%dom_vol, evaluate_error_core, calc_covariance_core )

    return
  contains 
!OCL SERIAL
    subroutine evaluate_error_core( base, &
      num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
      numsol_mean_lc, exactsol_mean_lc                     )
      implicit none
      class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
      real(RP), intent(inout) :: num_error_l1_lc(base%var_num)
      real(RP), intent(inout) :: num_error_l2_lc(base%var_num)
      real(RP), intent(inout) :: num_error_linf_lc(base%var_num)
      real(RP), intent(inout) :: numsol_mean_lc(base%var_num) 
      real(RP), intent(inout) :: exactsol_mean_lc(base%var_num) 

      integer :: lcdomid      
      class(LocalMesh2D), pointer :: lcmesh
      real(RP), allocatable :: q(:,:,:)
      real(RP), allocatable :: qexact(:,:,:)
      real(RP), allocatable :: qexact_intrp(:,:,:)
     !---------------------------------------------------------------------------

      do lcdomid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(lcdomid)
        allocate( q(lcmesh%refElem2D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact(lcmesh%refElem2D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num) )
  
        call set_data_lc( this, q, qexact, qexact_intrp,       & ! (in)
          lcmesh, lcmesh%refElem2D, base%epos_intrp            ) ! (out)
  
        call this%Evaluate_error_lc( &
          num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
          numsol_mean_lc, exactsol_mean_lc,                    &
          q, qexact, qexact_intrp, lcmesh, lcmesh%refElem2D    )
  
        deallocate( q, qexact, qexact_intrp )
      end do

      return
    end subroutine evaluate_error_core

!OCL SERIAL
    subroutine calc_covariance_core( base, &
      cov_numsol_numsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc, &
      numsol_mean, exactsol_mean                      )
      implicit none
      class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
      real(RP), intent(inout) :: cov_numsol_numsol_lc(base%var_num)
      real(RP), intent(inout) :: cov_numsol_exactsol_lc(base%var_num)
      real(RP), intent(inout) :: cov_exactsol_exactsol_lc(base%var_num)
      real(RP), intent(in) :: numsol_mean(base%var_num) 
      real(RP), intent(in) :: exactsol_mean(base%var_num) 

      integer :: lcdomid      
      class(LocalMesh2D), pointer :: lcmesh
      real(RP), allocatable :: q(:,:,:)
      real(RP), allocatable :: qexact(:,:,:)
      real(RP), allocatable :: qexact_intrp(:,:,:)
     !---------------------------------------------------------------------------

      do lcdomid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(lcdomid)
        allocate( q(lcmesh%refElem2D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact(lcmesh%refElem2D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num) )
  
        call set_data_lc( this, q, qexact, qexact_intrp, & ! (in)
          lcmesh, lcmesh%refElem2D, base%epos_intrp      ) ! (out)
  
        call this%Evaluate_covariance_lc( &
          cov_numsol_exactsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc,     &
          q, qexact, qexact_intrp, numsol_mean, exactsol_mean, lcmesh, lcmesh%refElem2D )
        
        deallocate( q, qexact, qexact_intrp )
      end do

      return
    end subroutine calc_covariance_core    
  end subroutine meshfield_analysis_numerror_2D_evaluate

!-- 3D --

  !-----------------------------------------------------------------------------
  !> Initialize a object to evaluate numerical error for 3D field
!OCL SERIAL
  subroutine meshfield_analysis_numerror_3D_Init( this, &
    porder_error_check, log_fname_base, log_step_interval, &
    refElem3D )
    implicit none
    class(MeshFieldAnalysisNumerror3D), intent(inout) :: this
    integer, intent(in) :: porder_error_check
    character(len=*), intent(in) :: log_fname_base
    integer, intent(in) :: log_step_interval
    type(HexahedralElement), intent(in) :: refElem3D
    !---------------------------------------------------------------------------

    call this%MeshFieldAnalysisNumerrorBase%Init( &
      porder_error_check, 3, refElem3D%Np, porder_error_check**3, &
      log_fname_base, log_step_interval                           )

    this%IntrpMat(:,:) = refElem3D%GenIntGaussLegendreIntrpMat( &
      this%PolyOrderErrorCheck,                                                         &  ! (in)
      this%intw_intrp, this%epos_intrp(:,1), this%epos_intrp(:,2), this%epos_intrp(:,3) )  ! (out)

    return
  end subroutine meshfield_analysis_numerror_3D_Init

!OCL SERIAL
  subroutine meshfield_analysis_numerror_3D_evaluate( &
    this, tstep, tsec, mesh, set_data_lc              )
    implicit none
    class(MeshFieldAnalysisNumerror3D), intent(inout) :: this
    integer, intent(in) :: tstep
    real(RP) :: tsec
    class(MeshBase3D), target, intent(in) :: mesh
    interface 
      subroutine set_data_lc( this_, q, qexact, qexact_intrp, lcmesh, elem3D, intrp_epos )
        import RP
        import ElementBase3D
        import LocalMesh3D
        import MeshFieldAnalysisNumerror3D
        class(MeshFieldAnalysisNumerror3D), intent(in) :: this_
        class(LocalMesh3D), intent(in) :: lcmesh
        class(ElementBase3D) :: elem3D
        real(RP), intent(out) :: q(elem3D%Np,lcmesh%Ne,this_%var_num)
        real(RP), intent(out) :: qexact(elem3D%Np,lcmesh%Ne,this_%var_num)
        real(RP), intent(out) :: qexact_intrp(this_%intrp_np,lcmesh%Ne,this_%var_num)
        real(RP), intent(in) :: intrp_epos(this_%intrp_np,this_%ndim)
      end subroutine set_data_lc
    end interface
    !---------------------------------------------------------------------------

    call this%Evaluate_base( tstep, tsec, mesh%dom_vol, evaluate_error_core, calc_covariance_core )

    return
  contains 
!OCL SERIAL
    subroutine evaluate_error_core( base, &
      num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
      numsol_mean_lc, exactsol_mean_lc                     )
      implicit none
      class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
      real(RP), intent(inout) :: num_error_l1_lc(base%var_num)
      real(RP), intent(inout) :: num_error_l2_lc(base%var_num)
      real(RP), intent(inout) :: num_error_linf_lc(base%var_num)
      real(RP), intent(inout) :: numsol_mean_lc(base%var_num) 
      real(RP), intent(inout) :: exactsol_mean_lc(base%var_num) 

      integer :: lcdomid      
      class(LocalMesh3D), pointer :: lcmesh
      real(RP), allocatable :: q(:,:,:)
      real(RP), allocatable :: qexact(:,:,:)
      real(RP), allocatable :: qexact_intrp(:,:,:)
     !---------------------------------------------------------------------------

      do lcdomid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(lcdomid)
        allocate( q(lcmesh%refElem3D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact(lcmesh%refElem3D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num) )
  
        call set_data_lc( this, q, qexact, qexact_intrp,       & ! (in)
          lcmesh, lcmesh%refElem3D, base%epos_intrp            ) ! (out)
  
        call this%Evaluate_error_lc( &
          num_error_l1_lc, num_error_l2_lc, num_error_linf_lc, &
          numsol_mean_lc, exactsol_mean_lc,                    &
          q, qexact, qexact_intrp, lcmesh, lcmesh%refElem3D    )
  
        deallocate( q, qexact, qexact_intrp )
      end do

      return
    end subroutine evaluate_error_core

!OCL SERIAL
    subroutine calc_covariance_core( base, &
      cov_numsol_numsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc, &
      numsol_mean, exactsol_mean                      )
      implicit none
      class(MeshFieldAnalysisNumerrorBase), intent(in) :: base
      real(RP), intent(inout) :: cov_numsol_numsol_lc(base%var_num)
      real(RP), intent(inout) :: cov_numsol_exactsol_lc(base%var_num)
      real(RP), intent(inout) :: cov_exactsol_exactsol_lc(base%var_num)
      real(RP), intent(in) :: numsol_mean(base%var_num) 
      real(RP), intent(in) :: exactsol_mean(base%var_num) 

      integer :: lcdomid
      class(LocalMesh3D), pointer :: lcmesh
      real(RP), allocatable :: q(:,:,:)
      real(RP), allocatable :: qexact(:,:,:)
      real(RP), allocatable :: qexact_intrp(:,:,:)
     !---------------------------------------------------------------------------

      do lcdomid=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(lcdomid)
        allocate( q(lcmesh%refElem3D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact(lcmesh%refElem3D%Np,lcmesh%Ne,base%var_num) )
        allocate( qexact_intrp(base%intrp_np,lcmesh%Ne,base%var_num) )
  
        call set_data_lc( this, q, qexact, qexact_intrp, & ! (in)
          lcmesh, lcmesh%refElem3D, base%epos_intrp      ) ! (out)
  
        call this%Evaluate_covariance_lc( &
          cov_numsol_exactsol_lc, cov_numsol_exactsol_lc, cov_exactsol_exactsol_lc,     &
          q, qexact, qexact_intrp, numsol_mean, exactsol_mean, lcmesh, lcmesh%refElem3D )
        
        deallocate( q, qexact, qexact_intrp )
      end do

      return
    end subroutine calc_covariance_core    
  end subroutine meshfield_analysis_numerror_3D_evaluate

!--- private

end module scale_meshfield_analysis_numerror