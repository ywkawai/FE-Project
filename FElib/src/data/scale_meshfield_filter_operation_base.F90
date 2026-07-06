!-------------------------------------------------------------------------------
!> module FElib / Data / Filter operation base
!!
!! @par Description
!!           This module provides a base class to apply filter operations to MeshField data
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_meshfield_filter_operation_base
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_element_modalfilter, only: ModalFilter
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfieldcomm_base, only: &
    MeshFieldContainer

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Base type to represent filter operation
  type, public :: MeshFieldFilterOperationBase
    integer :: operator_type
    integer :: hHaloSize

    integer :: Nnode_h1D_reconst

    !- Reconstruction matrices    
    real(RP), allocatable :: FilterMat_h1D(:,:)

    !- Reconstruction matrices
    real(RP), allocatable :: Minv_Ml_tr(:,:)
    real(RP), allocatable :: Minv_Mc_tr(:,:)
    real(RP), allocatable :: Minv_Mr_tr(:,:)
    real(RP), allocatable :: IntrpMat(:,:)

    !- Reconstruction matrices (2)
    real(RP), allocatable :: Ml_tr(:,:)
    real(RP), allocatable :: Mc_tr(:,:)
    real(RP), allocatable :: Mr_tr(:,:)

    !- Modal Filter
    type(ModalFilter) :: m_filter

  contains
    procedure, private :: Prepair_FilterMat => MeshFieldFilterOperationBase_prepair_filter_matrix
    procedure, private :: Prepair_ReconstMat => MeshFieldFilterOperationBase_prepair_reconstruct_matrix
    procedure, private :: Prepair_ReconstMat2 => MeshFieldFilterOperationBase_prepair_reconstruct2_matrix
  end type MeshFieldFilterOperationBase

  public :: MeshFieldFilterOperationBase_Init
  public :: MeshFieldFilterOperationBase_Final

  public :: MeshFieldFilterOperationBase_apply_filter1d_x
  public :: MeshFieldFilterOperationBase_apply_reconst1d_x
  public :: MeshFieldFilterOperationBase_apply_reconst1d_x_2
  public :: MeshFieldFilterOperationBase_apply_filter1d_y
  public :: MeshFieldFilterOperationBase_apply_reconst1d_y
  public :: MeshFieldFilterOperationBase_apply_reconst1d_y_2

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  private :: calc_filter_kenrnel
 

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  integer, parameter, public :: FILTER_OPTRTYPE_CONVFILTER   = 1
  integer, parameter, public :: FILTER_OPTRTYPE_RECONSTRUCT  = 2
  integer, parameter, public :: FILTER_OPTRTYPE_RECONSTRUCT2 = 3
  integer, parameter, public :: FILTER_OPTRTYPE_MODALFILTER  = 4

contains
!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_Init( this, &
    Nnode_h1D, &
    FilterOptrType, FilterShape, FilterWidthFac, Nnode_h1D_reconst )
    implicit none
    class(MeshFieldFilterOperationBase), intent(inout) :: this
    integer, intent(in) :: Nnode_h1D
    character(*), intent(in) :: FilterOptrType
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: FilterWidthFac
    integer, intent(in) :: Nnode_h1D_reconst
    !------------------------------

    select case(FilterOptrType)
    case ('ConvolFilter')
      call this%Prepair_FilterMat( Nnode_h1D, FilterShape, FilterWidthFac )
      this%operator_type = FILTER_OPTRTYPE_CONVFILTER
    case ('Reconstruction')
      call this%Prepair_ReconstMat( Nnode_h1D, Nnode_h1D_reconst )
      this%operator_type = FILTER_OPTRTYPE_RECONSTRUCT
      this%Nnode_h1D_reconst = Nnode_h1D_reconst
    case ('Reconstruction2')
      call this%Prepair_ReconstMat2( Nnode_h1D, Nnode_h1D_reconst )
      this%operator_type = FILTER_OPTRTYPE_RECONSTRUCT2
      this%Nnode_h1D_reconst = Nnode_h1D_reconst
    case default
      LOG_INFO('MeshFieldFilterOperationBase_Init',*) "Unsupported filter operation is specified. Check!", trim(FilterOptrType)
      call PRC_abort
    end select
    return
  end subroutine MeshFieldFilterOperationBase_Init

  subroutine MeshFieldFilterOperationBase_Final( this )
    implicit none
    class(MeshFieldFilterOperationBase), intent(inout) :: this
    !------------------------------

    if ( allocated(this%FilterMat_h1D) ) then
      deallocate( this%FilterMat_h1D )
    end if
    if ( allocated(this%Minv_Ml_tr) ) then
      deallocate( this%Minv_Ml_tr, this%Minv_Mc_tr, this%Minv_Mr_tr, this%IntrpMat )
    end if
    if ( allocated(this%Ml_tr) ) then
      deallocate( this%Ml_tr, this%Mc_tr, this%Mr_tr )
    end if

    return
  end subroutine MeshFieldFilterOperationBase_Final

!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_apply_filter1d_x( q, q0, Filter1D, Npx, Npy, Npz, Ne, NeA, Nnode_h1D )
    implicit none
    integer, intent(in) :: Npx, Npy, Npz, Ne
    integer, intent(in) :: NeA
    integer, intent(in) :: Nnode_h1D
    real(RP), intent(out) :: q(Npx,Npy,Npz,NeA)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,-Npy+1:2*Npy,Npz,Ne)
    real(RP), intent(in) :: Filter1D(Nnode_h1D,-Npx+1:2*Npx)

    integer :: ke, pz, py, px, p
    real(RP) :: tmp
    !-------------------------------------

    !$omp parallel do private(ke,px,py,pz,p,tmp) collapse(2)
    do ke=1, Ne
    do pz=1, Npz
      do py=1, Npy
      do px=1, Npx
        tmp = 0.0_RP
        do p=-Npx+1, Npx+Npx
          tmp = tmp + Filter1D(px,p) * q0(p,py,pz,ke)
        end do
        q(px,py,pz,ke) = tmp
      end do
      end do
    end do
    end do
    return
  end subroutine MeshFieldFilterOperationBase_apply_filter1d_x

!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_apply_reconst1d_x( q, q0, Minv_Ml_tr, Minv_Mc_tr, Minv_Mr_tr, IntrpMat, &
    Npx, Npy, Npz, Ne, Npx_reconst, lmesh )
    implicit none
    class(LocalMeshBase), intent(in) :: lmesh
    integer, intent(in) :: Npx, Npy, Npz, Ne
    integer, intent(in) :: Npx_reconst
    real(RP), intent(out) :: q(Npx,Npy,Npz,lmesh%NeA)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,-Npy+1:2*Npy,Npz,Ne)
    real(RP), intent(in) :: Minv_Ml_tr(Npx,Npx_reconst)
    real(RP), intent(in) :: Minv_Mc_tr(Npx,Npx_reconst)
    real(RP), intent(in) :: Minv_Mr_tr(Npx,Npx_reconst)
    real(RP), intent(in) :: IntrpMat(Npx,Npx_reconst)

    integer :: ke, pz, py, px

    integer :: i,j

    real(RP) :: tmp(Npx_reconst)
    real(RP) :: s
    !-------------------------------------

    !$omp parallel do private(ke,px,py,pz,i,j, tmp,s) collapse(2)
    do ke=1, Ne
    do pz=1, Npz
    do py=1, Npy
      do i=1, Npx_reconst
        s = 0.0_RP
        do j=1, Npx
          s = s &
            + Minv_Ml_tr(j,i) * q0(j-Npx,py,pz,ke) &
            + Minv_Mc_tr(j,i) * q0(    j,py,pz,ke) &
            + Minv_Mr_tr(j,i) * q0(j+Npx,py,pz,ke)
        end do
        tmp(i) = s
      end do
      q(:,py,pz,ke) = matmul(IntrpMat, tmp)
    end do
    end do
    end do
    return
  end subroutine MeshFieldFilterOperationBase_apply_reconst1d_x

!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_apply_reconst1d_x_2( q, q0, Minv_Ml_tr, Minv_Mc_tr, Minv_Mr_tr, &
    Npx, Npy, Npz, Ne, Npx_reconst, lmesh )
    implicit none
    class(LocalMeshBase), intent(in) :: lmesh
    integer, intent(in) :: Npx, Npy, Npz, Ne
    integer, intent(in) :: Npx_reconst
    real(RP), intent(out) :: q(Npx,Npy,Npz,lmesh%NeA)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,-Npy+1:2*Npy,Npz,Ne)
    real(RP), intent(in) :: Minv_Ml_tr(Npx,Npx)
    real(RP), intent(in) :: Minv_Mc_tr(Npx,Npx)
    real(RP), intent(in) :: Minv_Mr_tr(Npx,Npx)

    integer :: ke, pz, py, px

    integer :: i,j

    real(RP) :: tmp(Npx)
    real(RP) :: s
    !-------------------------------------

    !$omp parallel do private(ke,px,py,pz,i,j, tmp,s) collapse(2)
    do ke=1, Ne
    do pz=1, Npz
    do py=1, Npy
      do i=1, Npx
        s = 0.0_RP
        do j=1, Npx
          s = s &
            + Minv_Ml_tr(j,i) * q0(j-Npx,py,pz,ke) &
            + Minv_Mc_tr(j,i) * q0(    j,py,pz,ke) &
            + Minv_Mr_tr(j,i) * q0(j+Npx,py,pz,ke)
        end do
        tmp(i) = s
      end do
      q(:,py,pz,ke) = tmp(:)
    end do
    end do
    end do
    return
  end subroutine MeshFieldFilterOperationBase_apply_reconst1d_x_2

!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_apply_filter1d_y( q, q0, Filter1D, Npx, Npy, Npz, Ne, NeA, Nnode_h1D )
    implicit none
    integer, intent(in) :: Npx, Npy, Npz, Ne
    integer, intent(in) :: NeA
    integer, intent(in) :: Nnode_h1D
    real(RP), intent(out) :: q(Npx,Npy,Npz,NeA)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,-Npy+1:2*Npy,Npz,Ne)
    real(RP), intent(in) :: Filter1D(Nnode_h1D,-Npy+1:2*Npy)

    integer :: ke, px, py, pz, p
    real(RP) :: tmp
    !-------------------------------------

    !$omp parallel do private(ke,px,py,pz,p,tmp) collapse(2)
    do ke=1, Ne
    do pz=1, Npz
      do py=1, Npy
      do px=1, Npx
        tmp = 0.0_RP
        do p=-Npy+1, Npy+Npy
          tmp = tmp + Filter1D(py,p) * q0(px,p,pz,ke)
        end do
        q(px,py,pz,ke) = tmp
      end do
      end do
    end do
    end do
    return
  end subroutine MeshFieldFilterOperationBase_apply_filter1d_y

  
!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_apply_reconst1d_y( q, q0, Minv_Ml_tr, Minv_Mc_tr, Minv_Mr_tr, IntrpMat, &
    Npx, Npy, Npz, Ne, Npy_reconst, lmesh )
    implicit none
    class(LocalMeshBase), intent(in) :: lmesh
    integer, intent(in) :: Npx, Npy, Npz, Ne
    integer, intent(in) :: Npy_reconst
    real(RP), intent(out) :: q(Npx,Npy,Npz,lmesh%NeA)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,-Npy+1:2*Npy,Npz,Ne)
    real(RP), intent(in) :: Minv_Ml_tr(Npy,Npy_reconst)
    real(RP), intent(in) :: Minv_Mc_tr(Npy,Npy_reconst)
    real(RP), intent(in) :: Minv_Mr_tr(Npy,Npy_reconst)
    real(RP), intent(in) :: IntrpMat(Npy,Npy_reconst)

    integer :: ke, pz, py, px

    integer :: i,j
    real(RP) :: tmp(Npy_reconst)
    real(RP) :: s
    real(RP) :: q0_L(Npy), q0_C(Npy), q0_R(Npy)
    !-------------------------------------

    !$omp parallel do private(ke,px,py,pz,i,j, tmp,s, &
    !$omp q0_L,q0_C,q0_R) collapse(2)
    do ke=1, Ne
    do pz=1, Npz
    do px=1, Npx
      do j=1, Npy
        q0_L(j) = q0(px,j-Npy,pz,ke)
        q0_C(j) = q0(px,j    ,pz,ke)
        q0_R(j) = q0(px,j+Npy,pz,ke)
      end do
      do i=1, Npy_reconst
        s = 0.0_RP
        do j=1, Npy
          s = s &
            + Minv_Ml_tr(j,i) * q0_L(j) &
            + Minv_Mc_tr(j,i) * q0_C(j) &
            + Minv_Mr_tr(j,i) * q0_R(j)
        end do
        tmp(i) = s
      end do
      q(px,:,pz,ke) = matmul(IntrpMat, tmp)
    end do
    end do
    end do
    return
  end subroutine MeshFieldFilterOperationBase_apply_reconst1d_y

!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_apply_reconst1d_y_2( q, q0, Minv_Ml_tr, Minv_Mc_tr, Minv_Mr_tr, &
    Npx, Npy, Npz, Ne, Npy_reconst, lmesh )
    implicit none
    class(LocalMeshBase), intent(in) :: lmesh
    integer, intent(in) :: Npx, Npy, Npz, Ne
    integer, intent(in) :: Npy_reconst
    real(RP), intent(out) :: q(Npx,Npy,Npz,lmesh%NeA)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,-Npy+1:2*Npy,Npz,Ne)
    real(RP), intent(in) :: Minv_Ml_tr(Npy,Npy)
    real(RP), intent(in) :: Minv_Mc_tr(Npy,Npy)
    real(RP), intent(in) :: Minv_Mr_tr(Npy,Npy)

    integer :: ke, pz, py, px

    integer :: i,j
    real(RP) :: tmp(Npy)
    real(RP) :: s
    real(RP) :: q0_L(Npy), q0_C(Npy), q0_R(Npy)
    !-------------------------------------

    !$omp parallel do private(ke,px,py,pz,i,j, tmp,s, &
    !$omp q0_L,q0_C,q0_R) collapse(2)
    do ke=1, Ne
    do pz=1, Npz
    do px=1, Npx
      do j=1, Npy
        q0_L(j) = q0(px,j-Npy,pz,ke)
        q0_C(j) = q0(px,j    ,pz,ke)
        q0_R(j) = q0(px,j+Npy,pz,ke)
      end do
      do i=1, Npy
        s = 0.0_RP
        do j=1, Npy
          s = s &
            + Minv_Ml_tr(j,i) * q0_L(j) &
            + Minv_Mc_tr(j,i) * q0_C(j) &
            + Minv_Mr_tr(j,i) * q0_R(j)
        end do
        tmp(i) = s
      end do
      q(px,:,pz,ke) = tmp(:)
    end do
    end do
    end do
    return
  end subroutine MeshFieldFilterOperationBase_apply_reconst1d_y_2

!-- Private routines -----------------

!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_prepair_filter_matrix( this, Nnode_h1D, FilterShape, FilterWidthFac )
    use scale_polynomial, only: &
      Polynomial_GenLagrangePoly, &
      Polynomial_GenGaussLegendrePt, &
      Polynomial_GenGaussLegendrePtIntWeight
    use scale_element_line, only: LineElement
    implicit none
    class(MeshFieldFilterOperationBase), intent(inout) :: this
    integer, intent(in) :: Nnode_h1D
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: FilterWidthFac

    integer :: p1, p2

    real(RP) :: filterW
    real(RP) :: filterW2
    type(LineElement) :: elem1D
    type(LineElement) :: elem1D_intrp

    real(RP), allocatable :: lag(:,:)
    real(RP), allocatable :: filter_func(:)

    integer :: NIntNode
    real(RP), allocatable :: r_int1D(:)
    real(RP), allocatable :: w_int1D(:)

    integer :: PolyOrder_h
    !--------------------------

    PolyOrder_h = Nnode_h1D - 1
    call elem1D%Init( PolyOrder_h, .false. )
    call elem1D_intrp%Init( min(2*PolyOrder_h, 11), .false. )
    
    ! LOG_INFO("prepair_filter_matrix Pos:",*) elem1D%x1(:)
    ! LOG_INFO("prepair_filter_matrix IntWeight:",*) elem1D%IntWeight_lgl(:)

    allocate( this%FilterMat_h1D(elem1D%Np,-elem1D%Np+1:elem1D%Np+elem1D%Np) )

    NIntNode = min(2*PolyOrder_h, 11)
    ! NGLnodes = elem1D_intrp%Np
    
    allocate( r_int1D(NIntNode), w_int1D(NIntNode) )
    r_int1D(:) = Polynomial_GenGaussLegendrePt( NIntNode )
    w_int1D(:) = Polynomial_GenGaussLegendrePtIntWeight( NIntNode )
    ! r_int1D(:) = elem1D_intrp%x1(:)
    ! w_int1D(:) = elem1D_intrp%IntWeight_lgl(:)

    allocate( filter_func(NIntNode) )

    allocate( lag(NIntNode,elem1D%PolyOrder+1) )
    lag(:,:) = Polynomial_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, r_int1D(:) )

    this%FilterMat_h1D(:,:) = 0.0_RP

    filterW = FilterWidthFac * 2.0_RP / real(elem1D%Np,kind=RP)
    this%hHaloSize = elem1D%Np

    ! LOG_INFO("prepair_filter_matrix FilterW:",*) filterW

    do p1=1, elem1D%Np
      !--
      call calc_filter_kenrnel( filter_func, &
        FilterShape, r_int1D(:) - elem1D%x1(p1), filterW, NIntNode )      
      do p2=1, elem1D%Np
        this%FilterMat_h1D(p1,p2) = sum( w_int1D(:) * lag(:,p2) * filter_func(:) )
      end do
      filterW2 = sum( w_int1D(:) * filter_func(:) )

      !--
      call calc_filter_kenrnel( filter_func, &
        FilterShape, r_int1D(:) - 2.0_RP - elem1D%x1(p1), filterW, NIntNode )
      
      ! this%FilterMat_h1D(p1,0) = sum( elem1D_intrp%IntWeight_lgl(:) * filter_func(:) )
      do p2=1, elem1D%Np
        this%FilterMat_h1D(p1,-elem1D%Np+p2) = sum( w_int1D(:) * lag(:,p2) * filter_func(:) )
      end do
      filterW2 = filterW2 + sum( w_int1D(:) * filter_func(:) )

      !--
      call calc_filter_kenrnel( filter_func, &
        FilterShape, r_int1D(:) + 2.0_RP - elem1D%x1(p1), filterW, NIntNode )
      
      ! this%FilterMat_h1D(p1,elem1D%Np+1) = sum( elem1D_intrp%IntWeight_lgl(:) * filter_func(:) )
      do p2=1, elem1D%Np
        this%FilterMat_h1D(p1,elem1D%Np+p2) = sum( w_int1D(:) * lag(:,p2) * filter_func(:) )
      end do
      filterW2 = filterW2 + sum( w_int1D(:) * filter_func(:) )
      
      ! LOG_INFO("prepair_filter_matrix",*) p1, ":", filterW2, ":", this%FilterMat_h1D(p1,:)
      this%FilterMat_h1D(p1,:) = this%FilterMat_h1D(p1,:) / filterW2 ! normalization
    end do

    call elem1D%Final()
    call elem1D_intrp%Final()
    return
  end subroutine MeshFieldFilterOperationBase_prepair_filter_matrix

!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_prepair_reconstruct_matrix( this, Nnode_h1D, Nnode_h1D_reconst )
    use scale_polynomial, only: &
      Polynomial_GenLagrangePoly, &
      Polynomial_GenGaussLegendrePt, &
      Polynomial_GenGaussLegendrePtIntWeight
    use scale_element_line, only: LineElement
    use scale_element_modalfilter, only: ModalFilter
    implicit none
    class(MeshFieldFilterOperationBase), intent(inout) :: this
    integer, intent(in) :: Nnode_h1D
    integer, intent(in) :: Nnode_h1D_reconst

    integer :: p1, p2

    type(LineElement) :: elem1D
    type(LineElement) :: elem1D_reconst

    real(RP), allocatable :: lagr_l(:,:)
    real(RP), allocatable :: lagr_c(:,:)
    real(RP), allocatable :: lagr_r(:,:)

    real(RP), allocatable :: rec_lagr_l(:,:)
    real(RP), allocatable :: rec_lagr_c(:,:)
    real(RP), allocatable :: rec_lagr_r(:,:)

    integer :: NIntNode
    real(RP), allocatable :: r_int1D(:)
    real(RP), allocatable :: w_int1D(:)

    real(RP) :: M_h1D_l(Nnode_h1D_reconst,Nnode_h1D)
    real(RP) :: M_h1D_c(Nnode_h1D_reconst,Nnode_h1D)
    real(RP) :: M_h1D_r(Nnode_h1D_reconst,Nnode_h1D)

    real(RP) :: Minv(Nnode_h1D_reconst,Nnode_h1D_reconst)
    real(RP) :: tmpM(Nnode_h1D_reconst,Nnode_h1D)

    real(RP), allocatable :: x_int(:)
    real(RP) :: x_c(Nnode_h1D)

    integer :: PolyOrder
    integer :: PolyOrder_reconst
    type(ModalFilter) :: modalFilter1D
    !--------------------------

    PolyOrder = Nnode_h1D - 1
    PolyOrder_reconst = Nnode_h1D_reconst - 1

    call elem1D%Init( PolyOrder, .false. )    
    call elem1D_reconst%Init( PolyOrder_reconst, .false. )    
    call modalFilter1D%Init(elem1D_reconst, 0.0_RP, 2E1_RP, 16)

    this%hHaloSize = elem1D%Np

    ! NIntNode = elem1D%Np 
    NIntNode = ceiling( 0.5_RP * ( PolyOrder + PolyOrder_reconst ) ) + 1

    allocate( r_int1D(NIntNode), w_int1D(NIntNode) )
    allocate( x_int(NIntNode) )
    r_int1D(:) = Polynomial_GenGaussLegendrePt( NIntNode )
    w_int1D(:) = Polynomial_GenGaussLegendrePtIntWeight( NIntNode )

    !-
    allocate( lagr_l(NIntNode,elem1D%PolyOrder+1), rec_lagr_l(NIntNode,PolyOrder_reconst+1) )
    allocate( lagr_c(NIntNode,elem1D%PolyOrder+1), rec_lagr_c(NIntNode,PolyOrder_reconst+1) )
    allocate( lagr_r(NIntNode,elem1D%PolyOrder+1), rec_lagr_r(NIntNode,PolyOrder_reconst+1) )
    
    !-
!    x_int(:) = -1.0_RP + 0.25_RP * ( 1.0_RP + r_int1D(:) )
    x_int(:) = -1.0_RP + (1.0_RP / 3.0_RP ) * ( 1.0_RP + r_int1D(:) )
    rec_lagr_l(:,:) = Polynomial_GenLagrangePoly( PolyOrder_reconst, elem1D_reconst%x1, x_int(:) )

!    x_int(:) = 0.0_RP + 0.5_RP * ( 1.0_RP + r_int1D(:) )
    x_int(:) = r_int1D(:)
    lagr_l(:,:) = Polynomial_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, x_int(:) )

    !-
!    x_int(:) = -0.5_RP + 0.5_RP * ( 1.0_RP + r_int1D(:) )
    x_int(:) = -1.0_RP/3.0_RP + (1.0_RP / 3.0_RP ) * ( 1.0_RP + r_int1D(:) )
    rec_lagr_c(:,:) = Polynomial_GenLagrangePoly( PolyOrder_reconst, elem1D_reconst%x1, x_int(:) )

    x_int(:) = r_int1D(:)
    lagr_c(:,:) = Polynomial_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, x_int(:) )

    !-
!    x_int(:) = 0.5_RP + 0.25_RP * ( 1.0_RP + r_int1D(:) )
    x_int(:) = +1.0_RP/3.0_RP + (1.0_RP / 3.0_RP ) * ( 1.0_RP + r_int1D(:) )
    rec_lagr_r(:,:) = Polynomial_GenLagrangePoly( PolyOrder_reconst, elem1D_reconst%x1, x_int(:) )

    ! x_int(:) = -1.0_RP + 0.5_RP * ( 1.0_RP + r_int1D(:) )
    x_int(:) = r_int1D(:)
    lagr_r(:,:) = Polynomial_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, x_int(:) )

    !$omp parallel do collapse(2)
    do p2=1, elem1D%Np
      do p1=1, elem1D_reconst%Np
        ! M_h1D_l(p1,p2) = 0.25_RP * sum( w_int1D(:) * lagr_l(:,p2) * rec_lagr_l(:,p1) )
        ! M_h1D_c(p1,p2) = 0.5_RP  * sum( w_int1D(:) * lagr_c(:,p2) * rec_lagr_c(:,p1) )
        ! M_h1D_r(p1,p2) = 0.25_RP * sum( w_int1D(:) * lagr_r(:,p2) * rec_lagr_r(:,p1) )
        M_h1D_l(p1,p2) = 1.0_RP/3.0_RP * sum( w_int1D(:) * lagr_l(:,p2) * rec_lagr_l(:,p1) )
        M_h1D_c(p1,p2) = 1.0_RP/3.0_RP * sum( w_int1D(:) * lagr_c(:,p2) * rec_lagr_c(:,p1) )
        M_h1D_r(p1,p2) = 1.0_RP/3.0_RP * sum( w_int1D(:) * lagr_r(:,p2) * rec_lagr_r(:,p1) )
      end do
    end do
    
    allocate( this%Minv_Ml_tr(elem1D%Np,Nnode_h1D_reconst) )
    allocate( this%Minv_Mc_tr(elem1D%Np,Nnode_h1D_reconst) )
    allocate( this%Minv_Mr_tr(elem1D%Np,Nnode_h1D_reconst) )
    allocate( this%IntrpMat(elem1D%Np,Nnode_h1D_reconst) )

    Minv(:,:) = elem1D_reconst%invM(:,:)
    ! Minv(:,:) = matmul( modalFilter1D%FilterMat, Minv )
    
    tmpM(:,:) = matmul(Minv, M_h1D_l)
    this%Minv_Ml_tr(:,:) = transpose(tmpM)

    tmpM(:,:) = matmul(Minv, M_h1D_c)
    this%Minv_Mc_tr(:,:) = transpose(tmpM)

    tmpM(:,:) = matmul(Minv, M_h1D_r)
    this%Minv_Mr_tr(:,:) = transpose(tmpM)

!    x_c(:) = -0.5_RP + 0.5_RP * ( 1.0_RP + elem1D%x1(:) )
    x_c(:) = - 1.0_RP/3.0_RP + (1.0_RP/3.0_RP) * ( 1.0_RP + elem1D%x1(:) )
    this%IntrpMat(:,:) = Polynomial_GenLagrangePoly( elem1D_reconst%PolyOrder, elem1D_reconst%x1, x_c(:) )

    !-
    call elem1D%Final()
    call elem1D_reconst%Final()
    call modalFilter1D%Final()
    return
  end subroutine MeshFieldFilterOperationBase_prepair_reconstruct_matrix

!OCL SERIAL
  subroutine MeshFieldFilterOperationBase_prepair_reconstruct2_matrix( this, Nnode_h1D, Nnode_h1D_reconst )
    use scale_const, only: &
      EPS => CONST_EPS
    use scale_polynomial, only: &
      Polynomial_GenLagrangePoly, &
      Polynomial_GenGaussLegendrePt, &
      Polynomial_GenGaussLegendrePtIntWeight
    use scale_element_line, only: LineElement
    use scale_element_modalfilter, only: ModalFilter
    implicit none
    class(MeshFieldFilterOperationBase), intent(inout) :: this
    integer, intent(in) :: Nnode_h1D
    integer, intent(in) :: Nnode_h1D_reconst

    integer :: p0, p1, p2
    real(RP) :: x0
    real(RP) :: xr_rec, xl_rec
    real(RP) :: xr, xl
    real(RP) :: coef_l, coef_c, coef_r

    type(LineElement) :: elem1D
    type(LineElement) :: elem1D_reconst

    real(RP), allocatable :: lagr_l(:,:)
    real(RP), allocatable :: lagr_c(:,:)
    real(RP), allocatable :: lagr_r(:,:)

    real(RP), allocatable :: rec_lagr_l(:,:)
    real(RP), allocatable :: rec_lagr_c(:,:)
    real(RP), allocatable :: rec_lagr_r(:,:)

    integer :: NIntNode
    real(RP), allocatable :: r_int1D(:)
    real(RP), allocatable :: w_int1D(:)

    real(RP) :: M_h1D_l(Nnode_h1D_reconst,Nnode_h1D)
    real(RP) :: M_h1D_c(Nnode_h1D_reconst,Nnode_h1D)
    real(RP) :: M_h1D_r(Nnode_h1D_reconst,Nnode_h1D)

    real(RP) :: Minv(Nnode_h1D_reconst,Nnode_h1D_reconst)

    real(RP) :: Minv_Ml(Nnode_h1D_reconst,Nnode_h1D,Nnode_h1D)
    real(RP) :: Minv_Mc(Nnode_h1D_reconst,Nnode_h1D,Nnode_h1D)
    real(RP) :: Minv_Mr(Nnode_h1D_reconst,Nnode_h1D,Nnode_h1D)

    real(RP) :: IntrpMat(1,Nnode_h1D_reconst)
    real(RP) :: tmpMat(1,Nnode_h1D)

    real(RP), allocatable :: x_int(:)
    real(RP) :: x_c(1)

    integer :: PolyOrder
    integer :: PolyOrder_reconst
    type(ModalFilter) :: modalFilter1D
    !--------------------------

    PolyOrder = Nnode_h1D - 1
    PolyOrder_reconst = Nnode_h1D_reconst - 1

    call elem1D%Init( PolyOrder, .false. )    
    call elem1D_reconst%Init( PolyOrder_reconst, .false. )    
    call modalFilter1D%Init(elem1D_reconst, 0.0_RP, 1D3, 16)

    this%hHaloSize = elem1D%Np

    ! NIntNode = elem1D%Np 
    NIntNode = ceiling( 0.5_RP * ( PolyOrder + PolyOrder_reconst ) ) + 1

    allocate( r_int1D(NIntNode), w_int1D(NIntNode) )
    allocate( x_int(NIntNode) )
    r_int1D(:) = Polynomial_GenGaussLegendrePt( NIntNode )
    w_int1D(:) = Polynomial_GenGaussLegendrePtIntWeight( NIntNode )

    !-
    allocate( lagr_l(NIntNode,elem1D%PolyOrder+1), rec_lagr_l(NIntNode,PolyOrder_reconst+1) )
    allocate( lagr_c(NIntNode,elem1D%PolyOrder+1), rec_lagr_c(NIntNode,PolyOrder_reconst+1) )
    allocate( lagr_r(NIntNode,elem1D%PolyOrder+1), rec_lagr_r(NIntNode,PolyOrder_reconst+1) )
    
    !-
    Minv(:,:) = elem1D_reconst%invM(:,:)
    ! Minv(:,:) = matmul( modalFilter1D%FilterMat, Minv )

    do p0=1, elem1D%Np
      x0 = elem1D%x1(p0)

      xl_rec = - 1.0_RP
      xr_rec = xl_rec + 0.5_RP * ( 1.0_RP - x0 )
      coef_l = 0.5_RP * ( xr_rec - xl_rec )
      x_int(:) = xl_rec + coef_l * ( 1.0_RP + r_int1D(:) )
      rec_lagr_l(:,:) = Polynomial_GenLagrangePoly( PolyOrder_reconst, elem1D_reconst%x1, x_int(:) )

      xl = x0
      xr = 1.0_RP
      x_int(:) = xl + 0.5_RP * ( xr - xl ) * ( 1.0_RP + r_int1D(:) )
      lagr_l(:,:) = Polynomial_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, x_int(:) )

      !-
      xl_rec = xr_rec
      xr_rec = xl_rec + 0.5_RP * 2.0_RP
      coef_c = 0.5_RP * ( xr_rec - xl_rec )
      x_int(:) = xl_rec + coef_c * ( 1.0_RP + r_int1D(:) )
      rec_lagr_c(:,:) = Polynomial_GenLagrangePoly( PolyOrder_reconst, elem1D_reconst%x1, x_int(:) )

      x_int(:) = r_int1D(:)
      lagr_c(:,:) = Polynomial_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, x_int(:) )

      !-
      xl_rec = xr_rec
      xr_rec = xl_rec + 0.5_RP * ( x0 + 1.0_RP )
      coef_r = 0.5_RP * ( xr_rec - xl_rec )
      x_int(:) = xl_rec + coef_r * ( 1.0_RP + r_int1D(:) )
      rec_lagr_r(:,:) = Polynomial_GenLagrangePoly( PolyOrder_reconst, elem1D_reconst%x1, x_int(:) )

      xl = -1.0_RP
      xr = x0
      x_int(:) = xl + 0.5_RP * ( xr - xl ) * ( 1.0_RP + r_int1D(:) )
      lagr_r(:,:) = Polynomial_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, x_int(:) )

      !$omp parallel do collapse(2)
      do p2=1, elem1D%Np
        do p1=1, elem1D_reconst%Np
          M_h1D_l(p1,p2) = coef_l * sum( w_int1D(:) * lagr_l(:,p2) * rec_lagr_l(:,p1) )
          M_h1D_c(p1,p2) = coef_c * sum( w_int1D(:) * lagr_c(:,p2) * rec_lagr_c(:,p1) )
          M_h1D_r(p1,p2) = coef_r * sum( w_int1D(:) * lagr_r(:,p2) * rec_lagr_r(:,p1) )
        end do
      end do

      Minv_Ml(:,:,p0) = matmul( Minv, M_h1D_l )
      Minv_Mc(:,:,p0) = matmul( Minv, M_h1D_c )
      Minv_Mr(:,:,p0) = matmul( Minv, M_h1D_r )
    end do

    x_c(:) = 0.0_RP
    IntrpMat(:,:) = Polynomial_GenLagrangePoly( elem1D_reconst%PolyOrder, elem1D_reconst%x1, x_c(:) )

    allocate( this%Ml_tr(elem1D%Np,elem1D%Np) )
    allocate( this%Mc_tr(elem1D%Np,elem1D%Np) )
    allocate( this%Mr_tr(elem1D%Np,elem1D%Np) )

    do p0=1, Nnode_h1D
      tmpMat(:,:) = matmul(IntrpMat, Minv_Ml(:,:,p0))
      this%Ml_tr(:,p0) = tmpMat(1,:)

      tmpMat(:,:) = matmul(IntrpMat, Minv_Mc(:,:,p0))
      this%Mc_tr(:,p0) = tmpMat(1,:)

      tmpMat(:,:) = matmul(IntrpMat, Minv_Mr(:,:,p0))
      this%Mr_tr(:,p0) = tmpMat(1,:)      
    end do

    !-
    call elem1D%Final()
    call elem1D_reconst%Final()
    call modalFilter1D%Final()
    return
  end subroutine MeshFieldFilterOperationBase_prepair_reconstruct2_matrix

!OCL SERIAL
  subroutine calc_filter_kenrnel( filter_kernel, &
      FilterShape, x, filter_width, Np )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: filter_kernel(Np)
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: x(Np)
    real(RP), intent(in) :: filter_width
    !------------------------------------------------

    select case(FilterShape)
    case( "GAUSSIAN")
      filter_kernel(:) = exp( - (x(:)/filter_width)**2 )
    case( "TOPHAT")
      filter_kernel(:) = 0.5_RP * ( sign(1.0_RP, x(:) + 0.5_RP * filter_width) - sign(1.0_RP, x(:) - 0.5_RP * filter_width) )
    case default
      LOG_ERROR("MeshFieldFilterOperation3D_calc_filter_kernel",*) "The specified FilterShape is not supported. Check!", FilterShape
      call PRC_abort
    end select
    return
  end subroutine calc_filter_kenrnel

end module scale_meshfield_filter_operation_base