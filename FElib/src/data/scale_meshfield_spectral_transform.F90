!> module FElib / Data / Utility
!!
!! @par Description
!!           A module for providing utility routines associated with spectral transformation for DG fields
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_meshfield_spectral_transform
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_const, only: &
    PI => CONST_PI
  use scale_prc, only: &
    PRC_abort
  
  use scale_polynominal, only: &
    Polynominal_GenLagrangePoly
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_element_line, only: &
      LineElement
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: &
    MeshField1D, MeshField2D, &
    MeshField1DList, MeshField2DList

  use scale_fast_fourier_transform, only: &
    FastFourierTransform1D
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  type, public :: MeshField_SpetralTransformBase
    integer :: eval_type_id
    integer :: ndim
    integer :: var_num

    real(RP), allocatable :: FFTIntrpMat(:,:)
    real(RP), allocatable :: FFT_xi(:)
    integer :: NsamplePtPerElem

    real(RP), allocatable :: IntIntrpMat(:,:)
    real(RP), allocatable :: IntWeight(:)
    real(RP), allocatable :: IntXi(:)
    integer :: NintGLpt
  end type MeshField_SpetralTransformBase

  !> A derived type for spectral transform of MeshField1D
  type, public, extends(MeshField_SpetralTransformBase) :: MeshField_SpetralTransform1D
    integer :: kall
    integer :: ks, ke

    real(RP) :: xmin_gl, xmax_gl
    real(RP) :: delx

    real(RP), allocatable :: k(:)
    real(RP), allocatable :: spectral_coef(:,:,:)

    type(FastFourierTransform1D) :: fft
  contains
    procedure :: Init => MeshField_SpetralTransform1D_Init
    procedure :: Final => MeshField_SpetralTransform1D_Final
    procedure :: Transform => MeshField_SpetralTransform1D_transform
  end type MeshField_SpetralTransform1D

  !> A derived type for spectral transform of MeshField2D
  type, public, extends(MeshField_SpetralTransformBase) :: MeshField_SpetralTransform2D
    integer :: kall
    integer :: ks, ke
    integer :: lall
    integer :: ls, le

    real(RP) :: xmin_gl, xmax_gl
    real(RP) :: delx
    real(RP) :: ymin_gl, ymax_gl
    real(RP) :: dely

    integer :: NprcX    
    integer :: NprcY
    integer :: NeGX, NeGY

    integer :: NsamplePtPerElem1D

    real(RP), allocatable :: k(:)
    real(RP), allocatable :: l(:)
    real(RP), allocatable :: spectral_coef(:,:,:,:)

    type(FastFourierTransform1D) :: fft_x
    type(FastFourierTransform1D) :: fft_y
  contains
    procedure :: Init => MeshField_SpetralTransform2D_Init
    procedure :: Final => MeshField_SpetralTransform2D_Final
    procedure :: Transform => MeshField_SpetralTransform2D_transform
  end type MeshField_SpetralTransform2D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: ST_EVALTYPE_SAMPLE_UNIFORM_PTS = 1
  integer, public, parameter :: ST_EVALTYPE_L2PROJECTION_1     = 2
  integer, public, parameter :: ST_EVALTYPE_L2PROJECTION_2     = 3

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
!- Base type

!OCL SERIAL
  subroutine MeshField_SpetralTransformBase_Init( this, eval_type, ndim, var_num, NintGLpt, NsamplePtPerElem )
    implicit none
    class(MeshField_SpetralTransformBase), intent(inout) :: this
    integer, intent(in) :: eval_type
    integer, intent(in) :: ndim
    integer, intent(in) :: var_num
    integer, intent(in), optional :: NintGLpt
    integer, intent(in), optional :: NsamplePtPerElem
    !--------------------------------------------
    
    this%eval_type_id = eval_type
    call check_eval_type_id( eval_type )

    this%ndim = ndim
    this%var_num = var_num

    if (     this%eval_type_id == ST_EVALTYPE_L2PROJECTION_1 &
        .or. this%eval_type_id == ST_EVALTYPE_L2PROJECTION_2 ) then
      if ( .not. present(NintGLpt) ) then
        LOG_INFO("MeshField_SpetralTransformBase_Init",*) "The order of Gaussian quadrature should be given. Check!"
        call PRC_abort
      end if
      this%NintGLPt = NintGLpt
    else
      this%NintGLPt = -1
    end if

    if ( this%eval_type_id == ST_EVALTYPE_SAMPLE_UNIFORM_PTS ) then
      if ( .not. present(NsamplePtPerElem) ) then
        LOG_INFO("MeshField_SpetralTransformBase_Init",*) "NsamplePtPerElem shouled be given. Check!"
      else
        this%NsamplePtPerElem = NsamplePtPerElem
      end if
    else
      this%NsamplePtPerElem = -1
    end if
    LOG_INFO("MeshField_SpetralTransformBase_Init",*) "NintGLPt=", this%NintGLpt, "NsamplePtPerElem=", this%NsamplePtPerElem

    return
  end subroutine MeshField_SpetralTransformBase_Init

  subroutine check_eval_type_id( eval_type ) 
    implicit none
    integer, intent(in) :: eval_type
    !-------------------------
    select case(eval_type)
!    case(ST_EVALTYPE_L2PROJECTION_1,ST_EVALTYPE_L2PROJECTION_2,ST_EVALTYPE_SAMPLE_UNIFORM_PTS)
    case( ST_EVALTYPE_L2PROJECTION_1, ST_EVALTYPE_SAMPLE_UNIFORM_PTS )
    case default
      LOG_INFO("MeshField_SpetralTransformBase_Init",*) "Unsupported evaluation type is specified. Check!"
      call PRC_abort
    end select
    return
  end subroutine check_eval_type_id

!OCL SERIAL
  subroutine MeshField_SpetralTransformBase_Final( this )
    implicit none
    class(MeshField_SpetralTransformBase), intent(inout) :: this
    !--------------------------------------------

    if ( this%NintGLpt > 0 ) &
      deallocate( this%IntIntrpMat, this%IntXi, this%IntWeight )
    
    if ( this%NsamplePtPerElem > 0 ) &
      deallocate( this%FFTIntrpMat, this%FFT_xi )

    return
  end subroutine MeshField_SpetralTransformBase_Final  

!- 1D

!OCL SERIAL
  subroutine MeshField_SpetralTransform1D_Init( this, eval_type, ks, ke, mesh1D, var_num, &
    GLQuadOrd, NsamplePtPerElem )
    implicit none
    class(MeshField_SpetralTransform1D), intent(inout) :: this
    integer, intent(in) :: eval_type
    integer, intent(in) :: ks, ke
    class(MeshBase1D), intent(in), target :: mesh1D
    integer, intent(in) :: var_num
    integer, optional :: GLQuadOrd
    integer, optional :: NsamplePtPerElem

    integer :: i
    integer :: Lx
    real(RP) :: dx

    class(LocalMesh1D), pointer :: lmesh
    class(ElementBase1D), pointer :: elem

    type(LineElement) :: elem_dummy
    !--------------------------------------------

    lmesh => mesh1D%lcmesh_list(1)
    elem => lmesh%refElem1D

    call MeshField_SpetralTransformBase_Init( this, eval_type, 1, var_num, &
      GLQuadOrd, NsamplePtPerElem )

    this%ks = ks; this%ke = ke
    this%kall = ke - ks + 1

    this%xmin_gl = mesh1D%xmin_gl
    this%xmax_gl = mesh1D%xmax_gl
    this%delx = 1.0_RP / real(mesh1D%NeG,kind=RP)

    allocate( this%k(ks:ke) )
    Lx = this%xmax_gl - this%xmin_gl
    do i=ks, ke
      this%k(i) = i * 2.0_RP * PI / Lx
    end do
    
    allocate( this%spectral_coef(ks:ke,2,var_num) )

    call elem_dummy%Init(elem%PolyOrder, .false.)

    if ( this%NintGLpt > 0 ) then
      allocate( this%IntIntrpMat(this%NintGLpt,elem%Np) )
      allocate( this%IntXi(this%NintGLpt), this%IntWeight(this%NintGLpt) )

      this%IntIntrpMat(:,:) = elem_dummy%GenIntGaussLegendreIntrpMat( GLQuadOrd, this%IntWeight, this%IntXi )
    end if

    if ( this%NsamplePtPerElem > 0 ) then
      allocate( this%FFTIntrpMat(this%NsamplePtPerElem,elem%Np) )
      allocate( this%FFT_Xi(this%NsamplePtPerElem) )

      dx = 2.0_RP / real(this%NsamplePtPerElem,kind=RP)
      do i=1, this%NsamplePtPerElem
        this%FFT_Xi(i) = -1.0_RP + real(i - 0.5_RP,kind=RP) * dx
      end do
      this%FFTIntrpMat(:,:) = Polynominal_GenLagrangePoly( elem_dummy%PolyOrder, elem_dummy%x1, this%FFT_xi )

      call this%fft%Init( NsamplePtPerElem * mesh1D%NeG )
    end if

    call elem_dummy%Final()

    return
  end subroutine MeshField_SpetralTransform1D_Init

!OCL SERIAL
  subroutine MeshField_SpetralTransform1D_Final( this )
    implicit none
    class(MeshField_SpetralTransform1D), intent(inout) :: this
    !--------------------------------------------

    call MeshField_SpetralTransformBase_Final( this )

    if ( this%NsamplePtPerElem > 0 ) then
      call this%fft%Final()
    end if

    deallocate( this%k, this%spectral_coef )
    
    return
  end subroutine MeshField_SpetralTransform1D_Final

!OCL SERIAL
  subroutine MeshField_SpetralTransform1D_transform( this, q_list, mesh_num )
    implicit none
    class(MeshField_SpetralTransform1D), intent(inout) :: this
    integer, intent(in) :: mesh_num
    type(MeshField1DList), target :: q_list(this%var_num,mesh_num)

    class(MeshBase1D), pointer :: mesh1D
    class(LocalMesh1D), pointer :: lmesh
    !--------------------------------------------

    mesh1D => q_list(1,1)%ptr%mesh
    lmesh => mesh1D%lcmesh_list(1)
    
    select case( this%eval_type_id )
    case( ST_EVALTYPE_SAMPLE_UNIFORM_PTS ) 
      call spectral_transform1D_DFT( this%spectral_coef, & ! (out)
        q_list, this%var_num, mesh_num,                                        & ! (in)
        this%ks, this%ke, lmesh%refElem1D%Np, lmesh%Ne, this%NsamplePtPerElem, & ! (in)
        this%fft, this%FFTIntrpMat, this%FFT_xi,                               & ! (in)
        this%xmin_gl, this%delx, this%xmax_gl - this%xmin_gl                   ) ! (in)
    case( ST_EVALTYPE_L2PROJECTION_1)
      call spectral_transform1D_L2projection( this%spectral_coef,    & ! (out)
        q_list, this%var_num, mesh_num,                              & ! (in)
        this%ks, this%ke, lmesh%refElem1D%Np, lmesh%Ne,              & ! (in)
        this%IntIntrpMat, this%IntXi, this%IntWeight, this%NintGLPt, & ! (in)
        this%xmin_gl, this%delx, this%xmax_gl - this%xmin_gl         ) ! (in)
    end select
    
    return
  end subroutine MeshField_SpetralTransform1D_transform

!- 2D

!OCL SERIAL
  subroutine MeshField_SpetralTransform2D_Init( this, eval_type, ks, ke, ls, le, mesh2D, var_num, &
    GLQuadOrd, NsamplePtPerElem1D )
    use scale_element_line, only: &
      LineElement
    use scale_mesh_rectdom2d, only: &
      MeshRectDom2D
    implicit none
    class(MeshField_SpetralTransform2D), intent(inout) :: this
    integer, intent(in) :: eval_type
    integer, intent(in) :: ks, ke
    integer, intent(in) :: ls, le
    class(MeshBase2D), intent(in), target :: mesh2D
    integer, intent(in) :: var_num
    integer, optional :: GLQuadOrd
    integer, optional :: NsamplePtPerElem1D

    integer :: i, j
    integer :: Lx, Ly
    real(RP) :: dx

    class(LocalMesh2D), pointer :: lmesh
    class(ElementBase2D), pointer :: elem

    type(LineElement) :: elem_dummy
    class(MeshBase2D), pointer :: ptr_mesh2D
    !--------------------------------------------

    lmesh => mesh2D%lcmesh_list(1)
    elem => lmesh%refElem2D

    call MeshField_SpetralTransformBase_Init( this, eval_type, 1, var_num, GLQuadOrd, NsamplePtPerElem1D )

    this%ks = ks; this%ke = ke
    this%kall = ke - ks + 1

    this%ls = ls; this%le = le
    this%lall = le - ls + 1

    select type(ptr_mesh2D => mesh2D)
    class is (MeshRectDom2D)
      this%xmin_gl = ptr_mesh2D%xmin_gl
      this%xmax_gl = ptr_mesh2D%xmax_gl
      this%delx = ( this%xmax_gl - this%xmin_gl ) / real(ptr_mesh2D%NeGX, kind=RP)

      this%ymin_gl = ptr_mesh2D%ymin_gl
      this%ymax_gl = ptr_mesh2D%ymax_gl
      this%dely = ( this%ymax_gl - this%ymin_gl ) / real(ptr_mesh2D%NeGY, kind=RP)

      this%NeGX = ptr_mesh2D%NeGX
      this%NeGY = ptr_mesh2D%NeGY
      this%NprcX = ptr_mesh2D%NprcX
      this%NprcY = ptr_mesh2D%NprcY
    class default
      LOG_INFO('MeshField_SpetralTransform2D_Init',*) 'Unexpected mesh type is given. Check!'
      call PRC_abort
    end select

    allocate( this%k(ks:ke) )
    Lx = this%xmax_gl - this%xmin_gl
    do i=ks, ke
      this%k(i) = i * 2.0_RP * PI / Lx
    end do

    allocate( this%l(ls:le) )
    Ly = this%ymax_gl - this%ymin_gl
    do j=ls, le
      this%l(j) = j * 2.0_RP * PI / Ly
    end do
    
    allocate( this%spectral_coef(ks:ke,ls:le,2,var_num) )

    call elem_dummy%Init(elem%PolyOrder, .false.)

    if ( this%NintGLpt > 0 ) then
      allocate( this%IntIntrpMat(this%NintGLpt,elem%Nfp) )
      allocate( this%IntXi(this%NintGLpt), this%IntWeight(this%NintGLpt) )

      this%IntIntrpMat(:,:) = elem_dummy%GenIntGaussLegendreIntrpMat( GLQuadOrd, this%IntWeight, this%IntXi )
    end if

    if ( this%NsamplePtPerElem > 0 ) then
      this%NsamplePtPerElem1D = NsamplePtPerElem1D
      allocate( this%FFTIntrpMat(this%NsamplePtPerElem1D,elem%Nfp) )
      allocate( this%FFT_Xi(this%NsamplePtPerElem1D) )

      dx = 2.0_RP / real(this%NsamplePtPerElem1D,kind=RP)
      do i=1, this%NsamplePtPerElem1D
        this%FFT_Xi(i) = -1.0_RP + real(i - 0.5_RP,kind=RP) * dx
      end do
      this%FFTIntrpMat(:,:) = Polynominal_GenLagrangePoly( elem_dummy%PolyOrder, elem_dummy%x1, this%FFT_xi )

      call this%fft_x%Init( this%NsamplePtPerElem1D * this%NeGX )
      call this%fft_y%Init( this%NsamplePtPerElem1D * this%NeGY )
    end if

    call elem_dummy%Final()
    return
  end subroutine MeshField_SpetralTransform2D_Init

!OCL SERIAL
  subroutine MeshField_SpetralTransform2D_Final( this )
    implicit none
    class(MeshField_SpetralTransform2D), intent(inout) :: this
    !--------------------------------------------

    call MeshField_SpetralTransformBase_Final( this )

    if ( this%NsamplePtPerElem > 0 ) then
      call this%fft_x%Final()
      call this%fft_y%Final()
    end if
    
    deallocate( this%k, this%l, this%spectral_coef )
    return
  end subroutine MeshField_SpetralTransform2D_Final

!OCL SERIAL
  subroutine MeshField_SpetralTransform2D_transform( this, q_list, mesh_num_x, mesh_num_y )
    implicit none
    class(MeshField_SpetralTransform2D), intent(inout) :: this
    integer, intent(in) :: mesh_num_x, mesh_num_y
    type(MeshField2DList), target :: q_list(this%var_num,mesh_num_x,mesh_num_y)

    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh2D), pointer :: lmesh
    !--------------------------------------------

    mesh2D => q_list(1,1,1)%ptr%mesh
    lmesh => mesh2D%lcmesh_list(1)
    
    select case( this%eval_type_id )
    case( ST_EVALTYPE_SAMPLE_UNIFORM_PTS )
      call spectral_transform2D_DFT( this%spectral_coef, &
        q_list, this%var_num, mesh_num_x, mesh_num_y,                                   & ! (in)
        this%ks, this%ke,  this%ls, this%le, lmesh%refElem2D%Nfp, lmesh%NeX, lmesh%NeY, & ! (in)
        this%NeGY, this%NprcY, this%NsamplePtPerElem1D,                                 & ! (in)
        this%fft_x, this%fft_y, this%FFTIntrpMat )
    case( ST_EVALTYPE_L2PROJECTION_1)
      call spectral_transform2D_L2projection( this%spectral_coef,                       & ! (out)
        q_list, this%var_num, mesh_num_x*mesh_num_y,                                    & ! (in)
        this%ks, this%ke,  this%ls, this%le, lmesh%refElem2D%Nfp, lmesh%NeX, lmesh%NeY, & ! (in)
        this%IntIntrpMat, this%IntXi, this%IntWeight, this%NintGLPt,                    & ! (in)
        this%xmin_gl, this%delx, this%xmax_gl - this%xmin_gl,                           & ! (in)
        this%ymin_gl, this%dely, this%ymax_gl - this%ymin_gl                            ) ! (in)
    end select
    
    return
  end subroutine MeshField_SpetralTransform2D_transform

!--- Private subroutines (1D)

!OCL SERIAL
  subroutine spectral_transform1D_DFT( spectral_coef, &
    q_list, var_num, mesh_num, ks, ke, Np, Ne, NsamplePerElem, &
    fft, FFTIntrpMat, FFT_xi, &
    xmin_gl, delx, Lx )

    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD    
    implicit none
    integer, intent(in) :: ks, ke
    integer, intent(in) :: Np
    integer, intent(in) :: Ne
    integer, intent(in) :: NsamplePerElem
    integer, intent(in) :: var_num
    integer, intent(in) :: mesh_num
    real(RP), intent(out) :: spectral_coef(ks:ke,2,var_num)
    type(MeshField1DList), target :: q_list(var_num,mesh_num)
    class(FastFourierTransform1D), intent(in) :: fft
    real(RP), intent(in) :: FFTIntrpMat(NsamplePerElem,Np)
    real(RP), intent(in) :: FFT_xi(NsamplePerElem)
    real(RP), intent(in) :: xmin_gl
    real(RP), intent(in) :: delx
    real(RP), intent(in) :: Lx

    real(RP) :: g_q(NsamplePerElem,Ne,mesh_num,var_num)
    real(RP) :: q_tmp(Np)
    complex(RP) :: s_q(NsamplePerElem*Ne*mesh_num,var_num)

    integer :: m
    integer :: v
    integer :: kel

    integer :: Nall
    integer :: kk

    real(RP) :: x_(NsamplePerElem)
    !-----------------------------------------------

    !$omp parallel do collapse(3) private(v,m,kel,q_tmp, x_)
    do v=1, var_num
    do m=1, mesh_num
    do kel=1, Ne
      q_tmp(:) = q_list(v,m)%ptr%local(1)%val(:,kel)
      g_q(:,kel,m,v) = matmul( FFTIntrpMat, q_tmp )
    end do
    end do
    end do

    !$omp parallel do
    do v=1, var_num
      call fft%Forward_real( g_q(:,:,:,v), s_q(:,v) )
    end do

    Nall = NSamplePerElem * Ne * mesh_num
    !$omp parallel do private(kk)
    do v=1, var_num
      do kk=1, Nall/2+1
        spectral_coef(kk-1,1,v) = real (s_q(kk,v))
        spectral_coef(kk-1,2,v) = aimag(s_q(kk,v))
      end do
      do kk=Nall/2+1, Nall
        spectral_coef(kk-1-Nall,1,v) = real (s_q(kk,v))
        spectral_coef(kk-1-Nall,2,v) = aimag(s_q(kk,v))
      end do
    end do
    return
  end subroutine spectral_transform1D_DFT
  
!OCL SERIAL
  subroutine spectral_transform1D_L2projection( spectral_coef, &
    q_list, var_num, mesh_num, ks, ke, Np, Ne, IntIntrpMat, IntXi, IntWeight, NintGLPt, &
    xmin_gl, delx, Lx )

    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD    
    implicit none
    integer, intent(in) :: ks, ke
    integer, intent(in) :: Np
    integer, intent(in) :: Ne
    integer, intent(in) :: var_num
    integer, intent(in) :: mesh_num
    class(MeshField1DList), target :: q_list(var_num,mesh_num)
    real(RP), intent(out) :: spectral_coef(ks:ke,2,var_num)
    integer, intent(in) :: NintGLpt
    real(RP), intent(in) :: IntIntrpMat(NintGLpt,Np)
    real(RP), intent(in) :: IntXi(NintGLpt)
    real(RP), intent(in) :: IntWeight(NintGLpt)
    real(RP), intent(in) :: xmin_gl
    real(RP), intent(in) :: delx
    real(RP), intent(in) :: Lx

    real(RP), allocatable  :: s_coef_lc(:,:,:)
    real(RP), allocatable  :: s_coef(:,:,:)
    real(RP), allocatable :: q_tmp(:,:,:)

    integer :: vec_size
    integer :: km

    class(MeshBase1D), pointer :: mesh1D
    class(LocalMesh1D), pointer :: lmesh
    integer :: ldom
    integer :: meshID
    class(LocalMeshField1D), pointer :: lcfield

    integer :: kel, v
    integer :: k, kk
    integer :: ierr
    !-----------------------------------------------

    vec_size = var_num
    km = ke

    allocate( s_coef_lc(vec_size,2,0:km) )
    allocate( s_coef(vec_size,2,0:km) )
    
    !-
    mesh1D => q_list(1,1)%ptr%mesh
    lmesh => mesh1D%lcmesh_list(1)

    !-

    allocate( q_tmp(Np,vec_size,Ne) )

    s_coef_lc(:,:,:) = 0.0_RP
    
    do meshID=1, mesh_num
      mesh1D => q_list(1,meshID)%ptr%mesh
      do ldom=1, mesh1D%LOCAL_MESH_NUM
        lmesh => mesh1D%lcmesh_list(ldom)
        lcfield => q_list(v,meshID)%ptr%local(ldom)
        !$omp parallel do collapse(2)
        do kel=lmesh%NeS, lmesh%NeE
        do v=1, vec_size
          q_tmp(:,v,kel) = lcfield%val(:,kel)
        end do
        end do      
        call spectral_transform1D_L2projection_lc( s_coef_lc, &
          q_tmp, 0, km, Np, Ne, vec_size, IntIntrpMat, IntXi, IntWeight, NintGLPt, &
          (lmesh%xmin - xmin_gl)/Lx-0.5_RP, delx/Lx  )
      end do
    end do

    ! global sum
    call MPI_AllReduce( s_coef_lc, s_coef, vec_size * (km+1) * 2, &
      MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr   )

    !-
    !$omp parallel do private(kk)
    do v=1, var_num
    do k=ks, ke
      kk = abs(k)
      spectral_coef(k,1,v) = s_coef(v,1,kk)
      spectral_coef(k,2,v) = s_coef(v,2,kk)
    end do
    end do

    return
  end subroutine spectral_transform1D_L2projection
!OCL SERIAL
  subroutine spectral_transform1D_L2projection_lc( spectral_coef, &
    q, ks, ke, Np, Ne, vec_size, IntIntrpMat, IntXi, IntWeight, NintGLPt, &
    xmin_lc, delx  )
    implicit none
    integer, intent(in) :: ks, ke
    integer, intent(in) :: Np
    integer, intent(in) :: Ne
    integer, intent(in) :: vec_size
    real(RP), intent(in) :: q(Np,vec_size,Ne)
    real(RP), intent(inout) :: spectral_coef(vec_size,2,ks:ke)
    integer, intent(in) :: NintGLpt
    real(RP), intent(in) :: IntIntrpMat(NintGLpt,Np)
    real(RP), intent(in) :: IntXi(NintGLpt)
    real(RP), intent(in) :: IntWeight(NintGLpt)
    real(RP), intent(in) :: xmin_lc
    real(RP), intent(in) :: delx

    integer :: k
    integer :: kel
    integer :: v

    real(RP) :: q_intrp(NintGLpt,vec_size,Ne)
    real(RP) :: int_X(NintGLpt,Ne)

    real(RP) :: phi(NintGLpt)
    real(RP) :: cos_kx(NintGLpt)
    real(RP) :: sin_kx(NintGLpt)

    real(RP) :: int_w(NintGLpt)
    !-------------------------------

    int_w(:) = 0.5_RP * IntWeight(:) * delx

    !$omp parallel private(kel,k,v,cos_kx,sin_kx,phi)
    !$omp do
    do kel=1, Ne
      int_X(:,kel) = 2.0_RP * PI * ( xmin_lc + delx * ( ( kel - 1.0_RP ) + 0.5_RP * ( 1.0_RP + intXi(:) ) ) )
      q_intrp(:,:,kel) = matmul(IntIntrpMat, q(:,:,kel))
    end do
    !$omp do
    do k=ks, ke
      do kel=1, Ne
        phi(:) = mod( k * int_X(:,kel), 2.0_RP * PI )
        cos_kx(:) = cos( phi )
        sin_kx(:) = sin( phi )

        do v=1, vec_size
          spectral_coef(v,1,k) = spectral_coef(v,1,k) + sum(int_w(:) * cos_kx(:) * q_intrp(:,v,kel))
          spectral_coef(v,2,k) = spectral_coef(v,2,k) - sum(int_w(:) * sin_kx(:) * q_intrp(:,v,kel))
        end do
      end do
    end do
    !$omp end parallel
    return
  end subroutine spectral_transform1D_L2projection_lc

!--- Private subroutines (2D)

!OCL SERIAL
  subroutine spectral_transform2D_DFT( spectral_coef, &
    q_list, var_num, mesh_num_x, mesh_num_y, ks, ke, ls, le, Np1D, NeX, NeY, NeGY, NprcY, NsamplePerElem1D, &
    fft_x, fft_y, FFTIntrpMat )

    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD    
    implicit none
    integer, intent(in) :: ks, ke
    integer, intent(in) :: ls, le
    integer, intent(in) :: Np1D
    integer, intent(in) :: NeX
    integer, intent(in) :: NeY, NeGY, NprcY
    integer, intent(in) :: NsamplePerElem1D
    integer, intent(in) :: var_num
    integer, intent(in) :: mesh_num_x
    integer, intent(in) :: mesh_num_y
    real(RP), intent(out) :: spectral_coef(ks:ke,ls:le,2,var_num)
    class(MeshField2DList), target :: q_list(var_num,mesh_num_x,mesh_num_y)
    class(FastFourierTransform1D), intent(in) :: fft_x
    class(FastFourierTransform1D), intent(in) :: fft_y
    real(RP), intent(in) :: FFTIntrpMat(NsamplePerElem1D,Np1D)

    real(RP) :: FFTIntrpMat_tr(Np1D,NsamplePerElem1D)

    real(RP) :: g_q(NsamplePerElem1D*NeX*mesh_num_x,NsamplePerElem1D*NeY*mesh_num_y*var_num)
    complex(RP) :: s_qx(NsamplePerElem1D*NeX*mesh_num_x,size(g_q,2))
    complex(RP) :: g_qy(NsamplePerElem1D*NeY*mesh_num_y*NprcY,var_num,NsamplePerElem1D*NeX*mesh_num_x/NprcY)
    complex(RP) :: s_q_lc(NsamplePerElem1D*NeGY,var_num,NsamplePerElem1D*NeX*mesh_num_x/NprcY)
    complex(RP) :: s_q(NsamplePerElem1D*NeGY,var_num,NsamplePerElem1D*NeX*mesh_num_x)

    integer :: v
    integer :: xdim
    integer :: prcy
    integer :: j

    integer :: Nall_x, Nall_y
    integer :: kk, ll, kk_os

    integer :: sendcount, recvcount
    complex(RP) :: sendbuf(NsamplePerElem1D*NeY*mesh_num_y*var_num,NsamplePerElem1D*NeX*mesh_num_x)
    complex(RP) :: recvbuf(NsamplePerElem1D*NeY*mesh_num_y,var_num,NsamplePerElem1D*NeX*mesh_num_x/NprcY,NprcY)
    integer :: nx, local_nx
    integer :: ny, local_ny
    integer :: ierr
    !-----------------------------------------------

    !-
    Nall_x = NsamplePerElem1D * NeX * mesh_num_x
    Nall_y = NsamplePerElem1D * NeGY

    !-
    FFTIntrpMat_tr(:,:) = transpose(FFTIntrpMat)
    call spectral_transform2D_DFT_sampling( g_q, &
      q_list, var_num, mesh_num_x, mesh_num_y, Np1D, NeX, NeY, NsamplePerElem1D, &
      FFTIntrpMat_tr )
   
    !-
    !$omp parallel do
    do v=1, size(g_q,2)
      call fft_x%Forward( g_q(:,v), s_qx(:,v) )
    end do

    !-
    local_nx = Nall_x
    nx = local_nx
    
    local_ny = NsamplePerElem1D * NeY * mesh_num_y * var_num
    ny = NsamplePerElem1D * NeGY * var_num
    
    sendbuf(:,:) = transpose(s_qx)
    sendcount = local_ny * (nx / NprcY)
    recvcount = sendcount
    call MPI_AlltoAll( sendbuf, sendcount, MPI_DOUBLE_COMPLEX, &
                       recvbuf, recvcount, MPI_DOUBLE_COMPLEX, PRC_LOCAL_COMM_WORLD, ierr )

    local_ny = NsamplePerElem1D * NeY * mesh_num_y                     
    do v=1, var_num
    do prcy=1, NprcY
    do xdim=1, NsamplePerElem1D*NeX*mesh_num_x/NprcY
    do j=1, local_ny
      g_qy(j+(prcy-1)*local_ny,v,xdim) = recvbuf(j,v,xdim,prcy)
    end do
    end do
    end do
    end do

    !-
    !$omp parallel do collapse(2)
    do v=1, var_num
    do xdim=1, nx/NprcY
      call fft_y%Forward( g_qy(:,v,xdim), s_q_lc(:,v,xdim) )
    end do
    end do

    sendcount = size(s_q_lc)
    recvcount = size(s_q_lc)
    call MPI_Allgather( s_q_lc, sendcount, MPI_DOUBLE_COMPLEX, &
                           s_q, recvcount, MPI_DOUBLE_COMPLEX, PRC_LOCAL_COMM_WORLD, ierr )

    !-
    !$omp parallel do private(v,kk,ll)
    do v=1, var_num
      do ll=1, Nall_y/2+1
        do kk=1, Nall_x/2+1
          spectral_coef(kk-1,ll-1,1,v) = real (s_q(ll,v,kk))
          spectral_coef(kk-1,ll-1,2,v) = aimag(s_q(ll,v,kk))
        end do
        do kk=Nall_x/2+1, Nall_x
          spectral_coef(kk-1-Nall_x,ll-1,1,v) = real (s_q(ll,v,kk))
          spectral_coef(kk-1-Nall_x,ll-1,2,v) = aimag(s_q(ll,v,kk))
        end do
      end do
      do ll=Nall_y/2+1, Nall_y
        do kk=1, Nall_x/2+1
          spectral_coef(kk-1,ll-1-Nall_y,1,v) = real (s_q(ll,v,kk))
          spectral_coef(kk-1,ll-1-Nall_y,2,v) = aimag(s_q(ll,v,kk))
        end do
        do kk=Nall_x/2+1, Nall_x
          spectral_coef(kk-1-Nall_x,ll-1-Nall_y,1,v) = real (s_q(ll,v,kk))
          spectral_coef(kk-1-Nall_x,ll-1-Nall_y,2,v) = aimag(s_q(ll,v,kk))
        end do
      end do
    end do
    return
  end subroutine spectral_transform2D_DFT

!OCL SERIAL
  subroutine spectral_transform2D_DFT_sampling( g_q, &
    q_list, var_num, mesh_num_x, mesh_num_y, Np1D, NeX, NeY, NsamplePerElem1D, &
    FFTIntrpMat_tr )

    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD    
    implicit none
    integer, intent(in) :: Np1D
    integer, intent(in) :: NeX, NeY
    integer, intent(in) :: NsamplePerElem1D
    integer, intent(in) :: var_num
    integer, intent(in) :: mesh_num_x
    integer, intent(in) :: mesh_num_y
    real(RP), intent(out) :: g_q(NsamplePerElem1D,NeX,mesh_num_x,NsamplePerElem1D,NeY,mesh_num_y,var_num)
    class(MeshField2DList), intent(in) :: q_list(var_num,mesh_num_x,mesh_num_y)
    real(RP), intent(in) :: FFTIntrpMat_tr(Np1D,NsamplePerElem1D)

    real(RP) :: q_tmp(Np1D,Np1D)
    real(RP) :: q_intrp1(NsamplePerElem1D,Np1D)
    real(RP) :: q_intrp2(NsamplePerElem1D,NsamplePerElem1D)

    integer :: mx, my, v
    integer :: kel_x, kel_y, kel
    integer :: p, p1, pp1, p2, pp2
    !-----------------------------------------------

    !$omp parallel do private(mx,my,v, kel_x,kel_y,kel, p,p1,pp1,p2,pp2, q_tmp,q_intrp1,q_intrp2) collapse(5)
    do v=1, var_num
    do my=1, mesh_num_y
    do mx=1, mesh_num_x
    do kel_y=1, NeY
    do kel_x=1, NeX
      kel = kel_x + (kel_y-1)*NeX
      do p2=1, Np1D
      do p1=1, Np1D
        p = p1 + (p2-1)*Np1D
        q_tmp(p1,p2) = q_list(v,mx,my)%ptr%local(1)%val(p,kel)
      end do
      end do

      q_intrp1(:,:) = 0.0_RP
      do p2=1, Np1D
      do p1=1, NsamplePerElem1D
      do pp1=1, Np1D
        q_intrp1(p1,p2) = q_intrp1(p1,p2) + FFTIntrpMat_tr(pp1,p1) * q_tmp(pp1,p2)
      end do
      end do
      end do      
      q_intrp2(:,:) = 0.0_RP
      do p2=1, NsamplePerElem1D
      do pp2=1, Np1D
      do p1=1, NsamplePerElem1D
        q_intrp2(p1,p2) = q_intrp2(p1,p2) + FFTIntrpMat_tr(pp2,p2) * q_intrp1(p1,pp2)
      end do
      end do
      end do      

      do p2=1, NsamplePerElem1D
        g_q(:,kel_x,mx,p2,kel_y,my,v) = q_intrp2(:,p2)
      end do
    end do
    end do
    end do
    end do
    end do
    return
  end subroutine spectral_transform2D_DFT_sampling

!OCL SERIAL
  subroutine spectral_transform2D_L2projection( spectral_coef, &
    q_list, var_num, mesh_num, ks, ke, ls, le, Np1D, NeX, NeY, &
    IntIntrpMat1D, IntXi1D, IntWeight1D, NintGLPt1D,           &
    xmin_gl, delx, Lx, ymin_gl, dely, Ly )

    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD
    implicit none
    integer, intent(in) :: ks, ke
    integer, intent(in) :: ls, le
    integer, intent(in) :: Np1D
    integer, intent(in) :: NeX, Ney
    integer, intent(in) :: var_num
    integer, intent(in) :: mesh_num
    type(MeshField2DList), target :: q_list(var_num,mesh_num)
    real(RP), intent(out) :: spectral_coef(ks:ke,ls:le,2,var_num)
    integer, intent(in) :: NintGLpt1D
    real(RP), intent(in) :: IntIntrpMat1D(NintGLpt1D,Np1D)
    real(RP), intent(in) :: IntXi1D(NintGLpt1D)
    real(RP), intent(in) :: IntWeight1D(NintGLpt1D)
    real(RP), intent(in) :: xmin_gl
    real(RP), intent(in) :: delx
    real(RP), intent(in) :: Lx
    real(RP), intent(in) :: ymin_gl
    real(RP), intent(in) :: dely
    real(RP), intent(in) :: Ly

    real(RP), allocatable  :: s_coef_lc(:,:,:,:)
    real(RP), allocatable  :: s_coef(:,:,:,:)
    real(RP), allocatable :: q_tmp(:,:,:)

    real(RP) :: IntIntrpMat1D_tr(Np1D,NintGLpt1D)

    integer :: vec_size
    integer :: km

    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh2D), pointer :: lmesh
    integer :: ldom
    integer :: meshID

    integer :: kel, v
    integer :: k, kk, l
    integer :: ierr
    !-----------------------------------------------

    vec_size = var_num
    km = ke

    allocate( s_coef_lc(vec_size,2,0:km,ls:le) )
    allocate( s_coef(vec_size,2,0:km,ls:le) )
    
    !-
    mesh2D => q_list(1,1)%ptr%mesh
    lmesh => mesh2D%lcmesh_list(1)

    IntIntrpMat1D_tr(:,:) = transpose(IntIntrpMat1D)
    !-

    allocate( q_tmp(Np1D**2,vec_size,NeX*NeY) )

    s_coef_lc(:,:,:,:) = 0.0_RP
    
    do meshID=1, mesh_num
      mesh2D => q_list(1,meshID)%ptr%mesh
      do ldom=1, mesh2D%LOCAL_MESH_NUM
        lmesh => mesh2D%lcmesh_list(ldom)
        
        !$omp parallel do collapse(2)
        do kel=lmesh%NeS, lmesh%NeE
        do v=1, vec_size
          q_tmp(:,v,kel) = q_list(v,meshID)%ptr%local(ldom)%val(:,kel)
        end do
        end do      
        call spectral_transform2D_L2projection_lc( s_coef_lc, &
          q_tmp, 0, km, ls, le, Np1D, NeX, NeY, vec_size,     &
          IntIntrpMat1D_tr, IntXi1D, IntWeight1D, NintGLpt1D, &
          (lmesh%xmin - xmin_gl)/Lx-0.5_RP, delx/Lx,          &
          (lmesh%ymin - ymin_gl)/Ly-0.5_RP, dely/Ly           )
      end do
    end do

    ! global sum
    call MPI_AllReduce( s_coef_lc, s_coef, vec_size * (km+1)*(le-ls+1) * 2, &
      MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr   )

    !-
    !$omp parallel do collapse(2) private(kk)
    do v=1, var_num
    do l=ls, le
      do k=ks, ke
        kk = abs(k)
        if ( k >= 0 ) then
          spectral_coef(k,l,1,v) = s_coef(v,1,kk,l)
          spectral_coef(k,l,2,v) = s_coef(v,2,kk,l)
        else
          spectral_coef(k,l,1,v) = s_coef(v,1,kk,-l)
          spectral_coef(k,l,2,v) = - s_coef(v,2,kk,-l)
        end if
      end do
    end do
    end do

    return
  end subroutine spectral_transform2D_L2projection
!OCL SERIAL
  subroutine spectral_transform2D_L2projection_lc( spectral_coef, &
    q, ks, ke, ls, le, Np1D, NeX, NeY, vec_size, IntIntrpMat1D_tr, IntXi1D, IntWeight1D, NintGLpt1D, &
    xmin_lc, delx, ymin_lc, dely  )
    implicit none
    integer, intent(in) :: ks, ke
    integer, intent(in) :: ls, le
    integer, intent(in) :: Np1D
    integer, intent(in) :: NeX
    integer, intent(in) :: NeY
    integer, intent(in) :: vec_size
    real(RP), intent(in) :: q(Np1D,Np1D,vec_size,NeX,NeY)
    real(RP), intent(inout) :: spectral_coef(vec_size,2,ks:ke,ls:le)
    integer, intent(in) :: NintGLpt1D
    real(RP), intent(in) :: IntIntrpMat1D_tr(Np1D,NintGLpt1D)
    real(RP), intent(in) :: IntXi1D(NintGLpt1D)
    real(RP), intent(in) :: IntWeight1D(NintGLpt1D)
    real(RP), intent(in) :: xmin_lc
    real(RP), intent(in) :: delx
    real(RP), intent(in) :: ymin_lc
    real(RP), intent(in) :: dely

    integer :: k, l
    integer :: kel_x, kel_y
    integer :: p, px, py
    integer :: v

    real(RP) :: q_intrp_tmp(NintGLpt1D,Np1D,vec_size)
    real(RP) :: q_intrp_tmp2(NintGLpt1D**2,vec_size)
    real(RP) :: q_intrp(NintGLpt1D**2,vec_size,NeX,NeY)
    real(RP) :: int_X(NintGLpt1D,NeX)
    real(RP) :: int_Y(NintGLpt1D,NeY)

    real(RP) :: phi(NintGLpt1D**2)
    real(RP) :: cos_phi(NintGLpt1D**2)
    real(RP) :: sin_phi(NintGLpt1D**2)

    real(RP) :: int_w(NintGLpt1D**2)
    !-------------------------------

    !$omp parallel private(kel_x,kel_y,k,l,v,p,px,py,cos_phi,sin_phi,phi,q_intrp_tmp,q_intrp_tmp2)
    !$omp do
    do py=1, NintGLpt1D
    do px=1, NintGLpt1D
      p = px + (py-1)*NintGLpt1D
      int_w(p) = 0.25_RP * delx * dely * IntWeight1D(px) * IntWeight1D(py)
    end do
    end do
    !$omp do
    do kel_x=1, NeX
      int_X(:,kel_x) = 2.0_RP * PI * ( xmin_lc + delx * ( ( kel_x - 1.0_RP ) + 0.5_RP * ( 1.0_RP + IntXi1D(:) ) ) )
    end do
    !$omp do
    do kel_y=1, NeY
      int_Y(:,kel_y) = 2.0_RP * PI * ( ymin_lc + dely * ( ( kel_y - 1.0_RP ) + 0.5_RP * ( 1.0_RP + IntXi1D(:) ) ) )
    end do
    !$omp do collapse(2)
    do kel_y=1, NeY
    do kel_x=1, NeX
      q_intrp_tmp(:,:,:) = 0.0_RP
      do v=1, vec_size
      do py=1, Np1D
      do px=1, NintGLpt1D
        q_intrp_tmp(px,py,v) = q_intrp_tmp(px,py,v) &
          + sum(IntIntrpMat1D_tr(:,px) * q(:,py,v,kel_x,kel_y))
      end do
      end do
      end do

      q_intrp_tmp2(:,:) = 0.0_RP
      do py=1, NintGLpt1D
      do v=1, vec_size
      do px=1, NintGLpt1D
        p = px + (py-1)*NintGLpt1D
        q_intrp_tmp2(p,v) = q_intrp_tmp2(p,v) &
          + sum(IntIntrpMat1D_tr(:,py) * q_intrp_tmp(px,:,v))
      end do
      end do
      end do
      q_intrp(:,:,kel_x,kel_y) = q_intrp_tmp2(:,:)
    end do
    end do
    !$omp do collapse(2)
    do l=ls, le
    do k=ks, ke
      do kel_y=1, NeY
      do kel_x=1, NeX
        do py=1, NintGLpt1D
        do px=1, NintGLpt1D
          p = px + (py-1)*NintGLpt1D
          phi(p) = mod( k * int_X(px,kel_x) + l * int_Y(py,kel_y), 2.0_RP * PI )
        end do
        end do
        cos_phi(:) = cos( phi(:) )
        sin_phi(:) = sin( phi(:) )

        do v=1, vec_size
          spectral_coef(v,1,k,l) = spectral_coef(v,1,k,l) + sum(int_w(:) * cos_phi(:) * q_intrp(:,v,kel_x,kel_y))
          spectral_coef(v,2,k,l) = spectral_coef(v,2,k,l) - sum(int_w(:) * sin_phi(:) * q_intrp(:,v,kel_x,kel_y))
        end do
      end do
      end do
    end do
    end do
    !$omp end parallel
    return
  end subroutine spectral_transform2D_L2projection_lc

end module scale_meshfield_spectral_transform
