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
  
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_base2d, only: MeshBase2D

  use scale_meshfield_base, only: &
    MeshField1D, MeshField2D, &
    MeshField1DList, MeshField2DList

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

    real(RP), allocatable :: IntIntrpMat(:,:)
    real(RP), allocatable :: IntWeight(:)
    real(RP), allocatable :: IntXi(:)
    integer :: NintGLpt
  end type MeshField_SpetralTransformBase

  type, public, extends(MeshField_SpetralTransformBase) :: MeshField_SpetralTransform1D
    integer :: kall
    integer :: ks, ke

    real(RP) :: xmin_gl, xmax_gl
    real(RP) :: delx

    integer :: NsamplePtPerElem

    real(RP), allocatable :: k(:)
    real(RP), allocatable :: spectral_coef(:,:,:)
  contains
    procedure :: Init => MeshField_SpetralTransform1D_Init
    procedure :: Final => MeshField_SpetralTransform1D_Final
    procedure :: Transform => MeshField_SpetralTransform1D_transform
  end type MeshField_SpetralTransform1D

  type, public, extends(MeshField_SpetralTransformBase) :: MeshField_SpetralTransform2D
    integer :: kall
    integer :: ks, ke
    integer :: lall
    integer :: ls, le

    real(RP) :: xmin_gl, xmax_gl
    real(RP) :: delx
    real(RP) :: ymin_gl, ymax_gl
    real(RP) :: dely

    integer :: NsamplePtPerElem1D

    real(RP), allocatable :: k(:)
    real(RP), allocatable :: l(:)
    real(RP), allocatable :: spectral_coef(:,:,:,:)
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
  subroutine MeshField_SpetralTransformBase_Init( this, eval_type, ndim, var_num, NintGLpt )
    implicit none
    class(MeshField_SpetralTransformBase), intent(inout) :: this
    integer, intent(in) :: eval_type
    integer, intent(in) :: ndim
    integer, intent(in) :: var_num
    integer, intent(in), optional :: NintGLpt
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
    return
  end subroutine MeshField_SpetralTransformBase_Init

  subroutine check_eval_type_id( eval_type ) 
    implicit none
    integer, intent(in) :: eval_type
    !-------------------------
    select case(eval_type)
!    case(ST_EVALTYPE_L2PROJECTION_1,ST_EVALTYPE_L2PROJECTION_2,ST_EVALTYPE_SAMPLE_UNIFORM_PTS)
    case(ST_EVALTYPE_L2PROJECTION_1)
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

    return
  end subroutine MeshField_SpetralTransformBase_Final  

!- 1D

!OCL SERIAL
  subroutine MeshField_SpetralTransform1D_Init( this, eval_type, ks, ke, mesh1D, var_num, &
    GLQuadOrd )
    use scale_element_line, only: &
      LineElement
    implicit none
    class(MeshField_SpetralTransform1D), intent(inout) :: this
    integer, intent(in) :: eval_type
    integer, intent(in) :: ks, ke
    class(MeshBase1D), intent(in), target :: mesh1D
    integer, intent(in) :: var_num
    integer, optional :: GLQuadOrd

    integer :: i
    integer :: Lx

    class(LocalMesh1D), pointer :: lmesh
    class(ElementBase1D), pointer :: elem

    type(LineElement) :: elem_dummy
    !--------------------------------------------

    lmesh => mesh1D%lcmesh_list(1)
    elem => lmesh%refElem1D

    call MeshField_SpetralTransformBase_Init( this, eval_type, 1, var_num, GLQuadOrd )

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

    if ( this%NintGLpt > 0 ) then
      allocate( this%IntIntrpMat(this%NintGLpt,elem%Np) )
      allocate( this%IntXi(this%NintGLpt), this%IntWeight(this%NintGLpt) )

      call elem_dummy%Init(elem%PolyOrder, .false.)
      this%IntIntrpMat(:,:) = elem_dummy%GenIntGaussLegendreIntrpMat( GLQuadOrd, this%IntWeight, this%IntXi )
      call elem_dummy%Final()
    end if

    return
  end subroutine MeshField_SpetralTransform1D_Init

!OCL SERIAL
  subroutine MeshField_SpetralTransform1D_Final( this )
    implicit none
    class(MeshField_SpetralTransform1D), intent(inout) :: this
    !--------------------------------------------

    call MeshField_SpetralTransformBase_Final( this )
    
    deallocate( this%k, this%spectral_coef )    
    return
  end subroutine MeshField_SpetralTransform1D_Final

!OCL SERIAL
  subroutine MeshField_SpetralTransform1D_transform( this, q_list, mesh_num )
    implicit none
    class(MeshField_SpetralTransform1D), intent(inout) :: this
    integer, intent(in) :: mesh_num
    class(MeshField1DList), target :: q_list(this%var_num,mesh_num)

    class(MeshBase1D), pointer :: mesh1D
    class(LocalMesh1D), pointer :: lmesh
    !--------------------------------------------

    mesh1D => q_list(1,1)%ptr%mesh
    lmesh => mesh1D%lcmesh_list(1)
    
    select case( this%eval_type_id )
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
    GLQuadOrd )
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

    integer :: i, j
    integer :: Lx, Ly

    class(LocalMesh2D), pointer :: lmesh
    class(ElementBase2D), pointer :: elem

    type(LineElement) :: elem_dummy
    class(MeshBase2D), pointer :: ptr_mesh2D
    !--------------------------------------------

    lmesh => mesh2D%lcmesh_list(1)
    elem => lmesh%refElem2D

    call MeshField_SpetralTransformBase_Init( this, eval_type, 1, var_num, GLQuadOrd )

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

    if ( this%NintGLpt > 0 ) then
      allocate( this%IntIntrpMat(this%NintGLpt,elem%Np) )
      allocate( this%IntXi(this%NintGLpt), this%IntWeight(this%NintGLpt) )

      call elem_dummy%Init(elem%PolyOrder, .false.)
      this%IntIntrpMat(:,:) = elem_dummy%GenIntGaussLegendreIntrpMat( GLQuadOrd, this%IntWeight, this%IntXi )
      call elem_dummy%Final()
    end if

    return
  end subroutine MeshField_SpetralTransform2D_Init

!OCL SERIAL
  subroutine MeshField_SpetralTransform2D_Final( this )
    implicit none
    class(MeshField_SpetralTransform2D), intent(inout) :: this
    !--------------------------------------------

    call MeshField_SpetralTransformBase_Final( this )
    
    deallocate( this%k, this%l, this%spectral_coef )
    return
  end subroutine MeshField_SpetralTransform2D_Final

!OCL SERIAL
  subroutine MeshField_SpetralTransform2D_transform( this, q_list, mesh_num )
    implicit none
    class(MeshField_SpetralTransform2D), intent(inout) :: this
    integer, intent(in) :: mesh_num
    class(MeshField2DList), target :: q_list(this%var_num,mesh_num)

    class(MeshBase2D), pointer :: mesh2D
    class(LocalMesh2D), pointer :: lmesh
    !--------------------------------------------

    mesh2D => q_list(1,1)%ptr%mesh
    lmesh => mesh2D%lcmesh_list(1)
    
    select case( this%eval_type_id )
    case( ST_EVALTYPE_L2PROJECTION_1)
      call spectral_transform2D_L2projection( this%spectral_coef,                       & ! (out)
        q_list, this%var_num, mesh_num,                                                 & ! (in)
        this%ks, this%ke,  this%ls, this%le, lmesh%refElem2D%Nfp, lmesh%NeX, lmesh%NeY, & ! (in)
        this%IntIntrpMat, this%IntXi, this%IntWeight, this%NintGLPt,                    & ! (in)
        this%xmin_gl, this%delx, this%xmax_gl - this%xmin_gl,                           & ! (in)
        this%ymin_gl, this%dely, this%ymax_gl - this%ymin_gl                            ) ! (in)
    end select
    
    return
  end subroutine MeshField_SpetralTransform2D_transform

!--- Private subroutines (1D)

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
    real(RP) :: Lx
    integer :: meshID

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
        !$omp parallel do collapse(2)
        do kel=lmesh%NeS, lmesh%NeE
        do v=1, vec_size
          q_tmp(:,v,kel) = q_list(v,meshID)%ptr%local(ldom)%val(:,kel)
        end do
        end do      
        call spectral_transform1D_L2projection_lc( s_coef_lc, &
          q_tmp, 0, km, Np, Ne, vec_size, IntIntrpMat, IntXi, IntWeight, NintGLPt, &
          (lmesh%xmin - xmin_gl)/Lx, delx  )
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
    class(MeshField2DList), target :: q_list(var_num,mesh_num)
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
          (lmesh%xmin - xmin_gl)/Lx, delx,                    &
          (lmesh%ymin - ymin_gl)/Ly, dely                     )
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
          + sum(IntIntrpMat1D_tr(:,px) * q_intrp_tmp(px,:,v))
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
