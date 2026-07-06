!-------------------------------------------------------------------------------
!> module FElib / Multigrid / Field set base
!!
!! @par Description
!!      Manage field sets for multigrid method
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_multigrid_fieldset_base

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_sparsemat, only: SparseMat

  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement

  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D

  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_meshfieldcomm_base, only: MeshFieldCommBase
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D
  use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
  use scale_meshfieldcomm_cubedspheredom3d, only: MeshFieldCommCubedSphereDom3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Base type for multigrid field set
  type :: MGFieldBase
    integer :: var_num   !< Number of variables
    integer :: level_id  !< Level ID in multigrid hierarchy

    integer :: aux_var_num_tot !< Total number of auxiliary variables

    class(MeshFieldCommBase), pointer :: var_comm_ptr !< Pointer to communicator for main variables
    class(MeshFieldCommBase), pointer :: aux_comm_ptr !< Pointer to communicator for auxiliary variables
  end type MGFieldBase

  !> Derived type for 2D multigrid field set
  type, extends(MGFieldBase), public :: MGFieldSet2D
    type(MeshField2D) :: dq
    type(MeshField2D) :: f
    type(MeshField2D) :: res
    type(MeshField2D), allocatable :: aux_var(:)

    type(SparseMat) :: Dx
    type(SparseMat) :: Dy
    type(SparseMat) :: Lift
  contains
    procedure :: Init => MGFieldSet2D_Init
    procedure :: Final => MGFieldSet2D_Final
  end type MGFieldSet2D

  !> Derived type for 3D multigrid field set
  type, extends(MGFieldBase), public :: MGFieldSet3D
    type(MeshField3D) :: dq
    type(MeshField3D) :: f
    type(MeshField3D) :: res
    type(MeshField3D), allocatable :: aux_var(:)

    type(SparseMat) :: Dx
    type(SparseMat) :: Dy
    type(SparseMat) :: Dz
    type(SparseMat) :: Lift
  contains
    procedure :: Init => MGFieldSet3D_Init
    procedure :: Final => MGFieldSet3D_Final
  end type MGFieldSet3D

contains

!- 2D
  !> Initialize an object to manage multigrid field set in 2D
!OCL SERIAL
  subroutine MGFieldSet2D_Init(this, mesh2D, aux_scalar_num, aux_vec_num, level_id)
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    implicit none
    class(MGFieldSet2D), intent(inout), target :: this
    class(MeshBase2D), intent(in) :: mesh2D
    integer, intent(in) :: aux_scalar_num
    integer, intent(in) :: aux_vec_num
    integer, intent(in) :: level_id

    integer :: iv

    type(MeshFieldCommRectDom2D), pointer :: var_comm2D
    type(MeshFieldCommRectDom2D), pointer :: aux_comm2D

    type(MeshFieldCommCubedSphereDom2D), pointer :: var_comm2D_cs
    type(MeshFieldCommCubedSphereDom2D), pointer :: aux_comm2D_cs
    !---------------------------------

    this%level_id = level_id
    
    call this%dq%Init( "dq", "1", mesh2D )    
    call this%f%Init( "f", "1", mesh2D )
    call this%res%Init( "res", "1", mesh2D )

    this%aux_var_num_tot = aux_scalar_num + 2*aux_vec_num
    if ( this%aux_var_num_tot > 0 ) then
      allocate( this%aux_var(this%aux_var_num_tot) )
      do iv=1, this%aux_var_num_tot
        call this%aux_var(iv)%Init( "aux_var", "1", mesh2D )
      end do
    end if
    
    select type(mesh2D)
    type is (MeshRectDom2D)
      allocate(var_comm2D)
      call var_comm2D%Init( 1, 0, 0, mesh2D )
      this%var_comm_ptr => var_comm2D
    type is (MeshCubedSphereDom2D)
      allocate(var_comm2D_cs)
      call var_comm2D_cs%Init( 1, 0, 0, mesh2D )
      this%var_comm_ptr => var_comm2D_cs
    end select

    if ( this%aux_var_num_tot > 0 ) then
      select type(mesh2D)
      type is (MeshRectDom2D)
        allocate( aux_comm2D )
        call aux_comm2D%Init( aux_scalar_num, aux_vec_num, 0, mesh2D )
        this%aux_comm_ptr => aux_comm2D
      type is (MeshCubedSphereDom2D)
        allocate( aux_comm2D_cs )
        call aux_comm2D_cs%Init( aux_scalar_num, aux_vec_num, 0, mesh2D )
        this%aux_comm_ptr => aux_comm2D_cs
      end select
    end if
  
    call this%Dx%Init( mesh2D%refElem2D%Dx1, storage_format='ELL')
    call this%Dy%Init( mesh2D%refElem2D%Dx2, storage_format='ELL')
    call this%Lift%Init( mesh2D%refElem2D%Lift, storage_format='ELL')
    return
  end subroutine MGFieldSet2D_Init

  !> Finalize an object to manage multigrid field set in 2D
!OCL SERIAL
  subroutine MGFieldSet2D_Final(this)
    implicit none
    class(MGFieldSet2D), intent(inout) :: this

    integer :: iv
    class(MeshFieldCommBase), pointer :: comm_ptr
    !---------------------------------

    call this%dq%Final()    
    call this%f%Final()    
    call this%res%Final()

    select type(comm_ptr => this%var_comm_ptr)
    type is (MeshFieldCommRectDom2D)
      call comm_ptr%Final()
    type is (MeshFieldCommCubedSphereDom2D)
      call comm_ptr%Final()
    end select
    deallocate(this%var_comm_ptr)

    if ( this%aux_var_num_tot > 0 ) then
      
      select type(comm_ptr => this%aux_comm_ptr)
      type is (MeshFieldCommRectDom2D)
        call comm_ptr%Final()
      type is (MeshFieldCommCubedSphereDom2D)
        call comm_ptr%Final()
      end select
      deallocate(this%aux_comm_ptr)

      do iv=1, this%aux_var_num_tot
        call this%aux_var(iv)%Final()
      end do
    end if

    call this%Dx%Final()
    call this%Dy%Final()

    return
  end subroutine MGFieldSet2D_Final

!- 3D
  !> Initialize an object to manage multigrid field set in 3D
!OCL SERIAL
  subroutine MGFieldSet3D_Init(this, mesh3D, aux_scalar_num, aux_hvec_num, level_id)
    use scale_mesh_cubedom3d, only: MeshCubeDom3D
    implicit none
    class(MGFieldSet3D), intent(inout), target :: this
    class(MeshBase3D), intent(in) :: mesh3D
    integer, intent(in) :: aux_scalar_num
    integer, intent(in) :: aux_hvec_num
    integer, intent(in) :: level_id

    integer :: iv

    type(MeshFieldCommCubeDom3D), pointer :: var_comm3D
    type(MeshFieldCommCubeDom3D), pointer :: aux_comm3D

    type(MeshFieldCommCubedSphereDom3D), pointer :: var_comm3D_cs
    type(MeshFieldCommCubedSphereDom3D), pointer :: aux_comm3D_cs
    !---------------------------------

    this%level_id = level_id
    
    call this%dq%Init( "dq", "1", mesh3D )    
    call this%f%Init( "f", "1", mesh3D )
    call this%res%Init( "res", "1", mesh3D )

    this%aux_var_num_tot = aux_scalar_num + 2*aux_hvec_num
    if ( this%aux_var_num_tot > 0 ) then
      allocate( this%aux_var(this%aux_var_num_tot) )
      do iv=1, this%aux_var_num_tot
        call this%aux_var(iv)%Init( "aux_var", "1", mesh3D )
      end do
    end if
    
    select type(mesh3D)
    type is (MeshCubeDom3D)
      allocate(var_comm3D)
      call var_comm3D%Init( 1, 0, 0, mesh3D )
      this%var_comm_ptr => var_comm3D
    type is (MeshCubedSphereDom3D)
      allocate(var_comm3D_cs)
      call var_comm3D_cs%Init( 1, 0, 0, mesh3D )
      this%var_comm_ptr => var_comm3D_cs
    end select

    !-
    if ( this%aux_var_num_tot > 0 ) then
      select type(mesh3D)
      type is (MeshCubeDom3D)
        allocate( aux_comm3D )
        call aux_comm3D%Init( aux_scalar_num, aux_hvec_num, 0, mesh3D )
        this%aux_comm_ptr => aux_comm3D
      type is (MeshCubedSphereDom3D)
        allocate( aux_comm3D_cs )
        call aux_comm3D_cs%Init( aux_scalar_num, aux_hvec_num, 0, mesh3D )
        this%aux_comm_ptr => aux_comm3D_cs
      end select
    end if

    !-
    call this%Dx%Init( mesh3D%refElem3D%Dx1, storage_format='ELL')
    call this%Dy%Init( mesh3D%refElem3D%Dx2, storage_format='ELL')
    call this%Dz%Init( mesh3D%refElem3D%Dx3, storage_format='ELL')
    call this%Lift%Init( mesh3D%refElem3D%Lift, storage_format='ELL')
    return
  end subroutine MGFieldSet3D_Init

  !> Finalize an object to manage multigrid field set in 3D
!OCL SERIAL
  subroutine MGFieldSet3D_Final(this)
    implicit none
    class(MGFieldSet3D), intent(inout) :: this

    integer :: iv
    class(MeshFieldCommBase), pointer :: comm_ptr
    !---------------------------------

    call this%dq%Final()    
    call this%f%Final()    
    call this%res%Final()

    select type(comm_ptr => this%var_comm_ptr)
    type is (MeshFieldCommCubeDom3D)
      call comm_ptr%Final()
    type is (MeshFieldCommCubedSphereDom3D)
      call comm_ptr%Final()
    end select
    deallocate(this%var_comm_ptr)

    if ( this%aux_var_num_tot > 0 ) then
      
      select type(comm_ptr => this%aux_comm_ptr)
      type is (MeshFieldCommCubeDom3D)
        call comm_ptr%Final()
      type is (MeshFieldCommCubedSphereDom3D)
        call comm_ptr%Final()
      end select
      deallocate(this%aux_comm_ptr)

      do iv=1, this%aux_var_num_tot
        call this%aux_var(iv)%Final()
      end do
    end if

    call this%Dx%Final()
    call this%Dy%Final()
    call this%Dz%Final()

    return
  end subroutine MGFieldSet3D_Final

end module scale_multigrid_fieldset_base