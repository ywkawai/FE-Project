!-------------------------------------------------------------------------------
!> module FElib / Mesh / 2D domain
!!
!! @par Description
!!      Manage mesh hierarchy of 2D domain for element-based methods
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

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D

  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, public :: MGFieldSet2D
    integer :: var_num
    type(MeshField2D) :: dq
    type(MeshField2D) :: f
    type(MeshField2D) :: res
    type(MeshField2D), allocatable :: aux_var(:)
    type(MeshField2D), allocatable :: aux_var_hvec(:)
    integer :: level_id

    type(MeshFieldCommRectDom2D) :: var_comm
    type(MeshFieldCommRectDom2D) :: aux_comm
    type(SparseMat) :: Dx
    type(SparseMat) :: Dy
  contains
    procedure :: Init => MGFieldSet2D_Init
    procedure :: Final => MGFieldSet2D_Final
  end type MGFieldSet2D

contains
!-
!OCL SERIAL
  subroutine MGFieldSet2D_Init(this, mesh2D, aux_var_num, aux_hvec_comp_num, level_id)
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    implicit none
    class(MGFieldSet2D), intent(inout) :: this
    class(MeshBase2D), intent(in) :: mesh2D
    integer, intent(in) :: aux_var_num
    integer, intent(in) :: aux_hvec_comp_num
    integer, intent(in) :: level_id

    integer :: iv
    !---------------------------------

    this%level_id = level_id
    call this%dq%Init( "dq", "1", mesh2D )    
    call this%f%Init( "f", "1", mesh2D )
    call this%res%Init( "res", "1", mesh2D )

    if ( aux_var_num > 0 ) then
      allocate( this%aux_var(aux_var_num) )
      do iv=1, aux_var_num
        call this%aux_var(iv)%Init( "aux_var", "1", mesh2D )
      end do
    end if

    if ( aux_hvec_comp_num > 0 ) then
      allocate( this%aux_var_hvec(aux_hvec_comp_num) )
      do iv=1, aux_hvec_comp_num
        call this%aux_var_hvec(iv)%Init( "aux_var_hvec", "1", mesh2D )
      end do
    end if
    
    select type(mesh2D)
    type is (MeshRectDom2D)
      call this%var_comm%Init( 1, 0, 0, mesh2D )
      if ( aux_var_num > 0 .or. aux_hvec_comp_num > 0 ) then
        call this%aux_comm%Init( aux_var_num, aux_hvec_comp_num, 0, mesh2D )
      end if
    end select

    call this%Dx%Init( mesh2D%refElem2D%Dx1, storage_format='ELL')
    call this%Dy%Init( mesh2D%refElem2D%Dx2, storage_format='ELL')
    return
  end subroutine MGFieldSet2D_Init

!OCL SERIAL
  subroutine MGFieldSet2D_Final(this)
    implicit none
    class(MGFieldSet2D), intent(inout) :: this

    integer :: iv
    !---------------------------------

    call this%dq%Final()    
    call this%f%Final()    
    call this%res%Final()
    call this%var_comm%Final()

    if ( allocated(this%aux_var) .or. allocated(this%aux_var_hvec) ) then
      call this%aux_comm%Final()
    end if

    if ( allocated(this%aux_var) ) then
      do iv=1, size(this%aux_var)
        call this%aux_var(iv)%Final()
      end do
      deallocate( this%aux_var)
    end if

    if ( allocated(this%aux_var_hvec) ) then
      do iv=1, size(this%aux_var_hvec)
        call this%aux_var_hvec(iv)%Final()
      end do
      deallocate( this%aux_var_hvec)
    end if

    return
  end subroutine MGFieldSet2D_Final

end module scale_multigrid_fieldset_base