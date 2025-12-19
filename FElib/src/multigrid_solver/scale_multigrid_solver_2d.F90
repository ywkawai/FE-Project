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
module scale_multigrid_solver_2d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D

  use scale_meshfield_base, only: MeshField2D

  use scale_mesh_hierarchy_base, only: &
    MESH_HIERARCHY_pMG_FINEST_LEVEL, &
    MESH_HIERARCHY_hMG_FINEST_LEVEL, &
    MESH_HIERARCHY_TYPE_pMG,         &
    MESH_HIERARCHY_TYPE_hMG
  use scale_mesh_hierarchy_2d, only: &
    MeshHierarchy2D, &
    MeshHierarchyLevel2D, &
    MeshHierarchyLocalMGData2D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, public :: MultiGridSolver2D
    type(MeshHierarchy2D), pointer :: mesh_hierarchy_ptr
  contains
    procedure :: Init => MultiGridSolver2D_Init
    procedure :: Final => MultiGridSolver2D_Final
    procedure :: Operate_pMG_restriction => MultiGridSolver2D_Operate_pMG_restriction
    procedure :: Operate_pMG_correction => MultiGridSolver2D_Operate_pMG_correction
    procedure :: Operate_hMG_restriction => MultiGridSolver2D_Operate_hMG_restriction
    procedure :: Operate_hMG_correction => MultiGridSolver2D_Operate_hMG_correction
  end type MultiGridSolver2D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
contains
!OCL SERIAL
  subroutine MultiGridSolver2D_Init( this, &
    mesh_hierarachy  )
    implicit none
    class(MultiGridSolver2D), intent(inout) :: this
    class(MeshHierarchy2D), intent(in), target :: mesh_hierarachy
    !-------------------------------------------------------------

    this%mesh_hierarchy_ptr => mesh_hierarachy
    return
  end subroutine MultiGridSolver2D_Init

!OCL SERIAL
  subroutine MultiGridSolver2D_Final(this)
    implicit none
    class(MultiGridSolver2D), intent(inout) :: this
    !-------------------------------------------------------------
    return
  end subroutine MultiGridSolver2D_Final

!OCL SERIAL
  subroutine MultiGridSolver2D_Operate_pMG_restriction( this, res_c, &
    res, p_lev )
    implicit none
    class(MultiGridSolver2D), intent(in), target :: this
    class(MeshField2D), intent(inout) :: res_c
    class(MeshField2D), intent(in) :: res
    integer, intent(in) :: p_lev

    class(MeshHierarchy2D), pointer :: hierarchy
    integer :: ldomID
    class(MeshBase2D), pointer :: mesh2D
    !------------------------------------------------

    if ( this%mesh_hierarchy_ptr%NUM_pMG_LEVEL <  p_lev+1 ) then
      call PRC_abort()
    end if

    hierarchy => this%mesh_hierarchy_ptr
    mesh2D => hierarchy%p_mesh_list(p_lev)%ptr
    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      call MultiGridSolver2D_pMG_operation( res_c%local(ldomID)%val,     &
        res%local(ldomID)%val,                                           &
        hierarchy%elem2D_list(p_lev), hierarchy%elem2D_list(p_lev+1),    &
        mesh2D%lcmesh_list(ldomID), hierarchy%p_level(p_lev)%pMat1D_f2c, &
        .false.  )
    end do
    return
  end subroutine MultiGridSolver2D_Operate_pMG_restriction

!OCL SERIAL
  subroutine MultiGridSolver2D_Operate_pMG_correction( this, dq, &
    cor_c, p_lev )
    implicit none
    class(MultiGridSolver2D), intent(in), target :: this
    class(MeshField2D), intent(inout) :: dq
    class(MeshField2D), intent(in) :: cor_c
    integer, intent(in) :: p_lev

    class(MeshHierarchy2D), pointer :: hierarchy

    integer :: ldomID
    class(MeshBase2D), pointer :: mesh2D
    !------------------------------------------------

    hierarchy => this%mesh_hierarchy_ptr
    mesh2D => hierarchy%p_mesh_list(p_lev)%ptr

    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      call MultiGridSolver2D_pMG_operation( dq%local(ldomID)%val,        &
        cor_c%local(ldomID)%val,                                         &
        hierarchy%elem2D_list(p_lev+1), hierarchy%elem2D_list(p_lev),    &
        mesh2D%lcmesh_list(ldomID), hierarchy%p_level(p_lev)%pMat1D_c2f, &
        .true. )
    end do
    return
  end subroutine MultiGridSolver2D_Operate_pMG_correction

!OCL SERIAL
  subroutine MultiGridSolver2D_Operate_hMG_restriction( this, res_c, &
    res, h_lev )
    use scale_meshfield_base, only: MeshField2D
    implicit none
    class(MultiGridSolver2D), intent(in), target :: this
    class(MeshField2D), intent(inout) :: res_c
    class(MeshField2D), intent(in) :: res
    integer, intent(in) :: h_lev

    class(MeshHierarchy2D), pointer :: hierarchy

    integer :: ldomID
    class(MeshBase2D), pointer :: mesh2D
    class(MeshBase2D), pointer :: mesh2D_c
    !------------------------------------------------

    hierarchy => this%mesh_hierarchy_ptr
    mesh2D => hierarchy%h_mesh_list(h_lev)%ptr
    mesh2D_c => hierarchy%h_mesh_list(h_lev+1)%ptr

    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      call MultiGridSolver2D_hMG_restriction_core( res_c%local(ldomID)%val,          &
        res%local(ldomID)%val,                                                       &
        hierarchy%h_level(h_lev)%mg_local(ldomID), mesh2D%lcmesh_list(ldomID),       &
        hierarchy%elem2D_list(hierarchy%NUM_pMG_LEVEL), mesh2D_c%lcmesh_list(ldomID) )
    end do
    return
  end subroutine MultiGridSolver2D_Operate_hMG_restriction

!OCL SERIAL
  subroutine MultiGridSolver2D_Operate_hMG_correction( this, dq, &
    cor_c, h_lev )
    use scale_meshfield_base, only: MeshField2D
    implicit none
    class(MultiGridSolver2D), intent(in), target :: this
    class(MeshField2D), intent(inout) :: dq
    class(MeshField2D), intent(in) :: cor_c
    integer, intent(in) :: h_lev

    integer :: ldomID
    class(MeshBase2D), pointer :: mesh2D
    class(MeshHierarchy2D), pointer :: hierarchy
    !------------------------------------------------

    hierarchy => this%mesh_hierarchy_ptr
    mesh2D => hierarchy%h_mesh_list(h_lev)%ptr

    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      call MultiGridSolver2D_hMG_correction_core( dq%local(ldomID)%val,     &
        cor_c%local(ldomID)%val, hierarchy%h_level(h_lev)%mg_local(ldomID), &
        mesh2D%lcmesh_list(ldomID), mesh2D%refElem2D                        )
    end do
    return
  end subroutine MultiGridSolver2D_Operate_hMG_correction

!-- private --------------------------------------------------------------


!OCL SERIAL
  subroutine MultiGridSolver2D_hMG_restriction_core( res_c_lc,      &
      res_lc, mg_local, lmesh, elem, lmesh_c )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh_c
    real(RP), intent(out) :: res_c_lc(elem%Np,lmesh_c%NeA)
    real(RP), intent(in) :: res_lc(elem%Np,lmesh%NeA)
    class(MeshHierarchyLocalMGData2D), intent(in) :: mg_local

    integer :: k, ke, ke_c
    real(RP) :: Ic2fT_lc(4,4)
    real(RP) :: tmp_c(elem%Np)
    real(RP) :: tmp2(elem%Np)
    !-----------------------------------------

    !$omp parallel private(ke, ke_c, tmp2, tmp_c, Ic2fT_lc)
    !$omp do
    do ke_c=lmesh_c%NeS, lmesh_c%NeE
      
      tmp_c(:) = 0.0_RP
      do k=1, 4
        ke = mg_local%If2c_emap(k,ke_c)

        tmp2(:) = matmul( elem%M, res_lc(:,ke) )
        Ic2fT_lc(:,:) = transpose( mg_local%Ic2f(:,:,ke) )
        tmp_c(:) = tmp_c(:) + matmul( Ic2fT_lc, tmp2(:) )
      end do

      res_c_lc(:,ke_c) = matmul(lmesh_c%refElem2D%invM, tmp_c(:))  * 0.25_RP        
    end do
    !$omp end parallel
    return
  end subroutine MultiGridSolver2D_hMG_restriction_core

!OCL SERIAL
  subroutine MultiGridSolver2D_hMG_correction_core( dq_lc,  &
    dq_c_lc, mg_local, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: dq_lc(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dq_c_lc(elem%Np,lmesh%NeA)
    class(MeshHierarchyLocalMGData2D), intent(in) :: mg_local

    integer :: k, ke, ke_c
    !-----------------------------------------

    !$omp parallel do private(ke, ke_c)
    do ke=lmesh%NeS, lmesh%NeE
      ke_c = mg_local%Ic2f_emap(ke)
      dq_lc(:,ke) = dq_lc(:,ke) + &
            matmul( mg_local%Ic2f(:,:,ke), dq_c_lc(:,ke_c) )
    end do
    return
  end subroutine MultiGridSolver2D_hMG_correction_core

!OCL SERIAL
  subroutine MultiGridSolver2D_pMG_operation( q_o, &
    q_i, elem2D_i, elem2D_o, lcmesh, pMat1D, is_added )
    implicit none
    class(ElementBase2D), intent(in) :: elem2D_i
    class(ElementBase2D), intent(in) :: elem2D_o
    class(LocalMesh2D), intent(in) :: lcmesh
    real(RP), intent(inout) :: q_o(elem2D_o%Nfp,elem2D_o%Nfp,lcmesh%NeA)
    real(RP), intent(in) :: q_i(elem2D_i%Nfp,elem2D_i%Nfp,lcmesh%NeA)
    real(RP), intent(in) :: pMat1D(elem2D_o%Nfp,elem2D_i%Nfp)
    logical, intent(in) :: is_added

    integer :: ke

    integer :: px, py
    integer :: pxx, pyy
    real(RP) :: tmp1
    real(RP) :: tmp2(elem2D_o%Nfp,elem2D_i%Nfp)
    real(RP) :: tmp3(elem2D_o%Nfp)

    real(RP) :: mat_tr(elem2D_i%Nfp,elem2D_o%Nfp)
    !-------------------------------------------

    mat_tr(:,:) = transpose(pMat1D)

    !$omp parallel do private(tmp1, tmp2, tmp3)
    do ke=lcmesh%NeS, lcmesh%NeE
      do py=1, elem2D_i%Nfp
      do pxx=1, elem2D_o%Nfp
        tmp1 = 0.0_RP
        do px=1, elem2D_i%Nfp
          tmp1 = tmp1 + mat_tr(px,pxx) * q_i(px,py,ke)
        end do
        tmp2(pxx,py) = tmp1
      end do
      end do

      if ( is_added ) then
        do pyy=1, elem2D_o%Nfp
          tmp3(:) = 0.0_RP
          do py=1, elem2D_i%Nfp
            do px=1, elem2D_o%Nfp
              tmp3(px) = tmp3(px) + mat_tr(py,pyy) * tmp2(px,py)
            end do
          end do
          q_o(:,pyy,ke) = q_o(:,pyy,ke)  + tmp3(:)
        end do
      else
        do pyy=1, elem2D_o%Nfp
          tmp3(:) = 0.0_RP
          do py=1, elem2D_i%Nfp
            do px=1, elem2D_o%Nfp
              tmp3(px) = tmp3(px) + mat_tr(py,pyy) * tmp2(px,py)
            end do
          end do
          q_o(:,pyy,ke) = tmp3(:)
        end do
      end if
    end do
    return
  end subroutine MultiGridSolver2D_pMG_operation
end module scale_multigrid_solver_2d
