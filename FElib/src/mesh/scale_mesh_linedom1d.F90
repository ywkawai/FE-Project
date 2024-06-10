!-------------------------------------------------------------------------------
!> module FElib / Mesh / 1D domain
!!
!! @par Description
!!      Mangage mesh data of 1D domain for element-based methods
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_mesh_linedom1d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_mesh_base1d, only: &
    MeshBase1D, MeshBase1D_Init, MeshBase1D_Final, &
    MeshBase1D_assignDomID, MeshBase1D_setGeometricInfo, MeshBase1D_setupLocalDom

  use scale_localmesh_1d, only: LocalMesh1D
  use scale_element_base, only: ElementBase1D
  use scale_element_line, only: LineElement

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(MeshBase1D), public :: MeshLineDom1D
  contains
    procedure :: Init => MeshLineDom1D_Init
    procedure :: Final => MeshLineDom1D_Final
    procedure :: Generate => MeshLineDom1D_generate
  end type MeshLineDom1D

  public :: MeshLineDom1D_Init, MeshLineDom1D_Final
  public :: MeshLineDom1D_generate
  
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
  subroutine MeshLineDom1D_Init(this,       & ! (inout)
    NeG,                                    & ! (in)
    dom_xmin, dom_xmax,                     & ! (in)
    refElem, NLocalMeshPerPrc,              & ! (in)
    nproc, myrank, FX )                       ! (in)
    
    implicit none

    class(MeshLineDom1D), intent(inout) :: this
    integer, intent(in) :: NeG
    real(RP), intent(in) :: dom_xmin
    real(RP), intent(in) :: dom_xmax   
    type(LineElement), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc
    integer, intent(in), optional :: nproc
    integer, intent(in), optional :: myrank
    real(RP), intent(in), optional :: FX(NeG+1)
    !-----------------------------------------------------------------------------

    call MeshBase1D_Init(this, &
      NeG,                                    &
      dom_xmin, dom_xmax,                     &
      refElem, NLocalMeshPerPrc,              &
      nproc, myrank, FX )

      this%dom_vol = (this%xmax_gl - this%xmin_gl)

  end subroutine MeshLineDom1D_Init

  subroutine MeshLineDom1D_Final( this ) ! (inout)
    implicit none
    class(MeshLineDom1D), intent(inout) :: this
    !-----------------------------------------------------------------------------
  
    call MeshBase1D_Final(this)

  end subroutine MeshLineDom1D_Final
  
  subroutine MeshLineDom1D_generate( this ) ! (inout)
    implicit none
    class(MeshLineDom1D), intent(inout), target :: this
            
    integer :: n
    integer :: p
    type(LocalMesh1D), pointer :: mesh

    integer :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)

    !integer :: TILE_NUM_PER_PANEL
    integer :: tileID
    
    !-----------------------------------------------------------------------------

   ! TILE_NUM_PER_PANEL = this%LOCAL_MESH_NUM_global / 1
    
    
    !--- Construct the connectivity of patches  (only master node)

    call MeshBase1D_assignDomID( this,    & ! (in)
        tileID_table, panelID_table,      & ! (out)
        pi_table )                          ! (out)

    !--- Setup local meshes managed by my process

    do n=1, this%LOCAL_MESH_NUM
      mesh => this%lcmesh_list(n)
      tileID = tileID_table(n, mesh%PRC_myrank+1)
      call MeshBase1D_setupLocalDom( mesh,              & ! (inout)
         tileID,  panelID_table(tileID),                & ! (in)
         pi_table(tileID), this%Nprc,                   & ! (in)
         this%xmin_gl, this%xmax_gl,                    & ! (in)
         this%NeG / this%Nprc, this%FX(:) )               ! (in)

      !---
      ! write(*,*) "** my_rank=", mesh%PRC_myrank
      ! write(*,*) " tileID:", mesh%tileID
      ! write(*,*) " pnlID:", mesh%panelID, "-- i (within a panel)=", pi_table(tileID)
      ! write(*,*) " local mesh:", n, "( total", this%LOCAL_MESH_NUM, ")"
      ! write(*,*) " panel_connect:", this%tilePanelID_globalMap(:,mesh%tileID)
      ! write(*,*) " tile_connect:", this%tileID_globalMap(:,mesh%tileID)
      ! write(*,*) " face_connect:", this%tileFaceID_globalMap(:,mesh%tileID)
      ! write(*,*) " domain size"
      ! write(*,*) "   NeX:", mesh%Ne
      ! write(*,*) "   [X]:",  mesh%xmin, mesh%xmax
    end do

    this%isGenerated = .true.

    return
  end subroutine MeshLineDom1D_generate
  
  !-------------------------------------------------
end module scale_mesh_linedom1d
