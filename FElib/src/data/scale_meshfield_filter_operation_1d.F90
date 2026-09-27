!-------------------------------------------------------------------------------
!> module FElib / Data / Filter operation 1D
!!
!! @par Description
!!           This module provides classes to apply filter operations to MeshField data
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_meshfield_filter_operation_1d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_element_base, only: ElementBase1D
  use scale_element_modalfilter, only: ModalFilter
  use scale_mesh_base1d, only: MeshBase1D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_1d, only: LocalMesh1D

  use scale_meshfield_base, only: &
    MeshField1D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_meshfieldcomm_1d, only: &
    MeshFieldComm1D
  use scale_meshfieldcomm_base, only: &
    MeshFieldContainer

  use scale_meshfield_filter_operation_base, only: &
    MeshFieldFilterOperationBase, &
    MeshFieldFilterOperationBase_Init, MeshFieldFilterOperationBase_Final, &
    FILTER_OPTRTYPE_CONVFILTER,                                            &
    FILTER_OPTRTYPE_RECONSTRUCT, FILTER_OPTRTYPE_RECONSTRUCT2,             &
    FILTER_OPTRTYPE_RECONSTRUCT2_GL, &
    FILTER_OPTRTYPE_INTERFACE_CORRECTION, &
    apply_filter1d_x => MeshFieldFilterOperationBase_apply_filter1d_x,       &
    apply_reconst1d_x => MeshFieldFilterOperationBase_apply_reconst1d_x,     &
    apply_reconst1d_x_2 => MeshFieldFilterOperationBase_apply_reconst1d_x_2, &
    apply_interface_correction1d_x => MeshFieldFilterOperationBase_apply_interface_correction1d_x
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Derived type to represent filter operation for 1D mesh field
  type, public, extends(MeshFieldFilterOperationBase) :: MeshFieldFilterOperation1D
    type(MeshFieldComm1D) :: vars_comm

    integer :: Nnode_GL
  contains
    procedure :: Init => MeshFieldFilterOperation1D_Init
    procedure :: Final => MeshFieldFilterOperation1D_Final
    procedure :: Apply => MeshFieldFilterOperation1D_apply_filter
    procedure :: Apply_GL_lc => MeshFieldFilterOperation1D_apply_filter_GL_lc
  end type MeshFieldFilterOperation1D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  private :: extract_tmp1D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  !> Initialize an object to represent filter operation for 1D mesh field
!OCL SERIAL
  subroutine MeshFieldFilterOperation1D_Init( this, &
    FilterOptrType,                  &
    FilterShape, FilterWidthFac,     &
    Nnode_1D_reconst,                &
    sfield_num,                      &
    mesh1D, Nnode_1D_GL, IF_r )
    implicit none
    class(MeshFieldFilterOperation1D), intent(inout) :: this
    class(MeshBase1D), intent(in) :: mesh1D
    character(*), intent(in) :: FilterOptrType
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: FilterWidthFac
    integer, intent(in) :: Nnode_1D_reconst
    integer, intent(in) :: sfield_num
    integer, intent(in), optional :: Nnode_1D_GL
    integer, intent(in), optional :: IF_r
    !------------------------------------------

    call MeshFieldFilterOperationBase_Init( this, mesh1D%refElem1D%Np, &
      FilterOptrType, FilterShape, FilterWidthFac, &
      Nnode_1D_reconst, Nnode_1D_GL, IF_r )
    
    if (present(Nnode_1D_GL)) then
      this%Nnode_GL = Nnode_1D_GL
    else
      this%Nnode_GL = mesh1D%refElem1D%Np
    end if

    select type(mesh1D)
    class is (MeshBase1D)
      call this%vars_comm%Init( sfield_num, 0, mesh1D, &
        haloSize_1D=this%hHaloSize )
    class default
      LOG_INFO('MeshFieldFilterOperation1D_Init',*) "Unsupported mesh type is specified. Check!"
      call PRC_abort
    end select
    return
  end subroutine MeshFieldFilterOperation1D_Init

  !> Finalize an object to represent filter operation for 1D mesh field
!OCL SERIAL
  subroutine MeshFieldFilterOperation1D_Final( this )
    implicit none
    class(MeshFieldFilterOperation1D), intent(inout) :: this
    !------------------------------------------
    call MeshFieldFilterOperationBase_Final( this )
    call this%vars_comm%Final()
    return
  end subroutine MeshFieldFilterOperation1D_Final

  !> Apply filter operation to 1D mesh field data
!OCL serial
  subroutine MeshFieldFilterOperation1D_apply_filter( this, q_list, &
    mesh1D )
    implicit none
    class(MeshFieldFilterOperation1D), intent(inout) :: this
    type(MeshField1D), intent(inout), target :: q_list(this%vars_comm%field_num_tot)
    class(MeshBase1D), intent(in), target :: mesh1D

    class(LocalMesh1D), pointer :: lmesh1D
    class(ElementBase1D), pointer :: elem1D
    real(RP), allocatable :: tmp1D(:,:,:)

    integer :: iv
    integer :: n
    type(MeshFieldContainer) :: comm_vars_list(this%vars_comm%field_num_tot)
    !-----------------------------------------------------

    do iv=1, this%vars_comm%field_num_tot
      comm_vars_list(iv)%field1D => q_list(iv)
    end do
    call this%vars_comm%Put(comm_vars_list, 1)
    call this%vars_comm%Exchange()
    call this%vars_comm%Get(comm_vars_list, 1)

    do n=1, mesh1D%LOCAL_MESH_NUM
      lmesh1D => mesh1D%lcmesh_list(n)
      elem1D => lmesh1D%refElem1D
      allocate( tmp1D(-this%hHaloSize+1:elem1D%Np+this%hHaloSize,0:2,lmesh1D%Ne) )

      do iv=1, this%vars_comm%field_num_tot
        call extract_tmp1D( tmp1D, &
          q_list(iv)%local(n)%val, q_list(iv)%local(n)%val, lmesh1D, elem1D, lmesh1D%VMapP, this%hHaloSize )
        
        select case (this%operator_type)
        case (FILTER_OPTRTYPE_CONVFILTER)
          call apply_filter1d_x( q_list(iv)%local(n)%val, &
            tmp1D, this%FilterMat_h1D, elem1D%Np, 1, 1, lmesh1D%Ne, lmesh1D%NeA, elem1D%Np )
        case (FILTER_OPTRTYPE_RECONSTRUCT)
          call apply_reconst1d_x( q_list(iv)%local(n)%val, &
            tmp1D, this%Minv_Ml_tr, this%Minv_Mc_tr, this%Minv_Mr_tr, this%IntrpMat, elem1D%Np, 1, 1, lmesh1D%Ne, this%Nnode_h1D_reconst, lmesh1D )
        case (FILTER_OPTRTYPE_RECONSTRUCT2)
          call apply_reconst1d_x_2( q_list(iv)%local(n)%val, &
            tmp1D, this%Ml_tr, this%Mc_tr, this%Mr_tr, elem1D%Np, 1, 1, lmesh1D%Ne, this%Nnode_h1D_reconst, lmesh1D )
        case (FILTER_OPTRTYPE_INTERFACE_CORRECTION)
          call apply_interface_correction1d_x( q_list(iv)%local(n)%val, &
            tmp1D, this%IF_gL, this%IF_gR, elem1D%Np, 1, 1, lmesh1D%Ne, lmesh1D )
        end select

      end do
      deallocate( tmp1D )
    end do    
    return
  end subroutine MeshFieldFilterOperation1D_apply_filter


  !> Apply filter operation to 1D mesh field data
!OCL serial
  subroutine MeshFieldFilterOperation1D_apply_filter_GL_lc( this, &
    q_GL,   &
    q_list, &
    lmesh, elem1D  )
    implicit none
    class(MeshFieldFilterOperation1D), intent(inout) :: this
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(out) :: q_GL(this%Nnode_GL,lmesh%Ne,this%vars_comm%field_num_tot)
    type(MeshField1D), intent(inout), target :: q_list(this%vars_comm%field_num_tot)

    real(RP), allocatable :: tmp1D(:,:,:)

    integer :: iv
    integer :: n
    type(MeshFieldContainer) :: comm_vars_list(this%vars_comm%field_num_tot)
    !-----------------------------------------------------

    do iv=1, this%vars_comm%field_num_tot
      comm_vars_list(iv)%field1D => q_list(iv)
    end do
    call this%vars_comm%Put(comm_vars_list, 1)
    call this%vars_comm%Exchange()
    call this%vars_comm%Get(comm_vars_list, 1)

    allocate( tmp1D(-this%hHaloSize+1:elem1D%Np+this%hHaloSize,0:2,lmesh%Ne) )

    n = lmesh%lcdomID
    do iv=1, this%vars_comm%field_num_tot
      call extract_tmp1D( tmp1D, &
        q_list(iv)%local(n)%val, q_list(iv)%local(n)%val, lmesh, elem1D, lmesh%VMapP, this%hHaloSize )
      
      select case (this%operator_type)
      case (FILTER_OPTRTYPE_RECONSTRUCT2_GL)
        call apply_reconst1d_x_2_GL( q_GL(:,:,iv), &
          tmp1D, this%Ml_tr, this%Mc_tr, this%Mr_tr, elem1D%Np, this%Nnode_GL, lmesh%Ne, this%Nnode_h1D_reconst, lmesh )
      end select

    end do
    return
  end subroutine MeshFieldFilterOperation1D_apply_filter_GL_lc
  
!- Private -----------------------

!OCL SERIAL
  subroutine apply_reconst1d_x_2_GL( q, q0, Minv_Ml_tr, Minv_Mc_tr, Minv_Mr_tr, &
    Npx, NpxGL, Ne, Npx_reconst, lmesh )
    implicit none
    class(LocalMeshBase), intent(in) :: lmesh
    integer, intent(in) :: Npx, NpxGL, Ne
    integer, intent(in) :: Npx_reconst
    real(RP), intent(out) :: q(NpxGL,lmesh%Ne)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,0:2,Ne)
    real(RP), intent(in) :: Minv_Ml_tr(Npx,NpxGL)
    real(RP), intent(in) :: Minv_Mc_tr(Npx,NpxGL)
    real(RP), intent(in) :: Minv_Mr_tr(Npx,NpxGL)

    integer :: ke

    integer :: i,j

    real(RP) :: tmp(NpxGL)
    real(RP) :: s
    !-------------------------------------

    !$omp parallel do private(ke,i,j, tmp,s)
    do ke=1, Ne
      do i=1, NpxGL
        s = 0.0_RP
        do j=1, Npx
          s = s &
            + Minv_Ml_tr(j,i) * q0(j-Npx,1,ke) &
            + Minv_Mc_tr(j,i) * q0(    j,1,ke) &
            + Minv_Mr_tr(j,i) * q0(j+Npx,1,ke)
        end do
        tmp(i) = s
      end do
      q(:,ke) = tmp(:)
    end do
    return
  end subroutine apply_reconst1d_x_2_GL

!OCL SERIAL
  subroutine extract_tmp1D( tmp1D, q0, q0_, lmesh1D, elem1D, vmapP, hHaloSize )
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh1D
    class(ElementBase1D), intent(in) :: elem1D
    integer, intent(in) :: hHaloSize
    real(RP), intent(out) :: tmp1D(-hHaloSize+1:elem1D%Np+hHaloSize,0:2,lmesh1D%Ne)
    real(RP), intent(in) :: q0(elem1D%Np*lmesh1D%NeA) 
    real(RP), intent(in) :: q0_(elem1D%Np,lmesh1D%NeA) 
    integer, intent(in) :: vmapP(elem1D%NfpTot,lmesh1D%Ne) 

    integer :: ke, i, f
    integer :: fso
    integer :: ph
    integer :: iP
    real(RP) :: halo_h(hHaloSize,2)
    !------------------------------

    do f=1, 2
      fso = lmesh1D%Ne * elem1D%Np + hHaloSize * ( f-1 )
      do ph=1, hHaloSize
        halo_h(ph,f) = q0(fso+ph)
      end do
    end do
    
    !$omp parallel private(ke, i, ph, iP)
    !$omp do
    do ke=1, lmesh1D%Ne
      tmp1D(1:elem1D%Np,1,ke) = q0_(:,ke)

      ! Face 2
      iP = vmapP(2,ke)  
      if ( iP <= lmesh1D%Ne * elem1D%Np ) then
        do ph=1, hHaloSize
          tmp1D(ph+elem1D%Np,1,ke) = q0(iP + (ph-1))
        end do
      else
        do ph=1, hHaloSize
          tmp1D(ph+elem1D%Np,1,ke) = halo_h(ph,2)
        end do
      end if           

      ! Face 1
      iP = vmapP(1,ke)  
      if ( iP <= lmesh1D%Ne * elem1D%Np ) then
        do ph=1, hHaloSize
          tmp1D(-ph+1,1,ke) = q0(iP - (ph-1))
        end do
      else
        do ph=1, hHaloSize
          tmp1D(-ph+1,1,ke) = halo_h(ph,1)
        end do
      end if         
    end do

    !$omp end parallel
    return
  end subroutine extract_tmp1D
end module scale_meshfield_filter_operation_1d