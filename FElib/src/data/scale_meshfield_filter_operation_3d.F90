!-------------------------------------------------------------------------------
!> module FElib / Data / Filter operation
!!
!! @par Description
!!           This module provides classes to apply filter operations to MeshField data
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_meshfield_filter_operation_3d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_modalfilter, only: ModalFilter
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D  

  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_mesh_rectdom2d, only: &
    MeshRectDom2D
  use scale_mesh_cubedom3d, only: &
    MeshCubeDom3D
  use scale_meshfieldcomm_cubedom3d, only: &
    MeshFieldCommCubeDom3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldContainer

  use scale_meshfield_filter_operation_base, only: &
    MeshFieldFilterOperationBase, &
    MeshFieldFilterOperationBase_Init, MeshFieldFilterOperationBase_Final, &
    FILTER_OPTRTYPE_CONVFILTER,                                            &
    FILTER_OPTRTYPE_RECONSTRUCT, FILTER_OPTRTYPE_RECONSTRUCT2,             &
    apply_filter1d_x => MeshFieldFilterOperationBase_apply_filter1d_x,       &
    apply_reconst1d_x => MeshFieldFilterOperationBase_apply_reconst1d_x,     &
    apply_reconst1d_x_2 => MeshFieldFilterOperationBase_apply_reconst1d_x_2, &
    apply_filter1d_y => MeshFieldFilterOperationBase_apply_filter1d_y,       &
    apply_reconst1d_y => MeshFieldFilterOperationBase_apply_reconst1d_y,     &
    apply_reconst1d_y_2 => MeshFieldFilterOperationBase_apply_reconst1d_y_2
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 

  !> Derived type to represent filter operation for 3D mesh field
  type, extends(MeshFieldFilterOperationBase), public :: MeshFieldFilterOperation3D
    type(MeshFieldCommCubeDom3D) :: vars_comm
    integer :: vHaloSize
  contains
    procedure :: Init => MeshFieldFilterOperation3D_Init
    procedure :: Final => MeshFieldFilterOperation3D_Final
    procedure :: Apply => MeshFieldFilterOperation3D_apply_filter
  end type MeshFieldFilterOperation3D


  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  private :: extract_tmp3D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains

  !> Initialize an object to represent filter operation for 3D mesh field
!OCL SERIAL
  subroutine MeshFieldFilterOperation3D_Init( this, &
    FilterOptrType,                  &
    FilterShape, FilterWidthFac,     &
    Nnode_h1D_reconst, &    
    sfield_num, hvfield_num, htensorfield_num, &
    mesh3D )
    implicit none
    class(MeshFieldFilterOperation3D), intent(inout) :: this
    class(MeshBase3D), intent(in) :: mesh3D
    character(*), intent(in) :: FilterOptrType
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: FilterWidthFac
    integer, intent(in) :: Nnode_h1D_reconst
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    !------------------------------------------

    call MeshFieldFilterOperationBase_Init( this, mesh3D%refElem3D%Nnode_h1D, &
      FilterOptrType, FilterShape, FilterWidthFac, &
      Nnode_h1D_reconst )

    this%vHaloSize = 1

    select type(mesh3D)
    class is (MeshCubeDom3D)
      call this%vars_comm%Init( sfield_num, hvfield_num, htensorfield_num, mesh3D, &
        haloSize_h1D=this%hHaloSize, haloSize_v=this%vHaloSize )
    class default
      LOG_INFO('MeshFieldFilterOperation3D_Init',*) "Unsupported mesh type is specified. Check!"
      call PRC_abort
    end select
    
    return
  end subroutine MeshFieldFilterOperation3D_Init

  !> Finalize an object to represent filter operation for 3D mesh field
!OCL SERIAL
  subroutine MeshFieldFilterOperation3D_Final( this )
    implicit none
    class(MeshFieldFilterOperation3D), intent(inout) :: this
    !------------------------------------------    
    call MeshFieldFilterOperationBase_Final( this )
    call this%vars_comm%Final()
    return
  end subroutine MeshFieldFilterOperation3D_Final

  !> Apply filter operation to 3D mesh field data
!OCL serial
  subroutine MeshFieldFilterOperation3D_apply_filter( this, q_list, &
    mesh3D )
    implicit none
    class(MeshFieldFilterOperation3D), intent(inout) :: this
    type(MeshField3D), intent(inout), target :: q_list(this%vars_comm%field_num_tot)
    class(MeshBase3D), intent(in), target :: mesh3D

    class(LocalMesh3D), pointer :: lmesh3D
    class(ElementBase3D), pointer :: elem3D
    real(RP), allocatable :: tmp3D(:,:,:,:)

    integer :: iv
    integer :: n
    type(MeshFieldContainer) :: comm_vars_list(this%vars_comm%field_num_tot)
    !-----------------------------------------------------

    do iv=1, this%vars_comm%field_num_tot
      comm_vars_list(iv)%field3d => q_list(iv)
    end do
    call this%vars_comm%Put(comm_vars_list, 1)
    call this%vars_comm%Exchange()
    call this%vars_comm%Get(comm_vars_list, 1)

    do n=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(n)
      elem3D => lmesh3D%refElem3D
      allocate( tmp3D(-this%hHaloSize+1:elem3D%Nnode_h1D+this%hHaloSize,-this%hHaloSize+1:elem3D%Nnode_h1D+this%hHaloSize,elem3D%Nnode_v,lmesh3D%Ne) )

      do iv=1, this%vars_comm%field_num_tot
        call extract_tmp3D( tmp3D, &
          q_list(iv)%local(n)%val, q_list(iv)%local(n)%val, lmesh3D, elem3D, lmesh3D%VMapP, this%hHaloSize, 0 )
        
        select case (this%operator_type)
        case (FILTER_OPTRTYPE_CONVFILTER)
          call apply_filter1d_x( q_list(iv)%local(n)%val, &
            tmp3D, this%FilterMat_h1D, elem3D%Nnode_h1D, elem3D%Nnode_h1D, elem3D%Nnode_v, lmesh3D%Ne, lmesh3D%NeA, elem3D%Nnode_h1D )
        case (FILTER_OPTRTYPE_RECONSTRUCT)
          call apply_reconst1d_x( q_list(iv)%local(n)%val, &
            tmp3D, this%Minv_Ml_tr, this%Minv_Mc_tr, this%Minv_Mr_tr, this%IntrpMat, elem3D%Nnode_h1D, elem3D%Nnode_h1D, elem3D%Nnode_v, lmesh3D%Ne, this%Nnode_h1D_reconst, lmesh3D )
        case (FILTER_OPTRTYPE_RECONSTRUCT2)
          call apply_reconst1d_x_2( q_list(iv)%local(n)%val, &
            tmp3D, this%Ml_tr, this%Mc_tr, this%Mr_tr, elem3D%Nnode_h1D, elem3D%Nnode_h1D, elem3D%Nnode_v, lmesh3D%Ne, this%Nnode_h1D_reconst, lmesh3D )
        end select

      end do
      deallocate( tmp3D )
    end do

    call this%vars_comm%Put(comm_vars_list, 1)
    call this%vars_comm%Exchange()
    call this%vars_comm%Get(comm_vars_list, 1)

    do n=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(n)
      elem3D => lmesh3D%refElem3D
      allocate( tmp3D(-this%hHaloSize+1:elem3D%Nnode_h1D+this%hHaloSize,-this%hHaloSize+1:elem3D%Nnode_h1D+this%hHaloSize,elem3D%Nnode_v,lmesh3D%Ne) )

      do iv=1, this%vars_comm%field_num_tot
        call extract_tmp3D( tmp3D, &
          q_list(iv)%local(n)%val, q_list(iv)%local(n)%val, lmesh3D, elem3D, lmesh3D%VMapP, this%hHaloSize, 1 )

        select case (this%operator_type)
        case (FILTER_OPTRTYPE_CONVFILTER)
          call apply_filter1d_y( q_list(iv)%local(n)%val, &
            tmp3D, this%FilterMat_h1D, elem3D%Nnode_h1D, elem3D%Nnode_h1D, elem3D%Nnode_v,lmesh3D%Ne, lmesh3D%NeA, elem3D%Nnode_h1D )
        case (FILTER_OPTRTYPE_RECONSTRUCT)
          call apply_reconst1d_y( q_list(iv)%local(n)%val, &
            tmp3D, this%Minv_Ml_tr, this%Minv_Mc_tr, this%Minv_Mr_tr, this%IntrpMat, elem3D%Nnode_h1D, elem3D%Nnode_h1D, elem3D%Nnode_v, lmesh3D%Ne, this%Nnode_h1D_reconst, lmesh3D )
        case (FILTER_OPTRTYPE_RECONSTRUCT2)
          call apply_reconst1d_y_2( q_list(iv)%local(n)%val, &
            tmp3D, this%Ml_tr, this%Mc_tr, this%Mr_tr, elem3D%Nnode_h1D, elem3D%Nnode_h1D, elem3D%Nnode_v, lmesh3D%Ne, this%Nnode_h1D_reconst, lmesh3D )
        end select
        
      end do
      deallocate( tmp3D )
    end do
    
    return
  end subroutine MeshFieldFilterOperation3D_apply_filter
  
!- Private -----------------------

!OCL SERIAL
  subroutine extract_tmp3D( tmp3D, q0, q0_, lmesh3D, elem3D, vmapP, hHaloSize, x_or_y )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    integer, intent(in) :: hHaloSize
    real(RP), intent(out) :: tmp3D(-hHaloSize+1:elem3D%Nnode_h1D+hHaloSize,-hHaloSize+1:elem3D%Nnode_h1D+hHaloSize,elem3D%Nnode_v,lmesh3D%Ne)
    real(RP), intent(in) :: q0(elem3D%Np*lmesh3D%NeA) 
    real(RP), intent(in) :: q0_(elem3D%Nnode_h1D,elem3D%Nnode_h1D,elem3D%Nnode_v,lmesh3D%NeA) 
    integer, intent(in) :: vmapP(elem3D%NfpTot,lmesh3D%Ne) 
    integer, intent(in) :: x_or_y !< x: 0, y: 1

    integer :: ke, ke_h, ke_x, ke_y, ke_z, px, py, ph, pz, i, f
    integer :: fso, fs, fe
    integer :: iP(elem3D%Nnode_h1D)
    real(RP) :: halo_h(elem3D%Nnode_h1D,elem3D%Nnode_v,hHaloSize,max(lmesh3D%NeX,lmesh3D%NeY),lmesh3D%NeZ,4)
    integer :: Neh_f(4)
    integer :: NehxNeZ_os(4)
    !------------------------------

    Neh_f(:) = (/ lmesh3D%NeX, lmesh3D%NeY, lmesh3D%NeX, lmesh3D%NeY /)
    NehxNeZ_os(1) = 0
    do f=2, 4
      NehxNeZ_os(f) = NehxNeZ_os(f-1) + Neh_f(f-1) * lmesh3D%NeZ 
    end do
    
    !$omp parallel private(ke, ke_x, ke_y, ke_h, ke_z, i, px, py, ph, pz, fso, fs, fe, iP, f)

    !$omp do collapse(2)
    do f=1, 4
      do ke_z=1, lmesh3D%NeZ
      do ke_h=1, Neh_f(f)
        fso = lmesh3D%Ne * elem3D%Np &
            + elem3D%Nfp_h*hHaloSize * ( NehxNeZ_os(f) + (ke_h-1) + (ke_z-1)*Neh_f(f) )

        do ph=1, hHaloSize
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D + (ph-1)*elem3D%Nfp_h
          fe = fso + elem3D%Nnode_h1D - 1
          halo_h(:,pz,ph,ke_h,ke_z,f) = q0(fs:fe)
        end do
        end do
      end do
      end do
    end do

    !$omp do collapse(3)
    do ke_z=1, lmesh3D%NeZ
    do ke_y=1, lmesh3D%NeY
    do ke_x=1, lmesh3D%NeX
      ke = ke_x + (ke_y-1) * lmesh3D%NeX + (ke_z-1) * lmesh3D%NeX * lmesh3D%NeY

      do pz=1, elem3D%Nnode_v
      do py=1, elem3D%Nnode_h1D
      do px=1, elem3D%Nnode_h1D
        tmp3D(px,py,pz,ke) = q0_(px,py,pz,ke)
      end do
      end do
      end do

      if ( x_or_y == 1 ) then ! y-direction
        ! Face 1
        fso = 0
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D; fe = fs + elem3D%Nnode_h1D - 1
          iP(:) = vmapP(fs:fe,ke)
          if ( iP(1) <= lmesh3D%Ne * elem3D%Np ) then
            do ph=1, hHaloSize
              tmp3D(1:elem3D%Nnode_h1D,-ph+1,pz,ke) = q0(iP(:) - (ph-1)*elem3D%Nnode_h1D)
            end do
          else
            do ph=1, hHaloSize
              tmp3D(1:elem3D%Nnode_h1D,-ph+1,pz,ke) = halo_h(:,pz,ph,ke_x,ke_z,1)
            end do
          end if
        enddo

        ! Face 3
        fso = 2 * elem3D%Nfp_h
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D; fe = fs + elem3D%Nnode_h1D - 1
          iP(:) = vmapP(fs:fe,ke)  
          if ( iP(1) <= lmesh3D%Ne * elem3D%Np ) then
            do ph=1, hHaloSize
              tmp3D(1:elem3D%Nnode_h1D,ph+elem3D%Nnode_h1D,pz,ke) = q0(iP(:) + (ph-1)*elem3D%Nnode_h1D)
            end do
          else
            do ph=1, hHaloSize
              tmp3D(1:elem3D%Nnode_h1D,ph+elem3D%Nnode_h1D,pz,ke) = halo_h(:,pz,ph,ke_x,ke_z,3)
            end do
          end if          
        enddo
      else if ( x_or_y == 0 ) then ! x-direction
        ! Face 2
        fso = elem3D%Nfp_h
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D; fe = fs + elem3D%Nnode_h1D - 1
          iP(:) = vmapP(fs:fe,ke)  
          if ( iP(1) <= lmesh3D%Ne * elem3D%Np ) then
            do ph=1, hHaloSize
              tmp3D(ph+elem3D%Nnode_h1D,1:elem3D%Nnode_h1D,pz,ke) = q0(iP(:) + (ph-1))
            end do
          else
            do ph=1, hHaloSize
              tmp3D(ph+elem3D%Nnode_h1D,1:elem3D%Nnode_h1D,pz,ke) = halo_h(:,pz,ph,ke_y,ke_z,2)
            end do
          end if           
        enddo

        ! Face 4
        fso = 3 * elem3D%Nfp_h
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D; fe = fs + elem3D%Nnode_h1D - 1
          iP(:) = vmapP(fs:fe,ke)  
          if ( iP(1) <= lmesh3D%Ne * elem3D%Np ) then
            do ph=1, hHaloSize
              tmp3D(-ph+1,1:elem3D%Nnode_h1D,pz,ke) = q0(iP(:) - (ph-1))
            end do
          else
            do ph=1, hHaloSize
              tmp3D(-ph+1,1:elem3D%Nnode_h1D,pz,ke) = halo_h(:,pz,ph,ke_y,ke_z,4)
            end do
          end if           
        enddo
      end if
    
    end do
    end do
    end do

    !$omp end parallel
    return
  end subroutine extract_tmp3D
  
end module scale_meshfield_filter_operation_3d