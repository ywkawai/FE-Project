!-------------------------------------------------------------------------------
!> module FElib / Data / Filter operation 2D
!!
!! @par Description
!!           This module provides classes to apply filter operations to MeshField data
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_meshfield_filter_operation_2d
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
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_2d, only: LocalMesh2D

  use scale_meshfield_base, only: &
    MeshField2D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_mesh_rectdom2d, only: &
    MeshRectDom2D
  use scale_meshfieldcomm_rectdom2d, only: &
    MeshFieldCommRectDom2D
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

  !> Derived type to represent filter operation for 2D mesh field
  type, public, extends(MeshFieldFilterOperationBase) :: MeshFieldFilterOperation2D
    type(MeshFieldCommRectDom2D) :: vars_comm
  contains
    procedure :: Init => MeshFieldFilterOperation2D_Init
    procedure :: Final => MeshFieldFilterOperation2D_Final
    procedure :: Apply => MeshFieldFilterOperation2D_apply_filter
  end type MeshFieldFilterOperation2D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private type & procedure
  !
  private :: extract_tmp2D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  !> Initialize an object to represent filter operation for 2D mesh field
!OCL SERIAL
  subroutine MeshFieldFilterOperation2D_Init( this, &
    FilterOptrType,                  &
    FilterShape, FilterWidthFac,     &
    Nnode_h1D_reconst,               &
    sfield_num, hvfield_num, htensorfield_num, &
    mesh2D )
    implicit none
    class(MeshFieldFilterOperation2D), intent(inout) :: this
    class(MeshBase2D), intent(in) :: mesh2D
    character(*), intent(in) :: FilterOptrType
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: FilterWidthFac
    integer, intent(in) :: Nnode_h1D_reconst
    integer, intent(in) :: sfield_num
    integer, intent(in) :: hvfield_num
    integer, intent(in) :: htensorfield_num
    !------------------------------------------

    call MeshFieldFilterOperationBase_Init( this, mesh2D%refElem2D%Nfp, &
      FilterOptrType, FilterShape, FilterWidthFac, &
      Nnode_h1D_reconst )

    select type(mesh2D)
    class is (MeshRectDom2D)
      call this%vars_comm%Init( sfield_num, hvfield_num, htensorfield_num, mesh2D, &
        haloSize_1D=this%hHaloSize )
    class default
      LOG_INFO('MeshFieldFilterOperation2D_Init',*) "Unsupported mesh type is specified. Check!"
      call PRC_abort
    end select
    return
  end subroutine MeshFieldFilterOperation2D_Init

  !> Finalize an object to represent filter operation for 2D mesh field
!OCL SERIAL
  subroutine MeshFieldFilterOperation2D_Final( this )
    implicit none
    class(MeshFieldFilterOperation2D), intent(inout) :: this
    !------------------------------------------
    call MeshFieldFilterOperationBase_Final( this )
    call this%vars_comm%Final()
    return
  end subroutine MeshFieldFilterOperation2D_Final

  !> Apply filter operation to 2D mesh field data
!OCL serial
  subroutine MeshFieldFilterOperation2D_apply_filter( this, q_list, &
    mesh2D )
    implicit none
    class(MeshFieldFilterOperation2D), intent(inout) :: this
    type(MeshField2D), intent(inout), target :: q_list(this%vars_comm%field_num_tot)
    class(MeshBase2D), intent(in), target :: mesh2D

    class(LocalMesh2D), pointer :: lmesh2D
    class(ElementBase2D), pointer :: elem2D
    real(RP), allocatable :: tmp2D(:,:,:)

    integer :: iv
    integer :: n
    type(MeshFieldContainer) :: comm_vars_list(this%vars_comm%field_num_tot)
    !-----------------------------------------------------

    do iv=1, this%vars_comm%field_num_tot
      comm_vars_list(iv)%field2d => q_list(iv)
    end do
    call this%vars_comm%Put(comm_vars_list, 1)
    call this%vars_comm%Exchange()
    call this%vars_comm%Get(comm_vars_list, 1)

    do n=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(n)
      elem2D => lmesh2D%refElem2D
      allocate( tmp2D(-this%hHaloSize+1:elem2D%Nfp+this%hHaloSize,-this%hHaloSize+1:elem2D%Nfp+this%hHaloSize,lmesh2D%Ne) )

      do iv=1, this%vars_comm%field_num_tot
        call extract_tmp2D( tmp2D, &
          q_list(iv)%local(n)%val, q_list(iv)%local(n)%val, lmesh2D, elem2D, lmesh2D%VMapP, this%hHaloSize, 0 )
        
        select case (this%operator_type)
        case (FILTER_OPTRTYPE_CONVFILTER)
          call apply_filter1d_x( q_list(iv)%local(n)%val, &
            tmp2D, this%FilterMat_h1D, elem2D%Nfp, elem2D%Nfp, 1, lmesh2D%Ne, lmesh2D%NeA, elem2D%Nfp )
        case (FILTER_OPTRTYPE_RECONSTRUCT)
          call apply_reconst1d_x( q_list(iv)%local(n)%val, &
            tmp2D, this%Minv_Ml_tr, this%Minv_Mc_tr, this%Minv_Mr_tr, this%IntrpMat, elem2D%Nfp, elem2D%Nfp, 1, lmesh2D%Ne, this%Nnode_h1D_reconst, lmesh2D )
        case (FILTER_OPTRTYPE_RECONSTRUCT2)
          call apply_reconst1d_x_2( q_list(iv)%local(n)%val, &
            tmp2D, this%Ml_tr, this%Mc_tr, this%Mr_tr, elem2D%Nfp, elem2D%Nfp, 1, lmesh2D%Ne, this%Nnode_h1D_reconst, lmesh2D )
        end select

      end do
      deallocate( tmp2D )
    end do

    call this%vars_comm%Put(comm_vars_list, 1)
    call this%vars_comm%Exchange()
    call this%vars_comm%Get(comm_vars_list, 1)

    do n=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(n)
      elem2D => lmesh2D%refElem2D
      allocate( tmp2D(-this%hHaloSize+1:elem2D%Nfp+this%hHaloSize,-this%hHaloSize+1:elem2D%Nfp+this%hHaloSize,lmesh2D%Ne) )

      do iv=1, this%vars_comm%field_num_tot
        call extract_tmp2D( tmp2D, &
          q_list(iv)%local(n)%val, q_list(iv)%local(n)%val, lmesh2D, elem2D, lmesh2D%VMapP, this%hHaloSize, 1 )

        select case (this%operator_type)
        case (FILTER_OPTRTYPE_CONVFILTER)
          call apply_filter1d_y( q_list(iv)%local(n)%val, &
            tmp2D, this%FilterMat_h1D, elem2D%Nfp, elem2D%Nfp, 1, lmesh2D%Ne, lmesh2D%NeA, elem2D%Nfp )
        case (FILTER_OPTRTYPE_RECONSTRUCT)
          call apply_reconst1d_y( q_list(iv)%local(n)%val, &
            tmp2D, this%Minv_Ml_tr, this%Minv_Mc_tr, this%Minv_Mr_tr, this%IntrpMat, elem2D%Nfp, elem2D%Nfp, 1, lmesh2D%Ne, this%Nnode_h1D_reconst, lmesh2D )
        case (FILTER_OPTRTYPE_RECONSTRUCT2)
          call apply_reconst1d_y_2( q_list(iv)%local(n)%val, &
            tmp2D, this%Ml_tr, this%Mc_tr, this%Mr_tr, elem2D%Nfp, elem2D%Nfp, 1, lmesh2D%Ne, this%Nnode_h1D_reconst, lmesh2D )
        end select
        
      end do
      deallocate( tmp2D )
    end do
    
    return
  end subroutine MeshFieldFilterOperation2D_apply_filter
  
!- Private -----------------------

!OCL SERIAL
  subroutine extract_tmp2D( tmp2D, q0, q0_, lmesh2D, elem2D, vmapP, hHaloSize, x_or_y )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    integer, intent(in) :: hHaloSize
    real(RP), intent(out) :: tmp2D(-hHaloSize+1:elem2D%Nfp+hHaloSize,-hHaloSize+1:elem2D%Nfp+hHaloSize,lmesh2D%Ne)
    real(RP), intent(in) :: q0(elem2D%Np*lmesh2D%NeA) 
    real(RP), intent(in) :: q0_(elem2D%Nfp,elem2D%Nfp,lmesh2D%NeA) 
    integer, intent(in) :: vmapP(elem2D%NfpTot,lmesh2D%Ne) 
    integer, intent(in) :: x_or_y !< x: 0, y: 1

    integer :: ke, ke_h, ke_x, ke_y, px, py, ph, i, f
    integer :: fso, fs, fe
    integer :: iP(elem2D%Nfp)
    real(RP) :: halo_h(elem2D%Nfp,hHaloSize,max(lmesh2D%NeX,lmesh2D%NeY),4)
    integer :: Neh_f(4)
    integer :: Neh_os(4)
    !------------------------------

    Neh_f(:) = (/ lmesh2D%NeX, lmesh2D%NeY, lmesh2D%NeX, lmesh2D%NeY /)
    Neh_os(1) = 0
    do f=2, 4
      Neh_os(f) = Neh_os(f-1) + Neh_f(f-1)
    end do
    
    !$omp parallel private(ke, ke_x, ke_y, ke_h, i, px, py, ph, fso, fs, fe, iP, f)

    !$omp do
    do f=1, 4
      do ke_h=1, Neh_f(f)
        fso = lmesh2D%Ne * elem2D%Np &
            + elem2D%Nfp*hHaloSize * ( Neh_os(f) + (ke_h-1) )

        do ph=1, hHaloSize
          fs = fso + 1 + (ph-1)*elem2D%Nfp
          fe = fso + elem2D%Nfp - 1
          halo_h(:,ph,ke_h,f) = q0(fs:fe)
        end do
      end do
    end do

    !$omp do collapse(2)
    do ke_y=1, lmesh2D%NeY
    do ke_x=1, lmesh2D%NeX
      ke = ke_x + (ke_y-1) * lmesh2D%NeX

      do py=1, elem2D%Nfp
      do px=1, elem2D%Nfp
        tmp2D(px,py,ke) = q0_(px,py,ke)
      end do
      end do

      if ( x_or_y == 1 ) then ! y-direction
        ! Face 1
        fso = 0
        fs = fso + 1; fe = fs + elem2D%Nfp - 1
        iP(:) = vmapP(fs:fe,ke)
        if ( iP(1) <= lmesh2D%Ne * elem2D%Np ) then
          do ph=1, hHaloSize
            tmp2D(1:elem2D%Nfp,-ph+1,ke) = q0(iP(:) - (ph-1)*elem2D%Nfp)
          end do
        else
          do ph=1, hHaloSize
            tmp2D(1:elem2D%Nfp,-ph+1,ke) = halo_h(:,ph,ke_x,1)
          end do
        end if

        ! Face 3
        fso = 2 * elem2D%Nfp
        fs = fso + 1; fe = fs + elem2D%Nfp - 1
        iP(:) = vmapP(fs:fe,ke)  
        if ( iP(1) <= lmesh2D%Ne * elem2D%Np ) then
          do ph=1, hHaloSize
            tmp2D(1:elem2D%Nfp,ph+elem2D%Nfp,ke) = q0(iP(:) + (ph-1)*elem2D%Nfp)
          end do
        else
          do ph=1, hHaloSize
            tmp2D(1:elem2D%Nfp,ph+elem2D%Nfp,ke) = halo_h(:,ph,ke_x,3)
          end do
        end if          
      else if ( x_or_y == 0 ) then ! x-direction
        ! Face 2
        fso = elem2D%Nfp
        fs = fso + 1; fe = fs + elem2D%Nfp - 1
        iP(:) = vmapP(fs:fe,ke)  
        if ( iP(1) <= lmesh2D%Ne * elem2D%Np ) then
          do ph=1, hHaloSize
            tmp2D(ph+elem2D%Nfp,1:elem2D%Nfp,ke) = q0(iP(:) + (ph-1))
          end do
        else
          do ph=1, hHaloSize
            tmp2D(ph+elem2D%Nfp,1:elem2D%Nfp,ke) = halo_h(:,ph,ke_y,2)
          end do
        end if           

        ! Face 4
        fso = 3 * elem2D%Nfp
        fs = fso + 1; fe = fs + elem2D%Nfp - 1
        iP(:) = vmapP(fs:fe,ke)  
        if ( iP(1) <= lmesh2D%Ne * elem2D%Np ) then
          do ph=1, hHaloSize
            tmp2D(-ph+1,1:elem2D%Nfp,ke) = q0(iP(:) - (ph-1))
          end do
        else
          do ph=1, hHaloSize
            tmp2D(-ph+1,1:elem2D%Nfp,ke) = halo_h(:,ph,ke_y,4)
          end do
        end if         
      end if
    end do
    end do

    !$omp end parallel
    return
  end subroutine extract_tmp2D
end module scale_meshfield_filter_operation_2d