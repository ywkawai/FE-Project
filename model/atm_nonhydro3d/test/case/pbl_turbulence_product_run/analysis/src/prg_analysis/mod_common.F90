#include "scaleFElib.h"
module mod_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_sparsemat, only: &
    sparsemat

  use scale_meshfield_base, only: MeshField1D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_1d, only: LocalMesh1D  
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_element_base, only: &
    ElementBase1D, ElementBase3D

  use scale_mesh_bndinfo, only: MeshBndInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: common_Init
  public :: common_Final
  public :: global_horizontal_mean
  public :: inquire_bound_flag
  
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  type(SparseMat), public :: Dx, Dy, Dz
  type(SparseMat), public :: Sx, Sy, Sz
  type(SparseMat), public :: Lift

  type(MeshBndInfo), allocatable :: VelBC_list(:)
  integer, allocatable:: velBC_ids(:)
  integer, parameter :: domBnd_South_ID = 1
  integer, parameter :: domBnd_East_ID  = 2
  integer, parameter :: domBnd_North_ID = 3
  integer, parameter :: domBnd_West_ID  = 4
  integer, parameter :: domBnd_Btm_ID   = 5
  integer, parameter :: domBnd_Top_ID   = 6
  integer, parameter :: DOM_BND_NUM     = 6  

contains
!OCL SERIAL
  subroutine common_Init( mesh )
    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID,     &
      BND_TYPE_PERIODIC_ID
        
    implicit none    
    class(MeshBase3D), intent(in) :: mesh

    integer :: n

    class(LocalMeshBase), pointer :: ptr_lcmesh
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: refElem3D
    !--------------------------------------------------------------------

    refElem3D => mesh%refElem3D
    call Dx%Init( refElem3D%Dx1, storage_format='ELL' )
    call Dy%Init( refElem3D%Dx2, storage_format='ELL' )
    call DZ%Init( refElem3D%Dx3, storage_format='ELL' )
    call Sx%Init( refElem3D%Sx1, storage_format='ELL' )
    call Sy%Init( refElem3D%Sx2, storage_format='ELL' )
    call Sz%Init( refElem3D%Sx3, storage_format='ELL' )
    call Lift%Init( refElem3D%Lift, storage_format='ELL' )


    !----
    allocate( velBC_ids(DOM_BND_NUM) )
    velBC_ids(domBnd_Btm_ID) = BND_TYPE_SLIP_ID
    velBC_ids(domBnd_Top_ID) = BND_TYPE_SLIP_ID
    velBC_ids(domBnd_North_ID) = BND_TYPE_PERIODIC_ID
    velBC_ids(domBnd_South_ID) = BND_TYPE_PERIODIC_ID
    velBC_ids(domBnd_East_ID) = BND_TYPE_PERIODIC_ID
    velBC_ids(domBnd_West_ID) = BND_TYPE_PERIODIC_ID

    allocate( velBC_list(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      call mesh%GetLocalMesh( n, ptr_lcmesh )
      select type (ptr_lcmesh)
      type is (LocalMesh3D)
        lcmesh3D => ptr_lcmesh
      end select

      call bnd_Init_lc( &
        VelBC_list(n),         & ! (inout)
        velBC_ids(:),          & ! (in)
        lcmesh3D%VMapB, mesh, lcmesh3D, lcmesh3D%refElem3D ) ! (in)
    end do    

    return
  end subroutine common_Init

!OCL SERIAL
  subroutine common_Final()
    implicit none
    !--------------------------------------------------------------------
    call Dx%Final(); call Dy%Final(); call Dz%Final()
    call Sx%Final(); call Sy%Final(); call Sz%Final()
    call Lift%Final()

    deallocate( velBC_ids, VelBC_list )

    return
  end subroutine common_Final

!OCL SERIAL
  subroutine global_horizontal_mean(field)
    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD    
    implicit none

    class(MeshField1D), intent(inout), target :: field
    real(RP), allocatable :: sum_tmp(:,:)
    real(RP), allocatable :: sum_out(:,:)
    integer :: ke_z

    class(LocalMesh1D), pointer :: lmesh1D
    class(ElementBase1D), pointer :: refElemV1D

    integer :: nproc
    integer :: ierr
    !-----------------------------------------

    lmesh1D => field%mesh%lcmesh_list(1)
    refElemV1D => lmesh1D%refElem1D

    allocate( sum_tmp(refElemV1D%Np,lmesh1D%Ne) )
    allocate( sum_out(refElemV1D%Np,lmesh1D%Ne) )
    do ke_z=1, lmesh1D%Ne
      sum_tmp(:,ke_z) = field%local(1)%val(:,ke_z)
    end do

    ! global sum
    call MPI_AllReduce( sum_tmp, sum_out, refElemV1D%Np * lmesh1D%Ne, &
      MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr       )
    
    call MPI_Comm_size( PRC_LOCAL_COMM_WORLD, nproc, ierr )

    do ke_z=1, lmesh1D%Ne
      field%local(1)%val(:,ke_z) = sum_out(:,ke_z) / dble( nproc )
    end do
      
    return
  end subroutine global_horizontal_mean

!OCL SERIAL
  subroutine inquire_bound_flag( & 
    is_bound,                                & ! (out)
    domID, vmapM, vmapP, vmapB, lmesh, elem  ) ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem    
    logical, intent(out) :: is_bound(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: domID
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: i, i_, iM, iP
    !-----------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      is_bound(i) = .false.

      if (i_ > 0) then  
        iM = vmapM(i)
        if ( VelBC_list(domID)%list(i_) == BND_TYPE_SLIP_ID ) then
          is_bound(i) = .true.
        else if ( VelBC_list(domID)%list(i_) == BND_TYPE_NOSLIP_ID ) then          
          is_bound(i) = .true.
        end if

      end if
    end do

    return
  end subroutine inquire_bound_flag

!-- private ----------------------

!OCL SERIAL
  subroutine bnd_Init_lc(   &
    velBCInfo,               &  ! (inout)
    velBCIDs,                &  ! (in)
    vmapB, mesh, lmesh, elem )  ! (in)

    use scale_mesh_bndinfo, only: BND_TYPE_NOSPEC_ID
    implicit none
    
    type(MeshBndInfo), intent(inout) :: velBCInfo
    integer, intent(in) :: velBCIDs(DOM_BND_NUM)
    class(MeshBase3D), intent(in) :: mesh
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: vmapB(:)

    integer :: tileID
    integer :: dom_bnd_sizes(DOM_BND_NUM)
    integer :: bnd_buf_size
    integer :: b, is_, ie_
    !-----------------------------------------------
    
    dom_bnd_sizes(:) = &
        elem%Nfp_h*lmesh%NeZ*(/ lmesh%NeX, lmesh%NeY, lmesh%NeX, lmesh%NeY, 0, 0 /) &
      + elem%Nfp_v*lmesh%NeX*lmesh%NeY*(/ 0, 0, 0, 0, 1, 1 /)
    bnd_buf_size = sum(dom_bnd_sizes)
    
    call velBCInfo%Init( bnd_buf_size )
    call velBCInfo%Set(1, bnd_buf_size, BND_TYPE_NOSPEC_ID)

    tileID = lmesh%tileID
    is_ = 1
    do b=1, DOM_BND_NUM
      ie_ = is_ + dom_bnd_sizes(b) - 1
      if ( mesh%tileID_globalMap(b,tileID) == tileID      &
           .and. mesh%tileFaceID_globalMap(b,tileID) == b ) then
        call velBCInfo%Set(is_, ie_, velBCIDs(b))
      end if
      is_ = ie_ + 1
    end do

    return
  end  subroutine bnd_Init_lc

end module mod_common