#include "scalelib.h"
program test_mesh_rectdom2d_quadrial
  use scale_precision
  use scale_prc
  use scale_io  
  use scale
  use scale_element_quadrilateral
  use scale_mesh_rectdom2d
  use scale_localmesh_2d

  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_equal
  implicit none

  integer, parameter :: NeGX = 4
  integer, parameter :: NeGY = 4
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin = -1.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP

  type(QuadrilateralElement) :: refElem
  integer, parameter :: PolyOrder = 2
  
  type(MeshRectDom2D) :: mesh
  integer :: n

  !-------------------------------------------------

  call init()
  do n=1, mesh%LOCAL_MESH_NUM
    call check_connectivity( mesh%lcmesh_list(n), refElem )
  end do  
  call final()

contains
  subroutine init()
    integer :: comm, myrank, nprocs
    logical :: ismaster

    call UNITTEST_init( "test_mesh_rectdom2d_quadrilateral",          &
      comm=comm, myrank=myrank, nprocs=nprocs, ismaster=ismaster )

    !------
    call refElem%Init(PolyOrder, .true.)

    call mesh%Init( &
      NeGX, NeGY,                             &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true.,                         &
      refElem, NLocalMeshPerPrc, 1, 1 )
    
    call mesh%Generate()
    
  end subroutine init

  subroutine final()
    
    call mesh%Final()
    call refElem%Final()
    
    call UNITTEST_final()

  end subroutine final

  subroutine check_connectivity( lcmesh, elem )
    type(LocalMesh2D), intent(in) :: lcmesh
    type(QuadrilateralElement), intent(in) :: elem 

    integer :: k, i, j, f, p
    integer :: EToE_ans(elem%Nfaces,lcmesh%Ne)
    integer :: EToF_ans(elem%Nfaces,lcmesh%Ne)
    integer :: VMapM_ans(elem%Nfp,elem%Nfaces,lcmesh%Ne)
    integer :: VMapP_ans(elem%Nfp,elem%Nfaces,lcmesh%Ne)
    integer :: Np, Nfp
    integer :: NeX, NeY

    !---------------------------

    !-- Set answer

    Np = elem%Np
    Nfp = elem%Nfp
    NeX = lcmesh%NeX
    NeY = lcmesh%NeY
    
    do j=1, NeY
    do i=1, NeX
      k = i + (j-1)*NeX
      EToE_ans(:,k) = (/ i + (j-2)*NeX, i+1 +(j-1)*NeX, i + j*NeX, i-1 +(j-1)*NeX/)
      EToF_ans(:,k) = (/ 3, 4, 1, 2 /)
      do p=1, elem%Nfp
        VMapM_ans(p,:,k) = (k-1)*Np &
          + (/ p, Nfp + (p-1)*Nfp, p + (Nfp-1)*Nfp, 1 + (p-1)*Nfp /)
        VMapP_ans(p,:,k) = (EToE_ans(:,k)-1)*Np &
          + (/ p + (Nfp-1)*Nfp, 1 + (p-1)*Nfp, p, Nfp + (p-1)*Nfp /)
      end do

      if (j==1) then
        EToE_ans(1,k) = k; EToF_ans(1,k) = 1
        do p=1, Nfp
         VMapP_ans(p,1,k) = lcmesh%Ne*Np + (p + (i-1)*Nfp)
        end do
      end if
      if (i==NeX) then
        EToE_ans(2,k) = k; EToF_ans(2,k) = 2
        do p=1, Nfp
         VMapP_ans(p,2,k) = lcmesh%Ne*Np + Nfp*NeX + (p + (j-1)*Nfp)
        end do
      end if
      if (j==NeY) then
        EToE_ans(3,k) = k; EToF_ans(3,k) = 3
        do p=1, Nfp
          VMapP_ans(p,3,k) = lcmesh%Ne*Np + Nfp*NeX + Nfp*NeY  + (p + (i-1)*Nfp)
        end do
      end if
      if (i==1) then
        EToE_ans(4,k) = k; EToF_ans(4,k) = 4
        do p=1, Nfp
          VMapP_ans(p,4,k) = lcmesh%Ne*Np + 2*Nfp*NeX + Nfp*NeY  + (p + (j-1)*Nfp)
        end do
      end if
    end do
    end do
    
    !-- Check the connectivity of 2D mesh. 

    write(*,*) "** my_rank=", lcmesh%PRC_myrank
    write(*,*) " tileID:", lcmesh%tileID
    !write(*,*) " pnlID:", lcmesh%panelID, "-- i (within a panel)=", pi_table(tileID)
    write(*,*) " local mesh:", n, "( total", mesh%LOCAL_MESH_NUM, ")"
    write(*,*) " panel_connect:", mesh%tilePanelID_globalMap(:,lcmesh%tileID)
    write(*,*) " tile_connect:", mesh%tileID_globalMap(:,lcmesh%tileID)
    write(*,*) " face_connect:", mesh%tileFaceID_globalMap(:,lcmesh%tileID)
    write(*,*) " domain size"
    write(*,*) "   NeX, NeY:", lcmesh%NeX, lcmesh%NeY
    write(*,*) "   [X]     :",  lcmesh%xmin, lcmesh%xmax   
    write(*,*) "   [Y]     :",  lcmesh%ymin, lcmesh%ymax   
    
    do k=1, lcmesh%Ne
      write(*,*) "k=", k
      call UNITTEST_assert_equal( lcmesh%EToE(k,:), EToE_ans(:,k), "EtoE", idx=k )
      call UNITTEST_assert_equal( lcmesh%EToF(k,:), EToF_ans(:,k), "EtoF", idx=k )
      call UNITTEST_assert_equal( lcmesh%VMapM(:,k), VMapM_ans(:,:,k), "VMapM", idx=k )
      call UNITTEST_assert_equal( lcmesh%VMapP(:,k), VMapP_ans(:,:,k), "VMapP", idx=k )
      write(*,*) "--------------------------"
    end do 
    write(*,*) "***********************************************************************"
  end subroutine check_connectivity


  !--------------------------
  
end program test_mesh_rectdom2d_quadrial
