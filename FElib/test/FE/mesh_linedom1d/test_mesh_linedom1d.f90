#include "scalelib.h"
program test_mesh_linedom1d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale
  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_equal
  implicit none

  integer, parameter :: NeGX = 4
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP

  type(LineElement) :: refElem
  integer, parameter :: PolyOrder = 2
  
  type(MeshLineDom1D) :: mesh
  integer :: n

  !--------------------------------------------------------

  call init()
  do n=1, mesh%LOCAL_MESH_NUM
    call check_connectivity( mesh%lcmesh_list(n) )
  end do
  call final()

contains
  subroutine init()
    integer :: comm, myrank, nprocs
    logical :: ismaster

    call UNITTEST_init( "test_mesh_linedom1d",          &
      comm=comm, myrank=myrank, nprocs=nprocs, ismaster=ismaster )

    !------
    call refElem%Init(PolyOrder, .true.)

    call mesh%Init( &
      NeGX,                              &
      dom_xmin, dom_xmax,                &
      refElem, NLocalMeshPerPrc )
    
    call mesh%Generate()
    
  end subroutine init

  subroutine final()
    
    call mesh%Final()
    call refElem%Final()
    
    call UNITTEST_final()

  end subroutine final

  subroutine check_connectivity( lcmesh )
    type(LocalMesh1D), intent(in) :: lcmesh

    integer :: k
    integer :: EToE_ans(2,NeGX)
    integer :: EToF_ans(2,NeGX)
    integer :: VMapM_ans(2,NeGX)
    integer :: VMapP_ans(2,NeGX)
    integer :: Np

    !---------------------------

    !-- Set answer

    Np = lcmesh%refElem%Np

    do k=lcmesh%NeS, lcmesh%NeE
      EToE_ans(:,k) = (/ k-1, k+1 /)
      EToF_ans(:,k) = (/ 2, 1 /)
      VMapM_ans(:,k) = (/ 1+(k-1)*Np, k*Np /)   
      VMapP_ans(:,k) = (/ (k-1)*Np, 1+k*Np /)
    end do
    
    EToE_ans(1,lcmesh%NeS)  = 1
    EToF_ans(1,lcmesh%NeS)  = 1
    VMapP_ans(1,lcmesh%NeS) = lcmesh%Ne*Np + 1
    
    EToE_ans(2,lcmesh%NeE) = lcmesh%NeE
    EToF_ans(2,lcmesh%NeE) = 2
    VMapP_ans(2,lcmesh%NeE) = lcmesh%Ne*Np + 2

    !-- Check the connectivity of 1D mesh. 

    write(*,*) "** my_rank=", lcmesh%PRC_myrank
    write(*,*) " tileID:", lcmesh%tileID
    write(*,*) " pnlID:", lcmesh%panelID!, "-- i (within a panel)=", pi_table(tileID)
    write(*,*) " local mesh:", n, "( total", mesh%LOCAL_MESH_NUM, ")"
    write(*,*) " panel_connect:", mesh%tilePanelID_globalMap(:,lcmesh%tileID)
    write(*,*) " tile_connect:", mesh%tileID_globalMap(:,lcmesh%tileID)
    write(*,*) " face_connect:", mesh%tileFaceID_globalMap(:,lcmesh%tileID)
    write(*,*) " domain size"
    write(*,*) "   NeX:", lcmesh%Ne
    write(*,*) "   [X]:",  lcmesh%xmin, lcmesh%xmax   

    do k=1, lcmesh%Ne
      write(*,*) "k=", k
      call UNITTEST_assert_equal( lcmesh%EToE(k,:), EToE_ans(:,k), "EtoE", idx=k )
      call UNITTEST_assert_equal( lcmesh%EToF(k,:), EToF_ans(:,k), "EtoF", idx=k )
      call UNITTEST_assert_equal( lcmesh%VMapM(:,k), VMapM_ans(:,k), "VMapM", idx=k )
      call UNITTEST_assert_equal( lcmesh%VMapP(:,k), VMapP_ans(:,k), "VMapP", idx=k )
      write(*,*) "--------------------------"
    end do 


    write(*,*) "***********************************************************************"
  end subroutine check_connectivity


end program test_mesh_linedom1d
