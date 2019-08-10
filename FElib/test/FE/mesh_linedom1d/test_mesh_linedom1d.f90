program test_mesh_linedom1d
  use scale_precision
  use scale_prc
  use scale_io  
  use scale
  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  implicit none

  integer, parameter :: NeGX = 4
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP

  type(LineElement) :: refElem
  integer, parameter :: PolyOrder = 2
  
  type(MeshLineDom1D) :: mesh
  integer :: n

  call init()
  do n=1, mesh%LOCAL_MESH_NUM
    call check_connectivity( mesh%lcmesh_list(n) )
  end do
  call final()

contains
  subroutine init()
    integer :: comm, myrank, nprocs
    logical :: ismaster

    call PRC_MPIstart( comm )
    
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
      call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call IO_setup( "test", allow_noconf = .true. )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
    
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
    
    call PRC_MPIfinish()

  end subroutine final

  subroutine check_connectivity( lcmesh )
    type(LocalMesh1D), intent(in) :: lcmesh

    integer :: k

    !---------------------------

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
      write(*,*) "EToE=", lcmesh%EToE(k,:)
      write(*,*) "EToF=", lcmesh%EToF(k,:)
      write(*,*) "VMapM=", lcmesh%VMapM(:,k)
      write(*,*) "VMapP=", lcmesh%VMapP(:,k)
      write(*,*) "--------------------------"
    end do 
    write(*,*) "***********************************************************************"
  end subroutine check_connectivity

end program test_mesh_linedom1d
