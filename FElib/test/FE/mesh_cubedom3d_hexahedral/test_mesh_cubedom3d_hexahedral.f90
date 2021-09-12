#include "scalelib.h"
program test_mesh_cubedom3d_hexahedral
  use scale_precision
  use scale_prc
  use scale_io  
  use scale
  use scale_element_hexahedral
  use scale_mesh_cubedom3d
  use scale_localmesh_3d

  implicit none

  integer, parameter :: NeGX = 3
  integer, parameter :: NeGY = 3
  integer, parameter :: NeGZ = 3
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin = -1.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP
  real(RP), parameter :: dom_zmin = -1.0_RP
  real(RP), parameter :: dom_zmax = +1.0_RP

  type(HexahedralElement) :: refElem
  integer, parameter :: PolyOrder_h = 2
  integer, parameter :: PolyOrder_v = 2
  
  type(MeshCubeDom3D) :: mesh
  integer :: lmeshID

  !-------------------------------------------------

  call init()
  do lmeshID=1, mesh%LOCAL_MESH_NUM
    call check_connectivity( mesh%lcmesh_list(lmeshID), refElem )
  end do  
  call final()

contains
  subroutine init()
    integer :: comm, myrank, nprocs
    logical :: ismaster

    !----------------------------------------------

    call PRC_MPIstart( comm )
    
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
      call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call IO_setup( "test", allow_noconf = .true. )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
    
    !------
    LOG_INFO("init",*)  'Initialize HexahedralElement ..'
    call refElem%Init(PolyOrder_h, PolyOrder_v, .true.)

    LOG_INFO("init",*)  'Initialzie MeshCubeDom3D ..'
    call mesh%Init( &
      NeGX, NeGY, NeGZ,                                           &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      .true., .true., .true.,                                     &
      refElem, NLocalMeshPerPrc, nprocs, 1 )
    
    LOG_INFO("init",*)  'Generate MeshCubeDom3D ..'    
    call mesh%Generate()

    LOG_INFO("init",*) 'Initialization has succeeded.'

    return
  end subroutine init

  subroutine final()
    LOG_INFO("init",*) 'Finalize ...'

    call mesh%Final()
    call refElem%Final()
    
    call PRC_MPIfinish()
    return
  end subroutine final

  subroutine check_connectivity( lcmesh, elem )
    type(LocalMesh3D), intent(in) :: lcmesh
    type(HexahedralElement), intent(in) :: elem 

    integer :: ke, i, j, k
    integer :: p, q, n
    integer :: EToE_ans(elem%Nfaces,lcmesh%Ne)
    integer :: EToF_ans(elem%Nfaces,lcmesh%Ne)
    integer :: VMapM_h_ans(elem%Nfp_h,elem%Nfaces_h,lcmesh%Ne)
    integer :: VMapP_h_ans(elem%Nfp_h,elem%Nfaces_h,lcmesh%Ne)
    integer :: VMapM_v_ans(elem%Nfp_v,elem%Nfaces_v,lcmesh%Ne)
    integer :: VMapP_v_ans(elem%Nfp_v,elem%Nfaces_v,lcmesh%Ne)
    integer :: Np, Nnode_h1D, Nnode_v
    integer :: Nfaces_h, Nfp_h, Nfp_v
    integer :: NeX, NeY, NeZ
    integer :: vs, ve
    !---------------------------

    !-- Set answer

    Np = elem%Np
    Nnode_h1D = elem%Nnode_h1D
    Nnode_v = elem%Nnode_v
    Nfaces_h = elem%Nfaces_h
    Nfp_h = elem%Nfp_h
    Nfp_v = elem%Nfp_v
    NeX = lcmesh%NeX
    NeY = lcmesh%NeY
    NeZ = lcmesh%NeZ
    
    LOG_INFO("check_connectivity",*) 'Set information about correct values..'

    do k=1, NeZ
    do j=1, NeY
    do i=1, NeX
      ke = i + (j-1)*NeX + (k-1)*NeX*NeY
      EToE_ans(:,ke) = (/ &
        elemID(lcmesh,i,j-1,k), elemID(lcmesh,i+1,j,k), elemID(lcmesh,i,j+1,k), elemID(lcmesh,i-1,j,k), &
        elemID(lcmesh,i,j,k-1), elemID(lcmesh,i,j,k+1) /)
      EToF_ans(:,ke) = (/ 3, 4, 1, 2, 6, 5 /)
      do q=1, Nnode_v
      do p=1, Nnode_h1D
        n = p + (q-1)*Nnode_h1D
        VMapM_h_ans(n,:,ke) = (ke-1)*Np &
          + (/ nodeID(elem,p,1,q), nodeID(elem,Nnode_h1D,p,q), nodeID(elem,p,Nnode_h1D,q), nodeID(elem,1,p,q)  /)
        VMapP_h_ans(n,:,ke) = (EToE_ans(1:Nfaces_h,ke)-1)*Np &
          + (/ nodeID(elem,p,Nnode_h1D,q), nodeID(elem,1,p,q), nodeID(elem,p,1,q), nodeID(elem,Nnode_h1D,p,q)  /)
      end do
      end do

      do q=1, Nnode_h1D
      do p=1, Nnode_h1D
        n = p + (q-1)*Nnode_h1D
        VMapM_v_ans(n,:,ke) = (ke-1)*Np &
          + (/ nodeID(elem,p,q,1), nodeID(elem,p,q,Nnode_v)  /)
        VMapP_v_ans(n,:,ke) = (EToE_ans(Nfaces_h+1:elem%Nfaces,ke)-1)*Np &
          + (/ nodeID(elem,p,q,Nnode_v), nodeID(elem,p,q,1)  /)
      end do
      end do

      if (j==1) then
        EToE_ans(1,ke) = ke; EToF_ans(1,ke) = 1
        do p=1, Nfp_h
         VMapP_h_ans(p,1,ke) = lcmesh%Ne*Np + (p + ((i-1) + (k-1)*NeX)*Nfp_h)
        end do
      end if
      if (i==NeX) then
        EToE_ans(2,ke) = ke; EToF_ans(2,ke) = 2
        do p=1, Nfp_h
         VMapP_h_ans(p,2,ke) = lcmesh%Ne*Np + Nfp_h*NeX*NeZ + (p + ((j-1) + (k-1)*NeY)*Nfp_h)
        end do
      end if
      if (j==NeY) then
        EToE_ans(3,ke) = ke; EToF_ans(3,ke) = 3
        do p=1, Nfp_h
          VMapP_h_ans(p,3,ke) = lcmesh%Ne*Np + Nfp_h*(NeX + NeY)*NeZ + (p + ((i-1) + (k-1)*NeX)*Nfp_h)
        end do
      end if
      if (i==1) then
        EToE_ans(4,ke) = ke; EToF_ans(4,ke) = 4
        do p=1, Nfp_h
          VMapP_h_ans(p,4,ke) = lcmesh%Ne*Np + Nfp_h*(2*NeX + NeY)*NeZ + (p + ((j-1) + (k-1)*NeY)*Nfp_h)
        end do
      end if
      if (k==1) then
        EToE_ans(5,ke) = ke; EToF_ans(5,ke) = 5
        do p=1, Nfp_v
          VMapP_v_ans(p,1,ke) = lcmesh%Ne*Np + 2*Nfp_h*(NeX + NeY)*NeZ + (p + ((i-1) + (j-1)*NeX)*Nfp_v)
        end do
      end if
      if (k==NeZ) then
        EToE_ans(6,ke) = ke; EToF_ans(6,ke) = 6
        do p=1, Nfp_v
          VMapP_v_ans(p,2,ke) = lcmesh%Ne*Np + 2*Nfp_h*(NeX + NeY)*NeZ + NeX*NeY*Nfp_v + (p + ((i-1) + (j-1)*NeX)*Nfp_v)
        end do
      end if
    end do
    end do
    end do
    
    !---
    LOG_INFO("check_connectivity",*) 'Check the connectivity of 3D mesh..'

    write(*,*) "** my_rank=", lcmesh%PRC_myrank
    write(*,*) " tileID:", lcmesh%tileID
    !write(*,*) " pnlID:", lcmesh%panelID, "-- i (within a panel)=", pi_table(tileID)
    write(*,*) " local mesh:", n, "( total", mesh%LOCAL_MESH_NUM, ")"
    write(*,*) " panel_connect:", mesh%tilePanelID_globalMap(:,lcmesh%tileID)
    write(*,*) " tile_connect:", mesh%tileID_globalMap(:,lcmesh%tileID)
    write(*,*) " face_connect:", mesh%tileFaceID_globalMap(:,lcmesh%tileID)
    write(*,*) " domain size"
    write(*,*) "   NeX, NeY, NeZ:", lcmesh%NeX, lcmesh%NeY, lcmesh%NeZ
    write(*,*) "   [X]     :",  lcmesh%xmin, lcmesh%xmax   
    write(*,*) "   [Y]     :",  lcmesh%ymin, lcmesh%ymax   
    write(*,*) "   [Z]     :",  lcmesh%zmin, lcmesh%zmax   
    
    do ke=1, lcmesh%Ne
      write(*,*) "k=", ke
      call assert(ke, lcmesh%EToE(ke,:), EToE_ans(:,ke), "EtoE", elem%Nfaces)
      call assert(ke, lcmesh%EToF(ke,:), EToF_ans(:,ke), "EtoF", elem%Nfaces)
      
      vs = 1
      ve = Nfp_h*elem%Nfaces_h
      call assert(ke, lcmesh%VMapM(vs:ve,ke), VMapM_h_ans(:,:,ke), "VMapM_h", Nfp_h*elem%Nfaces_h)
      call assert(ke, lcmesh%VMapP(vs:ve,ke), VMapP_h_ans(:,:,ke), "VMapP_h", Nfp_h*elem%Nfaces_h)
     
      vs = ve+1
      ve = elem%NfpTot
      call assert(ke, lcmesh%VMapM(vs:ve,ke), VMapM_v_ans(:,:,ke), "VMapM_v", Nfp_v*elem%Nfaces_v)
      call assert(ke, lcmesh%VMapP(vs:ve,ke), VMapP_v_ans(:,:,ke), "VMapP_v", Nfp_v*elem%Nfaces_v)
      write(*,*) "--------------------------"
    end do 
    write(*,*) "***********************************************************************"
    
    return
  end subroutine check_connectivity

  subroutine assert(k, vals, ans, name, val_size)
    integer, intent(in) :: val_size
    integer, intent(in) :: k
    integer, intent(in) :: vals(val_size)
    integer, intent(in) :: ans(val_size)
    character(*), intent(in) :: name

    real(RP), parameter :: EPS = 1.0E-15_RP

    !--------------------------------------

    write(*,*) trim(name), "=", vals(:)
    if ( sum(dble(vals(:) - ans(:))**2) > EPS ) then
      LOG_ERROR('chechk_connectivity',*) 'The value of '//trim(name)//' is unexcepted!', &
        ' k=', k, ": val=", vals(:), " ans=", ans(:)
      call PRC_abort
    end if    

    return
  end subroutine assert

  function elemID(lmesh, i, j, k) result(eid)
    type(LocalMesh3D), intent(in) :: lmesh
    integer, intent(in) :: i, j, k
    integer :: eid
    !--------------------------------------
    eid = i + (j-1)*lmesh%NeX + (k-1)*lmesh%NeX*lmesh%NeY
    return
  end function elemID

  function nodeID(elem, i, j, k) result(nid)
    type(HexahedralElement), intent(in) :: elem
    integer, intent(in) :: i, j, k
    integer :: nid
    !--------------------------------------
    nid = i + (j-1)*elem%Nnode_h1D + (k-1)*elem%Nnode_h1D**2
    return
  end function nodeID

  !--------------------------
  
end program test_mesh_cubedom3d_hexahedral
