#include "scalelib.h"
program test_mesh_rectdom2d_quadrial
  use scale_precision
  use scale_prc
  use scale_io  
  use scale
  use scale_element_quadrilateral
  use scale_mesh_rectdom2d
  use scale_localmesh_2d

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
      NeGX, NeGY,                             &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true.,                         &
      refElem, NLocalMeshPerPrc )
    
    call mesh%Generate()
    
  end subroutine init

  subroutine final()
    
    call mesh%Final()
    call refElem%Final()
    
    call PRC_MPIfinish()

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

    !---------------------------

    !-- Set answer

    Np = elem%Np
    Nfp = elem%Nfp

    do j=1, NeGX
    do i=1, NeGX
      k = i + (j-1)*NeGX
      EToE_ans(:,k) = (/ i + (j-2)*NeGX, i+1 +(j-1)*NeGX, i + j*NeGX, i-1 +(j-1)*NeGX/)
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
      if (i==NeGX) then
        EToE_ans(2,k) = k; EToF_ans(2,k) = 2
        do p=1, Nfp
         VMapP_ans(p,2,k) = lcmesh%Ne*Np + Nfp*NeGX + (p + (j-1)*Nfp)
        end do
      end if
      if (j==NeGY) then
        EToE_ans(3,k) = k; EToF_ans(3,k) = 3
        do p=1, Nfp
          VMapP_ans(p,3,k) = lcmesh%Ne*Np + Nfp*NeGX + Nfp*NeGY  + (p + (i-1)*Nfp)
        end do
      end if
      if (i==1) then
        EToE_ans(4,k) = k; EToF_ans(4,k) = 4
        do p=1, Nfp
          VMapP_ans(p,4,k) = lcmesh%Ne*Np + 2*Nfp*NeGX + Nfp*NeGY  + (p + (j-1)*Nfp)
        end do
      end if
    end do
    end do
    

    !-- Check the connectivity of 2D mesh. 

    write(*,*) "** my_rank=", lcmesh%PRC_myrank
    write(*,*) " tileID:", lcmesh%tileID
    write(*,*) " pnlID:", lcmesh%panelID!, "-- i (within a panel)=", pi_table(tileID)
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
      call assert(k, lcmesh%EToE(k,:), EToE_ans(:,k), "EtoE", elem%Nfaces)
      call assert(k, lcmesh%EToF(k,:), EToF_ans(:,k), "EtoF", elem%Nfaces)
      call assert(k, lcmesh%VMapM(:,k), VMapM_ans(:,:,k), "VMapM", Nfp*elem%Nfaces)
      call assert(k, lcmesh%VMapP(:,k), VMapP_ans(:,:,k), "VMapP", Nfp*elem%Nfaces)
      write(*,*) "--------------------------"
    end do 
    write(*,*) "***********************************************************************"
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
    if ( sum((vals(:) - ans(:))**2) > EPS ) then
      LOG_ERROR('chechk_connectivity',*) 'The value of '//trim(name)//' is unexcepted!', &
        ' k=', k, ": val=", vals(:), " ans=", ans(:)
      call PRC_abort
    end if    
  end subroutine assert

  !--------------------------
  
end program test_mesh_rectdom2d_quadrial
