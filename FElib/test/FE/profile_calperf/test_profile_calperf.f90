#include "scalelib.h"
program test_profile_calperf
  use scale_precision
  use scale_prc
  use scale_prof
  use scale_io  
  use scale_file_history
  use scale_const, only: &
    const_setup,         &
    PI => CONST_PI
  use scale

  use scale_sparsemat
  use scale_element_base
  use scale_element_hexahedral
  use scale_localmesh_3d
  use scale_mesh_cubedom3d
  
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D

  implicit none

  integer, parameter :: NprcX = 1
  integer, parameter :: NeX   = 8
  integer, parameter :: NprcY = 4
  integer, parameter :: NeY   = 8
  integer, parameter :: NeZ   = 8
  integer, parameter :: NLocalMeshPerPrc = 1

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin = -1.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP
  real(RP), parameter :: dom_zmin = -1.0_RP
  real(RP), parameter :: dom_zmax = +1.0_RP

  type(HexahedralElement) :: refElem
  integer, parameter :: PolyOrder_h = 7
  integer, parameter :: PolyOrder_v = 7
  type(SparseMat) :: Dx_csr, Dy_csr, Dz_csr, Lift_csr  
  type(SparseMat) :: Dx_ell, Dy_ell, Dz_ell, Lift_ell
  
  type(MeshCubeDom3D), target :: mesh

  type(MeshField3D), target :: q(5)

  type(MeshFieldCommCubeDom3D) :: fields_comm

  integer :: n
  type(LocalMesh3D), pointer :: lcmesh  
  !----------------------------------------------

  call init()

  !----

  call PROF_setprefx('MAIN')

  do n=1, 20
    LOG_PROGRESS('(i5,a)') n, "* Profile data communication.."
    call perf_comm()
    LOG_PROGRESS('(i5,a)') n, "* Profile SpMV operation .."
    call perf_sparsemat( q(1)%local(1)%val, &
      Dx_csr, Dy_csr, Dz_csr, Lift_csr,     &
      mesh%lcmesh_list(1), refElem, '_CSR'  )

    call perf_sparsemat( q(1)%local(1)%val, &
      Dx_ell, Dy_ell, Dz_ell, Lift_ell,     &
      mesh%lcmesh_list(1), refElem, '_ELL'  )

    LOG_PROGRESS('(i5,a)') n, "* Profile SpMV operation (OMP).."

    call perf_sparsemat_omp( q(1)%local(1)%val,  &
      Dx_csr, Dy_csr, Dz_csr, Lift_csr,          &
      mesh%lcmesh_list(1), refElem, '_CSR'       )

    call perf_sparsemat_omp( q(1)%local(1)%val,  &
      Dx_ell, Dy_ell, Dz_ell, Lift_ell,          &
      mesh%lcmesh_list(1), refElem, '_ELL'       )      

    LOG_PROGRESS('(i5,a)') n, "* Profile SpMV operation (MultiVar1).."

    call perf_sparsemat_multivar1( &
      q(1)%local(1)%val, q(2)%local(1)%val, &
      q(3)%local(1)%val, q(4)%local(1)%val, &
      q(5)%local(1)%val,                    &
      Dx_csr, Dy_csr, Dz_csr, Lift_csr,     &
      mesh%lcmesh_list(1), refElem, '_CSR'  )

    call perf_sparsemat_multivar1( &
      q(1)%local(1)%val, q(2)%local(1)%val, &
      q(3)%local(1)%val, q(4)%local(1)%val, &
      q(5)%local(1)%val,                    &
      Dx_ell, Dy_ell, Dz_ell, Lift_ell,     &
      mesh%lcmesh_list(1), refElem, '_ELL'  )

    LOG_PROGRESS('(i5,a)') n, "* Profile SpMV operation (MultiVar2).."

    call perf_sparsemat_multivar2( &
      q(1)%local(1)%val, q(2)%local(1)%val, &
      q(3)%local(1)%val, q(4)%local(1)%val, &
      q(5)%local(1)%val,                    &
      Dx_csr, Dy_csr, Dz_csr, Lift_csr,     &
      mesh%lcmesh_list(1), refElem, '_CSR'  )

    call perf_sparsemat_multivar2( &
      q(1)%local(1)%val, q(2)%local(1)%val, &
      q(3)%local(1)%val, q(4)%local(1)%val, &
      q(5)%local(1)%val,                    &
      Dx_ell, Dy_ell, Dz_ell, Lift_ell,     &
      mesh%lcmesh_list(1), refElem, '_ELL'  )

    if( IO_L ) call flush(IO_FID_LOG)
  end do

  call final()

contains
!OCL SERIAL
  subroutine perf_sparsemat( &
    q_, Dx, Dy, Dz, Lift, lmesh, elem, mat_suffix  )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh  
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)
    type(SparseMat) :: Dx, Dy, Dz, Lift
    character(*), intent(in) :: mat_suffix

    integer :: ke
    real(RP) :: Fx(elem%Np)
    real(RP) :: Fy(elem%Np)
    real(RP) :: Fz(elem%Np)
    real(RP) :: LiftDelFlux(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    !-----------------------------------

    call PROF_rapstart('Dx'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke), Fx)
    end do 
    call PROF_rapend('Dx'//trim(mat_suffix), 0)

    call PROF_rapstart('Dy'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dy, q_(:,ke), Fy)
    end do 
    call PROF_rapend('Dy'//trim(mat_suffix), 0)

    call PROF_rapstart('Dz'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dz, q_(:,ke), Fz)
    end do 
    call PROF_rapend('Dz'//trim(mat_suffix), 0)

    del_flux(:,:) = 1.0_RP
    call PROF_rapstart('Lift'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Lift, del_flux(:,ke), LiftDelFlux)
    end do 
    call PROF_rapend('Lift'//trim(mat_suffix), 0)

    call PROF_rapstart('merge'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke), Fx)
      call sparsemat_matmul(Dy, q_(:,ke), Fy)
      call sparsemat_matmul(Dz, q_(:,ke), Fz)
      call sparsemat_matmul(Lift, del_flux(:,ke), LiftDelFlux)
    end do 
    call PROF_rapend('merge'//trim(mat_suffix), 0)

    return
  end subroutine perf_sparsemat

!OCL SERIAL
  subroutine perf_sparsemat_multivar1( &
    q1_, q2_, q3_, q4_, q5_, Dx, Dy, Dz, Lift, lmesh, elem, mat_suffix  )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh  
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(in) :: q1_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q2_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q3_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q4_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q5_(elem%Np,lmesh%NeA)
    type(SparseMat) :: Dx, Dy, Dz, Lift
    character(*), intent(in) :: mat_suffix

    integer :: ke
    real(RP) :: Fx(elem%Np)
    real(RP) :: Fy(elem%Np)
    real(RP) :: Fz(elem%Np)
    real(RP) :: LiftDelFlux(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,5)
    !-----------------------------------

    call PROF_rapstart('multivar1_Dx'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q1_(:,ke), Fx)
      call sparsemat_matmul(Dx, q2_(:,ke), Fx)
      call sparsemat_matmul(Dx, q3_(:,ke), Fx)
      call sparsemat_matmul(Dx, q4_(:,ke), Fx)
      call sparsemat_matmul(Dx, q5_(:,ke), Fx)
    end do 
    call PROF_rapend('multivar1_Dx'//trim(mat_suffix), 0)

    call PROF_rapstart('multivar1_Dy'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dy, q1_(:,ke), Fy)
      call sparsemat_matmul(Dy, q2_(:,ke), Fy)
      call sparsemat_matmul(Dy, q3_(:,ke), Fy)
      call sparsemat_matmul(Dy, q4_(:,ke), Fy)
      call sparsemat_matmul(Dy, q5_(:,ke), Fy)
    end do 
    call PROF_rapend('multivar1_Dy'//trim(mat_suffix), 0)

    call PROF_rapstart('multivar1_Dz'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dz, q1_(:,ke), Fz)
      call sparsemat_matmul(Dz, q2_(:,ke), Fz)
      call sparsemat_matmul(Dz, q3_(:,ke), Fz)
      call sparsemat_matmul(Dz, q4_(:,ke), Fz)
      call sparsemat_matmul(Dz, q5_(:,ke), Fz)
    end do 
    call PROF_rapend('multivar1_Dz'//trim(mat_suffix), 0)

    del_flux(:,:,:) = 1.0_RP
    call PROF_rapstart('multivar1_Lift'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Lift, del_flux(:,ke,1), LiftDelFlux)
      call sparsemat_matmul(Lift, del_flux(:,ke,2), LiftDelFlux)
      call sparsemat_matmul(Lift, del_flux(:,ke,3), LiftDelFlux)
      call sparsemat_matmul(Lift, del_flux(:,ke,4), LiftDelFlux)
      call sparsemat_matmul(Lift, del_flux(:,ke,5), LiftDelFlux)
    end do 
    call PROF_rapend('multivar1_Lift'//trim(mat_suffix), 0)

    call PROF_rapstart('multivar1_merge'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q1_(:,ke), Fx)
      call sparsemat_matmul(Dy, q1_(:,ke), Fy)
      call sparsemat_matmul(Dz, q1_(:,ke), Fz)
      call sparsemat_matmul(Lift, del_flux(:,ke,1), LiftDelFlux)

      call sparsemat_matmul(Dx, q2_(:,ke), Fx)
      call sparsemat_matmul(Dy, q2_(:,ke), Fy)
      call sparsemat_matmul(Dz, q2_(:,ke), Fz)
      call sparsemat_matmul(Lift, del_flux(:,ke,2), LiftDelFlux)

      call sparsemat_matmul(Dx, q3_(:,ke), Fx)
      call sparsemat_matmul(Dy, q3_(:,ke), Fy)
      call sparsemat_matmul(Dz, q3_(:,ke), Fz)
      call sparsemat_matmul(Lift, del_flux(:,ke,3), LiftDelFlux)

      call sparsemat_matmul(Dx, q4_(:,ke), Fx)
      call sparsemat_matmul(Dy, q4_(:,ke), Fy)
      call sparsemat_matmul(Dz, q4_(:,ke), Fz)
      call sparsemat_matmul(Lift, del_flux(:,ke,4), LiftDelFlux)

      call sparsemat_matmul(Dx, q5_(:,ke), Fx)
      call sparsemat_matmul(Dy, q5_(:,ke), Fy)
      call sparsemat_matmul(Dz, q5_(:,ke), Fz)
      call sparsemat_matmul(Lift, del_flux(:,ke,5), LiftDelFlux)
    end do 
    call PROF_rapend('multivar1_merge'//trim(mat_suffix), 0)

    return
  end subroutine perf_sparsemat_multivar1

!OCL SERIAL
  subroutine perf_sparsemat_multivar2( &
    q1_, q2_, q3_, q4_, q5_, Dx, Dy, Dz, Lift, lmesh, elem, mat_suffix  )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh  
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(in) :: q1_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q2_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q3_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q4_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q5_(elem%Np,lmesh%NeA)
    type(SparseMat) :: Dx, Dy, Dz, Lift
    character(*), intent(in) :: mat_suffix

    integer :: ke, p
    real(RP) :: Fx(5,elem%Np)
    real(RP) :: Fy(5,elem%Np)
    real(RP) :: Fz(5,elem%Np)
    real(RP) :: LiftDelFlux(5,elem%Np)
    real(RP) :: del_flux(5,elem%NfpTot,lmesh%Ne)
    real(RP) :: Flx(5,elem%Np)
    !-----------------------------------

    call PROF_rapstart('multivar2_Dx'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      do p=1, elem%Np
        Flx(1,p) = q1_(p,ke)
        Flx(2,p) = q2_(p,ke)
        Flx(3,p) = q3_(p,ke)
        Flx(4,p) = q4_(p,ke)
        Flx(5,p) = q5_(p,ke)
      end do
      call sparsemat_matmul(Dx, Flx, Fx)
    end do 
    call PROF_rapend('multivar2_Dx'//trim(mat_suffix), 0)

    call PROF_rapstart('multivar2_Dy'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      do p=1, elem%Np
        Flx(1,p) = q1_(p,ke)
        Flx(2,p) = q2_(p,ke)
        Flx(3,p) = q3_(p,ke)
        Flx(4,p) = q4_(p,ke)
        Flx(5,p) = q5_(p,ke)
      end do
      call sparsemat_matmul(Dy, Flx, Fy)
    end do 
    call PROF_rapend('multivar2_Dy'//trim(mat_suffix), 0)

    call PROF_rapstart('multivar2_Dz'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      do p=1, elem%Np
        Flx(1,p) = q1_(p,ke)
        Flx(2,p) = q2_(p,ke)
        Flx(3,p) = q3_(p,ke)
        Flx(4,p) = q4_(p,ke)
        Flx(5,p) = q5_(p,ke)
      end do
      call sparsemat_matmul(Dz, Flx, Fz)
    end do 
    call PROF_rapend('multivar2_Dz'//trim(mat_suffix), 0)

    del_flux(:,:,:) = 1.0_RP
    call PROF_rapstart('multivar2_Lift'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Lift, del_flux(:,:,ke), LiftDelFlux)
    end do 
    call PROF_rapend('multivar2_Lift'//trim(mat_suffix), 0)

    call PROF_rapstart('multivar2_merge'//trim(mat_suffix), 0)
    do ke=lmesh%NeS, lmesh%NeE
      do p=1, elem%Np
        Flx(1,p) = q1_(p,ke)
        Flx(2,p) = q2_(p,ke)
        Flx(3,p) = q3_(p,ke)
        Flx(4,p) = q4_(p,ke)
        Flx(5,p) = q5_(p,ke)
      end do
      call sparsemat_matmul(Dx, Flx, Fx)

      do p=1, elem%Np
        Flx(1,p) = q1_(p,ke)
        Flx(2,p) = q2_(p,ke)
        Flx(3,p) = q3_(p,ke)
        Flx(4,p) = q4_(p,ke)
        Flx(5,p) = q5_(p,ke)
      end do
      call sparsemat_matmul(Dy, Flx, Fy)

      do p=1, elem%Np
        Flx(1,p) = q1_(p,ke)
        Flx(2,p) = q2_(p,ke)
        Flx(3,p) = q3_(p,ke)
        Flx(4,p) = q4_(p,ke)
        Flx(5,p) = q5_(p,ke)
      end do
      call sparsemat_matmul(Dz, Flx, Fz)

      call sparsemat_matmul(Lift, del_flux(:,:,ke), LiftDelFlux)
    end do 
    call PROF_rapend('multivar2_merge'//trim(mat_suffix), 0)

    return
  end subroutine perf_sparsemat_multivar2

!OCL SERIAL
  subroutine perf_sparsemat_omp( &
      q_, Dx, Dy, Dz, Lift, lmesh, elem, mat_suffix  )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh  
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    character(*), intent(in) :: mat_suffix

    integer :: ke
    real(RP) :: Fx(elem%Np)
    real(RP) :: Fy(elem%Np)
    real(RP) :: Fz(elem%Np)
    real(RP) :: LiftDelFlux(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    !-----------------------------------

    call PROF_rapstart('Dx_omp'//trim(mat_suffix), 0)
    !$omp parallel do private(Fx)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke), Fx)
    end do 
    call PROF_rapend('Dx_omp'//trim(mat_suffix), 0)

    call PROF_rapstart('Dy_omp'//trim(mat_suffix), 0)
    !$omp parallel do private(Fy)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dy, q_(:,ke), Fy)
    end do 
    call PROF_rapend('Dy_omp'//trim(mat_suffix), 0)

    call PROF_rapstart('Dz_omp'//trim(mat_suffix), 0)
    !$omp parallel do private(Fz)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dz, q_(:,ke), Fz)
    end do 
    call PROF_rapend('Dz_omp'//trim(mat_suffix), 0)

    del_flux(:,:) = 1.0_RP
    call PROF_rapstart('Lift_omp'//trim(mat_suffix), 0)
    !$omp parallel do private(LiftDelFlux)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Lift, del_flux(:,ke), LiftDelFlux)
    end do 
    call PROF_rapend('Lift_omp'//trim(mat_suffix), 0)

    call PROF_rapstart('merge_omp'//trim(mat_suffix), 0)
    !$omp parallel do private(Fx, Fy, Fz, LiftDelFlux)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,ke), Fx)
      call sparsemat_matmul(Dy, q_(:,ke), Fy)
      call sparsemat_matmul(Dz, q_(:,ke), Fz)
      call sparsemat_matmul(Lift, del_flux(:,ke), LiftDelFlux)
    end do 
    call PROF_rapend('merge_omp'//trim(mat_suffix), 0)

    return
  end subroutine perf_sparsemat_omp

!------------------------------------------------
  subroutine init()
    implicit none
    integer :: comm, myrank, nprocs
    logical :: ismaster
    
    integer :: ke, p
    integer :: v
    !---------

    call PRC_MPIstart( comm )
    
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
        
    ! setup scale_io
    call IO_setup( "test_profile_calperf", "test.conf" )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
    
    ! setup PROF
    call PROF_setup

    call PROF_setprefx('INIT')
    call PROF_rapstart('Initialize', 0)

    ! setup constants
    call CONST_setup

    !------
    call refElem%Init(PolyOrder_h, PolyOrder_v, .false.)

    call Dx_csr%Init( refElem%Dx1, storage_format='CSR' )
    call Dy_csr%Init( refElem%Dx2, storage_format='CSR' )
    call Dz_csr%Init( refElem%Dx3, storage_format='CSR' )
    call Lift_csr%Init( refElem%Lift, storage_format='CSR' )
    call Lift_csr%Print()

    call Dx_ell%Init( refElem%Dx1, storage_format='ELL' )
    call Dy_ell%Init( refElem%Dx2, storage_format='ELL' )
    call Dz_ell%Init( refElem%Dx3, storage_format='ELL' )
    call Lift_ell%Init( refElem%Lift, storage_format='ELL' )
    call Lift_ell%Print()

    call mesh%Init( &
      NprcX*NeX, NprcY*NeY, NeZ,                                  &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax, &
      .true., .true., .true.,                                     &
      refElem, NLocalMeshPerPrc, NprcX, NprcY )
    
    call mesh%Generate()

    !---
    do v=1, 5
      call q(v)%Init( "q", "1", mesh )
    end do
    call fields_comm%Init(5, 0, 0, mesh)

    !---
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)

      do v=1, 5
      do ke=lcmesh%NeS, lcmesh%NeE
      do p=1, refElem%Np
        q(v)%local(n)%val(p,ke) = sin( PI * lcmesh%pos_en(p,ke,1) ) &
                                * sin( PI * lcmesh%pos_en(p,ke,2) ) &
                                * sin( PI * lcmesh%pos_en(p,ke,3) )                  
      end do
      end do
      end do
    end do
    !---

    call PROF_rapend('Initialize', 0)

    return
  end subroutine init

  subroutine final()
    implicit none

    integer :: v
    !------------------------------------------------------

    call PROF_setprefx('FIN')
    call PROF_rapstart('All', 1)

    do v=1, 5
      call q(v)%Final()
    end do
    call fields_comm%Final()

    call mesh%Final()
    call Dx_csr%Final(); call Dy_csr%Final(); call Dz_csr%Final()
    call Lift_csr%Final()    
    call Dx_ell%Final(); call Dy_ell%Final(); call Dz_ell%Final()
    call Lift_ell%Final()    
    call refElem%Final()
    
    call PROF_rapend  ('All', 1)
    call PROF_rapreport

    call PRC_MPIfinish()

    return
  end subroutine final

  subroutine perf_comm()
    implicit none
    type(MeshFieldContainer) :: field_list(5)
    integer :: v
    !---------------------------

    call PROF_rapstart('Communication', 0)    
    do v=1, 5
      field_list(v)%field3d => q(v)
    end do
    call PROF_rapstart('Communication_put', 1)    
    call fields_comm%Put(field_list, 1)
    call PROF_rapend('Communication_put', 1)    

    call PROF_rapstart('Communication_exchange', 1)    
    call fields_comm%Exchange()
    call PROF_rapend('Communication_exchange', 1)    

    call PROF_rapstart('Communication_get', 1)    
    call fields_comm%Get(field_list, 1)
    call PROF_rapend('Communication_get', 1)
    call PROF_rapend('Communication', 0)    

    return
  end subroutine perf_comm

end program test_profile_calperf