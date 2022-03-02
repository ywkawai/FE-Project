#include "scalelib.h"
program prg_pbltb_analysis
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort

  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D

  use scale_element_line, only: LineElement
  use scale_localmesh_1d, only: LocalMesh1D
  use scale_mesh_linedom1d, only: MeshLineDom1D

  use scale_localmesh_2d, only: LocalMesh2D

  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D

  use scale_sparsemat, only: SparseMat
  
  use scale_meshfield_base, only: &
    MeshField1D, MeshField3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldContainer
  use scale_meshfieldcomm_cubedom3d, only: &
    MeshFieldCommCubeDom3D
  
  use scale_file_base_meshfield, only: &
    FILE_base_meshfield
  
  use mod_common, only: &
    global_horizontal_mean
  
  use mod_diag_tb, only: &
    diag_tb_process
  
  use scale_mesh_base3d, only: &
    DIMTYPEID_XYZT => MeshBase3D_DIMTYPEID_XYZT
  
  implicit none

   ! MPI parameters
  integer                 :: nprocs                      ! number of processes               (execution)
  integer                 :: myrank                      ! my rank                           (execution)
  logical                 :: ismaster                    ! master process?                   (execution)

  type(HexahedralElement) :: refElem3D
  type(MeshCubeDom3D), target :: mesh3D

  type(LineElement) :: refElemV1D
  type(MeshLineDom1D) :: meshV1D

  integer :: istep

  type(MeshField3D) :: U
  type(MeshField3D) :: V
  type(MeshField3D) :: W
  type(MeshField3D), target :: PT
  type(MeshField3D), target :: PRES

  type(MeshField3D), target :: DDENS
  type(MeshField3D), target :: MOMX
  type(MeshField3D), target :: MOMY
  type(MeshField3D), target :: MOMZ
  type(MeshField3D), target :: DRHOT
  type(MeshField3D), target :: PRES_hyd
  type(MeshField3D), target :: DENS_hyd

  type(MeshFieldCommCubeDom3D) :: vars_comm
  type(MeshFieldContainer) :: comm_vars_ptr(9)

  type(MeshField1D) :: DENS_V1D
  type(MeshField1D) :: RHOT_V1D
  type(MeshField1D) :: MOMZ_V1D
  type(MeshField1D) :: PT_V1D
  type(MeshField1D) :: W_V1D
  type(MeshField1D) :: W_PRIM3_V1D
  type(MeshField1D) :: HEAT_EDDYFLX_V1D
  type(MeshField1D) :: MOMZ_EDDYFLX_V1D
  integer, parameter :: DIAGVID_DENS = 1
  integer, parameter :: DIAGVID_RHOT = 2
  integer, parameter :: DIAGVID_MOMZ = 3
  integer, parameter :: DIAGVID_PT   = 4
  integer, parameter :: DIAGVID_W    = 5
  integer, parameter :: DIAGVID_W_PRIM3      = 6
  integer, parameter :: DIAGVID_HEAT_EDDYFLX = 7
  integer, parameter :: DIAGVID_MOMZ_EDDYFLX = 8


  type(FILE_base_meshfield) :: in_file
  type(FILE_base_meshfield) :: in_bs_file
  character(len=H_LONG) :: in_filebase
  character(len=H_LONG) :: in_bs_filebase
  integer :: num_step       = 1
  real(RP) :: start_time0   = 12600.0_RP 
  real(RP) :: output_tintrv = 60.0_RP


  type(FILE_base_meshfield) :: out_file_V1D
  character(*), parameter :: dtype = 'REAL8'

  integer :: n
  class(LocalMesh3D), pointer :: lcmesh3D
  real(RP) :: harea

  real(RP) :: start_time, end_time

  !-----------------------------------------------------------------------------


  call initialize()

  !-
  do istep=1, num_step
    LOG_PROGRESS(*) "istep=", istep
    start_time = start_time0 + (istep-2)*output_tintrv
    end_time   = start_time0 + (istep-1)*output_tintrv

    if (istep==1) then
      call in_bs_file%Read_Var(DIMTYPEID_XYZT, "DENS_hyd", DENS_hyd, step=1)
      call in_bs_file%Read_Var(DIMTYPEID_XYZT, "PRES_hyd", PRES_hyd, step=1)
    end if
    call in_file%Read_Var(DIMTYPEID_XYZT, "DDENS", DDENS, step=istep)
    call in_file%Read_Var(DIMTYPEID_XYZT, "U", U, step=istep)
    call in_file%Read_Var(DIMTYPEID_XYZT, "V", V, step=istep)
    call in_file%Read_Var(DIMTYPEID_XYZT, "W", W, step=istep)
    call in_file%Read_Var(DIMTYPEID_XYZT, "PT", PT, step=istep)

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)
      call set_prgvars( &
        MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,         & ! (out)
        PRES%local(n)%val,                                                                   & ! (out)
        DDENS%local(n)%val, U%local(n)%val, V%local(n)%val, W%local(n)%val, PT%local(n)%val, & ! (in)
        DENS_hyd%local(n)%val, PRES_hyd%local(n)%val, lcmesh3D, lcmesh3D%refElem3D           ) ! (in)
    end do
    call vars_comm%Put( comm_vars_ptr, 1 )
    call vars_comm%Exchange()
    call vars_comm%Get( comm_vars_ptr, 1 )
    
    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)
      call analyze_lc_1( &
        DENS_V1D%local(1)%val, RHOT_V1D%local(1)%val, MOMZ_V1D%local(1)%val,    &
        DDENS%local(n)%val, PT%local(n)%val, W%local(n)%val,                    &
        DENS_hyd%local(n)%val,                                                  &
        n, lcmesh3D, refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D, &
        meshV1D%lcmesh_list(1), refElemV1D )
    end do

    call global_horizontal_mean( DENS_V1D )
    call global_horizontal_mean( RHOT_V1D )
    call global_horizontal_mean( MOMZ_V1D )

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)
      call analyze_lc_2( &
        PT_V1D%local(1)%val, W_V1D%local(1)%val,                                           &
        W_PRIM3_V1D%local(1)%val,                                                          &
        HEAT_EDDYFLX_V1D%local(1)%val, MOMZ_EDDYFLX_V1D%local(1)%val,                      &
        DDENS%local(n)%val, PT%local(n)%val, W%local(n)%val,                               &
        DENS_V1D%local(1)%val, RHOT_V1D%local(1)%val, MOMZ_V1D%local(1)%val,               &
        DENS_hyd%local(n)%val,                                                             &
        n, lcmesh3D, refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D,            &
        meshV1D%lcmesh_list(1), refElemV1D )
    end do

    call global_horizontal_mean( W_PRIM3_V1D )
    call global_horizontal_mean( HEAT_EDDYFLX_V1D )
    call global_horizontal_mean( MOMZ_EDDYFLX_V1D )

    call diag_tb_process( istep, start_time, end_time, &
      DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
      DENS_hyd, PRES_hyd, PRES, PT, DENS_V1D,          &
      mesh3D, meshV1D                                  )

    if (ismaster) then
      call out_file_V1D%Write_var1D( DIAGVID_DENS, DENS_V1D, start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_RHOT, RHOT_V1D, start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_PT, PT_V1D, start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_W , W_V1D , start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_W_PRIM3 , W_PRIM3_V1D , start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_HEAT_EDDYFLX, HEAT_EDDYFLX_V1D, start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_MOMZ_EDDYFLX, MOMZ_EDDYFLX_V1D, start_time, end_time )
    end if
  end do
  !-
  call Finalize()

contains
!OCL SERIAL
  subroutine set_prgvars( MOMX_, MOMY_, MOMZ_, DRHOT_, PRES_,  &
    DDENS_, U_, V_, W_, PT_, DENS_hyd_, PRES_hyd_,             &
    lmesh, elem )

    use scale_const, only: &
      PRES00 => CONST_PRE00, &
      Rdry => CONST_Rdry,    &
      CPdry => CONST_CPdry,  &
      CVdry => CONST_CVdry
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: DRHOT_(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: PRES_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: U_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: V_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: W_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PT_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd_(elem%Np,lmesh%NeA)

    integer :: ke
    real(RP) :: DENS(elem%Np), RHOT_hyd(elem%Np)
    !--------------------------------------------------------

    !$omp parallel do private(DENS, RHOT_hyd)
    do ke=lmesh%NeS, lmesh%NeE
      DENS(:) = DENS_hyd_(:,ke) + DDENS_(:,ke)
      RHOT_hyd(:) = PRES00/Rdry * (PRES_hyd_(:,ke)/PRES00)**(CVdry/CPdry)

      MOMX_(:,ke) = DENS(:) * U_(:,ke)
      MOMY_(:,ke) = DENS(:) * V_(:,ke)
      MOMZ_(:,ke) = DENS(:) * W_(:,ke)
      DRHOT_(:,ke) = DENS(:) * PT_(:,ke) - RHOT_hyd(:)
      PRES_(:,ke) = PRES00 * ( Rdry / PRES00 * ( RHOT_hyd(:) + DRHOT_(:,ke) ) )**(CPdry/CVdry)
    end do

    return
  end subroutine set_prgvars

!OCL SERIAL
  subroutine analyze_lc_1( DENS_V1D_, RHOT_V1D_, MOMZ_V1D_,  &
    DDENS_, PT_, W_,                                         &
    DENS_hyd_,                                               &
    domID, lmesh, elem, lmeshH2D, elemH2D, lmeshV1D, elemV1D )

    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmeshH2D
    class(ElementBase2D), intent(in) :: elemH2D
    class(LocalMesh1D), intent(in) :: lmeshV1D
    class(ElementBase1D), intent(in) :: elemV1D
    real(RP), intent(inout) :: DENS_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(inout) :: RHOT_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(inout) :: MOMZ_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(in)    :: DDENS_   (elem%Np,lmesh%NeA)
    real(RP), intent(in)    :: PT_      (elem%Np,lmesh%NeA)
    real(RP), intent(in)    :: W_       (elem%Np,lmesh%NeA)
    real(RP), intent(in)    :: DENS_hyd_(elem%Np,lmesh%NeA)
    integer, intent(in) :: domID

    integer :: ke_x, ke_y, ke_z
    integer :: ke, ke2D
    integer :: p

    real(RP) :: DENS_lc(elem%Np)
    integer :: hSliceID(elem%Nnode_h1D**2)
    real(RP) :: int_w(elem%Nnode_h1D**2)
    !---------------------------------------------------------------------

    if (domID==1) then
      !$omp parallel do
      do ke_z=1, lmesh%NeZ
        DENS_V1D_(:,ke_z) = 0.0_RP
        RHOT_V1D_(:,ke_z) = 0.0_RP
        MOMZ_V1D_(:,ke_z) = 0.0_RP
      end do
      harea = 0.0_RP
    end if

    !$omp parallel do private( ke_x, ke_y, p, ke2D, ke, hSliceID, &
    !$omp int_w, DENS_lc )
    do ke_z=1, lmesh%NeZ
      do ke_y=1, lmesh%NeY
      do ke_x=1, lmesh%NeX
        ke2D = ke_x + (ke_y-1)*lmesh%NeX 
        ke = ke2D + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_lc(:) = DENS_hyd_(:,ke) + DDENS_(:,ke)

        do p=1, elem%Nnode_v
          hSliceID(:) = elem%Hslice(:,p)
          int_w(:) = elemH2D%IntWeight_lgl(:) * lmeshH2D%J(:,ke2D)

          DENS_V1D_(p,ke_z) = DENS_V1D_(p,ke_z) &
                            + sum( int_w(:) * DENS_lc(hSliceID(:)) )
          RHOT_V1D_(p,ke_z) = RHOT_V1D_(p,ke_z) &
                            + sum( int_w(:) * DENS_lc(hSliceID(:)) * PT_(hSliceID(:),ke) )
          MOMZ_V1D_(p,ke_z) = MOMZ_V1D_(p,ke_z) &
                            + sum( int_w(:) * DENS_lc(hSliceID(:)) * W_ (hSliceID(:),ke) )
        end do
      end do    
      end do  
    end do
    hSliceID(:) = elem%Hslice(:,1)
    do ke_y=1, lmesh%NeY
    do ke_x=1, lmesh%NeX
      ke2D = ke_x + (ke_y-1)*lmesh%NeX
      harea = harea + sum( elemH2D%IntWeight_lgl(:) * lmeshH2D%J(:,ke2D) )
    end do    
    end do  

    if ( domID == mesh3D%LOCAL_MESH_NUM ) then
      !$omp parallel do
      do ke_z=1, lmesh%NeZ
        DENS_V1D_(:,ke_z) = DENS_V1D_(:,ke_z) / harea
        RHOT_V1D_(:,ke_z) = RHOT_V1D_(:,ke_z) / harea
        MOMZ_V1D_(:,ke_z) = MOMZ_V1D_(:,ke_z) / harea
      end do
    end if

    return
  end subroutine  analyze_lc_1

!OCL SERIAL
  subroutine analyze_lc_2( &
    PT_V1D_, W_V1D_,                                         &
    W_PRIM3_V1D_, HEAT_EDDYFLX_V1D_, MOMZ_EDDYFLX_V1D_,      &
    DDENS_, PT_, W_,                                         &
    DENS_V1D_, RHOT_V1D_, MOMZ_V1D_,                         &    
    DENS_hyd_,                                               &
    domID, lmesh, elem, lmeshH2D, elemH2D, lmeshV1D, elemV1D )

    use scale_const, only: &
      CpDry => CONST_CPdry
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmeshH2D
    class(ElementBase2D), intent(in) :: elemH2D
    class(LocalMesh1D), intent(in) :: lmeshV1D
    class(ElementBase1D), intent(in) :: elemV1D
    real(RP), intent(inout) :: PT_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(inout) :: W_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(inout) :: W_PRIM3_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(inout) :: HEAT_EDDYFLX_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(inout) :: MOMZ_EDDYFLX_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(in) :: DDENS_   (elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PT_      (elem%Np,lmesh%NeA)
    real(RP), intent(in) :: W_       (elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(in) :: RHOT_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(in) :: MOMZ_V1D_(elemV1D%Np,lmeshV1D%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem%Np,lmesh%NeA)
    integer, intent(in) :: domID

    integer :: ke_x, ke_y, ke_z
    integer :: ke, ke2D
    integer :: p

    real(RP) :: DENS_lc(elem%Np)
    integer :: hSliceID(elem%Nnode_h1D**2)
    real(RP) :: int_w(elem%Nnode_h1D**2)

    real(RP) :: W_hm_V1D
    !---------------------------------------------------------------------

    if (domID==1) then
      !$omp parallel do
      do ke_z=1, lmesh%NeZ
        PT_V1D_(:,ke_z) = 0.0_RP
        W_V1D_ (:,ke_z) = 0.0_RP
        W_PRIM3_V1D_(:,ke_z) = 0.0_RP
        HEAT_EDDYFLX_V1D_(:,ke_z) = 0.0_RP
        MOMZ_EDDYFLX_V1D_(:,ke_z) = 0.0_RP
      end do
      harea = 0.0_RP
    end if

    !$omp parallel do private( ke_x, ke_y, p, ke2D, ke, hSliceID, &
    !$omp int_w, DENS_lc, W_hm_V1D )
    do ke_z=1, lmesh%NeZ
      do ke_y=1, lmesh%NeY
      do ke_x=1, lmesh%NeX
        ke2D = ke_x + (ke_y-1)*lmesh%NeX 
        ke = ke2D + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_lc(:) = DENS_hyd_(:,ke) + DDENS_(:,ke)

        do p=1, elem%Nnode_v
          hSliceID(:) = elem%Hslice(:,p)
          int_w(:) = elemH2D%IntWeight_lgl(:) * lmeshH2D%J(:,ke2D)

          W_hm_V1D = MOMZ_V1D_(p,ke_z) / DENS_V1D_(p,ke_z)

          HEAT_EDDYFLX_V1D_(p,ke_z) = HEAT_EDDYFLX_V1D_(p,ke_z) &
                            + sum( int_w(:) * DENS_lc(hSliceID(:))                              &
                              * ( W_ (hSliceID(:),ke) - W_hm_V1D                              ) &
                              * ( PT_(hSliceID(:),ke) - RHOT_V1D_(p,ke_z) / DENS_V1D_(p,ke_z) ) )
          
          MOMZ_EDDYFLX_V1D_(p,ke_z) = MOMZ_EDDYFLX_V1D_(p,ke_z) &
                            + sum( int_w(:) * DENS_lc(hSliceID(:))                              &
                              * ( W_ (hSliceID(:),ke) - W_hm_V1D                              ) &
                              * ( W_ (hSliceID(:),ke) - W_hm_V1D                              ) )
          
          W_PRIM3_V1D_(p,ke_z) = W_PRIM3_V1D_(p,ke_z) &
                              + sum( int_w(:) * DENS_lc(hSliceID(:))                            &
                              * ( W_ (hSliceID(:),ke) - W_hm_V1D                              ) &
                              * ( W_ (hSliceID(:),ke) - W_hm_V1D                              ) &
                              * ( W_ (hSliceID(:),ke) - W_hm_V1D                              ) )
  
        end do
      end do    
      end do  
    end do
    hSliceID(:) = elem%Hslice(:,1)
    do ke_y=1, lmesh%NeY
    do ke_x=1, lmesh%NeX
      ke2D = ke_x + (ke_y-1)*lmesh%NeX
      harea = harea + sum( elemH2D%IntWeight_lgl(:) * lmeshH2D%J(:,ke2D) )
    end do    
    end do  

    if ( domID == mesh3D%LOCAL_MESH_NUM ) then
      !$omp parallel do
      do ke_z=1, lmesh%NeZ
        PT_V1D_(:,ke_z) = RHOT_V1D_(:,ke_z) / DENS_V1D_(:,ke_z)
        W_V1D_ (:,ke_z) = MOMZ_V1D_(:,ke_z) / DENS_V1D_(:,ke_z)
        HEAT_EDDYFLX_V1D_(:,ke_z) = CpDry * HEAT_EDDYFLX_V1D_(:,ke_z) / harea
        MOMZ_EDDYFLX_V1D_(:,ke_z) = MOMZ_EDDYFLX_V1D_(:,ke_z) / DENS_V1D_(:,ke_z) / harea
        W_PRIM3_V1D_(:,ke_z) = W_PRIM3_V1D_(:,ke_z) / DENS_V1D_(:,ke_z) / harea
      end do
    end if

    return
  end subroutine  analyze_lc_2

  !----------------------------------------------------------------------------

!OCL SERIAL
  subroutine initialize()
    use scale_prc, only: &
      PRC_MPIstart,        &
      PRC_SINGLECOM_setup, &
      PRC_ERRHANDLER_setup
    use scale_const, only: &
      CONST_setup
    use scale_file, only: &
      FILE_setup    
    use mod_diag_tb, only: &
      diag_tb_Init
    use mod_common, only: &
      common_Init        
    use scale_mesh_base1d, only: &
      DIMTYPEID_ZT => MeshBase1D_DIMTYPEID_XT
    implicit none

    integer :: ierr
    integer :: comm                     ! communicator   (execution)
    logical :: fileexist
    character(len=H_LONG) :: cnf_fname  ! config file for launcher

    character(len=H_LONG) :: out_filebase_tb
    character(len=H_LONG) :: out_filebase_V1D
  
    !-
    namelist / PARAM_PBLTB_ANALYSIS / &
      in_filebase,                    &
      in_bs_filebase,                 &
      out_filebase_V1D,               &
      out_filebase_tb,                &
      start_time0,                    &
      output_tintrv,                  &
      num_step

    !-
    integer :: k
    logical :: is_spec_FZ    
    integer, parameter :: FZ_nmax = 1000
    real(RP) :: FZ(FZ_nmax)
  
    integer :: NprcX = 24
    integer :: NprcY = 8
    integer :: NeX   = 1
    integer :: NeY   = 3
    integer :: NeZ   = 12
    integer, parameter :: NLocalMeshPerPrc = 1
    
    integer :: PolyOrder_h = 7
    integer :: PolyOrder_v = 7
    logical  :: LumpedMassMatFlag = .false.

    real(RP) :: dom_xmin = 0.0_RP
    real(RP) :: dom_xmax = 9.6E3_RP
    real(RP) :: dom_ymin = 0.0_RP
    real(RP) :: dom_ymax = 9.6E3_RP
    real(RP) :: dom_zmin = 0.0_RP
    real(RP) :: dom_zmax = 3.0E3_RP    

    namelist / PARAM_ATMOS_MESH / &
      dom_xmin, dom_xmax,                          &
      dom_ymin, dom_ymax,                          &
      dom_zmin, dom_zmax,                          &
      FZ,                                          &
      NeX, NeY, NeZ,                               &
      PolyOrder_h, PolyOrder_v, LumpedMassMatFlag, &
      NprcX, NprcY
        
    !-----------------------------------------------------------------------

      ! start MPI
    call PRC_MPIstart( comm ) ! [OUT]

    ! setup MPI communicator
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
                              nprocs,  & ! [OUT]
                              myrank,  & ! [OUT]
                              ismaster ) ! [OUT]
  
    call PRC_ERRHANDLER_setup( use_fpm = .false., & ! [IN]
                               master  = .false.  ) ! [IN]

    cnf_fname = IO_ARG_getfname( ismaster )

    ! setup standard I/O
    call IO_setup( "PBLTB_ANALYSIS", cnf_fname )

    ! setup Log
    call IO_LOG_setup( myrank, ismaster )
    ! setup constants
    call CONST_setup
    ! setup fie I/O
    call FILE_setup( myrank )
    
    LOG_NEWLINE
    LOG_INFO("PBLTB_ANALYSIS",*) 'Setup'

    !-
    in_filebase = "./history"
    out_filebase_tb  = "./analysis/history_tb"
    out_filebase_V1D = "./analysis/history_v1D"

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_PBLTB_ANALYSIS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("PBLTB_ANALYSIS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("PBLTB_ANALYSIS_setup",*) 'Not appropriate names in namelist PARAM_PBLTB_ANALYSIS. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_PBLTB_ANALYSIS)

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_MESH,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("PBLTB_ANALYSIS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("PBLTB_ANALYSIS_setup",*) 'Not appropriate names in namelist PARAM_ATM_MESH. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_MESH)

    !- Mesh
    LOG_INFO("PBLTB_ANALYSIS",*) 'Setup mesh..'

    call refElem3D%Init( PolyOrder_h, PolyOrder_v, LumpedMassMatFlag )
    call refElemV1D%Init( PolyOrder_v, LumpedMassMatFlag )

    is_spec_FZ = .true.
    do k=1, NeZ+1
      if (FZ(k) < 0.0_RP) then
        is_spec_FZ = .false.
      end if
    end do

    if (is_spec_FZ) then
      call mesh3D%Init( NprcX*NeX, NprcY*NeY, NeZ, &
        dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,         &
        .true., .true., .false., refElem3D, NLocalMeshPerPrc, NprcX, NprcY, &
        FZ=FZ )
    else
      call mesh3D%Init( NprcX*NeX, NprcY*NeY, NeZ, &
        dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,         &
        .true., .true., .false., refElem3D, NLocalMeshPerPrc, NprcX, NprcY  )
    end if
    call mesh3D%Generate()

    call meshV1D%Init( NeZ, dom_zmin, dom_zmax, refElemV1D, 1, &
                       nproc=1, myrank=0, FX=FZ(1:NeZ+1)       )
    call meshV1D%Generate()

    ! Data
    LOG_INFO("PBLTB_ANALYSIS",*) 'Setup variables..'

    call U%Init("U", "m/s", mesh3D)
    call V%Init("V", "m/s", mesh3D)
    call W%Init("W", "m/s", mesh3D)
    call PT%Init("PT", "K", mesh3D)
    call PRES%Init("PRES", "Pa", mesh3D)

    call DDENS%Init("DDENS", "kg/m3", mesh3D)
    call MOMX%Init("MOMX", "kg/m3.m/s", mesh3D)
    call MOMY%Init("MOMY", "kg/m3.m/s", mesh3D)
    call MOMZ%Init("MOMZ", "kg/m3.m/s", mesh3D)
    call DRHOT%Init("DRHOT", "K", mesh3D)

    call PRES_hyd%Init("PRES_hyd", "Pa", mesh3D)
    call DENS_hyd%Init("DENS_hyd", "kg/m3", mesh3D)

    call DENS_V1D%Init("DENS", "kg/m3", meshV1D)
    call RHOT_V1D%Init("RHOT", "kg/m3.K", meshV1D)
    call MOMZ_V1D%Init("MOMZ", "kg/m3.m/s", meshV1D)
    call PT_V1D%Init("PT", "K", meshV1D)
    call W_V1D %Init("W", "m/s", meshV1D)
    call W_PRIM3_V1D%Init("W_PRIM3", "m/s", meshV1D)
    call HEAT_EDDYFLX_V1D%Init("HEAT_EDDYFLX", "J.m-2.s-1", meshV1D)
    call MOMZ_EDDYFLX_V1D%Init("MOMZ_EDDYFLX", "m/s.kg.m-2.s-1", meshV1D)

    ! Input & Output   

    LOG_INFO("PBLTB_ANALYSIS",*) 'Open input file ..'

    call in_file%Init(5, mesh3D=mesh3D)
    call in_file%Open(in_filebase)

    call in_bs_file%Init(2, mesh3D=mesh3D)
    call in_bs_file%Open(in_bs_filebase)

    if (ismaster) then
      LOG_INFO("PBLTB_ANALYSIS",*) 'Setup output file for vertical 1D data..'

      call out_file_V1D%Init(9, mesh1D=meshV1D)
      call out_file_V1D%Create( out_filebase_V1D, 'PBL turbulence analysis', &
        dtype, fileexist, myrank=myrank )
      call out_file_V1D%Def_Var(DENS_V1D, "horizontal averaged density", DIAGVID_DENS, DIMTYPEID_ZT, dtype, &
        timeinv=output_tintrv )
      call out_file_V1D%Def_Var(RHOT_V1D, "horizontal averaged density*potential temperature", DIAGVID_RHOT, DIMTYPEID_ZT, dtype, &
        timeinv=output_tintrv )
      call out_file_V1D%Def_Var(PT_V1D, "horizontal averaged PT", DIAGVID_PT, DIMTYPEID_ZT, dtype, &
        timeinv=output_tintrv )
      call out_file_V1D%Def_Var( W_V1D,  "horizontal averaged W",  DIAGVID_W, DIMTYPEID_ZT, dtype, &
        timeinv=output_tintrv )
      call out_file_V1D%Def_Var( W_PRIM3_V1D,  "horizontal averaged W_PRIM3",  DIAGVID_W_PRIM3, DIMTYPEID_ZT, dtype, &
        timeinv=output_tintrv )
      call out_file_V1D%Def_Var(HEAT_EDDYFLX_V1D, "horizontal averaged eddy heat flux", DIAGVID_HEAT_EDDYFLX, DIMTYPEID_ZT, dtype, &
        timeinv=output_tintrv )
      call out_file_V1D%Def_Var(MOMZ_EDDYFLX_V1D, "horizontal averaged eddy momentum flux in z-direction", &
        DIAGVID_MOMZ_EDDYFLX, DIMTYPEID_ZT, dtype, timeinv=output_tintrv )      
      call out_file_V1D%End_def()
    end if

    ! Communication
    call vars_comm%Init( 9, 0, mesh3D )
    comm_vars_ptr(1)%field3d => DDENS
    comm_vars_ptr(2)%field3d => MOMX
    comm_vars_ptr(3)%field3d => MOMY
    comm_vars_ptr(4)%field3d => MOMZ
    comm_vars_ptr(5)%field3d => DRHOT
    comm_vars_ptr(6)%field3d => DENS_hyd
    comm_vars_ptr(7)%field3d => PRES_hyd
    comm_vars_ptr(8)%field3d => PT
    comm_vars_ptr(9)%field3d => PRES

    ! Diagnose
    call common_Init( mesh3D )
    call diag_tb_Init( mesh3D, meshV1D, &
      out_filebase_tb, dtype, output_tintrv, &
      myrank, ismaster, NLocalMeshPerPrc )


    LOG_INFO("PBLTB_ANALYSIS",*) 'Setup has been finished.'

    return
  end subroutine initialize

!OCL SERIAL
  subroutine finalize()
    use scale_prc, only: &
      PRC_mpibarrier, &
      PRC_MPIfinish
    use scale_file, only: &
      FILE_close_all

    use mod_diag_tb, only: &
      diag_tb_Final
    use mod_common, only: &
      common_Final
    implicit none
    !---------------------------------------------

    ! Diagnose
    call diag_tb_Final()
    call common_Final()


    ! IO
    call in_file%Close()
    call in_file%Final()

    call in_bs_file%Close()
    call in_bs_file%Final()

    if (ismaster) then
      call out_file_V1D%Close()
      call out_file_V1D%Final()
    end if

    ! Communication
    call vars_comm%Final()

    ! Field
    call U%Final()
    call V%Final()
    call W%Final()
    call PT%Final()
    call PRES%Final()

    call DDENS%Final()
    call MOMX%Final()
    call MOMY%Final()
    call MOMZ%Final()
    call DRHOT%Final()

    call PRES_hyd%Final()
    call DENS_hyd%Final()

    call DENS_V1D%Final()
    call RHOT_V1D%Final()
    call MOMZ_V1D%Final()
    call PT_V1D%Final()
    call W_V1D%Final()
    call W_PRIM3_V1D%Final()
    call HEAT_EDDYFLX_V1D%Final()
    call MOMZ_EDDYFLX_V1D%Final()

    ! Mesh
    call mesh3D%Final()
    call meshV1D%Final()
    ! Element
    call refElem3D%Final()
    call refElemV1D%Final()

    ! stop MPI
    call PRC_mpibarrier
    call PRC_MPIfinish 

    return
  end subroutine finalize

end program prg_pbltb_analysis