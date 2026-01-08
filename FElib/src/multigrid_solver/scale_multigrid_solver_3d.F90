!-------------------------------------------------------------------------------
!> module FElib / Multigrid / Solver 3D
!!
!! @par Description
!!      Manage a multigrid solver for 3D domain
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_multigrid_solver_3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_sparsemat, only: SparseMat

  use scale_element_base, only: ElementBase3D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement

  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D

  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D

  use scale_mesh_hierarchy_base, only: &
    pMG_FINEST_LEVEL => MESH_HIERARCHY_pMG_FINEST_LEVEL, &
    hMG_FINEST_LEVEL => MESH_HIERARCHY_hMG_FINEST_LEVEL, &
    MESH_HIERARCHY_TYPE_pMG,                             &
    MESH_HIERARCHY_TYPE_hMG
  use scale_mesh_hierarchy_3d, only: &
    MeshHierarchy3D,      &
    MeshHierarchyLevel3D, &
    MeshHierarchyLocalMGData3D
  
  use scale_multigrid_solver_base, only: MultiGridSolverBase
  use scale_multigrid_smoother_base, only: &
    MGSmootherBase3D, MGSmoother_PRE_ID, MGSmoother_POST_ID
  use scale_multigrid_fieldset_base, only: MGFieldSet3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  !> Derived type for multigrid solver in 3D domain
  type, extends(MultiGridSolverBase), public :: MultiGridSolver3D
    type(MeshHierarchy3D), pointer :: mesh_hierarchy_ptr  !< Pointer to an object to manage 3D mesh hierarchy
    class(MGSmootherBase3D), pointer :: mg_smoother_ptr   !< Pointer to an object to manage 3D MG smoother

    type(MGFieldSet3D), allocatable :: fields_h(:)        !< Array of objects to manage MG field sets for h-MG
    type(MGFieldSet3D), allocatable :: fields_p(:)        !< Array of objects to manage MG field sets for p-MG
  contains
    procedure :: Init => MultiGridSolver3D_Init
    procedure :: Final => MultiGridSolver3D_Final
    procedure :: Solve => MultiGridSolver3D_solve
    procedure :: do_Vcycle => MultiGridSolver3D_do_Vcycle
    !-
    procedure :: do_hMG_Vcycle => MultiGridSolver3D_do_hMG_Vcycle
    procedure :: Operate_pMG_restriction => MultiGridSolver3D_Operate_pMG_restriction
    procedure :: Operate_pMG_correction => MultiGridSolver3D_Operate_pMG_correction
    procedure :: Operate_hMG_restriction => MultiGridSolver3D_Operate_hMG_restriction
    procedure :: Operate_hMG_correction => MultiGridSolver3D_Operate_hMG_correction
  end type MultiGridSolver3D

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
contains
  !> Initialize an object for multigrid solver in 3D domain
!OCL SERIAL
  subroutine MultiGridSolver3D_Init( this, &
    mesh_hierarchy, mg_smoother,    &
    aux_var_num, aux_vec_num        )
    use scale_multigrid_solver_base, only: MultiGridSolverBase_Init
    implicit none
    class(MultiGridSolver3D), intent(inout) :: this
    class(MeshHierarchy3D), intent(in), target :: mesh_hierarchy
    class(MGSmootherBase3D), intent(in), target :: mg_smoother
    integer, intent(in) :: aux_var_num
    integer, intent(in) :: aux_vec_num

    integer :: lev_p
    integer :: lev_h
    !-------------------------------------------------------------

    call MultiGridSolverBase_Init( this, mesh_hierarchy )

    this%mesh_hierarchy_ptr => mesh_hierarchy
    this%mg_smoother_ptr => mg_smoother

    !- Prepare p-MG field sets
    allocate( this%fields_p( mesh_hierarchy%NUM_pMG_LEVEL ) )
    do lev_p=1, mesh_hierarchy%NUM_pMG_LEVEL
      call this%fields_p(lev_p)%Init( mesh_hierarchy%p_mesh_list(lev_p)%ptr, &
        aux_var_num, aux_vec_num, lev_p )
    end do

    !- Prepare h-MG field sets
    allocate( this%fields_h( mesh_hierarchy%NUM_hMG_LEVEL ) )
    do lev_h=1, mesh_hierarchy%NUM_hMG_LEVEL
      call this%fields_h(lev_h)%Init( mesh_hierarchy%h_mesh_list(lev_h)%ptr, &
        aux_var_num, aux_vec_num, lev_h )
    end do
    return
  end subroutine MultiGridSolver3D_Init

  !> Finalize an object for multigrid solver in 3D domain
!OCL SERIAL
  subroutine MultiGridSolver3D_Final(this)
    use scale_multigrid_solver_base, only: MultiGridSolverBase_Final
    implicit none
    class(MultiGridSolver3D), intent(inout) :: this

    integer :: lev_p
    integer :: lev_h
    !-------------------------------------------------------------

    call MultiGridSolverBase_Final( this )

    if (allocated(this%fields_p) ) then
      do lev_p=1, this%mesh_hierarchy_ptr%NUM_pMG_LEVEL
        call this%fields_p(lev_p)%Final()
      end do
    end if
    deallocate( this%fields_p )

    if (allocated(this%fields_h) ) then
      do lev_h=1, this%mesh_hierarchy_ptr%NUM_hMG_LEVEL
        call this%fields_h(lev_h)%Final()
      end do
    end if
    deallocate( this%fields_h )

    return
  end subroutine MultiGridSolver3D_Final

  !> Solve a linear system using multigrid method in 3D domain
!OCL SERIAL
  subroutine MultiGridSolver3D_solve(this, q, &
    f )
    implicit none
    class(MultiGridSolver3D), intent(inout) :: this

    class(MeshField3D), intent(inout), target :: q
    class(MeshField3D), intent(inout) :: f

    integer :: vcyc_itr

    class(LocalMesh3D), pointer :: lmesh
    integer :: ldomID
    integer :: ke
    !-------------------------------------------------------------

    do ldomID=1, q%mesh%LOCAL_MESH_NUM
      lmesh => q%mesh%lcmesh_list(ldomID)
      do ke=lmesh%NeS, lmesh%NeE
        this%fields_p(pMG_FINEST_LEVEL)%dq%local(ldomID)%val(:,ke) = q%local(ldomID)%val(:,ke)
      end do
    end do

    this%current_p_lev = pMG_FINEST_LEVEL-1
    this%current_h_lev = hMG_FINEST_LEVEL-1

    !-
    do vcyc_itr=1, this%vcyc_num_max
      LOG_INFO("MultiGridSolver3D_solve",*) "V-cycle iteration:", vcyc_itr

      call this%do_Vcycle( pMG_FINEST_LEVEL, f, vcyc_itr )
      if ( this%Is_converged(this%mg_smoother_ptr) ) then
        LOG_INFO("MultiGridSolver3D_solve",*) "V-cycle converged: vcyc_itr=", vcyc_itr
        exit
      end if
    end do

    !-
    do ldomID=1, q%mesh%LOCAL_MESH_NUM
      lmesh => q%mesh%lcmesh_list(ldomID)
      do ke=lmesh%NeS, lmesh%NeE
        q%local(ldomID)%val(:,ke) = this%fields_p(pMG_FINEST_LEVEL)%dq%local(ldomID)%val(:,ke)
      end do
    end do
    return
  end subroutine MultiGridSolver3D_solve

  !> Do a V-cycle in 3D domain
!OCL SERIAL
  recursive subroutine MultiGridSolver3D_do_Vcycle( this, &
    mg_level, f_in, vcyc_itr )
    implicit none
    class(MultiGridSolver3D), intent(inout), target :: this
    integer, intent(in) :: mg_level
    type(MeshField3D), intent(inout) :: f_in
    integer, intent(in) :: vcyc_itr

    logical :: invoke_hMG

    class(MGFieldSet3D), pointer :: fs_p
    class(MeshHierarchy3D), pointer :: mesh_hierarchy
    !---------------------------------------------------------------------

    this%current_p_lev = this%current_p_lev + 1    

    mesh_hierarchy => this%mesh_hierarchy_ptr
    fs_p => this%fields_p(mg_level)

    LOG_INFO("MultiGridSolver3D_do_Vcycle",*) "Start: p_level=", this%current_p_lev

    !- Pre-relaxation
    call this%mg_smoother_ptr%Do_smoothing( fs_p%dq, fs_p%res,                        &
      f_in, fs_p%aux_var, fs_p%var_comm_ptr, fs_p%aux_comm_ptr,                       &
      fs_p%Dx, fs_p%Dy, fs_p%Dz, fs_p%Lift, mesh_hierarchy%p_mesh_list(mg_level)%ptr, & 
      this%current_p_lev, this%current_h_lev, MGSmoother_PRE_ID )

    if ( vcyc_itr == 1 .and. this%current_p_lev == pMG_FINEST_LEVEL ) then
      call this%mg_smoother_ptr%Get_initial_residual_statistics( &
        this%history_residual_l2_initial, this%history_residual_max_initial )
    end if
    call this%mg_smoother_ptr%Output_residual_history()

    ! if ( vcyc_itr == 1 ) call Output_tmp_data(this, fs_p, f_in, vcyc_itr, "_pre")

    if ( mg_level == mesh_hierarchy%NUM_pMG_LEVEL .and. mesh_hierarchy%NUM_hMG_LEVEL == 0 ) then
      ! if ( vcyc_itr == 1 ) call Output_tmp_data(this, fs_p, f_in, vcyc_itr, "_post")
      LOG_INFO("MultiGridSolver3D_do_Vcycle",*) "End: p_level=", this%current_p_lev
      this%current_p_lev = this%current_p_lev - 1 
      return
    end if

    !-
    invoke_hMG = ( mg_level+1 >= mesh_hierarchy%NUM_pMG_LEVEL .and. mesh_hierarchy%NUM_hMG_LEVEL > 0 )

    if ( invoke_hMG ) then
      !- Restriction
      call this%Operate_pMG_restriction( this%fields_h(hMG_FINEST_LEVEL)%f, &
        this%fields_p(mg_level)%res, mg_level )

      !- Advance node in the V-cycle
      call this%do_hMG_Vcycle( hMG_FINEST_LEVEL, this%fields_h(hMG_FINEST_LEVEL)%f )
      
      !- Correction
      call this%Operate_pMG_correction( this%fields_p(mg_level)%dq, &
        this%fields_h(hMG_FINEST_LEVEL)%dq, mg_level )
    else
      !- Restriction
      call this%Operate_pMG_restriction( this%fields_p(mg_level+1)%f, &
        this%fields_p(mg_level)%res, mg_level )
      
      !- Advance node in the V-cycle
      call this%do_Vcycle( mg_level+1, this%fields_p(mg_level+1)%f, vcyc_itr )

      !- Correction
      call this%Operate_pMG_correction( this%fields_p(mg_level)%dq, &
        this%fields_p(mg_level+1)%dq, mg_level )
    end if

    !- Post-relaxation
    call this%mg_smoother_ptr%Do_smoothing( fs_p%dq, fs_p%res,                        &
      f_in, fs_p%aux_var, fs_p%var_comm_ptr, fs_p%aux_comm_ptr,                       &
      fs_p%Dx, fs_p%Dy, fs_p%Dz, fs_p%Lift, mesh_hierarchy%p_mesh_list(mg_level)%ptr, &
      this%current_p_lev, this%current_h_lev, MGSmoother_POST_ID )

    ! if ( vcyc_itr == 1 ) call Output_tmp_data(this, fs_p, f_in, vcyc_itr, "_post")

    call this%mg_smoother_ptr%Output_residual_history()
    LOG_INFO("MultiGridSolver3D_do_Vcycle",*) "End: p_level=", this%current_p_lev

    this%current_p_lev = this%current_p_lev - 1
    return
  end subroutine MultiGridSolver3D_do_Vcycle

!OCL SERIAL
  recursive subroutine MultiGridSolver3D_do_hMG_Vcycle( this, mg_level, f_in )
    implicit none
    class(MultiGridSolver3D), intent(inout), target :: this
    integer, intent(in) :: mg_level
    type(MeshField3D), intent(inout) :: f_in

    class(MeshHierarchy3D), pointer :: mesh_hierarchy
    class(MGFieldSet3D), pointer :: fs_h
    !----------------------------------------------

    mesh_hierarchy => this%mesh_hierarchy_ptr
    fs_h => this%fields_h(mg_level)

    this%current_h_lev = this%current_h_lev + 1    

    LOG_INFO("MultiGridSolver3D_do_hMG_Vcycle",*) "Start: h_level=", this%current_h_lev

    if ( mg_level == mesh_hierarchy%NUM_hMG_LEVEL ) then
      ! It should be replaced by a direct solver in the future
      call this%mg_smoother_ptr%Do_smoothing( fs_h%dq, fs_h%res,                        &
        f_in, fs_h%aux_var, fs_h%var_comm_ptr, fs_h%aux_comm_ptr,                       &
        fs_h%Dx, fs_h%Dy, fs_h%Dz, fs_h%Lift, mesh_hierarchy%h_mesh_list(mg_level)%ptr, &
        this%current_p_lev, this%current_h_lev, MGSmoother_PRE_ID                       )

      call this%mg_smoother_ptr%Output_residual_history()

      LOG_INFO("MultiGridSolver3D_do_hMG_Vcycle",*) "End: h_level=", this%current_h_lev
      this%current_h_lev = this%current_h_lev - 1
      return
    end if

    !- Pre-relaxation
    call this%mg_smoother_ptr%Do_smoothing( fs_h%dq, fs_h%res,                        &
      f_in, fs_h%aux_var, fs_h%var_comm_ptr, fs_h%aux_comm_ptr,                       &
      fs_h%Dx, fs_h%Dy, fs_h%Dz, fs_h%Lift, mesh_hierarchy%h_mesh_list(mg_level)%ptr, &
      this%current_p_lev, this%current_h_lev, MGSmoother_PRE_ID                       )

    !- Restriction
    call this%Operate_hMG_restriction( this%fields_h(mg_level+1)%f, &
      this%fields_h(mg_level)%res, mg_level )

    !- Advance node in the V-cycle
    call this%do_hMG_Vcycle( mg_level+1, this%fields_h(mg_level+1)%f )

    !- Correction
    call this%Operate_hMG_correction( this%fields_h(mg_level)%dq, &
      this%fields_h(mg_level+1)%dq, mg_level )

    !- Post-relaxation
    call this%mg_smoother_ptr%Do_smoothing( fs_h%dq, fs_h%res,                        &
      f_in, fs_h%aux_var, fs_h%var_comm_ptr, fs_h%aux_comm_ptr,                       &
      fs_h%Dx, fs_h%Dy, fs_h%Dz, fs_h%Lift, mesh_hierarchy%h_mesh_list(mg_level)%ptr, &
      this%current_p_lev, this%current_h_lev, MGSmoother_POST_ID                      )

    LOG_INFO("MultiGridSolver3D_do_hMG_Vcycle",*) "End: h_level=", this%current_h_lev
    this%current_h_lev = this%current_h_lev - 1
    return
  end subroutine MultiGridSolver3D_do_hMG_Vcycle

!--
!OCL SERIAL
  subroutine MultiGridSolver3D_Operate_pMG_restriction( this, res_c, &
    res, p_lev )
    implicit none
    class(MultiGridSolver3D), intent(in), target :: this
    class(MeshField3D), intent(inout) :: res_c
    class(MeshField3D), intent(in) :: res
    integer, intent(in) :: p_lev

    class(MeshHierarchy3D), pointer :: hierarchy
    integer :: ldomID
    class(MeshBase3D), pointer :: mesh3D
    !------------------------------------------------

    if ( this%mesh_hierarchy_ptr%NUM_pMG_LEVEL <  p_lev+1 ) then
      call PRC_abort()
    end if

    hierarchy => this%mesh_hierarchy_ptr
    mesh3D => hierarchy%p_mesh_list(p_lev)%ptr
    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      call MultiGridSolver3D_pMG_operation( res_c%local(ldomID)%val,     &
        res%local(ldomID)%val,                                           &
        hierarchy%elem3D_list(p_lev), hierarchy%elem3D_list(p_lev+1),    &
        mesh3D%lcmesh_list(ldomID), hierarchy%p_level(p_lev)%pMat1D_f2c, &
        .false.  )
    end do
    return
  end subroutine MultiGridSolver3D_Operate_pMG_restriction

!OCL SERIAL
  subroutine MultiGridSolver3D_Operate_pMG_correction( this, dq, &
    cor_c, p_lev )
    implicit none
    class(MultiGridSolver3D), intent(in), target :: this
    class(MeshField3D), intent(inout) :: dq
    class(MeshField3D), intent(in) :: cor_c
    integer, intent(in) :: p_lev

    class(MeshHierarchy3D), pointer :: hierarchy

    integer :: ldomID
    class(MeshBase3D), pointer :: mesh3D
    !------------------------------------------------

    hierarchy => this%mesh_hierarchy_ptr
    mesh3D => hierarchy%p_mesh_list(p_lev)%ptr

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      call MultiGridSolver3D_pMG_operation( dq%local(ldomID)%val,        &
        cor_c%local(ldomID)%val,                                         &
        hierarchy%elem3D_list(p_lev+1), hierarchy%elem3D_list(p_lev),    &
        mesh3D%lcmesh_list(ldomID), hierarchy%p_level(p_lev)%pMat1D_c2f, &
        .true. )
    end do
    return
  end subroutine MultiGridSolver3D_Operate_pMG_correction

!OCL SERIAL
  subroutine MultiGridSolver3D_Operate_hMG_restriction( this, res_c, &
    res, h_lev )
    use scale_meshfield_base, only: MeshField3D
    implicit none
    class(MultiGridSolver3D), intent(in), target :: this
    class(MeshField3D), intent(inout) :: res_c
    class(MeshField3D), intent(in) :: res
    integer, intent(in) :: h_lev

    class(MeshHierarchy3D), pointer :: hierarchy

    integer :: ldomID
    class(MeshBase3D), pointer :: mesh3D
    class(MeshBase3D), pointer :: mesh3D_c
    !------------------------------------------------

    hierarchy => this%mesh_hierarchy_ptr
    mesh3D => hierarchy%h_mesh_list(h_lev)%ptr
    mesh3D_c => hierarchy%h_mesh_list(h_lev+1)%ptr

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      call MultiGridSolver3D_hMG_restriction_core( res_c%local(ldomID)%val,          &
        res%local(ldomID)%val,                                                       &
        hierarchy%h_level(h_lev)%mg_local(ldomID), mesh3D%lcmesh_list(ldomID),       &
        hierarchy%elem3D_list(hierarchy%NUM_pMG_LEVEL), mesh3D_c%lcmesh_list(ldomID) )
    end do
    return
  end subroutine MultiGridSolver3D_Operate_hMG_restriction

!OCL SERIAL
  subroutine MultiGridSolver3D_Operate_hMG_correction( this, dq, &
    cor_c, h_lev )
    use scale_meshfield_base, only: MeshField3D
    implicit none
    class(MultiGridSolver3D), intent(in), target :: this
    class(MeshField3D), intent(inout) :: dq
    class(MeshField3D), intent(in) :: cor_c
    integer, intent(in) :: h_lev

    integer :: ldomID
    class(MeshBase3D), pointer :: mesh3D
    class(MeshHierarchy3D), pointer :: hierarchy
    !------------------------------------------------

    hierarchy => this%mesh_hierarchy_ptr
    mesh3D => hierarchy%h_mesh_list(h_lev)%ptr

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      call MultiGridSolver3D_hMG_correction_core( dq%local(ldomID)%val,     &
        cor_c%local(ldomID)%val, hierarchy%h_level(h_lev)%mg_local(ldomID), &
        mesh3D%lcmesh_list(ldomID), mesh3D%refElem3D                        )
    end do
    return
  end subroutine MultiGridSolver3D_Operate_hMG_correction

!-- private --------------------------------------------------------------

!OCL SERIAL
  subroutine Output_tmp_data( this, fs, fin, vcyc_itr, postfix )
    use scale_prc, only: PRC_myrank
    use scale_file_base_meshfield, only: FILE_base_meshfield
    use scale_mesh_base3d, only: MESHBASE3D_DIMTYPEID_XYZ
    implicit none
    class(MultiGridSolver3D), intent(in) :: this
    class(MGFieldSet3D), intent(in), target :: fs
    class(MeshField3D), intent(in) :: fin
    integer, intent(in) :: vcyc_itr
    character(len=*), intent(in) :: postfix

    type(FILE_base_meshfield) :: file
    class(MeshBase3D), pointer :: mesh
    character(len=H_MID) :: fname
    logical :: fileexisted

    integer, parameter :: DQ_VID = 1
    integer, parameter :: RES_VID = 2
    integer, parameter :: FIN_VID = 3
    !-----------------------------------------------------

    select type (mesh => fs%dq%mesh)
    type is (MeshCubeDom3D)
      call file%Init(2, mesh3D=mesh)
    end select

    write(fname,'(a,i2.2,a,i2.2,a)') "tmp_plev", this%current_p_lev, "_vcyc", vcyc_itr, trim(postfix)
    call file%Create( fname, "MG", "REAL8", fileexisted, myrank=PRC_myrank )
    call file%Def_Var( "dq", "1", "dq", DQ_VID, MESHBASE3D_DIMTYPEID_XYZ, "REAL8")
    call file%Def_Var( "res", "1", "res", RES_VID, MESHBASE3D_DIMTYPEID_XYZ, "REAL8")
    call file%Def_Var( "fin", "1", "fin", FIN_VID, MESHBASE3D_DIMTYPEID_XYZ, "REAL8")
    call file%End_def()

    call file%Write_var3D( DQ_VID, fs%dq, 0.0_RP, 1.0_RP)
    call file%Write_var3D( RES_VID, fs%res, 0.0_RP, 1.0_RP)
    call file%Write_var3D( FIN_VID, fin, 0.0_RP, 1.0_RP)
    call file%Close()
    call file%Final()
    return
  end subroutine Output_tmp_data

!OCL SERIAL
  subroutine MultiGridSolver3D_hMG_restriction_core( res_c_lc,      &
      res_lc, mg_local, lmesh, elem, lmesh_c )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh3D), intent(in) :: lmesh_c
    real(RP), intent(out) :: res_c_lc(elem%Np,lmesh_c%NeA)
    real(RP), intent(in) :: res_lc(elem%Np,lmesh%NeA)
    class(MeshHierarchyLocalMGData3D), intent(in) :: mg_local

    integer :: k, ke, ke_c
    real(RP) :: Ic2fT_lc(8,8)
    real(RP) :: tmp_c(elem%Np)
    real(RP) :: tmp2(elem%Np)
    !-----------------------------------------

    !$omp parallel private(ke, ke_c, tmp2, tmp_c, Ic2fT_lc)
    !$omp do
    do ke_c=lmesh_c%NeS, lmesh_c%NeE
      
      tmp_c(:) = 0.0_RP
      do k=1, 4
        ke = mg_local%If2c_emap(k,ke_c)

        tmp2(:) = matmul( elem%M, res_lc(:,ke) )
        Ic2fT_lc(:,:) = transpose( mg_local%Ic2f(:,:,ke) )
        tmp_c(:) = tmp_c(:) + matmul( Ic2fT_lc, tmp2(:) )
      end do

      res_c_lc(:,ke_c) = matmul(lmesh_c%refElem3D%invM, tmp_c(:))  * 0.25_RP        
    end do
    !$omp end parallel
    return
  end subroutine MultiGridSolver3D_hMG_restriction_core

!OCL SERIAL
  subroutine MultiGridSolver3D_hMG_correction_core( dq_lc,  &
    dq_c_lc, mg_local, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: dq_lc(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dq_c_lc(elem%Np,lmesh%NeA)
    class(MeshHierarchyLocalMGData3D), intent(in) :: mg_local

    integer :: k, ke, ke_c
    !-----------------------------------------

    !$omp parallel do private(ke, ke_c)
    do ke=lmesh%NeS, lmesh%NeE
      ke_c = mg_local%Ic2f_emap(ke)
      dq_lc(:,ke) = dq_lc(:,ke) + &
            matmul( mg_local%Ic2f(:,:,ke), dq_c_lc(:,ke_c) )
    end do
    return
  end subroutine MultiGridSolver3D_hMG_correction_core

!OCL SERIAL
  subroutine MultiGridSolver3D_pMG_operation( q_o, &
    q_i, elem3D_i, elem3D_o, lcmesh, pMat1D, is_added )
    implicit none
    class(ElementBase3D), intent(in) :: elem3D_i
    class(ElementBase3D), intent(in) :: elem3D_o
    class(LocalMesh3D), intent(in) :: lcmesh
    real(RP), intent(inout) :: q_o(elem3D_o%Nnode_h1D,elem3D_o%Nnode_h1D,elem3D_o%Nnode_v,lcmesh%NeA)
    real(RP), intent(in) :: q_i(elem3D_i%Nnode_h1D,elem3D_i%Nnode_h1D,elem3D_i%Nnode_v,lcmesh%NeA)
    real(RP), intent(in) :: pMat1D(elem3D_o%Nnode_h1D,elem3D_i%Nnode_h1D)
    logical, intent(in) :: is_added

    integer :: ke

    integer :: px, py, pz
    integer :: pxx, pyy, pzz
    real(RP) :: tmp1
    real(RP) :: tmp2(elem3D_o%Nnode_h1D,elem3D_i%Nnode_h1D)
    real(RP) :: tmp3(elem3D_o%Nnode_h1D,elem3D_o%Nnode_h1D,elem3D_i%Nnode_v)
    real(RP) :: tmp4(elem3D_o%Nnode_h1D)
    real(RP) :: tmp_h(elem3D_o%Nnode_h1D,elem3D_o%Nnode_h1D)

    real(RP) :: mat_tr(elem3D_i%Nnode_h1D,elem3D_o%Nnode_h1D)
    !-------------------------------------------

    mat_tr(:,:) = transpose(pMat1D)

    !$omp parallel do private(tmp1, tmp2, tmp3, tmp4, tmp_h)
    do ke=lcmesh%NeS, lcmesh%NeE

      do pz=1, elem3D_i%Nnode_v
        do py=1, elem3D_i%Nnode_h1D
        do pxx=1, elem3D_o%Nnode_h1D
          tmp1 = 0.0_RP
          do px=1, elem3D_i%Nnode_h1D
            tmp1 = tmp1 + mat_tr(px,pxx) * q_i(px,py,pz,ke)
          end do
          tmp2(pxx,py) = tmp1
        end do
        end do

        do pyy=1, elem3D_o%Nnode_h1D
          tmp4(:) = 0.0_RP
          do py=1, elem3D_i%Nnode_h1D
            do px=1, elem3D_o%Nnode_h1D
              tmp4(px) = tmp4(px) + mat_tr(py,pyy) * tmp2(px,py)
            end do
          end do
          tmp3(:,pyy,pz) = tmp4(:)
        end do
      end do

      if ( is_added ) then
        do pzz=1, elem3D_o%Nnode_v
          tmp_h(:,:) = 0.0_RP
          do pz=1, elem3D_i%Nnode_v
            do py=1, elem3D_o%Nnode_h1D
            do px=1, elem3D_o%Nnode_h1D
              tmp_h(px,py) = tmp_h(px,py) + mat_tr(pz,pzz) * tmp3(px,py,pz)
            end do
            end do
          end do
          q_o(:,:,pzz,ke) = q_o(:,:,pzz,ke) + tmp_h(:,:)
        end do
      else
        do pzz=1, elem3D_o%Nnode_v
          tmp_h(:,:) = 0.0_RP
          do pz=1, elem3D_i%Nnode_v
            do py=1, elem3D_o%Nnode_h1D
            do px=1, elem3D_o%Nnode_h1D
              tmp_h(px,py) = tmp_h(px,py) + mat_tr(pz,pzz) * tmp3(px,py,pz)
            end do
            end do
          end do
          q_o(:,:,pzz,ke) = tmp_h(:,:)
        end do
      end if
    end do
    return
  end subroutine MultiGridSolver3D_pMG_operation
end module scale_multigrid_solver_3d
