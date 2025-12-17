#include "scaleFElib.h"
module mod_poisson2d_mg
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale_prc
  use scale_io

  use scale_sparsemat, only: SparseMat
  
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D
  
  use scale_mesh_hierarchy_base, only: &
    pMG_FINEST_LEVEL => MESH_HIERARCHY_pMG_FINEST_LEVEL, &
    hMG_FINEST_LEVEL => MESH_HIERARCHY_hMG_FINEST_LEVEL

  use scale_mesh_hierarchy_2d, only: MeshHierarchy2D

  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D
  
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: Poisson2d_mg_Init
  public :: Poisson2d_mg_Final
  public :: Poisson2d_MG_solve
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

  type, public :: FieldSet
    type(MeshField2D) :: dq
    type(MeshField2D) :: f
    type(MeshField2D) :: res
    type(MeshField2D) :: qx
    type(MeshField2D) :: qy
    integer :: level_id

    type(MeshFieldCommRectDom2D) :: var_comm
    type(MeshFieldCommRectDom2D) :: aux_comm
    type(SparseMat) :: Dx
    type(SparseMat) :: Dy
  contains
    procedure :: Init => FieldSet_Init
    procedure :: Final => FieldSet_Final
  end type FieldSet

  type(FieldSet), allocatable, target :: fields_h(:)
  type(FieldSet), allocatable, target :: fields_p(:)

  integer, allocatable :: pMG_ITR_NUM_list(:)
  integer, allocatable :: hMG_ITR_NUM_list(:)

  type(MeshHierarchy2D) :: mesh_hierarchy

  real(RP) :: gtau

contains
  !> Initialization
  subroutine Poisson2d_mg_Init( mesh )
    use mod_poisson2d_smoother, only: &
      Poisson2d_smoother_Init    
    implicit none
    class(MeshBase2D), intent(in) :: mesh
    
    integer :: lev_h
    integer :: lev_p

    integer, parameter :: pMG_LV_LIST_MAX = 16
    integer :: P_LEVEL_NUM
    integer :: P_LEVEL_LIST(pMG_LV_LIST_MAX)
    integer :: P_ITR_NUM_LIST(pMG_LV_LIST_MAX)


    integer, parameter :: hMG_LV_LIST_MAX = 16
    integer :: H_LEVEL_NUM
    integer :: NeGX_list(hMG_LV_LIST_MAX)
    integer :: NeGY_list(hMG_LV_LIST_MAX)
    integer :: H_ITR_NUM_LIST(hMG_LV_LIST_MAX)

    namelist /PARAM_Poisson2D_MG/ &
      H_LEVEL_NUM, NeGX_list, NeGY_list, H_ITR_NUM_LIST, &
      P_LEVEL_NUM, P_LEVEL_LIST, P_ITR_NUM_LIST

    integer :: ierr
    !---------------------------------------------------------------------------


    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_Poisson2D_MG,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("Poisson2D_MG_Init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("Poisson2D_MG_Init",*) 'Not appropriate names in namelist PARAM_Poisson2D_MG. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_Poisson2D_MG)

    !- Setup an object to manage the mesh hierarchy

    LOG_INFO("Poisson2D_MG_Init",*) "Setup mesh hierarchy"

    call mesh_hierarchy%Init( mesh, &
      P_LEVEL_LIST(1:P_LEVEL_NUM), P_LEVEL_NUM,                       &
      NeGX_list(1:H_LEVEL_NUM), NeGY_list(1:H_LEVEL_NUM), H_LEVEL_NUM )

    !- Setup objects to manage the field set
    
    LOG_INFO("Poisson2D_MG_Init",*) "Setup field sets"
      
    allocate( fields_h(mesh_hierarchy%NUM_hMG_LEVEL) )
    do lev_h=1, mesh_hierarchy%NUM_hMG_LEVEL
      call fields_h(lev_h)%Init( mesh_hierarchy%h_mesh_list(lev_h)%ptr, lev_h )
    end do

    allocate( fields_p(mesh_hierarchy%NUM_pMG_LEVEL) )
    do lev_p=1, mesh_hierarchy%NUM_pMG_LEVEL
      call fields_p(lev_p)%Init( mesh_hierarchy%p_mesh_list(lev_p)%ptr, lev_p )
    end do

    !-
    allocate( pMG_ITR_NUM_list(P_LEVEL_NUM) )
    pMG_ITR_NUM_list(:) = P_ITR_NUM_LIST(1:P_LEVEL_NUM)

    allocate( hMG_ITR_NUM_list(H_LEVEL_NUM) )
    hMG_ITR_NUM_list(:) = H_ITR_NUM_LIST(1:H_LEVEL_NUM)

    !-
    LOG_INFO("Poisson2D_MG_Init",*) "Setup smoother"
    call Poisson2d_smoother_Init( mesh )

    !-

    return
  end subroutine Poisson2d_mg_Init

  !> Finalization
  subroutine Poisson2d_mg_Final()
    use mod_poisson2d_smoother, only: Poisson2d_smoother_Final
    implicit none

    integer :: lev_h
    integer :: lev_p
    !---------------------------------------------------------------------------

    !-
    call Poisson2d_smoother_Final()

    !-
    do lev_h=1, size(fields_h)
      call fields_h(lev_h)%Final()
    end do

    do lev_p=1, size(fields_p)
      call fields_p(lev_p)%Final()
    end do

    !-
    call mesh_hierarchy%Final()
    return
  end subroutine Poisson2d_mg_Final

!OCL SERIAL
  subroutine Poisson2d_MG_solve( q, &
    f )
    implicit none
    class(MeshField2D), intent(inout), target :: q
    class(MeshField2D), intent(inout) :: f

    logical :: cal_res_flag
    integer :: itr

    class(LocalMesh2D), pointer :: lmesh
    integer :: ldomID
    integer :: ke
    !-------------------------------------

    !-
    do ldomID=1, q%mesh%LOCAL_MESH_NUM
      lmesh => q%mesh%lcmesh_list(ldomID)
      do ke=lmesh%NeS, lmesh%NeE
        fields_p(pMG_FINEST_LEVEL)%dq%local(ldomID)%val(:,ke) = q%local(ldomID)%val(:,ke)
      end do
    end do

    !-
    do itr=1, 40
      cal_res_flag = ( mod(itr,10) == 0 .or. itr == 1 )
      call Poisson2d_mg_Vcycle( pMG_FINEST_LEVEL, f )
    end do

    !-
    do ldomID=1, q%mesh%LOCAL_MESH_NUM
      lmesh => q%mesh%lcmesh_list(ldomID)
      do ke=lmesh%NeS, lmesh%NeE
        q%local(ldomID)%val(:,ke) = fields_p(pMG_FINEST_LEVEL)%dq%local(ldomID)%val(:,ke)
      end do
    end do
    return
  end subroutine Poisson2d_MG_solve

!OCL SERIAL
  recursive subroutine Poisson2d_mg_Vcycle( mg_level, f_in )
    use mod_poisson2d_smoother, only: Poisson2d_smoother_advance_itr_1step
    implicit none
    integer, intent(in) :: mg_level
    type(MeshField2D), intent(inout) :: f_in

    real(RP) :: itr_res_eps
    integer :: m

    logical :: invoke_hMG

    integer :: itr_num
    logical :: cal_res_flag

    type(FieldSet), pointer :: fs_p

    logical :: is_mg_top
    logical :: zero_initial_guess

    class(LocalMesh2D), pointer :: lmesh2D
    integer :: ldomID
    integer :: ke
    !----------------------------------------------

    invoke_hMG = ( mg_level+1 >= mesh_hierarchy%NUM_pMG_LEVEL .and. mesh_hierarchy%NUM_hMG_LEVEL > 0 )
    is_mg_top = ( mg_level == 1 )

    itr_num = pMG_ITR_NUM_list(mg_level)

    fs_p => fields_p(mg_level)

    LOG_INFO("Poisson2d_mg_Vcycle",*) "mg_level=", mg_level

    !- Pre-relaxation
    do m=1, itr_num
      cal_res_flag = (m == 1 .or. m==itr_num .or. mod(m,5)==0)
      zero_initial_guess = ( mg_level /= 1 .and. m==1 )

      call Poisson2d_smoother_advance_itr_1step( fs_p%dq, fs_p%res, &
        f_in, fs_p%qx, fs_p%qy, m, fs_p%var_comm, fs_p%aux_comm,    &
        fs_p%Dx, fs_p%Dy, mesh_hierarchy%p_mesh_list(mg_level)%ptr, &
        cal_res_flag, zero_initial_guess, is_mg_top )
    end do
    if ( mg_level == mesh_hierarchy%NUM_pMG_LEVEL .and. mesh_hierarchy%NUM_hMG_LEVEL == 0 ) then
      LOG_INFO("Poisson2d_mg_Vcycle",*) "End: mg_level=", mg_level
      return
    end if

    !-
    if ( invoke_hMG ) then
      !- Restriction
      call mesh_hierarchy%Operate_pMG_restriction( fields_h(hMG_FINEST_LEVEL)%f, &
        fields_p(mg_level)%res, mg_level )

      !- Advance node in the V-cycle
      call Poisson2d_hMG_Vcycle( hMG_FINEST_LEVEL, fields_h(hMG_FINEST_LEVEL)%f )
      
      !- Correction
      call mesh_hierarchy%Operate_pMG_correction( fields_p(mg_level)%dq, &
        fields_h(hMG_FINEST_LEVEL)%dq, mg_level )
      
      ! do ldomID=1, fs_p%f%mesh%LOCAL_MESH_NUM
      !   lmesh2D => fs_p%f%mesh%lcmesh_list(ldomID)
      !   do ke=lmesh2D%NeS, lmesh2D%NeE
      !     fields_p(mg_level)%f%local(ldomID)%val(:,ke) = fields_h(hMG_FINEST_LEVEL)%f%local(ldomID)%val(:,ke)
      !   end do
      ! end do
    else
      !- Restriction
      call mesh_hierarchy%Operate_pMG_restriction( fields_p(mg_level+1)%f, &
        fields_p(mg_level)%res, mg_level )
      
      !- Advance node in the V-cycle
      call Poisson2d_mg_Vcycle( mg_level+1, fields_p(mg_level+1)%f )

      !- Correction
      ! LOG_INFO("Poisson2d_mg_Vcycle",*) "mg_level=", mg_level, "correction"
      call mesh_hierarchy%Operate_pMG_correction( fields_p(mg_level)%dq, &
        fields_p(mg_level+1)%dq, mg_level )
    end if

    !- Post-relaxation
    do m=1, itr_num
      cal_res_flag = (m == 1 .or. m==ITR_NUM .or. mod(m,5)==0)

      call Poisson2d_smoother_advance_itr_1step( fs_p%dq, fs_p%res, &
        f_in, fs_p%qx, fs_p%qy, m, fs_p%var_comm, fs_p%aux_comm,    &
        fs_p%Dx, fs_p%Dy, mesh_hierarchy%p_mesh_list(mg_level)%ptr, &
        cal_res_flag, .false., is_mg_top )
    end do

    LOG_INFO("Poisson2d_mg_Vcycle",*) "End: mg_level=", mg_level

    return
  end subroutine Poisson2d_mg_Vcycle

!OCL SERIAL
  recursive subroutine Poisson2d_hMG_Vcycle( mg_level, f_in )
    use mod_poisson2d_smoother, only: Poisson2d_smoother_advance_itr_1step
    implicit none
    integer, intent(in) :: mg_level
    type(MeshField2D), intent(inout) :: f_in

    real(RP) :: itr_res_eps
    integer :: m

    integer :: itr_num
    logical :: cal_res_flag

    type(FieldSet), pointer :: fs_h
    logical :: zero_initial_guess
    !----------------------------------------------
    itr_num = hMG_ITR_NUM_list(mg_level)

    fs_h => fields_h(mg_level)

    LOG_INFO("Poisson2d_hMG_Vcycle",*) "mg_level=", mg_level

    if ( mg_level == mesh_hierarchy%NUM_hMG_LEVEL ) then
      LOG_INFO("Poisson2d_hMG_Vcycle",*) "Check f", f_in%local(1)%val(:,1)

      ! Direct solver
      do m=1, itr_num
        zero_initial_guess = ( m==1 )
        cal_res_flag = (m == 1 .or. m==itr_num .or. mod(m,5)==0)

        call Poisson2d_smoother_advance_itr_1step( fs_h%dq, fs_h%res, &
          f_in, fs_h%qx, fs_h%qy, m, fs_h%var_comm, fs_h%aux_comm,    &
          fs_h%Dx, fs_h%Dy, mesh_hierarchy%h_mesh_list(mg_level)%ptr, &
          cal_res_flag, zero_initial_guess, .false. )
      end do      
      return
    end if

    !- Pre-relaxation
    do m=1, itr_num
      cal_res_flag = (m == 1 .or. m==ITR_NUM .or. mod(m,5)==0)
      zero_initial_guess = ( m==1 )

      call Poisson2d_smoother_advance_itr_1step( fs_h%dq, fs_h%res, &
        f_in, fs_h%qx, fs_h%qy, m, fs_h%var_comm, fs_h%aux_comm,    &
        fs_h%Dx, fs_h%Dy, mesh_hierarchy%h_mesh_list(mg_level)%ptr, &
        cal_res_flag, zero_initial_guess, .false. )
    end do

    !-
    !- Restriction
    call mesh_hierarchy%Operate_hMG_restriction( fields_h(mg_level+1)%f, &
      fields_h(mg_level)%res, mg_level )

    !- Advance node in the V-cycle
    call Poisson2d_hMG_Vcycle( mg_level+1, fields_h(mg_level+1)%f )

    !- Correction
    ! LOG_INFO("Poisson2d_mg_Vcycle",*) "mg_level=", mg_level, "correction"
    call mesh_hierarchy%Operate_hMG_correction( fields_h(mg_level)%dq, &
      fields_h(mg_level+1)%dq, mg_level )

    !- Post-relaxation
    do m=1, itr_num
      cal_res_flag = (m == 1 .or. m==ITR_NUM .or. mod(m,5)==0)

      call Poisson2d_smoother_advance_itr_1step( fs_h%dq, fs_h%res, &
        f_in, fs_h%qx, fs_h%qy, m, fs_h%var_comm, fs_h%aux_comm,    &
        fs_h%Dx, fs_h%Dy, mesh_hierarchy%h_mesh_list(mg_level)%ptr, &
        cal_res_flag, .false., .false. )
    end do

    LOG_INFO("Poisson2d_hMG_Vcycle",*) "End: mg_level=", mg_level

    return
  end subroutine Poisson2d_hMG_Vcycle

!--------------

  subroutine FieldSet_Init(this, mesh2D, level_id)
    use scale_mesh_rectdom2d, only: MeshRectDom2D
    implicit none
    class(FieldSet), intent(inout) :: this
    class(MeshBase2D), intent(in) :: mesh2D
    integer, intent(in) :: level_id
    !---------------------------------

    this%level_id = level_id
    call this%dq%Init( "dq", "1", mesh2D )    
    call this%f%Init( "f", "1", mesh2D )
    call this%res%Init( "res", "1", mesh2D )
    call this%qx%Init( "qx", "1", mesh2D )
    call this%qy%Init( "qy", "1", mesh2D )

    select type(mesh2D)
    type is (MeshRectDom2D)
      call this%var_comm%Init( 1, 0, 0, mesh2D )
      call this%aux_comm%Init( 0, 1, 0, mesh2D )
    end select

    call this%Dx%Init( mesh2D%refElem2D%Dx1, storage_format='ELL')
    call this%Dy%Init( mesh2D%refElem2D%Dx2, storage_format='ELL')
    return
  end subroutine FieldSet_Init

  subroutine FieldSet_Final(this)
    implicit none
    class(FieldSet), intent(inout) :: this
    !---------------------------------
    call this%var_comm%Final()
    call this%aux_comm%Final()

    call this%dq%Final()    
    call this%f%Final()    
    call this%res%Final()
    call this%qx%Final()
    call this%qy%Final()
    return
  end subroutine FieldSet_Final

end module mod_poisson2d_mg
