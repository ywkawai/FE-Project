!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_interp_vcoord
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
    PRC_myrank, PRC_abort
  use scale_const, only: &
    PI => CONST_PI,          &
    RPlanet => CONST_radius, &
    UNDEF => CONST_UNDEF,    &
    EPS   => CONST_EPS  
  use scale_element_base, only: ElementBase2D, ElementBase3D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_mesh_topography, only: &
    MeshTopography
  use scale_file_base_meshfield, only: &
    FILE_base_meshfield
  use scale_meshfield_base, only: MeshField3D

  use mod_regrid_mesh_base, only: &
    regrid_mesh_base,                    &
    REGRID_MESHTYPE_CUBEDSPHERE3D_ID,    &
    REGRID_MESHTYPE_LONLAT3D_ID,         &
    REGRID_MESHTYPE_STRUCTURED3D_ID,     &
    INMESH_ID => REGRID_MESH_BASE_IN_ID, &
    OUTMESH_ID => REGRID_MESH_BASE_OUT_ID
  
  use mod_regrid_nodemap, only: regrid_nodemap
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  type :: vinterpcoef_info
    real(RP), allocatable :: xi2v_coef(:,:)
    integer , allocatable :: xi2v_idx_k (:,:,:)
    integer , allocatable :: xi2v_idx_p (:,:,:)    
  end type

  type, public :: regrid_interp_vcoord
    integer :: vintrp_typeid
    type(MeshTopography) :: topography

    type(vinterpcoef_info), allocatable :: vintrp_info(:)
    type(MeshField3D) :: pres
    type(HexahedralElement) :: elem   

    class(regrid_mesh_base), pointer :: out_mesh_ptr
    type(regrid_mesh_base) :: out_mesh

    type(MeshCubeDom3D) :: mesh3D
    type(MeshCubedSphereDom3D) :: csmesh3D

    logical :: extrapolate

    type(MeshField3D), public :: vintrp_var3D    
  contains
    procedure :: Init => regrid_interp_vcoord_Init
    procedure :: Final => regrid_interp_vcoord_Final
    procedure :: Update_weight => regrid_interp_vcoord_update_weight
    procedure :: Interpolate => regrid_interp_vcoord_interpolate
  end type

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: REGRID_VCOORD_MODEL_ID  = 1
  integer, public, parameter :: REGRID_VCOORD_HEIGHT_ID = 2
  integer, public, parameter :: REGRID_VCOORD_PRESS_ID  = 3
  integer, public, parameter :: REGRID_VCOORD_SIGMA_ID  = 4
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains

!OCL SERIAL
  subroutine regrid_interp_vcoord_Init( this,         &
    out_mesh, nodeMap_list )

    use scale_mesh_base2d, only: &
      MF2D_XY => MeshBase2D_DIMTYPEID_XY   
    use scale_mesh_base3d, only: &
      MeshBase3D_DIMTYPEID_X, MeshBase3D_DIMTYPEID_Y, MeshBase3D_DIMTYPEID_Z,  &
      MeshBase3D_DIMTYPEID_XYZ, MeshBase3D_DIMTYPEID_XYZT      
    use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D
    use scale_meshfieldcomm_cubedspheredom2d, only: MeshFieldCommCubedSphereDom2D
    use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D
    use scale_meshfieldcomm_cubedspheredom3d, only: MeshFieldCommCubedSphereDom3D
    use scale_meshutil_vcoord, only: MESH_VCOORD_TERRAIN_FOLLOWING_ID

    implicit none

    class(regrid_interp_vcoord), intent(inout), target :: this
    class(regrid_mesh_base), intent(inout), target :: out_mesh
    type(regrid_nodemap), intent(in), target :: nodeMap_list(:)

    character(len=H_SHORT) :: vintrp_name
    character(len=H_LONG ) :: in_topofile_basename = ''       ! Basename of the input  file 
    character(len=H_SHORT) :: topo_varname         = 'topo'  
    character(len=H_SHORT) :: pres_varname         = 'PRES'  

    integer, parameter :: FZ_nmax  = 1000 
    integer :: out_NeZ             = -1    
    integer :: out_PolyOrder_v     = 1    
    real(RP) :: out_dom_vmin
    real(RP) :: out_dom_vmax       
    real(RP) :: out_FZ(FZ_nmax)
    logical :: extrapolate

    namelist / PARAM_REGRID_VCOORD / &
      vintrp_name,          &
      in_topofile_basename, &
      topo_varname,         &
      pres_varname,         &
      out_NeZ,              &  
      out_PolyOrder_v,      &   
      out_dom_vmin,         &
      out_dom_vmax,         &             
      out_Fz,               &
      extrapolate

    integer :: ierr

    type(FILE_base_meshfield) :: topo_file
    type(MeshFieldCommRectDom2D) :: comm2d
    type(MeshFieldCommCubeDom3D) :: comm3d
    type(MeshFieldCommCubedSphereDom2D) :: comm2d_cs
    type(MeshFieldCommCubedSphereDom3D) :: comm3d_cs

    integer :: n
    type(regrid_nodemap), pointer :: nmap
    type(LocalMesh3D), pointer :: lcmesh
    type(ElementBase3D), pointer :: elem

    class(MeshBase2D), pointer :: ptr_mesh2D
    class(regrid_mesh_base), pointer :: in_mesh

    class(MeshBase3D), pointer :: vintrp_ptr_outmesh3D
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("regrid_vcoord_Init",*) 'Setup'

    vintrp_name = 'MODEL'
    out_FZ(:)   = UNDEF
    extrapolate = .false.

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_REGRID_VCOORD,iostat=ierr)
    if ( ierr < 0 ) then !--- missing
        LOG_INFO("regrid_vcoord",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("regrid_vcoord",*) 'Not appropriate names in namelist PARAM_REGRID_VCOORD. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_REGRID_VCOORD)
    
    select case( vintrp_name )
    case( 'MODEL' )
      this%vintrp_typeid= REGRID_VCOORD_MODEL_ID
      this%out_mesh_ptr => out_mesh
    case( 'HEIGHT' )
      this%vintrp_typeid = REGRID_VCOORD_HEIGHT_ID   
    case( 'PRESSURE' )
      this%vintrp_typeid = REGRID_VCOORD_PRESS_ID
    case( 'SIGMA' )
      this%vintrp_typeid = REGRID_VCOORD_SIGMA_ID      
    case default
      LOG_ERROR("regrid_vcoord_Init",*) 'Not appropriate vintrp_type. Check!', vintrp_name
      call PRC_abort
    end select

    !- Vertical interpolation  (MODEL)

    do n=1, out_mesh%ptr_mesh3D%LOCAL_MESH_NUM
      lcmesh => out_mesh%ptr_mesh3D%lcmesh_list(n)
      !---
      nmap => nodeMap_list(n)
      in_mesh => nmap%in_mesh_list(1)
      call regrid_vcoord_search_pos_model( this, &
        nmap%elem_k(:,:), nmap%elem_z(:,:),                               & ! (out)
        lcmesh, lcmesh%refElem3D, nmap%local_domID(:,:), nmap%lcprc(:,:), & ! (in)
        nmap%elem_i(:,:), nmap%elem_j(:,:),                               & ! (in)
        nmap%in_mesh_list, in_mesh%NeZ, in_mesh%NeX, in_mesh%NeY          ) ! (in)
    end do

    this%extrapolate = extrapolate

    !- Prepair topography data

    call out_mesh%ptr_mesh3D%GetMesh2D( ptr_mesh2D )

    call this%topography%Init( topo_varname, ptr_mesh2D )
    if ( in_topofile_basename /= '' ) then
      select type( ptr_mesh2D )
      class is (MeshRectDom2D)
        call topo_file%Init( 1, mesh2D=ptr_mesh2D )
      class is (MeshCubedSphereDom2D)
        call topo_file%Init( 1, meshCubedSphere2D=ptr_mesh2D )
      end select
      call topo_file%Open( in_topofile_basename, PRC_myrank )
      call topo_file%Read_Var( MF2D_XY, topo_varname, this%topography%topo )
      call topo_file%Close()
      call topo_file%Final()
    end if

    !-- Setup vertical coordinate 

    select type( ptr_mesh2D )
    class is (MeshRectDom2D)
      call comm2d%Init( 1, 0, 0, ptr_mesh2D )
      call comm3d%Init( 1, 1, 0, out_mesh%mesh3D )
      call this%topography%SetVCoordinate( out_mesh%ptr_mesh3D,  & ! (inout)
        MESH_VCOORD_TERRAIN_FOLLOWING_ID, out_mesh%dom_zmax,     & ! (in)
        comm3d, comm2d )                                           ! (in)
      call comm2d%Final()
      call comm3d%Final()
    class is (MeshCubedSphereDom2D)
      call comm2d_cs%Init( 1, 0, 0, ptr_mesh2D )
      call comm3d_cs%Init( 1, 1, 0, out_mesh%csmesh3D )
      call this%topography%SetVCoordinate( out_mesh%ptr_mesh3D,  & ! (inout)
        MESH_VCOORD_TERRAIN_FOLLOWING_ID, out_mesh%dom_zmax,     & ! (in)
        comm3d_cs, comm2d_cs )                                     ! (in)
      call comm2d_cs%Final()
      call comm3d_cs%Final()
    end select

    !--

    if ( this%vintrp_typeid /= REGRID_VCOORD_MODEL_ID ) then  
  
      call this%elem%Init( out_mesh%elem3D%PolyOrder_h, out_PolyOrder_v, .true. )

      call this%out_mesh%Init( OUTMESH_ID, out_mesh%mesh_type_id, .false. )
      this%out_mesh_ptr => this%out_mesh

      select type( ptr_mesh2D )
      class is (MeshRectDom2D)

        call this%mesh3D%Init( out_mesh%NeGX, out_mesh%NeGY, out_NeZ,                 &
          out_mesh%dom_xmin, out_mesh%dom_xmax, out_mesh%dom_ymin, out_mesh%dom_ymax, &
          out_dom_vmin, out_dom_vmax,                                                 &
          out_mesh%isPeriodicX, out_mesh%isPeriodicY, out_mesh%isPeriodicZ,           &
          this%elem, out_mesh%NLocalMeshPerPRC,                                       &
          out_mesh%NprcX, out_mesh%NprcY,                                             &
          FZ=out_FZ(1:out_NeZ+1) ) 

        if ( out_mesh%mesh_type_id == REGRID_MESHTYPE_LONLAT3D_ID ) then
          call this%mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_X, 'lon', 'degree_east', 'longitude' )
          call this%mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_Y, 'lat', 'degree_north', 'latitude' )
          call this%mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_XYZ, 'lonlatv', 'degree', 'longitude,latitude,altitude' )
          call this%mesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_XYZT, 'lonlatzv', 'degree', 'longitude,latitude,altitude' )    
        end if

        this%out_mesh_ptr%ptr_mesh3D => this%mesh3D

      class is (MeshCubedSphereDom2D)

        call this%csmesh3D%Init( out_mesh%NeGX, out_mesh%NeGY, out_NeZ, RPlanet,      &
          out_dom_vmin, out_dom_vmax,                                                 &
          this%elem, out_mesh%NLocalMeshPerPRC,                                       &
          nproc=out_mesh%Nprc, FZ=out_FZ(1:out_NeZ+1) ) 

          this%out_mesh_ptr%ptr_mesh3D => this%csmesh3D

      end select
            
      vintrp_ptr_outmesh3D => this%out_mesh_ptr%ptr_mesh3D
      call vintrp_ptr_outmesh3D%Generate()

      if ( this%vintrp_typeid == REGRID_VCOORD_PRESS_ID ) then
        call vintrp_ptr_outmesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_Z, 'p', 'Pa', 'altitude (preesure coordinate)', &
          positive_down = .true. )
        call this%pres%Init( trim(pres_varname), "Pa", out_mesh%ptr_mesh3D )
      else if ( this%vintrp_typeid == REGRID_VCOORD_SIGMA_ID ) then
        call vintrp_ptr_outmesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_Z, 'sig', '1', 'altitude (sigma coordinate)', &
          positive_down = .true. )
        call this%pres%Init( trim(pres_varname), "Pa", out_mesh%ptr_mesh3D )
      else if ( this%vintrp_typeid == REGRID_VCOORD_HEIGHT_ID ) then
        call vintrp_ptr_outmesh3D%SetDimInfo( MeshBase3D_DIMTYPEID_Z, 'z', 'm', 'altitude (height coordinate)' )        
      end if   

      !------------------------

      allocate( this%vintrp_info(vintrp_ptr_outmesh3D%LOCAL_MESH_NUM) ) 

      do n=1, vintrp_ptr_outmesh3D%LOCAL_MESH_NUM        
        lcmesh => vintrp_ptr_outmesh3D%lcmesh_list(n)
        elem => lcmesh%refElem3D 
        allocate( this%vintrp_info(n)%xi2v_coef(elem%Np,lcmesh%Ne) )
        allocate( this%vintrp_info(n)%xi2v_idx_k(elem%Np,lcmesh%Ne,2) ) 
        allocate( this%vintrp_info(n)%xi2v_idx_p(elem%Np,lcmesh%Ne,2) )        
      end do
    
    end if

    call this%vintrp_var3D%Init( 'vintrp_var3D', '', this%out_mesh_ptr%ptr_mesh3D )

    return
  end subroutine regrid_interp_vcoord_Init

!OCL SERIAL  
  subroutine regrid_interp_vcoord_Final( this )
    implicit none
    class(regrid_interp_vcoord), intent(inout) :: this

    integer :: n
    class(MeshBase3D), pointer :: mesh3D
    !-------------------------------------------

    call this%vintrp_var3D%Final()
    call this%topography%Final()

    if ( this%vintrp_typeid == REGRID_VCOORD_PRESS_ID ) &
      call this%pres%Final()
    
    if ( allocated(this%vintrp_info) ) then
      do n=1, size(this%vintrp_info)
        deallocate( this%vintrp_info(n)%xi2v_coef )
        deallocate( this%vintrp_info(n)%xi2v_idx_k ) 
        deallocate( this%vintrp_info(n)%xi2v_idx_p )               
      end do
      deallocate( this%vintrp_info )      
    end if

    if (this%vintrp_typeid /= REGRID_VCOORD_MODEL_ID ) then
      call this%out_mesh%Final()
      call this%elem%Final()
    end if

    nullify( this%out_mesh_ptr )

    return
  end subroutine regrid_interp_vcoord_Final

!OCL SERIAL  
  subroutine regrid_interp_vcoord_update_weight( this, &
    istep, out_mesh, nodeMap_list,                     &
    GP_flag, out_mesh_GP, GPMat                        )

    use mod_regrid_interp_field, only: &
      regrid_interp_field_Interpolate
    implicit none
    class(regrid_interp_vcoord), intent(inout) :: this
    integer, intent(in) :: istep    
    class(regrid_mesh_base), intent(in), target :: out_mesh
    type(regrid_nodemap), intent(in) :: nodeMap_list(:)
    logical, intent(in) :: GP_flag
    class(regrid_mesh_base), intent(in), target :: out_mesh_GP
    real(RP), intent(in) :: GPMat(:,:)

    integer :: n
    class(MeshBase3D), pointer :: mesh3D
    !---------------------------------------------------------------

    select case( this%vintrp_typeid )
    case ( REGRID_VCOORD_MODEL_ID )
      return      
    case ( REGRID_VCOORD_PRESS_ID, REGRID_VCOORD_SIGMA_ID )
      call regrid_interp_field_Interpolate( this%pres,     &
        istep, this%pres%varname, out_mesh, nodeMap_list,  &
        GP_flag, out_mesh_GP, GPMat                        )  
    end select
    
    mesh3D => out_mesh%ptr_mesh3D
    do n=1, mesh3D%LOCAL_MESH_NUM
      call interp_vcoord_update_weight_core( this, n,           &
        mesh3D%lcmesh_list(n), mesh3D%refElem3D,                &
        this%out_mesh_ptr%ptr_mesh3D%lcmesh_list(n), this%elem  )
    end do

    return
  end subroutine regrid_interp_vcoord_update_weight

!OCL SERIAL  
  subroutine regrid_interp_vcoord_interpolate( this, istep, mesh_ref, field_ref )
    implicit none
    class(regrid_interp_vcoord), intent(inout) :: this
    integer, intent(in) :: istep    
    class(MeshBase3D), intent(in) :: mesh_ref
    class(MeshField3D), intent(in) :: field_ref

    integer :: n
    !---------------------------------------------------------------

    if ( this%vintrp_typeid == REGRID_VCOORD_MODEL_ID ) return

    do n=1, this%out_mesh_ptr%ptr_mesh3D%LOCAL_MESH_NUM
      call interp_vcoord_interpolate_core( this%vintrp_info(n), n,   &
        this%vintrp_var3D%local(n)%val, field_ref%local(n)%val,      &
        mesh_ref%lcmesh_list(n), mesh_ref%lcmesh_list(n)%refElem3D,  &
        this%out_mesh_ptr%ptr_mesh3D%lcmesh_list(n), this%elem       )
    end do
       
    return
  end subroutine regrid_interp_vcoord_interpolate

!- private -------------------------------------

!OCL SERIAL
  subroutine interp_vcoord_interpolate_core(  &
    vintrp_info, domID, var_intrp, var_ref, lcmesh_ref, elem_ref, lcmesh, elem  )

    implicit none
    type(vinterpcoef_info), intent(inout) :: vintrp_info    
    integer, intent(in) :: domID
    class(LocalMesh3D), intent(in) :: lcmesh_ref
    class(ElementBase3D), intent(in) :: elem_ref
    class(LocalMesh3D), intent(in) :: lcmesh   
    class(ElementBase3D), intent(in) :: elem    
    real(RP), intent(out) :: var_intrp(elem    %Np,lcmesh    %NeA)
    real(RP), intent(in) ::  var_ref  (elem_ref%Np,lcmesh_ref%NeA)

    integer :: ke, p
    integer :: ke_r, p_r, kke_r, pp_r
    real(RP) :: coef
    !--------------------------------------------------

    !$omp parallel do collapse(2) private(ke, p, ke_r, p_r, kke_r, pp_r, coef) 
    do ke=lcmesh%NeS, lcmesh%NeE
    do p=1, elem%Np
      ke_r  = vintrp_info%xi2v_idx_k(p,ke,1)
      kke_r = vintrp_info%xi2v_idx_k(p,ke,2)
      p_r   = vintrp_info%xi2v_idx_p(p,ke,1)
      pp_r  = vintrp_info%xi2v_idx_p(p,ke,2)
      coef = vintrp_info%xi2v_coef(p,ke)
      var_intrp(p,ke) = UNDEF
      if ( ke_r  > 0 ) var_intrp(p,ke) = coef * var_ref(p_r,ke_r) 
      if ( kke_r > 0 ) var_intrp(p,ke) = var_intrp(p,ke) + (1.0_RP - coef) * var_ref(pp_r,kke_r)
    end do
    end do

    return
  end subroutine interp_vcoord_interpolate_core

!OCL SERIAL
  subroutine interp_vcoord_update_weight_core(       &
    this, domID, lcmesh_ref, elem_ref, lcmesh, elem  )

    implicit none

    class(regrid_interp_vcoord), intent(inout) :: this    
    integer, intent(in) :: domID    
    class(LocalMesh3D), intent(in), target :: lcmesh_ref
    class(ElementBase3D), intent(in) :: elem_ref
    class(LocalMesh3D), intent(in) :: lcmesh   
    class(ElementBase3D), intent(in) :: elem

    integer :: p, ke
    integer :: p_ref, pp_ref, ke_ref, kke_ref
    integer :: ke_xy, ke_z, ke_zz
    integer :: p_xy, p_z, p_zz
    integer :: ke_sfc_ref, ke_top_ref
    integer :: p_sfc_ref, p_top_ref
    real(RP) :: height_ref(elem_ref%Np,lcmesh_ref%Ne)
    real(RP) :: height    (elem%Np    ,lcmesh%Ne    )

    integer :: indx_k(elem%Np,2)
    integer :: indx_p(elem%Np,2)
    real(RP) :: vfact(elem%Np)

    class(LocalMesh2D), pointer :: lcmesh2D
    real(RP), allocatable :: SfcPres(:,:)

    !-------------------------------------------------------------
    
    if ( this%vintrp_typeid == REGRID_VCOORD_PRESS_ID ) then
      !$omp parallel
      !$omp do
      do ke=lcmesh_ref%NeS, lcmesh_ref%NeE
       height_ref(:,ke) = - log( this%pres%local(domID)%val(:,ke) )
      end do
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
       height(:,ke) = - log( lcmesh%pos_en(:,ke,3) )
      end do      
      !$omp end parallel
    else if ( this%vintrp_typeid == REGRID_VCOORD_SIGMA_ID ) then
      lcmesh2D => lcmesh_ref%lcmesh2D
      allocate( SfcPres(lcmesh2D%refElem2D%Np,lcmesh2D%NeA) )
      !$omp parallel private(ke, ke_xy)
      !$omp do
      do ke_xy=lcmesh2D%NeS, lcmesh2D%NeE
        SfcPres(:,ke_xy) = this%pres%local(domID)%val(elem_ref%Hslice(:,1),ke_xy)
      end do
      !$omp do
      do ke=lcmesh_ref%NeS, lcmesh_ref%NeE
        ke_xy = lcmesh_ref%EMap3Dto2D(ke)
        height_ref(:,ke) = - log( this%pres%local(domID)%val(:,ke) )      &
                           + log( SfcPres(elem_ref%IndexH2Dto3D(:),ke_xy) )
      end do
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
        ke_xy = lcmesh%EMap3Dto2D(ke)
        height(:,ke) = - log( lcmesh%pos_en(:,ke,3) )
      end do      
      !$omp end parallel      
    else if ( this%vintrp_typeid == REGRID_VCOORD_HEIGHT_ID ) then
      !$omp parallel
      !$omp do
      do ke=lcmesh_ref%NeS, lcmesh_ref%NeE
       height_ref(:,ke) = lcmesh_ref%zlev(:,ke)
      end do
      !$omp do
      do ke=lcmesh%NeS, lcmesh%NeE
       height(:,ke) = lcmesh%pos_en(:,ke,3)
      end do      
      !$omp end parallel
    end if

    !$omp parallel do collapse(2) private( &
    !$omp ke_z, ke_xy, ke, ke_sfc_ref, ke_top_ref,    &
    !$omp indx_k, indx_p, vfact,                      &
    !$omp p_z, p_xy, p, p_sfc_ref, p_top_ref,         &
    !$omp ke_zz, p_zz, ke_ref, p_ref, kke_ref, pp_ref )
    do ke_z=1, lcmesh%NeZ    
    do ke_xy=1, lcmesh%NeX * lcmesh%NeY
      ke = ke_xy + (ke_z - 1) * lcmesh%NeX * lcmesh%NeY
      ke_sfc_ref = ke_xy
      ke_top_ref = ke_xy + (lcmesh_ref%NeZ - 1) * lcmesh%NeX * lcmesh%NeY

      indx_k(:,:) = -1
      indx_p(:,:) = -1
      vfact(:)    = UNDEF     
      do p_z=1, elem%Nnode_v      
      do p_xy=1, elem%Nnode_h1D**2
        p = p_xy + (p_z - 1) * elem%Nnode_h1D**2
        p_sfc_ref = p_xy
        p_top_ref = p_xy + (elem_ref%Nnode_v -1) * elem%Nnode_h1D**2

        if ( height(p,ke) < height_ref(p_sfc_ref,ke_sfc_ref) - EPS ) then
          if ( this%extrapolate ) then
            indx_k(p,1) = ke_sfc_ref; indx_p(p,1) = p_sfc_ref
            vfact(p) = 1.0_RP
          end if
        else if ( height(p,ke) < height_ref(p_sfc_ref,ke_sfc_ref) ) then
          indx_k(p,1) = ke_sfc_ref; indx_p(p,1) = p_sfc_ref
          vfact(p) = 1.0_RP          
        else if ( height(p,ke) > height_ref(p_top_ref,ke_top_ref) + EPS ) then
          if ( this%extrapolate ) then
            indx_k(p,1) = ke_top_ref; indx_p(p,1) = p_top_ref
            vfact(p) = 1.0_RP          
          end if
        else if ( height(p,ke) >= height_ref(p_top_ref,ke_top_ref) ) then 
          indx_k(p,1) = ke_top_ref; indx_p(p,1) = p_top_ref
          vfact(p) = 1.0_RP          
        else  

          search : do ke_zz = 1, lcmesh_ref%NeZ          
          do p_zz = 1, elem_ref%Nnode_v              
            ke_ref = ke_xy + (ke_zz-1) * lcmesh%NeX * lcmesh%NeY
            p_ref = p_xy + (p_zz-1) * elem%Nnode_h1D**2
            if ( p_zz < elem_ref%Nnode_v ) then
              pp_ref = p_xy + p_zz * elem%Nnode_h1D**2
              kke_ref = ke_ref
            else
              pp_ref = p_xy + elem%Nnode_h1D**2
              kke_ref = ke_xy + ke_zz * lcmesh%NeX * lcmesh%NeY
            end if
            if (      height(p,ke) >= height_ref(p_ref,ke_ref)  &
                .and. height(p,ke) < height_ref(pp_ref,kke_ref) ) then
              
              indx_k(p,1) = ke_ref; indx_k(p,2) = kke_ref;
              indx_p(p,1) =  p_ref; indx_p(p,2) =  pp_ref;
              vfact(p) = ( height_ref(pp_ref,kke_ref) - height    (p,ke)       ) &
                       / ( height_ref(pp_ref,kke_ref) - height_ref(p_ref,ke_ref) )
              
              exit search
            end if
          end do
          end do search

        end if
      end do
      end do

      this%vintrp_info(domID)%xi2v_idx_k(:,ke,1) = indx_k(:,1)
      this%vintrp_info(domID)%xi2v_idx_k(:,ke,2) = indx_k(:,2)  
      this%vintrp_info(domID)%xi2v_idx_p(:,ke,1) = indx_p(:,1)
      this%vintrp_info(domID)%xi2v_idx_p(:,ke,2) = indx_p(:,2) 
      this%vintrp_info(domID)%xi2v_coef(:,ke) = vfact(:)         
    end do
    end do

    return    
  end subroutine interp_vcoord_update_weight_core

!OCL SERIAL
  subroutine regrid_vcoord_search_pos_model( this, elem_k, elem_out, &
    lcmesh, elem3D, local_domID, lcprc, elem_i, elem_j,              &
    in_mesh_list, in_NeGZ, in_NeX, in_NeY                            )
 
    implicit none

    type(regrid_interp_vcoord), intent(in) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    type(ElementBase3D), intent(in) :: elem3D
    integer, intent(out) :: elem_k(elem3D%Np,lcmesh%NeX*lcmesh%NeY*lcmesh%NeZ)
    real(RP), intent(out) :: elem_out(elem3D%Np,lcmesh%NeX*lcmesh%NeY*lcmesh%NeZ)
    integer, intent(in) :: local_domID(elem3D%Nnode_h1D**2,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: lcprc(elem3D%Nnode_h1D**2,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: elem_i(elem3D%Nnode_h1D**2,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: elem_j(elem3D%Nnode_h1D**2,lcmesh%NeX*lcmesh%NeY)
    type(regrid_mesh_base), intent(in), target :: in_mesh_list(:)
    integer, intent(in) :: in_NeGZ
    integer, intent(in) :: in_NeX
    integer, intent(in) :: in_NeY
   
    integer :: in_n
    integer :: in_prc
    integer :: ke, p
    integer :: ke_h, p_h
    integer :: ke_z, ke_z2
    integer :: p_z

    type(LocalMesh3D), pointer :: in_lcmesh    
    real(RP) :: in_Z0, in_Z1
    integer :: in_ke3D  
    !-------------------------------------------
        
    !$omp parallel private( &
    !$omp ke_h, p_h, in_lcmesh, in_prc, in_n,            &
    !$omp ke_z, ke_z2, p_z, ke, p, in_ke3D, in_Z0, in_Z1 )

    !$omp workshare
    elem_k(:,:) = -1
    !$omp end workshare
    
    !$omp do collapse(2)
    do ke_h=1, lcmesh%NeX * lcmesh%NeY
    do p_h=1, elem3D%Nnode_h1D**2

      in_n   = local_domID(p_h,ke_h)
      in_prc = lcprc(p_h,ke_h)
      
      if ( in_n > 0 .and. in_prc > 0 .and.                 &
           elem_i(p_h,ke_h) > 0 .and. elem_j(p_h,ke_h) > 0 ) then

        in_lcmesh => in_mesh_list(in_prc)%ptr_mesh3D%lcmesh_list(in_n)

        do ke_z=1, lcmesh%NeZ
        do p_z=1, elem3D%Nnode_v
          ke = ke_h + (ke_z-1)*lcmesh%NeX*lcmesh%NeY
          p = p_h + (p_z-1)*elem3D%Nnode_h1D**2
          do ke_z2=1, in_NeGZ
            in_ke3D = elem_i(p_h,ke_h) + (elem_j(p_h,ke_h) - 1) * in_NeX &
                    + (ke_z2 - 1) * in_NeX * in_NeY
            in_Z0 = in_lcmesh%pos_ev(in_lcmesh%EToV(in_ke3D,1),3)
            in_Z1 = in_lcmesh%pos_ev(in_lcmesh%EToV(in_ke3D,5),3)
            if ( in_Z0 <= lcmesh%pos_en(p,ke,3) .and. lcmesh%pos_en(p,ke,3) <= in_Z1 ) then
              elem_k  (p,ke) = ke_z2
              elem_out(p,ke) = lcmesh%pos_en(p,ke,3)
              exit
            end if
          end do
        end do
        end do
      end if

    end do
    end do
   !$omp end do
   !$omp end parallel

    return
  end subroutine regrid_vcoord_search_pos_model

end module mod_regrid_interp_vcoord