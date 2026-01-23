#include "scaleFElib.h"
module mod_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_const, only: &
    CpDry => CONST_CPdry

  use scale_sparsemat, only: &
    sparsemat

  use scale_meshfield_base, only: MeshField1D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_1d, only: LocalMesh1D  
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_element_base, only: &
    ElementBase1D, ElementBase3D

  use scale_mesh_bndinfo, only: MeshBndInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: common_Init
  public :: common_Final
  !----
  public :: global_horizontal_mean
  !----
  public :: get_reconstructed_flux

  public :: get_del_flux_cent  
  public :: get_del_flux_cent_heat
  !----
  public :: calc_stab_bndflx
  
  !----
  public :: inquire_bound_flag
  
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  type(SparseMat), public :: Dx, Dy, Dz
  type(SparseMat), public :: Sx, Sy, Sz
  type(SparseMat), public :: Lift

  type(MeshBndInfo), allocatable :: VelBC_list(:)
  integer, allocatable:: velBC_ids(:)
  type(MeshBndInfo), allocatable :: thermBC_list(:)
  integer, allocatable:: thermBC_ids(:)

  integer, parameter :: domBnd_South_ID = 1
  integer, parameter :: domBnd_East_ID  = 2
  integer, parameter :: domBnd_North_ID = 3
  integer, parameter :: domBnd_West_ID  = 4
  integer, parameter :: domBnd_Btm_ID   = 5
  integer, parameter :: domBnd_Top_ID   = 6
  integer, parameter :: DOM_BND_NUM     = 6  

  real(RP), private :: BTM_FIXED_HEAT_FLUX  = 15.88_RP
  real(RP), private :: TOP_FIXED_HEAT_FLUX  = 15.88_RP  ! ztop=1.6 km

  real(RP) :: StabCoef_bnd = 0.0_RP

  integer, parameter :: BC_HEAT_FIXED_TEMP    = 1
  integer, parameter :: BC_MOM_NoSlip         = 1
  integer, parameter :: BC_HEAT_FIXED_FLUX    = 2
  integer, parameter :: BC_MOM_CONST_BULKCOEF = 2

  integer :: BTM_BC_HEAT_TYPEID
  integer :: TOP_BC_HEAT_TYPEID

  real(RP), allocatable :: gLB(:), gRB(:)
  integer, allocatable :: vmapM_Z(:,:)
  integer, allocatable :: vmapP_Z(:,:)
  integer, allocatable :: nz_Z(:,:)

  private :: common_gen_vmapZ

contains
!OCL SERIAL
  subroutine common_Init( mesh )
    use scale_polynominal, only: &
      Polynominal_GenLegendrePoly
    use scale_mesh_bndinfo, only: &
      BND_TYPE_NOSLIP_ID,   &
      BND_TYPE_SLIP_ID,     &
      BND_TYPE_PERIODIC_ID, &
      BND_TYPE_FIXVAL_ID
        
    implicit none    
    class(MeshBase3D), target, intent(in) :: mesh

    integer :: n
    integer :: p2D, pZ1D, p

    class(LocalMeshBase), pointer :: ptr_lcmesh
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: refElem3D

    integer :: Nnode_v
    real(RP), allocatable :: nodepos_v(:)
    real(RP), allocatable :: L_pp1(:,:)

    character(H_SHORT) :: BTM_BC_TYPE_HEAT
    character(H_SHORT) :: TOP_BC_TYPE_HEAT

    namelist / PARAM_RBCONV_ANALYSIS_BC / &
      BTM_BC_TYPE_HEAT,   &
      TOP_BC_TYPE_HEAT,   &
      BTM_FIXED_HEAT_FLUX, &
      TOP_FIXED_HEAT_FLUX, &
      StabCoef_bnd

    integer :: ierr
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("common_Init",*) 'Setup'

    refElem3D => mesh%refElem3D
    call Dx%Init( refElem3D%Dx1, storage_format='ELL' )
    call Dy%Init( refElem3D%Dx2, storage_format='ELL' )
    call DZ%Init( refElem3D%Dx3, storage_format='ELL' )
    call Sx%Init( refElem3D%Sx1, storage_format='ELL' )
    call Sy%Init( refElem3D%Sx2, storage_format='ELL' )
    call Sz%Init( refElem3D%Sx3, storage_format='ELL' )
    call Lift%Init( refElem3D%Lift, storage_format='ELL' )

    !----
    lcmesh3D => mesh%lcmesh_list(1)

    Nnode_v = refElem3D%Nnode_v
    allocate( gLB(refElem3D%Np), gRB(refElem3D%Np) )
    allocate( nodepos_v(Nnode_v), L_pp1(Nnode_v,refElem3D%PolyOrder_v+2) )

    nodepos_v(:) = refElem3D%x3(refElem3D%Colmask(:,1))
    L_pp1(:,:) = Polynominal_GenLegendrePoly( refElem3D%PolyOrder_v+1, nodepos_v(:) )
    
    do pZ1D=1, Nnode_v
    do p2D=1, refElem3D%Nnode_h1D**2
      p = p2D + (pZ1D-1)*refElem3D%Nnode_h1D**2
      gLB(p) = 0.5_RP * real( (-1)**mod(refElem3D%PolyOrder_v,2), RP) &
             * ( L_pp1(pZ1D,Nnode_v) - L_pp1(pZ1D,Nnode_v+1) )
      gRB(p) = 0.5_RP                                        &
             * ( L_pp1(pZ1D,Nnode_v) + L_pp1(pZ1D,Nnode_v+1) )
    end do
    end do

    !----
    allocate( velBC_ids(DOM_BND_NUM) )
    velBC_ids(domBnd_Btm_ID) = BND_TYPE_NOSLIP_ID
    velBC_ids(domBnd_Top_ID) = BND_TYPE_NOSLIP_ID
    velBC_ids(domBnd_North_ID) = BND_TYPE_PERIODIC_ID
    velBC_ids(domBnd_South_ID) = BND_TYPE_PERIODIC_ID
    velBC_ids(domBnd_East_ID) = BND_TYPE_PERIODIC_ID
    velBC_ids(domBnd_West_ID) = BND_TYPE_PERIODIC_ID

    allocate( thermBC_ids(DOM_BND_NUM) )
    thermBC_ids(domBnd_Btm_ID) = BND_TYPE_FIXVAL_ID
    thermBC_ids(domBnd_Top_ID) = BND_TYPE_FIXVAL_ID
    thermBC_ids(domBnd_North_ID) = BND_TYPE_PERIODIC_ID
    thermBC_ids(domBnd_South_ID) = BND_TYPE_PERIODIC_ID
    thermBC_ids(domBnd_East_ID) = BND_TYPE_PERIODIC_ID
    thermBC_ids(domBnd_West_ID) = BND_TYPE_PERIODIC_ID

    allocate( velBC_list(mesh%LOCAL_MESH_NUM) )
    allocate( thermBC_list(mesh%LOCAL_MESH_NUM) )

    do n=1, mesh%LOCAL_MESH_NUM
      call mesh%GetLocalMesh( n, ptr_lcmesh )
      select type (ptr_lcmesh)
      type is (LocalMesh3D)
        lcmesh3D => ptr_lcmesh
      end select

      call bnd_Init_lc( &
        VelBC_list(n),         & ! (inout)
        velBC_ids(:),          & ! (in)
        thermBC_list(n),         & ! (inout)
        thermBC_ids(:),          & ! (in)
        lcmesh3D%VMapB, mesh, lcmesh3D, lcmesh3D%refElem3D ) ! (in)
    end do    

    !--- read namelist
    BTM_BC_TYPE_HEAT = "FixedTemp"
    TOP_BC_TYPE_HEAT = "FixedTemp"

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_RBCONV_ANALYSIS_BC ,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("diag_tb_Init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("diag_tb_Init",*) 'Not appropriate names in namelist PARAM_RBCONV_ANALYSIS_BC . Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_RBCONV_ANALYSIS_BC)

    BTM_BC_HEAT_TYPEID = get_bctype_heat_id(BTM_BC_TYPE_HEAT)
    TOP_BC_HEAT_TYPEID = get_bctype_heat_id(TOP_BC_TYPE_HEAT)

    ! Flux reconstruction
    allocate( vmapM_Z(refElem3D%NfpTot,lcmesh3D%NeZ), vmapP_Z(refElem3D%NfpTot,lcmesh3D%NeZ) )
    allocate( nz_Z(refElem3D%NfpTot,lcmesh3D%NeZ) )
    call common_gen_vmapZ( vmapM_Z, vmapP_Z, nz_Z, &
                  lcmesh3D, refElem3D              )

    return
  end subroutine common_Init

!OCL SERIAL
  subroutine common_Final()
    implicit none
    !--------------------------------------------------------------------
    call Dx%Final(); call Dy%Final(); call Dz%Final()
    call Sx%Final(); call Sy%Final(); call Sz%Final()
    call Lift%Final()

    deallocate( velBC_ids, VelBC_list )

    return
  end subroutine common_Final

!OCL SERIAL
  subroutine global_horizontal_mean(field)
    use mpi
    use scale_prc, only: &
      PRC_LOCAL_COMM_WORLD    
    implicit none

    class(MeshField1D), intent(inout), target :: field
    real(RP), allocatable :: sum_tmp(:,:)
    real(RP), allocatable :: sum_out(:,:)
    integer :: ke_z

    class(LocalMesh1D), pointer :: lmesh1D
    class(ElementBase1D), pointer :: refElemV1D

    integer :: nproc
    integer :: ierr
    !-----------------------------------------

    lmesh1D => field%mesh%lcmesh_list(1)
    refElemV1D => lmesh1D%refElem1D

    allocate( sum_tmp(refElemV1D%Np,lmesh1D%Ne) )
    allocate( sum_out(refElemV1D%Np,lmesh1D%Ne) )
    do ke_z=1, lmesh1D%Ne
      sum_tmp(:,ke_z) = field%local(1)%val(:,ke_z)
    end do

    ! global sum
    call MPI_AllReduce( sum_tmp, sum_out, refElemV1D%Np * lmesh1D%Ne, &
      MPI_DOUBLE_PRECISION, MPI_SUM, PRC_LOCAL_COMM_WORLD, ierr       )
    
    call MPI_Comm_size( PRC_LOCAL_COMM_WORLD, nproc, ierr )

    do ke_z=1, lmesh1D%Ne
      field%local(1)%val(:,ke_z) = sum_out(:,ke_z) / dble( nproc )
    end do
      
    return
  end subroutine global_horizontal_mean

!OCL SERIAL  
  subroutine calc_stab_bndflx( &
    SGS_HEAT_EDDYFLX,                                  &
    PT,  DENS_hyd, DDENS,     &
    lcmesh, elem, lcmesh2D, elem2D )

    use scale_localmesh_2d, only: LocalMesh2D
    use scale_element_base, only: ElementBase2D  
    use scale_sparsemat, only: &
      SparseMat, &
      sparsemat_matmul
    implicit none

    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: SGS_HEAT_EDDYFLX(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: PT(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem%Np,lcmesh%NeA)

    integer :: ke, ke2D, ke_top, ke_btm
    integer :: ij

    logical :: cal_tend_flag

    real(RP) :: Kh_top(elem%Np,lcmesh%Ne2D)
    real(RP) :: Kh_btm(elem%Np,lcmesh%Ne2D)
    real(RP) :: DzPT(elem%Np,lcmesh%NeA)
    real(RP) :: Fz(elem%Np)
    real(RP) :: lz
    !---------------------------------------------


    lz = 1.5_RP * ( lcmesh%zmax - lcmesh%zmin ) / real( lcmesh%NeZ * elem%Nnode_v, kind=RP )

    !$omp parallel workshare
    SGS_HEAT_EDDYFLX(:,:,:) = 0.0_RP
    !$omp end parallel workshare
    
    !$omp parallel do private(ke2D, ke_top, ke_btm, Fz)
    do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
      ke_btm = ke2D
      ke_top = ke2D + (lcmesh%NeZ-1)*lcmesh2D%Ne
      Kh_btm(:,ke2D) = StabCoef_bnd * exp( - (lcmesh%zlev(:,ke_btm) / lz)**2 )
      Kh_top(:,ke2D) = StabCoef_bnd * exp( - ( ( lcmesh%zmax - lcmesh%zlev(:,ke_top) ) / lz)**2 )

      call sparsemat_matmul( Dz, PT(:,ke_btm), Fz)
      DzPT(:,ke_btm) = lcmesh%Escale(:,ke_btm,3,3) * Fz(:)
      SGS_HEAT_EDDYFLX(:,1,ke2D) = ( DENS_hyd(:,ke_btm) + DDENS(:,ke_btm) ) * Kh_btm(:,ke2D) * DzPT(:,ke_btm)

      call sparsemat_matmul( Dz, PT(:,ke_top), Fz)
      DzPT(:,ke_top) = lcmesh%Escale(:,ke_top,3,3) * Fz(:)
      SGS_HEAT_EDDYFLX(:,lcmesh%NeZ,ke2D) = ( DENS_hyd(:,ke_top) + DDENS(:,ke_top) ) * Kh_top(:,ke2D) * DzPT(:,ke_top)
    end do      

    return
  end subroutine calc_stab_bndflx

!OCL SERIAL
  subroutine get_reconstructed_flux( reconst_flux, &
    flux, del_flux, lcmesh, elem                   )

    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem    
    real(RP), intent(out) :: reconst_flux(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: flux(elem%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: del_flux(elem%NfpTot,lcmesh%Ne)

    integer :: ke, ke2D, ke_z
    integer :: p, p2D, p_z
    integer :: ps_L, pe_L
    integer :: ps_R, pe_R

    real(RP) :: lift_L(elem%Np), lift_R(elem%Np)
    !-----------------------------------------------

    ps_L = elem%Nfp_h * elem%Nfaces_h + 1
    pe_L = ps_L + elem%Nfp_v - 1
    ps_R = pe_L + 1
    pe_R = ps_R + elem%Nfp_v - 1

    !$omp parallel do collapse(2) &
    !$omp private(ke, p, p2D, p_z, lift_L, lift_R)
    do ke2D=1,lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D

      do p_z=1, elem%Nnode_v
      do p2D=1, elem%Nfp_v
        p = p2D + (p_z-1)*elem%Nfp_v
        lift_L(p) = del_flux(ps_L+p2D-1,ke) 
        lift_R(p) = del_flux(ps_R+p2D-1,ke) 
      end do
      end do
      reconst_flux(:,ke) = flux(:,ke_z,ke2D)     &
        + lift_L(:) * gLB(:) + lift_R(:) * gRB(:) 
    end do
    end do

    return
  end subroutine get_reconstructed_flux

!OCL SERIAL
  subroutine get_del_flux_cent_heat( del_flux,    &
    flux, lcmesh, elem )

    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot,lcmesh%Ne)
    real(RP), intent(in) :: flux(elem%Np*lcmesh%NeZ,lcmesh%Ne2D)

    real(RP) :: sfc_flux(elem%Nfp_v,lcmesh%Ne2D)
    real(RP) :: top_flux(elem%Nfp_v,lcmesh%Ne2D)
    !--------------------------------

    sfc_flux(:,:) = - BTM_FIXED_HEAT_FLUX / CpDry
    top_flux(:,:) = - TOP_FIXED_HEAT_FLUX / CpDry
    call get_del_flux_cent( del_flux,    &
      flux, lcmesh, elem, TOP_BC_HEAT_TYPEID, BTM_BC_HEAT_TYPEID,  &
      sfc_flux, top_flux )

    return
  end subroutine get_del_flux_cent_heat

!OCL SERIAL
  subroutine get_del_flux_cent( del_flux,    &
    flux, lcmesh, elem, bndtype_top_id, bndtype_btm_id,  &
    sfc_flux, top_flux )

    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem
    real(RP), intent(out) :: del_flux(elem%NfpTot,lcmesh%Ne)
    real(RP), intent(in) :: flux(elem%Np*lcmesh%NeZ,lcmesh%Ne2D)
    integer, intent(in) :: bndtype_top_id ! mirror: 1, specified flux: 2
    integer, intent(in) :: bndtype_btm_id ! mirror: 1, specified flux: 2
    real(RP), optional, intent(in) :: sfc_flux(elem%Nfp_v,lcmesh%Ne2D)
    real(RP), optional, intent(in) :: top_flux(elem%Nfp_v,lcmesh%Ne2D)

    integer :: ke2D, ke_z, ke
    integer :: p2D, p
    integer :: iM(elem%NfpTot), iP(elem%NfpTot)
    !-----------------------------------------

    !$omp parallel do collapse(2) private(ke, iM, iP)
    do ke_z=1, lcmesh%NeZ
    do ke2D=1, lcmesh%Ne2D
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      iM(:) = vmapM_Z(:,ke_z); iP(:) = vmapP_Z(:,ke_z)

      del_flux(:,ke) = &
        0.5_RP * ( flux(iP(:),ke2D) - flux(iM(:),ke2D) )
    end do
    end do

    !omp parallel do private(Ke2D, ke, iM, p2D, p)
    do ke2D=1, lcmesh%Ne2D
      !---
      ke = ke2D
      iM(:) = vmapM_Z(:,1)
      if (bndtype_btm_id == 1) then
        do p2D=1, elem%Nfp_v
          p = elem%Nfp_h * elem%Nfaces_h + p2D
          del_flux(p,ke) = 0.0_RP
        end do
      else if (bndtype_btm_id == 2) then
        do p2D=1, elem%Nfp_v
          p = elem%Nfp_h * elem%Nfaces_h + p2D
          del_flux(p,ke) = sfc_flux(p2D,ke2D) - flux(iM(p),ke2D)
        end do
      end if

      !--
      ke = ke2D + (lcmesh%NeZ-1)*lcmesh%Ne2D
      iM(:) = vmapM_Z(:,lcmesh%NeZ)
      if (bndtype_top_id == 1) then
        do p2D=1, elem%Nfp_v
          p = elem%Nfp_h * elem%Nfaces_h + elem%Nfp_v + p2D
          del_flux(p,ke) = 0.0_RP
        end do
      else if (bndtype_top_id == 2) then
        do p2D=1, elem%Nfp_v
          p = elem%Nfp_h * elem%Nfaces_h + elem%Nfp_v + p2D
          del_flux(p,ke) = top_flux(p2D,ke2D) - flux(iM(p),ke2D)
        end do
      end if
    end do

    return
  end subroutine get_del_flux_cent

!OCL SERIAL
  subroutine inquire_bound_flag( & 
    is_bound,                                & ! (out)
    domID, vmapM, vmapP, vmapB, lmesh, elem  ) ! (in)

    use scale_mesh_bndinfo, only: &
      BND_TYPE_SLIP_ID, BND_TYPE_NOSLIP_ID, &
      BND_TYPE_ADIABAT_ID
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem    
    logical, intent(out) :: is_bound(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: domID
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapB(:)

    integer :: i, i_, iM, iP
    !-----------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iP = vmapP(i)
      i_ = iP - elem%Np*lmesh%NeE
      is_bound(i) = .false.

      if (i_ > 0) then  
        iM = vmapM(i)
        if ( VelBC_list(domID)%list(i_) == BND_TYPE_SLIP_ID ) then
          is_bound(i) = .true.
        else if ( VelBC_list(domID)%list(i_) == BND_TYPE_NOSLIP_ID ) then          
          is_bound(i) = .true.
        end if

      end if
    end do

    return
  end subroutine inquire_bound_flag
    
!-- private ----------------------

!OCL SERIAL
  subroutine bnd_Init_lc(   &
    velBCInfo,               &  ! (inout)
    velBCIDs,                &  ! (in)
    thermBCInfo,               &  ! (inout)
    thermBCIDs,                &  ! (in)
    vmapB, mesh, lmesh, elem )  ! (in)

    use scale_mesh_bndinfo, only: BND_TYPE_NOSPEC_ID
    implicit none
    
    type(MeshBndInfo), intent(inout) :: velBCInfo
    integer, intent(in) :: velBCIDs(DOM_BND_NUM)
    type(MeshBndInfo), intent(inout) :: thermBCInfo
    integer, intent(in) :: thermBCIDs(DOM_BND_NUM)
    class(MeshBase3D), intent(in) :: mesh
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(in) :: vmapB(:)

    integer :: tileID
    integer :: dom_bnd_sizes(DOM_BND_NUM)
    integer :: bnd_buf_size
    integer :: b, is_, ie_
    !-----------------------------------------------
    
    dom_bnd_sizes(:) = &
        elem%Nfp_h*lmesh%NeZ*(/ lmesh%NeX, lmesh%NeY, lmesh%NeX, lmesh%NeY, 0, 0 /) &
      + elem%Nfp_v*lmesh%NeX*lmesh%NeY*(/ 0, 0, 0, 0, 1, 1 /)
    bnd_buf_size = sum(dom_bnd_sizes)
    
    call velBCInfo%Init( bnd_buf_size )
    call velBCInfo%Set(1, bnd_buf_size, BND_TYPE_NOSPEC_ID)

    call thermBCInfo%Init( bnd_buf_size )
    call thermBCInfo%Set(1, bnd_buf_size, BND_TYPE_NOSPEC_ID)

    tileID = lmesh%tileID
    is_ = 1
    do b=1, DOM_BND_NUM
      ie_ = is_ + dom_bnd_sizes(b) - 1
      if ( mesh%tileID_globalMap(b,tileID) == tileID      &
           .and. mesh%tileFaceID_globalMap(b,tileID) == b ) then
        call velBCInfo%Set(is_, ie_, velBCIDs(b))
        call thermBCInfo%Set(is_, ie_, thermBCIDs(b))
      end if
      is_ = ie_ + 1
    end do

    return
  end  subroutine bnd_Init_lc

!OCL SERIAL
  subroutine common_gen_vmapZ( &
    vmapM, vmapP, nz,          & ! (out)
    lmesh, elem   ) ! (in)
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(elementbase3D), intent(in) :: elem
    integer, intent(out) :: vmapM(elem%NfpTot,lmesh%NeZ)
    integer, intent(out) :: vmapP(elem%NfpTot,lmesh%NeZ)    
    integer, intent(out) :: nz(elem%NfpTot,lmesh%NeZ)    

    integer :: ke_z
    integer :: f
    integer :: vs, ve
    !-------------------------------------------------------

    !$omp parallel private(f, vs, ve)
    !$omp do
    do ke_z=1, lmesh%NeZ
      do f=1, elem%Nfaces_h
        vs = 1 + (f-1)*elem%Nfp_h
        ve = vs + elem%Nfp_h - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_h(:,f) + (ke_z-1)*elem%Np
      end do
      do f=1, elem%Nfaces_v
        vs = elem%Nfp_h*elem%Nfaces_h + 1 + (f-1)*elem%Nfp_v
        ve = vs + elem%Nfp_v - 1
        vmapM(vs:ve,ke_z) = elem%Fmask_v(:,f) + (ke_z-1)*elem%Np
      end do
      vmapP(:,ke_z) = vmapM(:,ke_z)
      nz(:,ke_z) = lmesh%normal_fn(:,1+(ke_z-1)*lmesh%Ne2D,3)
    end do
    !$omp do
    do ke_z=1, lmesh%NeZ
      vs = elem%Nfp_h*elem%Nfaces_h + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z > 1) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,2) + (ke_z-2)*elem%Np

      vs = elem%Nfp_h*elem%Nfaces_h + elem%Nfp_v + 1
      ve = vs + elem%Nfp_v - 1
      if (ke_z < lmesh%NeZ) &
        vmapP(vs:ve,ke_z) = elem%Fmask_v(:,1) + ke_z*elem%Np
    end do
    !$omp end parallel

    return
  end subroutine common_gen_vmapZ

  function get_bctype_heat_id(bctype) result(type_id)
    character(*), intent(in) :: bctype
    integer :: type_id
    !------------------
    select case( trim(bctype) )
    case('FixedFlux')
      type_id = BC_HEAT_FIXED_FLUX
    case('FixedTemp')
      type_id = BC_HEAT_FIXED_TEMP
    case default
      LOG_ERROR("RBconv_get_bctype_heat_id",*) 'Not appropriate boundary condition. Check!', trim(bctype)
      call PRC_abort
    end select
    return
  end function get_bctype_heat_id  

end module mod_common