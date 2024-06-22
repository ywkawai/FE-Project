#include "scaleFElib.h"
module mod_diag_tb
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,    &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00

  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase1D, ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_1d, only: LocalMesh1D  
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base1d, only: MeshBase1D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField1D, MeshField3D

  use scale_file_base_meshfield, only: &
    FILE_base_meshfield
  use scale_mesh_base3d, only: &
    DIMTYPEID_XYZT => MeshBase3D_DIMTYPEID_XYZT
  use scale_mesh_base1d, only: &
    DIMTYPEID_ZT => MeshBase1D_DIMTYPEID_XT
  
  use scale_atm_phy_tb_dgm_smg, only: &
    atm_phy_tb_dgm_smg_Init,    &
    atm_phy_tb_dgm_smg_Final,   &
    atm_phy_tb_dgm_smg_cal_grad
  
  use mod_common, only: &
    global_horizontal_mean,     &
    Dx, Dy, Dz, Sx, Sy, Sz, Lift

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: diag_tb_Init
  public :: diag_tb_Final
  public :: diag_tb_process
  
  type(MeshField3D), target :: T11, T12, T13, T21, T22, T23, T31, T32, T33
  type(MeshField3D), target :: DF1, DF2, DF3
  type(MeshField3D), target :: Nu, Kh
  type(MeshField3D), target :: TKE

  type(MeshField1D) :: SGS_HEAT_EDDYFLX_V1D
  type(MeshField1D) :: SGS_MOMZ_EDDYFLX_V1D
  type(MeshField1D) :: Nu_V1D
  type(MeshField1D) :: Kh_V1D

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  type(FILE_base_meshfield) :: out_file

  integer, parameter :: DIAGVID_TB_NU   = 1
  integer, parameter :: DIAGVID_TB_KH   = 2
  integer, parameter :: DIAGVID_TB_Ri   = 3
  integer, parameter :: DIAGVID_TB_T33  = 4
  integer, parameter :: DIAGVID_TB_DF3  = 5


  type(FILE_base_meshfield) :: out_file_V1D
  integer, parameter :: DIAGVID_TB_SGS_HEAT_EDDYFLUX_V1D   = 1
  integer, parameter :: DIAGVID_TB_SGS_MOMZ_EDDYFLUX_V1D   = 2
  integer, parameter :: DIAGVID_TB_Nu_V1D            = 3
  integer, parameter :: DIAGVID_TB_Kh_V1D            = 4

  logical :: ismaster
  integer :: NLocalMeshPerPrc


  abstract interface
    subroutine atm_phy_tb_cal_grad( &
      T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (out)
      DF1, DF2, DF3,                                              & ! (out)
      TKE, Nu, Kh,                                                & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,    & ! (in)
      PRES, PT,                                                   & ! (in)
      Dx, Dy, Dz, Sx, Sy, Sz, Lift, lmesh, elem, lmesh2D, elem2D, & ! (in)
      is_bound                                                    ) ! (in)
      import RP
      import LocalMesh3D
      import LocalMesh2D
      import ElementBase3D
      import ElementBase2D
      import SparseMat
      implicit none

      class(LocalMesh3D), intent(in) :: lmesh
      class(ElementBase3D), intent(in) :: elem
      class(LocalMesh2D), intent(in) :: lmesh2D
      class(ElementBase2D), intent(in) :: elem2D
      real(RP), intent(out) :: T11(elem%Np,lmesh%NeA), T12(elem%Np,lmesh%NeA), T13(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: T21(elem%Np,lmesh%NeA), T22(elem%Np,lmesh%NeA), T23(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: T31(elem%Np,lmesh%NeA), T32(elem%Np,lmesh%NeA), T33(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: DF1(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: DF2(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: DF3(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: TKE(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: Nu(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: Kh(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PT(elem%Np,lmesh%NeA)
      type(SparseMat), intent(in) :: Dx, Dy, Dz
      type(SparseMat), intent(in) :: Sx, Sy, Sz
      type(SparseMat), intent(in) :: Lift
      logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
    end subroutine atm_phy_tb_cal_grad
  end interface
  procedure (atm_phy_tb_cal_grad), pointer :: tbsolver_cal_grad => null()

  abstract interface    
    subroutine atm_phy_tb_final()
      implicit none
    end subroutine atm_phy_tb_final
  end interface 
  procedure (atm_phy_tb_final), pointer :: tbsolver_final => null()

contains
!OCL SERIAL
  subroutine diag_tb_Init( tb_scheme, mesh, meshV1D, &
    out_filebase_tb, dtype, out_tintrv,   &
    myrank, is_master_, NLocalMeshPerPrc_ )

    use scale_atm_phy_tb_dgm_dns, only: &
      atm_phy_tb_dgm_dns_Init, &
      atm_phy_tb_dgm_dns_cal_grad, &
      atm_phy_tb_dgm_dns_Final
    use scale_atm_phy_tb_dgm_smg, only: &
      atm_phy_tb_dgm_smg_Init, &
      atm_phy_tb_dgm_smg_cal_grad, &
      atm_phy_tb_dgm_smg_Final
    implicit none    
    character(len=*), intent(in) :: tb_scheme
    class(MeshCubeDom3D), intent(in) :: mesh
    class(MeshBase1D), intent(in) :: meshV1D
    character(len=*), intent(in) :: out_filebase_tb
    character(len=*), intent(in) :: dtype
    real(RP), intent(in) :: out_tintrv
    integer, intent(in) :: myrank
    logical, intent(in) :: is_master_
    integer, intent(in) :: NLocalMeshPerPrc_

    logical fileexist
    !--------------------------------------------------------------------

    select case(trim(tb_scheme))
    case ("DNS")
      call atm_phy_tb_dgm_dns_Init( mesh )
      tbsolver_cal_grad => atm_phy_tb_dgm_dns_cal_grad
      tbsolver_final => atm_phy_tb_dgm_dns_Final
    case ("SMAGORINSKY")
      call atm_phy_tb_dgm_smg_Init( mesh )
      tbsolver_cal_grad => atm_phy_tb_dgm_smg_cal_grad
      tbsolver_final => atm_phy_tb_dgm_smg_Final
    case default
      LOG_ERROR("diag_tb_Init",*) 'Invalid tb_scheme. Check!'
      call PRC_abort      
    end select

    call T11%Init("T11", "kg.m-3.m2.s-2", mesh)
    call T12%Init("T12", "kg.m-3.m2.s-2", mesh)
    call T13%Init("T13", "kg.m-3.m2.s-2", mesh)

    call T21%Init("T21", "kg.m-3.m2.s-2", mesh)
    call T22%Init("T22", "kg.m-3.m2.s-2", mesh)
    call T23%Init("T23", "kg.m-3.m2.s-2", mesh)

    call T31%Init("T31", "kg.m-3.m2.s-2", mesh)
    call T32%Init("T32", "kg.m-3.m2.s-2", mesh)
    call T33%Init("T33", "kg.m-3.m2.s-2", mesh)

    call TKE%Init("TKE", "m2/s2", mesh)

    call DF1%Init("DF1", "kg.m-3.m2.s-1.K", mesh)
    call DF2%Init("DF2", "kg.m-3.m2.s-1.K", mesh)
    call DF3%Init("DF3", "kg.m-3.m2.s-1.K", mesh)

    call Nu%Init( "Nu", "m2/s", mesh )
    call Kh%Init( "Kh", "m2/s", mesh )

    call SGS_HEAT_EDDYFLX_V1D%Init("SGS_HEAT_EDDYFLX", "J.m-2.s-1", meshV1D)
    call SGS_MOMZ_EDDYFLX_V1D%Init("SGS_MOMZ_EDDYFLX", "m/s.kg.m-2.s-1", meshV1D)
    call Nu_V1D%Init( "Nu", "m2/s", meshV1D )
    call Kh_V1D%Init( "Kh", "m2/s", meshV1D )

    !
    ! call out_file%Init(12, mesh3D=mesh)
    ! call out_file%Create( out_filebase_tb, "turbulence process", &
    !  dtype, fileexist, myrank=myrank )
    ! call out_file%Def_Var( Nu, "Eddy viscosity", DIAGVID_TB_NU, DIMTYPEID_XYZT, &
    !   dtype, timeinv=out_tintrv )
    ! call out_file%Def_Var( Kh, "Eddy diffusivity", DIAGVID_TB_KH, DIMTYPEID_XYZT, &
    !   dtype, timeinv=out_tintrv )
    ! call out_file%Def_Var( TKE, "TKE", DIAGVID_TB_TKEsgs, DIMTYPEID_XYZT, &
    !   dtype, timeinv=out_tintrv )
    ! call out_file%Def_Var( dPTdz, "DF3", DIAGVID_TB_dPTdz, DIMTYPEID_XYZT, &
    !   dtype, timeinv=out_tintrv )

    ! call out_file%End_def()

    ismaster = is_master_
    NLocalMeshPerPrc = NLocalMeshPerPrc_
    if (ismaster) then
      call out_file_V1D%Init(4, mesh1D=meshV1D)
      call out_file_V1D%Create( trim(out_filebase_tb)//"_V1D", 'PBL turbulence analysis', &
        dtype, fileexist, myrank=myrank )
      call out_file_V1D%Def_Var(SGS_HEAT_EDDYFLX_V1D, "horizontal averaged SGS eddy heat flux", &
        DIAGVID_TB_SGS_HEAT_EDDYFLUX_V1D, DIMTYPEID_ZT, dtype, timeinv=out_tintrv )
      call out_file_V1D%Def_Var(SGS_MOMZ_EDDYFLX_V1D, "horizontal averaged SGS eddy momentum flux in z-direction", &
        DIAGVID_TB_SGS_MOMZ_EDDYFLUX_V1D, DIMTYPEID_ZT, dtype, timeinv=out_tintrv )      
      call out_file_V1D%Def_Var(Nu_V1D, "Horizontally averaged eddy viscosity", &
        DIAGVID_TB_Nu_V1D, DIMTYPEID_ZT, dtype, timeinv=out_tintrv )      
      call out_file_V1D%Def_Var(Kh_V1D, "Horizontally averaged eddy diffusivity", &
        DIAGVID_TB_Kh_V1D, DIMTYPEID_ZT, dtype, timeinv=out_tintrv )      
      call out_file_V1D%End_def()
    end if

    return
  end subroutine diag_tb_Init

!OCL SERIAL
  subroutine diag_tb_process( &
    istep, start_time, end_time,     &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,  &
    DENS_hyd, PRES_hyd, PRES, PT,    &
    DENS_v1D,                        &
    dyn_bnd, mesh, meshV1D         )
  
    use scale_atm_phy_tb_dgm_smg, only: &
      atm_phy_tb_dgm_smg_cal_grad
    use scale_atm_dyn_dgm_bnd, only: &
      AtmDynBnd
    use mod_common, only: &
      inquire_bound_flag
    implicit none

    integer, intent(in) :: istep
    real(RP), intent(in) :: start_time, end_time
    class(MeshField3D), intent(inout) :: DDENS
    class(MeshField3D), intent(inout) :: MOMX
    class(MeshField3D), intent(inout) :: MOMY
    class(MeshField3D), intent(inout) :: MOMZ
    class(MeshField3D), intent(inout) :: DRHOT
    class(MeshField3D), intent(in) :: DENS_hyd
    class(MeshField3D), intent(in) :: PRES_hyd
    class(MeshField3D), intent(inout) :: PRES
    class(MeshField3D), intent(inout) :: PT
    class(MeshField1D), intent(in) :: DENS_v1D
    class(AtmDynBnd), intent(in) :: dyn_bnd
    class(MeshCubeDom3D), intent(in), target :: mesh
    class(MeshBase1D), intent(in), target :: meshV1D

    integer :: n
    class(LocalMesh3D), pointer :: lmesh
    class(LocalMesh1D), pointer :: lmeshV1D

    logical, allocatable :: is_bound(:,:)

    real(RP), allocatable :: Rtot(:,:)
    real(RP), allocatable :: CPtot(:,:)
    !--------------------------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(n)

      allocate( is_bound(lmesh%refElem%NfpTot,lmesh%Ne) )
      call inquire_bound_flag( is_bound,              &
        n, lmesh%VMapM, lmesh%VMapP, lmesh%VMapB,     &
        lmesh, lmesh%refElem3D                        )

      allocate( Rtot(lmesh%refElem3D%Np,lmesh%NeA) )
      allocate( CPtot(lmesh%refElem3D%Np,lmesh%NeA) )
      !$omp parallel workshare
      Rtot(:,:) = Rdry
      CPtot(:,:) = CPdry
      !$omp end parallel workshare

      call dyn_bnd%ApplyBC_Grad_TBVARS_lc( n, &
        DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val,         &
        MOMZ%local(n)%val,            & 
        PT%local(n)%val, PRES%local(n)%val,                               &
        DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                                & 
        Rtot, CPtot, lmesh%Gsqrt, lmesh%GsqrtH, lmesh%GIJ(:,:,1,1), lmesh%GIJ(:,:,1,2), lmesh%GIJ(:,:,2,2), &
        lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),  &
        lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), &
        lmesh%VMapM, lmesh%VMapP, lmesh%VMapB, lmesh, lmesh%refElem3D,          &
        lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D )

      call tbsolver_cal_grad( &
        T11%local(n)%val, T12%local(n)%val, T13%local(n)%val,             &
        T21%local(n)%val, T22%local(n)%val, T23%local(n)%val,             &
        T31%local(n)%val, T32%local(n)%val, T33%local(n)%val,             &
        DF1%local(n)%val, DF2%local(n)%val, DF3%local(n)%val,             &
        TKE%local(n)%val,                                                 &
        Nu%local(n)%val, Kh%local(n)%val,                                 &
        DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val,         &
        MOMZ%local(n)%val,  DRHOT%local(n)%val,                           &
        DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                     & 
        PRES%local(n)%val, PT%local(n)%val,                               &
        Dx, Dy, Dz, Sx, Sy, Sz, Lift,                                     &
        lmesh, lmesh%refElem3D, lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D, &
        is_bound )    
      
      deallocate( is_bound, Rtot, CPtot )
    end do
    
    do n=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(n)
      lmeshV1D => meshV1D%lcmesh_list(1)
      call cal_auxvars( &
        SGS_HEAT_EDDYFLX_V1D%local(n)%val,                       &
        SGS_MOMZ_EDDYFLX_V1D%local(n)%val,                       &
        Nu_V1D%local(n)%val,                                     &
        Kh_V1D%local(n)%val,                                     &
        T33%local(n)%val, DF3%local(n)%val,                      &
        Nu%local(n)%val, Kh%local(n)%val,                        &
        DDENS%local(n)%val, PT%local(n)%val,                     &
        DENS_hyd%local(n)%val,                                   &
        DENS_v1D%local(n)%val,                                   &
        n, lmesh, lmesh%refElem3D, lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D, &
        lmeshV1D, lmeshV1D%refElem1D )
    end do

    call global_horizontal_mean( SGS_HEAT_EDDYFLX_V1D )
    call global_horizontal_mean( SGS_MOMZ_EDDYFLX_V1D )

    ! call out_file%Write_var3D( DIAGVID_TB_NU, Nu, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_KH, Kh, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_Ri, Ri, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_Sabs, Sabs, start_time, end_time )

    ! call out_file%Write_var3D( DIAGVID_TB_S11, S11, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_S12, S12, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_S22, S22, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_S23, S23, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_S31, S31, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_S33, S33, start_time, end_time )
    ! call out_file%Write_var3D( DIAGVID_TB_TKEsgs, TKE, start_time, end_time )

    ! call out_file%Write_var3D( DIAGVID_TB_dPTdz, dPTdz, start_time, end_time )

    if (ismaster) then
      call out_file_V1D%Write_var1D( DIAGVID_TB_SGS_HEAT_EDDYFLUX_V1D, &
        SGS_HEAT_EDDYFLX_V1D, start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_TB_SGS_MOMZ_EDDYFLUX_V1D, &
        SGS_MOMZ_EDDYFLX_V1D, start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_TB_Nu_V1D, &
        Nu_V1D, start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_TB_Kh_V1D, &
        Kh_V1D, start_time, end_time )
    end if

    return
  end subroutine diag_tb_process

!OCL SERIAL
  subroutine diag_tb_Final()
    implicit none   
    !--------------------------------------------------------------------

    call tbsolver_final()

    call TKE%Final()

    call Nu%Final()
    call Kh%Final()

    call SGS_HEAT_EDDYFLX_V1D%Final()
    call SGS_MOMZ_EDDYFLX_V1D%Final()

    ! call out_file%Close()
    ! call out_file%Final()
    if (ismaster) then
      call out_file_V1D%Close()
      call out_file_V1D%Final()
    end if

    return
  end subroutine diag_tb_Final

!- private --------------------------------------
!OCL SERIAL
  subroutine cal_auxvars( &
    SGS_HEAT_EDDYFLX_V1D_,              &
    SGS_MOMZ_EDDYFLX_V1D_,              &
    Nu_V1D_, Kh_V1D_,                   &
    T33_, DF3_, Nu_, Kh_,                              & 
    DDENS_, PT_, DENS_hyd_,                            &
    DENS_V1D_,                                         &
    domID, lmesh, elem3D, lmeshH2D, elemH2D, lmeshV1D, elem1D )

    use scale_const, only: &
      EPS => CONST_EPS,  &
      Grav => CONST_GRAV, &
      CpDry => CONST_CPdry
    use mod_common, only: &
      get_del_flux_cent,     &
      get_reconstructed_flux
    
    implicit none   
    integer, intent(in) :: domID
    class(LocalMesh3D), intent(in), target :: lmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh2D), intent(in) :: lmeshH2D
    class(ElementBase2D), intent(in) :: elemH2D
    class(LocalMesh1D), intent(in) :: lmeshV1D
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(out) :: SGS_HEAT_EDDYFLX_V1D_(elem1D%Np,lmeshV1D%NeA)
    real(RP), intent(out) :: SGS_MOMZ_EDDYFLX_V1D_(elem1D%Np,lmeshV1D%NeA)
    real(RP), intent(out) :: Nu_V1D_(elem1D%Np,lmeshV1D%NeA)
    real(RP), intent(out) :: Kh_V1D_(elem1D%Np,lmeshV1D%NeA)
    real(RP), intent(in) :: T33_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: DF3_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: Nu_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: Kh_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: PT_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_v1D_(elem1D%Np,lmeshV1D%NeA)

    class(LocalMesh2D), pointer :: lmesh2D

    integer :: ke, ke2D
    integer :: ke_x, ke_y, ke_z
    integer :: p
    real(RP) :: tmp(elem3D%Np)
    real(RP) :: harea

    real(RP) :: DENS_lc(elem3D%Np)
    integer :: hSliceID(elem3D%Nnode_h1D**2)
    real(RP) :: int_w(elem3D%Nnode_h1D**2)

    real(RP) :: SGS_HEAT_EDDYFLX_z(elem3D%Np,lmesh%NeZ,lmesh%Ne2D)
    real(RP) :: SGS_HEAT_EDDYFLX_reconst(elem3D%Np,lmesh%NeA)
    real(RP) :: del_flux_sgs_heat_eddyflux(elem3D%NfpTot,lmesh%Ne)

    real(RP) :: SGS_MOMZ_EDDYFLX_z(elem3D%Np,lmesh%NeZ,lmesh%Ne2D)
    real(RP) :: SGS_MOMZ_EDDYFLX_reconst(elem3D%Np,lmesh%NeA)
    real(RP) :: del_flux_sgs_momz_eddyflux(elem3D%NfpTot,lmesh%Ne)

    real(RP) :: sfc_heat_flux(elem3D%Nfp_v,lmesh%Ne2D)
    real(RP) :: sfc_momz_flux(elem3D%Nfp_v,lmesh%Ne2D)
    !--------------------------------------------------------------------

    lmesh2D => lmesh%lcmesh2D

    if (domID==1) then
      !$omp parallel do
      do ke_z=1, lmesh%NeZ
        SGS_HEAT_EDDYFLX_V1D_(:,ke_z) = 0.0_RP
        SGS_MOMZ_EDDYFLX_V1D_(:,ke_z) = 0.0_RP
        Nu_V1D_(:,ke_z) = 0.0_RP
        Kh_V1D_(:,ke_z) = 0.0_RP
      end do
      harea = 0.0_RP
    end if

    !$omp parallel private(ke, ke2D, tmp)
    ! !$omp do 
    ! do ke=lmesh%NeS, lmesh%NeE
    !   tmp(:) = 2.0_RP * ( S11_(:,ke)**2 + S22_(:,ke)**2 + S33_(:,ke)**2 ) &
    !          + 4.0_RP * ( S31_(:,ke)**2 + S12_(:,ke)**2 + S23_(:,ke)**2 )
    !   Sabs_(:,ke) = sqrt(tmp(:))

    !   Ri_(:,ke) = Grav / PT_(:,ke) * dPTdz_(:,ke) / max( Sabs_(:,ke)**2, EPS )
    ! end do
    !$omp do
    do ke2D=lmesh2D%NeS, lmesh2D%NeE
      sfc_heat_flux(:,ke2D) = 0.0_RP
      sfc_momz_flux(:,ke2D) = 0.0_RP
    end do
    !$omp end parallel
    !$omp parallel do private( ke_x, ke_y, p, ke2D, ke, DENS_lc )
    do ke_z=1, lmesh%NeZ
      do ke_y=1, lmesh%NeY
      do ke_x=1, lmesh%NeX
        ke2D = ke_x + (ke_y-1)*lmesh%NeX 
        ke = ke2D + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_lc(:) = DENS_hyd_(:,ke) + DDENS_(:,ke)

        SGS_HEAT_EDDYFLX_z(:,ke_z,ke2D) = DF3_(:,ke)
        SGS_MOMZ_EDDYFLX_z(:,ke_z,ke2D) = T33_(:,ke)        
      end do    
      end do  
    end do    

    call get_del_flux_cent( del_flux_sgs_heat_eddyflux, &
      SGS_HEAT_EDDYFLX_z, lmesh, elem3D, 1, 1  )
    call get_reconstructed_flux( SGS_HEAT_EDDYFLX_reconst, & 
      SGS_HEAT_EDDYFLX_z, del_flux_sgs_heat_eddyflux, lmesh, elem3D    )

    call get_del_flux_cent( del_flux_sgs_momz_eddyflux, &
      SGS_MOMZ_EDDYFLX_z, lmesh, elem3D, 1, 1  )
    call get_reconstructed_flux( SGS_MOMZ_EDDYFLX_reconst, & 
      SGS_MOMZ_EDDYFLX_z, del_flux_sgs_momz_eddyflux, lmesh, elem3D    )

    !$omp parallel do private( ke_x, ke_y, p, ke2D, ke, hSliceID, &
    !$omp int_w, DENS_lc )
    do ke_z=1, lmesh%NeZ
      do ke_y=1, lmesh%NeY
      do ke_x=1, lmesh%NeX
        ke2D = ke_x + (ke_y-1)*lmesh%NeX 
        ke = ke2D + (ke_z-1)*lmesh%NeX*lmesh%NeY
        DENS_lc(:) = DENS_hyd_(:,ke) + DDENS_(:,ke)

        do p=1, elem3D%Nnode_v
          hSliceID(:) = elem3D%Hslice(:,p)
          int_w(:) = elemH2D%IntWeight_lgl(:) * lmeshH2D%J(:,ke2D)

          SGS_HEAT_EDDYFLX_V1D_(p,ke_z) = SGS_HEAT_EDDYFLX_V1D_(p,ke_z)       &
                            + sum( int_w(:) * SGS_HEAT_EDDYFLX_reconst(hSliceID(:),ke)    )
          
          SGS_MOMZ_EDDYFLX_V1D_(p,ke_z) = SGS_MOMZ_EDDYFLX_V1D_(p,ke_z) &
                            + sum( int_w(:) * SGS_MOMZ_EDDYFLX_reconst(hSliceID(:),ke)    )

          Nu_V1D_(p,ke_z) = Nu_V1D_(p,ke_z) &
                            + sum( int_w(:) * Nu_(hSliceID(:),ke)    )
          Kh_V1D_(p,ke_z) = Kh_V1D_(p,ke_z) &
                            + sum( int_w(:) * Kh_(hSliceID(:),ke)    )
        end do
      end do    
      end do  
    end do    
    hSliceID(:) = elem3D%Hslice(:,1)
    do ke_y=1, lmesh%NeY
    do ke_x=1, lmesh%NeX
      ke2D = ke_x + (ke_y-1)*lmesh%NeX
      harea = harea + sum( elemH2D%IntWeight_lgl(:) * lmeshH2D%J(:,ke2D) )
    end do    
    end do  

    if ( domID == NLocalMeshPerPrc ) then
      !$omp parallel do
      do ke_z=1, lmesh%NeZ
        SGS_HEAT_EDDYFLX_V1D_(:,ke_z) = CpDry * SGS_HEAT_EDDYFLX_V1D_(:,ke_z) / harea
        SGS_MOMZ_EDDYFLX_V1D_(:,ke_z) = SGS_MOMZ_EDDYFLX_V1D_(:,ke_z) / DENS_v1D_(:,ke_z) / harea
        Nu_V1D_(:,ke_z) = Nu_V1D_(:,ke_z) / harea
        Kh_V1D_(:,ke_z) = Kh_V1D_(:,ke_z) / harea
      end do
    end if

    return
  end subroutine cal_auxvars
end module mod_diag_tb
