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
  
  type(MeshField3D), target :: S11, S12, S22, S23, S31, S33
  type(MeshField3D), target :: TKE
  type(MeshField3D), target :: dPTdx, dPTdy, dPTdz
  type(MeshField3D), target :: Nu
  type(MeshField3D), target :: Kh
  type(MeshField3D), target :: Sabs
  type(MeshField3D), target :: Ri

  type(MeshField1D) :: SGS_HEAT_EDDYFLX_V1D
  type(MeshField1D) :: SGS_MOMZ_EDDYFLX_V1D

  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  type(FILE_base_meshfield) :: out_file

  integer, parameter :: DIAGVID_TB_NU   = 1
  integer, parameter :: DIAGVID_TB_KH   = 2
  integer, parameter :: DIAGVID_TB_Ri   = 3
  integer, parameter :: DIAGVID_TB_Sabs = 4
  integer, parameter :: DIAGVID_TB_S11  = 5
  integer, parameter :: DIAGVID_TB_S12  = 6
  integer, parameter :: DIAGVID_TB_S22  = 7
  integer, parameter :: DIAGVID_TB_S23  = 8
  integer, parameter :: DIAGVID_TB_S31  = 9
  integer, parameter :: DIAGVID_TB_S33  = 10
  integer, parameter :: DIAGVID_TB_TKEsgs = 11
  integer, parameter :: DIAGVID_TB_dPTdz  = 12

  type(FILE_base_meshfield) :: out_file_V1D
  integer, parameter :: DIAGVID_TB_SGS_HEAT_EDDYFLUX_V1D   = 1
  integer, parameter :: DIAGVID_TB_SGS_MOMZ_EDDYFLUX_V1D   = 2

  logical :: ismaster
  integer :: NLocalMeshPerPrc


contains
!OCL SERIAL
  subroutine diag_tb_Init( mesh, meshV1D, &
    out_filebase_tb, dtype, out_tintrv,   &
    myrank, is_master_, NLocalMeshPerPrc_ )

    use scale_atm_phy_tb_dgm_smg, only: &
      atm_phy_tb_dgm_smg_Init
    implicit none    
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

    call atm_phy_tb_dgm_smg_Init( mesh )

    call S11%Init("S11", "s-1", mesh)
    call S12%Init("S12", "s-1", mesh)
    call S22%Init("S22", "s-1", mesh)
    call S23%Init("S23", "s-1", mesh)
    call S31%Init("S31", "s-1", mesh)
    call S33%Init("S33", "s-1", mesh)
    call TKE%Init("TKE", "m2/s2", mesh)

    call dPTdx%Init("dPTdx", "K/m", mesh)
    call dPTdy%Init("dPTdy", "K/m", mesh)
    call dPTdz%Init("dPTdz", "K/m", mesh)

    call Nu%Init( "Nu", "m2/s", mesh )
    call Kh%Init( "Kh", "m2/s", mesh )
    call Ri%Init( "Ri", "1", mesh )

    call Sabs%Init("Sabs", "s-1", mesh)

    call SGS_HEAT_EDDYFLX_V1D%Init("SGS_HEAT_EDDYFLX", "J.m-2.s-1", meshV1D)
    call SGS_MOMZ_EDDYFLX_V1D%Init("SGS_MOMZ_EDDYFLX", "m/s.kg.m-2.s-1", meshV1D)

    !
    call out_file%Init(12, mesh3D=mesh)
    call out_file%Create( out_filebase_tb, "turbulence process", &
     dtype, fileexist, myrank=myrank )
    call out_file%Def_Var( Nu, "Eddy viscosity", DIAGVID_TB_NU, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( Kh, "Eddy diffusivity", DIAGVID_TB_KH, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( Ri, "Richardson number", DIAGVID_TB_Ri, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( Sabs, "Sabs", DIAGVID_TB_Sabs, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( S11, "S11", DIAGVID_TB_S11, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( S12, "S12", DIAGVID_TB_S12, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( S22, "S22", DIAGVID_TB_S22, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( S23, "S23", DIAGVID_TB_S23, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( S31, "S31", DIAGVID_TB_S31, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( S33, "S33", DIAGVID_TB_S33, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( TKE, "TKE", DIAGVID_TB_TKEsgs, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )
    call out_file%Def_Var( dPTdz, "dPTdz", DIAGVID_TB_dPTdz, DIMTYPEID_XYZT, &
      dtype, timeinv=out_tintrv )

    call out_file%End_def()

    ismaster = is_master_
    NLocalMeshPerPrc = NLocalMeshPerPrc_
    if (ismaster) then
      call out_file_V1D%Init(2, mesh1D=meshV1D)
      call out_file_V1D%Create( trim(out_filebase_tb)//"_V1D", 'PBL turbulence analysis', &
        dtype, fileexist, myrank=myrank )
      call out_file_V1D%Def_Var(SGS_HEAT_EDDYFLX_V1D, "horizontal averaged SGS eddy heat flux", &
        DIAGVID_TB_SGS_HEAT_EDDYFLUX_V1D, DIMTYPEID_ZT, dtype, timeinv=out_tintrv )
      call out_file_V1D%Def_Var(SGS_MOMZ_EDDYFLX_V1D, "horizontal averaged SGS eddy momentum flux in z-direction", &
        DIAGVID_TB_SGS_MOMZ_EDDYFLUX_V1D, DIMTYPEID_ZT, dtype, timeinv=out_tintrv )      
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
    mesh, meshV1D         )
  
    use scale_atm_phy_tb_dgm_smg, only: &
      atm_phy_tb_dgm_smg_cal_grad
    use mod_common, only: &
      inquire_bound_flag
    implicit none

    integer, intent(in) :: istep
    real(RP), intent(in) :: start_time, end_time
    class(MeshField3D), intent(in) :: DDENS
    class(MeshField3D), intent(in) :: MOMX
    class(MeshField3D), intent(in) :: MOMY
    class(MeshField3D), intent(in) :: MOMZ
    class(MeshField3D), intent(in) :: DRHOT
    class(MeshField3D), intent(in) :: DENS_hyd
    class(MeshField3D), intent(in) :: PRES_hyd
    class(MeshField3D), intent(in) :: PRES
    class(MeshField3D), intent(in) :: PT
    class(MeshField1D), intent(in) :: DENS_v1D
    class(MeshCubeDom3D), intent(in), target :: mesh
    class(MeshBase1D), intent(in), target :: meshV1D

    integer :: n
    class(LocalMesh3D), pointer :: lmesh
    class(LocalMesh1D), pointer :: lmeshV1D

    logical, allocatable :: is_bound(:,:)
    !--------------------------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(n)

      allocate( is_bound(lmesh%refElem%NfpTot,lmesh%Ne) )
      call inquire_bound_flag( is_bound,              &
        n, lmesh%VMapM, lmesh%VMapP, lmesh%VMapB,     &
        lmesh, lmesh%refElem3D                        )

      call atm_phy_tb_dgm_smg_cal_grad( &
        S11%local(n)%val, S12%local(n)%val, S22%local(n)%val,             &
        S23%local(n)%val, S31%local(n)%val, S33%local(n)%val,             &
        TKE%local(n)%val,                                                 &
        dPTdx%local(n)%val, dPTdy%local(n)%val, dPTdz%local(n)%val,       &
        Nu%local(n)%val, Kh%local(n)%val,                                 &
        DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val,         &
        MOMZ%local(n)%val,  DRHOT%local(n)%val,                           &
        DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                     & 
        PRES%local(n)%val, PT%local(n)%val,                               &
        Dx, Dy, Dz, Sx, Sy, Sz, Lift,                                     &
        lmesh, lmesh%refElem3D, lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D, &
        is_bound )    
      
      deallocate( is_bound )
    end do
    
    do n=1, mesh%LOCAL_MESH_NUM
      lmesh => mesh%lcmesh_list(n)
      lmeshV1D => meshV1D%lcmesh_list(1)
      call cal_auxvars( &
        Ri%local(n)%val, Sabs%local(n)%val,                      &
        SGS_HEAT_EDDYFLX_V1D%local(n)%val,                       &
        SGS_MOMZ_EDDYFLX_V1D%local(n)%val,                       &
        S11%local(n)%val, S12%local(n)%val, S22%local(n)%val,    &
        S23%local(n)%val, S31%local(n)%val, S33%local(n)%val,    &
        TKE%local(n)%val, dPTdz%local(n)%val,                                      &
        Nu%local(n)%val, Kh%local(n)%val,                        &
        DDENS%local(n)%val, PT%local(n)%val,                     &
        DENS_hyd%local(n)%val,                                   &
        DENS_v1D%local(n)%val,                                   &
        n, lmesh, lmesh%refElem3D, lmesh%lcmesh2D, lmesh%lcmesh2D%refElem2D, &
        lmeshV1D, lmeshV1D%refElem1D )
    end do

    call global_horizontal_mean( SGS_HEAT_EDDYFLX_V1D )
    call global_horizontal_mean( SGS_MOMZ_EDDYFLX_V1D )

    call out_file%Write_var3D( DIAGVID_TB_NU, Nu, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_KH, Kh, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_Ri, Ri, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_Sabs, Sabs, start_time, end_time )

    call out_file%Write_var3D( DIAGVID_TB_S11, S11, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_S12, S12, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_S22, S22, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_S23, S23, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_S31, S31, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_S33, S33, start_time, end_time )
    call out_file%Write_var3D( DIAGVID_TB_TKEsgs, TKE, start_time, end_time )

    call out_file%Write_var3D( DIAGVID_TB_dPTdz, dPTdz, start_time, end_time )

    if (ismaster) then
      call out_file_V1D%Write_var1D( DIAGVID_TB_SGS_HEAT_EDDYFLUX_V1D, &
        SGS_HEAT_EDDYFLX_V1D, start_time, end_time )
      call out_file_V1D%Write_var1D( DIAGVID_TB_SGS_MOMZ_EDDYFLUX_V1D, &
        SGS_MOMZ_EDDYFLX_V1D, start_time, end_time )
    end if

    return
  end subroutine diag_tb_process

!OCL SERIAL
  subroutine diag_tb_Final()
    use scale_atm_phy_tb_dgm_smg, only: &
      atm_phy_tb_dgm_smg_Final
    implicit none   
    !--------------------------------------------------------------------

    call atm_phy_tb_dgm_smg_Final()

    call S11%Final(); call S12%Final()
    call S22%Final(); call S23%Final()
    call S31%Final(); call S33%Final()
    call Sabs%Final()
    call TKE%Final()
    
    call dPTdx%Final(); call dPTdy%Final(); call dPTdz%Final()

    call Nu%Final()
    call Kh%Final()
    call Ri%Final()

    call SGS_HEAT_EDDYFLX_V1D%Final()
    call SGS_MOMZ_EDDYFLX_V1D%Final()

    call out_file%Close()
    call out_file%Final()
    if (ismaster) then
      call out_file_V1D%Close()
      call out_file_V1D%Final()
    end if

    return
  end subroutine diag_tb_Final

!- private --------------------------------------
!OCL SERIAL
  subroutine cal_auxvars( &
    Ri_, Sabs_,                         &
    SGS_HEAT_EDDYFLX_V1D_,              &
    SGS_MOMZ_EDDYFLX_V1D_,              &
    S11_, S12_, S22_, S23_, S31_, S33_, &
    TKE_, dPTdz_, Nu_, Kh_,                            & 
    DDENS_, PT_, DENS_hyd_,                            &
    DENS_V1D_,                                         &
    domID, lmesh, elem3D, lmeshH2D, elemH2D, lmeshV1D, elem1D )

    use scale_const, only: &
      EPS => CONST_EPS,  &
      Grav => CONST_GRAV, &
      CpDry => CONST_CPdry
    
    implicit none   
    integer, intent(in) :: domID
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh2D), intent(in) :: lmeshH2D
    class(ElementBase2D), intent(in) :: elemH2D
    class(LocalMesh1D), intent(in) :: lmeshV1D
    class(ElementBase1D), intent(in) :: elem1D
    real(RP), intent(out) :: Ri_(elem3D%Np,lmesh%NeA)
    real(RP), intent(out) :: Sabs_(elem3D%Np,lmesh%NeA)
    real(RP), intent(out) :: SGS_HEAT_EDDYFLX_V1D_(elem1D%Np,lmeshV1D%NeA)
    real(RP), intent(out) :: SGS_MOMZ_EDDYFLX_V1D_(elem1D%Np,lmeshV1D%NeA)
    real(RP), intent(in) :: S11_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: S12_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: S22_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: S23_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: S31_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: S33_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: TKE_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: dPTdz_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: Nu_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: Kh_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: PT_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem3D%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_v1D_(elem1D%Np,lmeshV1D%NeA)

    integer :: ke, ke2D
    integer :: ke_x, ke_y, ke_z
    integer :: p
    real(RP) :: tmp(elem3D%Np)
    real(RP) :: harea

    real(RP) :: DENS_lc(elem3D%Np)
    integer :: hSliceID(elem3D%Nnode_h1D**2)
    real(RP) :: int_w(elem3D%Nnode_h1D**2)
    !--------------------------------------------------------------------

    if (domID==1) then
      !$omp parallel do
      do ke_z=1, lmesh%NeZ
        SGS_HEAT_EDDYFLX_V1D_(:,ke_z) = 0.0_RP
        SGS_MOMZ_EDDYFLX_V1D_(:,ke_z) = 0.0_RP
      end do
      harea = 0.0_RP
    end if

    !$omp parallel do private(tmp)
    do ke=lmesh%NeS, lmesh%NeE
      tmp(:) = 2.0_RP * ( S11_(:,ke)**2 + S22_(:,ke)**2 + S33_(:,ke)**2 ) &
             + 4.0_RP * ( S31_(:,ke)**2 + S12_(:,ke)**2 + S23_(:,ke)**2 )
      Sabs_(:,ke) = sqrt(tmp(:))

      Ri_(:,ke) = Grav / PT_(:,ke) * dPTdz_(:,ke) / max( Sabs_(:,ke)**2, EPS )
    end do
    
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
                            + sum( int_w(:) * DENS_lc(hSliceID(:))    &
                              * Kh_(hSliceID(:),ke) * dPTdz_(hSliceID(:),ke)    )
          
          SGS_MOMZ_EDDYFLX_V1D_(p,ke_z) = SGS_MOMZ_EDDYFLX_V1D_(p,ke_z) &
                            + sum( int_w(:) * DENS_lc(hSliceID(:))      &
                              * ( 2.0_RP * Nu_(hSliceID(:),ke) &
                                 * ( 2.0_RP * S33_(hSliceID(:),ke) - S11_(hSliceID(:),ke) - S22_(hSliceID(:),ke) ) / 3.0_RP &
                                 - TKE_(hSliceID(:),ke) * 2.0_RP/3.0_RP ) )

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
      end do
    end if

    return
  end subroutine cal_auxvars
end module mod_diag_tb
