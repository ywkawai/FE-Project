!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!      Explicit numerical diffusion for Atmospheric dynamical process. 
!!      For the discretization, the local DGM (e.g., Cockburn and Shu, 1998) is used. 
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_numdiff
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo
  use scale_model_mesh_manager, only: ModelMesh3D
    
  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    PRGVAR_NUM, &
    DENS_VID => PRGVAR_DDENS_ID, THERM_VID => PRGVAR_THERM_ID,&
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID, DENSHYD_VID => AUXVAR_DENSHYDRO_ID
  
  use scale_atm_dyn_dgm_bnd, only: AtmDynBnd

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !

  type, public :: AtmDyn_Nonhydro3D_Numdiff
    integer ::  ND_LAPLACIAN_NUM
    real(RP) :: ND_COEF_h
    real(RP) :: ND_COEF_v
    real(RP) :: dtsec

    type(MeshField3D), allocatable :: NUMDIFF_FLUX_VARS3D(:)
    type(ModelVarManager) :: NUMDIFF_FLUX_manager
    integer :: NUMDIFF_FLUX_commid

    type(MeshField3D), allocatable :: NUMDIFF_TEND_VARS3D(:)
    type(ModelVarManager) :: NUMDIFF_TEND_manager
    integer :: NUMDIFF_TEND_commid

  contains
    procedure :: Init => atm_dyn_dgm_nonhydro3d_numdiff_Init
    procedure :: Final => atm_dyn_dgm_nonhydro3d_numdiff_Final
    procedure :: Apply => atm_dyn_dgm_nonhydro3d_numdiff_Apply
  end type AtmDyn_Nonhydro3D_Numdiff
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: NUMDIFF_FLUX_NUM   = 3
  integer, public, parameter :: NUMDIFFFLX_X_ID    = 1
  integer, public, parameter :: NUMDIFFFLX_Y_ID    = 2
  integer, public, parameter :: NUMDIFFFLX_Z_ID    = 3

  type(VariableInfo), public :: ATMOS_DYN_NUMDIFF_FLUX_VINFO(NUMDIFF_FLUX_NUM)
  DATA ATMOS_DYN_NUMDIFF_FLUX_VINFO / &
    VariableInfo( NUMDIFFFLX_X_ID, 'DIFFFLX_X', 'flux in x-direction',  &
                  '?.m/s',  3, 'XYZ',  ''                                           ),   & 
    VariableInfo( NUMDIFFFLX_Y_ID, 'DIFFFLX_Y', 'flux in y-direction',  &
                  '?.m/s',  3, 'XYZ',  ''                                           ),   &
    VariableInfo( NUMDIFFFLX_Z_ID, 'DIFFFLX_Z', 'flux in z-direction',  &
                  '?.m/s',  3, 'XYZ',  ''                                           )    /
    
  !-
  integer, public, parameter :: NUMDIFF_TEND_NUM   = 2
  integer, public, parameter :: NUMDIFF_LAPLAH_ID  = 1
  integer, public, parameter :: NUMDIFF_LAPLAV_ID  = 2
  
  type(VariableInfo), public :: ATMOS_DYN_NUMDIFF_TEND_VINFO(NUMDIFF_TEND_NUM)
  DATA ATMOS_DYN_NUMDIFF_TEND_VINFO / &
    VariableInfo( NUMDIFF_LAPLAH_ID, 'NUMDIFF_LAPLAH', 'tendency due to nundiff',  &
                  '?/s',  3, 'XYZ',  ''                                                   ), &
    VariableInfo( NUMDIFF_LAPLAV_ID, 'NUMDIFF_LAPLAV', 'tendency due to nundiff',  &
                  '?/s',  3, 'XYZ',  ''                                                   )  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------

  private :: numdiff_tend 
  private :: numdiff_cal_laplacian 
  private :: numdiff_cal_flx

  private :: cal_del_flux_lap
  private :: cal_del_flux_lap_with_coef
  private :: cal_del_gradDiffVar

contains
!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_numdiff_Init( this, model_mesh3D, dtsec )
    use scale_prc, only: PRC_abort
    implicit none
    class(AtmDyn_Nonhydro3D_Numdiff), intent(inout) :: this
    class(ModelMesh3D), intent(inout), target :: model_mesh3D
    real(RP), intent(in) :: dtsec

    class(MeshBase3D), pointer :: mesh3D
    !--------------------------------------------

    integer ::  ND_LAPLACIAN_NUM = 1
    real(RP) :: ND_COEF_h        = 0.0_RP
    real(RP) :: ND_COEF_v        = 0.0_RP

    namelist /PARAM_ATMOS_DYN_NUMDIFF/ &
      ND_LAPLACIAN_NUM,                &
      ND_COEF_h, ND_COEF_v

    integer :: ierr

    integer :: v
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_NUMDIFF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup_numdiff",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup_numdiff",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_NUMDIFF. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN_NUMDIFF)

    this%ND_LAPLACIAN_NUM = ND_LAPLACIAN_NUM
    this%ND_COEF_H = ND_COEF_h
    this%ND_COEF_v = ND_COEF_v

    this%dtsec = dtsec

    !------------
    mesh3D => model_mesh3D%ptr_mesh

    call this%NUMDIFF_FLUX_manager%Init()
    allocate( this%NUMDIFF_FLUX_VARS3D(NUMDIFF_FLUX_NUM) )

    do v = 1, NUMDIFF_FLUX_NUM
      call this%NUMDIFF_FLUX_manager%Regist(  &
        ATMOS_DYN_NUMDIFF_FLUX_VINFO(v), mesh3D, & ! (in) 
        this%NUMDIFF_FLUX_VARS3D(v),             & ! (inout) 
        .false., fill_zero=.true.                ) ! (in)
    end do

    call model_mesh3D%Create_communicator( &
      NUMDIFF_FLUX_NUM, 0, 0,          & ! (in)
      this%NUMDIFF_FLUX_manager,       & ! (inout)
      this%NUMDIFF_FLUX_VARS3D(:),     & ! (in)
      this%NUMDIFF_FLUX_commid         ) ! (out)

    !-
    call this%NUMDIFF_TEND_manager%Init()
    allocate( this%NUMDIFF_TEND_VARS3D(NUMDIFF_TEND_NUM) )

    do v = 1, NUMDIFF_TEND_NUM
      call this%NUMDIFF_TEND_manager%Regist( &
        ATMOS_DYN_NUMDIFF_TEND_VINFO(v), mesh3D, & ! (in) 
        this%NUMDIFF_TEND_VARS3D(v),             & ! (inout)
        .false., fill_zero=.true.                ) ! (in)      
    end do

    call model_mesh3D%Create_communicator( &
      NUMDIFF_TEND_NUM, 0, 0,          & ! (in)
      this%NUMDIFF_TEND_manager,       & ! (inout)
      this%NUMDIFF_TEND_VARS3D(:),     & ! (in)
      this%NUMDIFF_TEND_commid         ) ! (out)

    return
  end subroutine atm_dyn_dgm_nonhydro3d_numdiff_Init

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_numdiff_Final( this )
    implicit none

    class(AtmDyn_Nonhydro3D_Numdiff), intent(inout) :: this
    !--------------------------------------------
    
    call this%NUMDIFF_FLUX_manager%Final()
    deallocate( this%NUMDIFF_FLUX_VARS3D )

    call this%NUMDIFF_TEND_manager%Final()
    deallocate( this%NUMDIFF_TEND_VARS3D )

    return
  end subroutine atm_dyn_dgm_nonhydro3d_numdiff_Final  

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_numdiff_Apply( this, &
    PROG_VARS, AUX_VARS, boundary_cond, &
    Dx, Dy, Dz, Lift, mesh )

    implicit none

    class(AtmDyn_Nonhydro3D_Numdiff), intent(inout) :: this
    class(ModelVarManager), intent(inout) :: PROG_VARS
    class(ModelVarManager), intent(inout) :: AUX_VARS
    class(AtmDynBnd), intent(in) :: boundary_cond
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    class(MeshBase3D), intent(in), target :: mesh

    class(MeshField3D), pointer :: DDENS, DENS_hyd
    !--------------------------------------------

    call PROG_VARS%Get3D(DENS_VID, DDENS)    
    call AUX_VARS%Get3D(DENSHYD_VID, DENS_hyd)

    call PROG_VARS%MeshFieldComm_Exchange()

    call apply_numfilter( this, THERM_VID, PROG_VARS, boundary_cond, DDENS, DENS_hyd, Dx, Dy, Dz, Lift, mesh )
    call apply_numfilter( this, MOMZ_VID, PROG_VARS, boundary_cond, DDENS, DENS_hyd, Dx, Dy, Dz, Lift, mesh )
    call apply_numfilter( this, MOMX_VID, PROG_VARS, boundary_cond, DDENS, DENS_hyd, Dx, Dy, Dz, Lift, mesh )
    call apply_numfilter( this, MOMY_VID, PROG_VARS, boundary_cond, DDENS, DENS_hyd, Dx, Dy, Dz, Lift, mesh )
    call apply_numfilter( this, DENS_VID, PROG_VARS, boundary_cond, DDENS, DENS_hyd, Dx, Dy, Dz, Lift, mesh )

    return
  end subroutine atm_dyn_dgm_nonhydro3d_numdiff_Apply

  !- private ------------------------------

!OCL SERIAL
  subroutine apply_numfilter( this, varid, prgvars_list, boundary_cond, DDENS, DENS_hyd, &
    Dx, Dy, Dz, Lift, mesh )

    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      PRGVAR_DDENS_ID
        
    implicit none
          
    class(AtmDyn_Nonhydro3D_Numdiff), intent(inout) :: this
    integer, intent(in) :: varid
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(AtmDynBnd), intent(in) :: boundary_cond
    class(MeshField3D), intent(in) :: DDENS
    class(MeshField3D), intent(in) :: DENS_hyd
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    class(MeshBase3D), intent(in), target :: mesh

    class(MeshField3D), pointer :: var
    class(MeshField3D), pointer :: ND_flx_x, ND_flx_y, ND_flx_z
    class(MeshField3D), pointer :: ND_lapla_h, ND_lapla_v

    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    integer :: ke

    integer :: nd_itr
    real(RP) :: nd_sign
    logical :: dens_weight_flag
    logical, allocatable :: is_bound(:,:)

    real(RP), allocatable :: tmp_tend(:,:)
    !-----------------------------------------

    nd_sign = (-1)**(mod(this%ND_LAPLACIAN_NUM+1,2))
    dens_weight_flag = (varid /= PRGVAR_DDENS_ID)

    call prgvars_list%Get3D(varid, var)
    call this%NUMDIFF_FLUX_manager%Get3D(NUMDIFFFLX_X_ID, ND_flx_x)
    call this%NUMDIFF_FLUX_manager%Get3D(NUMDIFFFLX_Y_ID, ND_flx_y)
    call this%NUMDIFF_FLUX_manager%Get3D(NUMDIFFFLX_Z_ID, ND_flx_z)

    call this%NUMDIFF_TEND_manager%Get3D(NUMDIFF_LAPLAH_ID, ND_lapla_h)
    call this%NUMDIFF_TEND_manager%Get3D(NUMDIFF_LAPLAV_ID, ND_lapla_v)

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)

      allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
      call boundary_cond%ApplyBC_numdiff_even_lc( var%local(n)%val, is_bound, varid, n, &
        lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),      &
        lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )
      
      call numdiff_cal_flx( &
        ND_flx_x%local(n)%val, ND_flx_y%local(n)%val, ND_flx_z%local(n)%val,             &
        var%local(n)%val, var%local(n)%val, DDENS%local(n)%val, DENS_hyd%local(n)%val,   &
        Dx, Dy, Dz, Lift, lcmesh, lcmesh%refElem3D, is_bound, dens_weight_flag           ) 

      deallocate( is_bound )
    end do

    !* Exchange halo data
    call this%NUMDIFF_FLUX_manager%MeshFieldComm_Exchange()

    do nd_itr=1, this%ND_LAPLACIAN_NUM-1
      do n = 1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
          
        allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
        call boundary_cond%ApplyBC_numdiff_odd_lc( &
          ND_flx_x%local(n)%val, ND_flx_y%local(n)%val, ND_flx_z%local(n)%val, is_bound, varid, n, &
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),               &
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )

        call numdiff_cal_laplacian( &
          ND_lapla_h%local(n)%val, ND_lapla_v%local(n)%val,                    &
          ND_flx_x%local(n)%val, ND_flx_y%local(n)%val, ND_flx_z%local(n)%val, &
          Dx, Dy, Dz, Lift, lcmesh, lcmesh%refElem3D, is_bound                 )
        
        deallocate( is_bound )
      end do
      !* Exchange halo data
      call this%NUMDIFF_TEND_manager%MeshFieldComm_Exchange()

      do n = 1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
          
        allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne) )
        call boundary_cond%ApplyBC_numdiff_even_lc( &
          ND_lapla_h%local(n)%val, is_bound, varid, n,                               &
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D         )
          
        call numdiff_cal_flx( &
          ND_flx_x%local(n)%val, ND_flx_y%local(n)%val, ND_flx_z%local(n)%val, &
          ND_lapla_h%local(n)%val, ND_lapla_v%local(n)%val,                    &
          DDENS%local(n)%val, DENS_hyd%local(n)%val,                           &
          Dx, Dy, Dz, Lift, lcmesh, lcmesh%refElem3D, is_bound, .false. ) 

        deallocate( is_bound )
      end do
      !* Exchange halo data
      call this%NUMDIFF_FLUX_manager%MeshFieldComm_Exchange()
    end do
    
    do n = 1, mesh%LOCAL_MESH_NUM

      allocate( is_bound(lcmesh%refElem%NfpTot,lcmesh%Ne), tmp_tend(lcmesh%refElem3D%Np,lcmesh%NeA) )

      call boundary_cond%ApplyBC_numdiff_odd_lc( &
        ND_flx_x%local(n)%val, ND_flx_y%local(n)%val, ND_flx_z%local(n)%val,       &
        is_bound, varid, n,                                                        &
        lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), &
        lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D )

      call numdiff_tend( tmp_tend(:,:),  &
        ND_flx_x%local(n)%val, ND_flx_y%local(n)%val, ND_flx_z%local(n)%val,   &
        DDENS%local(n)%val, DENS_hyd%local(n)%val,                             &
        nd_sign * this%ND_COEF_H, nd_sign * this%ND_COEF_V,                    &
        Dx, Dy, Dz, Lift, lcmesh, lcmesh%refElem3D, is_bound, dens_weight_flag ) 

      !$omp parallel do
      do ke=lcmesh%NeS, lcmesh%NeE
        var%local(n)%val(:,ke) = var%local(n)%val(:,ke) + this%dtsec * tmp_tend(:,ke)
      end do

      deallocate( is_bound, tmp_tend )
    end do

    return
  end subroutine apply_numfilter

!OCL SERIAL  
  subroutine numdiff_tend( &
    tend_,                                                     & ! (out)
    GxV_, GyV_, GzV_,                                          & ! (in)
    DDENS_, DENS_hyd, diffcoef_h, diffcoef_v,                  & ! (in)
    Dx, Dy, Dz, Lift, lmesh, elem, is_bound, mul_dens_flag     )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    real(RP), intent(inout)  :: tend_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GyV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: diffcoef_h
    real(RP), intent(in)  :: diffcoef_v
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
    logical, intent(in) :: mul_dens_flag

    real(RP) ::  del_flux(elem%NfpTot,lmesh%Ne)
    real(RP) :: coef_h(elem%Np), coef_v(elem%Np)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    
    integer :: ke

    !--------------------------------------------
    
    call cal_del_flux_lap_with_coef( del_flux,                                & ! (out)
      GxV_, GyV_, GzV_,                                                       & ! (in)
      DDENS_, DENS_hyd, diffcoef_h, diffcoef_v,                               & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound, mul_dens_flag )                                    ! (in)
        
    !$omp parallel do private( &
    !$omp coef_h, coef_v, Fx, Fy, Fz, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE

      if (mul_dens_flag) then
        coef_h(:) = diffcoef_h * (DENS_hyd(:,ke) + DDENS_(:,ke))
        coef_v(:) = diffcoef_v * (DENS_hyd(:,ke) + DDENS_(:,ke))
      else
        coef_h(:) = diffcoef_h
        coef_v(:) = diffcoef_v
      end if
      call sparsemat_matmul(Dx, coef_h(:) * GxV_(:,ke), Fx)
      call sparsemat_matmul(Dy, coef_h(:) * GyV_(:,ke), Fy)
      call sparsemat_matmul(Dz, coef_v(:) * GzV_(:,ke) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke), LiftDelFlx)

      tend_(:,ke) =   ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)      &
          + lmesh%Escale(:,ke,2,2) * Fy(:)      &
          + lmesh%Escale(:,ke,3,3) * Fz(:)      &
          + LiftDelFlx(:)                       )
    end do

    return
  end subroutine numdiff_tend

!OCL SERIAL  
  subroutine cal_del_flux_lap_with_coef( del_flux, & ! (out)
    GxV_, GyV_, GzV_,                              & ! (in)
    DDENS_, DENS_hyd, coef_h, coef_v,              & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem,         & ! (in)
    is_bound, mul_dens_flag )                        ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  GxV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GyV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzV_(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in)  :: coef_h
    real(RP), intent(in)  :: coef_v
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: mul_dens_flag
    
    integer :: i, iP, iM
    real(RP) :: weightP, weightM
    !------------------------------------------------------------------------

    !$omp parallel do                        &
    !$omp private( iM, iP, weightP, weightM  ) 
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      if (mul_dens_flag) then
        weightM = 0.5_RP * (DENS_hyd(iM) + DDENS_(iM))
        weightP = 0.5_RP * (DENS_hyd(iP) + DDENS_(iP))
      else
        weightM = 0.5_RP
        weightP = 0.5_RP
      end if

      if (is_bound(i)) then
        del_flux(i) = &
            coef_h * ( weightP * GxV_(iP) - weightM * GxV_(iM) ) * nx(i) &
          + coef_h * ( weightP * GyV_(iP) - weightM * GyV_(iM) ) * ny(i) &
          + coef_v * ( weightP * GzV_(iP) - weightM * GzV_(iM) ) * nz(i) 
      else
        del_flux(i) = &
            ( 1.0_RP + sign(1.0_RP,nx(i)) ) * coef_h * ( weightP * GxV_(iP) - weightM * GxV_(iM) ) * nx(i) &
          + ( 1.0_RP + sign(1.0_RP,ny(i)) ) * coef_h * ( weightP * GyV_(iP) - weightM * GyV_(iM) ) * ny(i) &
          + ( 1.0_RP + sign(1.0_RP,nz(i)) ) * coef_v * ( weightP * GzV_(iP) - weightM * GzV_(iM) ) * nz(i) 
      end if
    end do

    return
  end subroutine cal_del_flux_lap_with_coef

  !--

!OCL SERIAL  
  subroutine numdiff_cal_laplacian( &
    lapla_h, lapla_v,                                  & ! (out)
    GxV_, GyV_, GzV_,                                  & ! (in)
    Dx, Dy, Dz, Lift, lmesh, elem, is_bound            ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    real(RP), intent(out) :: lapla_h(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: lapla_v(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GxV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GyV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: GzV_(elem%Np,lmesh%NeA)
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)

    real(RP) ::  del_flux_h(elem%NfpTot,lmesh%Ne)
    real(RP) ::  del_flux_v(elem%NfpTot,lmesh%Ne)
    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    
    integer :: ke
    !--------------------------------------------
    
    call cal_del_flux_lap( del_flux_h, del_flux_v,                            & ! (out)
      GxV_, GyV_, GzV_,                                                       & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound )                                                   ! (in)
        
    !$omp parallel do private( &
    !$omp Fx, Fy, Fz, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE

      call sparsemat_matmul(Dx, GxV_(:,ke), Fx)
      call sparsemat_matmul(Dy, GyV_(:,ke), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux_h(:,ke), LiftDelFlx)

      lapla_h(:,ke) = ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)      &
          + lmesh%Escale(:,ke,2,2) * Fy(:)      &
          + LiftDelFlx(:)                       &
          )
      
      call sparsemat_matmul(Dz, GzV_(:,ke) ,Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke)*del_flux_v(:,ke), LiftDelFlx)
      lapla_v(:,ke) = (  &
            lmesh%Escale(:,ke,3,3) * Fz(:)      &
          + LiftDelFlx(:)                       &
          )
    end do

    return
  end subroutine numdiff_cal_laplacian 

!OCL SERIAL
  subroutine cal_del_flux_lap( del_flux_h, del_flux_v, & ! (out)
    GxV_, GyV_, GzV_,                                  & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound    ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux_h(elem%NfpTot*lmesh%Ne)
    real(RP), intent(out) ::  del_flux_v(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  GxV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GyV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  GzV_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    !------------------------------------------------------------------------

    !$omp parallel do         &
    !$omp private( iM, iP ) 
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
          
      if (is_bound(i)) then
        del_flux_h(i) = 0.5_RP * ( ( GxV_(iP) - GxV_(iM) ) * nx(i) &
                                 + ( GyV_(iP) - GyV_(iM) ) * ny(i) )
        del_flux_v(i) = 0.5_RP * ( GzV_(iP) - GzV_(iM) ) * nz(i)
      else
        del_flux_h(i) = 0.5_RP * ( ( 1.0_RP + sign(1.0_RP,nx(i)) ) * ( GxV_(iP) - GxV_(iM) ) * nx(i) &
                                 + ( 1.0_RP + sign(1.0_RP,ny(i)) ) * ( GyV_(iP) - GyV_(iM) ) * ny(i) )
        del_flux_v(i) = 0.5_RP * ( 1.0_RP + sign(1.0_RP,nz(i)) ) * ( GzV_(iP) - GzV_(iM) ) * nz(i)
      end if
    end do

    return
  end subroutine cal_del_flux_lap

  !-------------------------------------------------------

!OCL SERIAL  
  subroutine numdiff_cal_flx( &
    GxV_, GyV_, GzV_,                                         & ! (out)
    Varh_, Varv_, DDENS_, DENS_hyd_,                          & ! (in)
    Dx, Dy, Dz, Lift, lmesh, elem, is_bound, divide_dens_flag ) ! (in)
    
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
    real(RP), intent(out)  :: GxV_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GyV_(elem%Np,lmesh%NeA)
    real(RP), intent(out)  :: GzV_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Varh_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: Varv_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem%Np,lmesh%NeA)
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne)
    logical, intent(in) :: divide_dens_flag

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: vh(elem%Np), vv(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,3)

    integer :: ke

    !------------------------------------------------------------------------------

    call cal_del_gradDiffVar( del_flux,                                       & ! (out)
      Varh_, Varv_, DDENS_, DENS_hyd_,                                        & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound, divide_dens_flag )                                 ! (in)

    !$omp parallel do private(Fx, Fy, Fz, LiftDelFlx, vh, vv)
    do ke=lmesh%NeS, lmesh%NeE

      if (divide_dens_flag) then
        vh(:) = Varh_(:,ke) / (DDENS_(:,ke) + DENS_hyd_(:,ke))
        vv(:) = Varv_(:,ke) / (DDENS_(:,ke) + DENS_hyd_(:,ke))
      else
        vh(:) = Varh_(:,ke)
        vv(:) = Varv_(:,ke)      
      end if

      call sparsemat_matmul(Dx, vh, Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,1), LiftDelFlx)
      GxV_(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dy, vh, Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,2), LiftDelFlx)
      GyV_(:,ke) = lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:)

      call sparsemat_matmul(Dz, vv, Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,3), LiftDelFlx)
      GzV_(:,ke) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)
    end do

    return
  end subroutine numdiff_cal_flx

!OCL SERIAL  
  subroutine cal_del_gradDiffVar( del_flux, & ! (out)
    VARh_, VARv_, DDENS_, DENS_hyd_,        & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem,  & ! (in)
    is_bound, divide_dens_flag              ) ! (in)
    
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(out) :: del_flux(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(in) :: Varh_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: Varv_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: divide_dens_flag
    
    integer :: i, iP, iM
    real(RP) :: delVarh, delVarv
    real(RP) :: weight_P, weight_M
    !------------------------------------------------------------------------

    !$omp parallel do private(iM, iP, delVarh, delVarv,  weight_P, weight_M)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      if (divide_dens_flag) then
        weight_P = 1.0_RP / (DDENS_(iP) + DENS_hyd_(iP))
        weight_M = 1.0_RP / (DDENS_(iM) + DENS_hyd_(iM))
      else
        weight_M = 1.0_RP
        weight_P = 1.0_RP
      end if

      delVarh = 0.5_RP * (Varh_(iP) * weight_P - Varh_(iM) * weight_M)
      delVarv = 0.5_RP * (Varv_(iP) * weight_P - Varv_(iM) * weight_M)
      if (is_bound(i)) then
        del_flux(i,1) = delVarh * nx(i)
        del_flux(i,2) = delVarh * ny(i)
        del_flux(i,3) = delVarv * nz(i)  
      else
        del_flux(i,1) = ( 1.0_RP - sign(1.0_RP,nx(i)) ) * delVarh * nx(i)
        del_flux(i,2) = ( 1.0_RP - sign(1.0_RP,ny(i)) ) * delVarh * ny(i)
        del_flux(i,3) = ( 1.0_RP - sign(1.0_RP,nz(i)) ) * delVarv * nz(i)  
      end if
    end do

    return
  end subroutine cal_del_gradDiffVar

  !------------------------------------------------
end module scale_atm_dyn_dgm_nonhydro3d_numdiff
