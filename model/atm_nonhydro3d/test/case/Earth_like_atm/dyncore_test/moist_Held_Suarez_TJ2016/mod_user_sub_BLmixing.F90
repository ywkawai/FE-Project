#include "scalelib.h"
module mod_user_sub_BLmixing

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort 
  use scale_const, only: &
    Rdry => CONST_Rdry,    &
    Rvap => CONST_Rvap,    &
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  & 
    LHV0 => CONST_LHV0,    &
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV, &
    OHM => CONST_OHM,   &
    RPlanet => CONST_RADIUS, &
    PI => CONST_PI
  use scale_tracer, only: &
    TRACER_inq_id

  use scale_sparsemat, only: SparseMat
  use scale_element_base, only: ElementBase1D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_element_modalfilter, only: ModalFilter

  use mod_atmos_vars, only: &
    AtmosVars, &
    AtmosVarsContainer, &
    AtmosVars_GetLocalMeshPrgVars,    &
    AtmosVars_GetLocalMeshPhyAuxVars, &
    AtmosVars_GetLocalMeshQTRC_Qv

  use mod_user_sub_Filter, only: &
    Filter

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  public :: USER_sub_BLmixing_Init
  public :: USER_sub_BLmixing_calc_tendency
  
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
  real(RP) :: BL_CE_MOM = 0.0_RP
  real(RP) :: BL_CE_PT = 0.0044_RP
  real(RP) :: BL_CE_QV = 0.0044_RP
  real(RP) :: BL_TAPER_Z = 600.0_RP
  real(RP) :: BL_SFC_STAB_VDIFFCoef = 0.0_RP

  logical :: APPLY_NewFilter
  type(Filter) :: newFilter

  type(MeshField3D) :: MOMX_tp_BL
  type(MeshField3D) :: MOMY_tp_BL
  type(MeshField3D) :: RHOT_tp_BL
  type(MeshField3D) :: RHOQ_tp_BL
  type(MeshField3D) :: VViscCoef_BL  
  type(MeshField3D) :: VDiffCoef_BL
  type(MeshField3D) :: HFLX_BL

contains
  subroutine USER_sub_BLmixing_Init( mesh3D )
    implicit none
    class(MeshBase3D), intent(in) :: mesh3D

    character(len=H_SHORT) :: FilterShape
    real(RP) :: FilterWidthFac

    namelist / PARAM_USER_BLmixing / &
      BL_CE_PT, BL_CE_QV, BL_CE_MOM,     &
      BL_SFC_STAB_VDIFFCoef, BL_TAPER_Z, &
      APPLY_NewFilter, FilterShape, FilterWidthFac

    integer :: ierr   
    !------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_BLImixing_setup",*) 'Setup'

    APPLY_NewFilter = .false.
    FilterShape     = "GAUSSIAN"
    FilterWidthFac  = 1.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER_BLmixing,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER_BLmixing. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER_BLmixing)

    !--
    call newFilter%Init( FilterShape, FilterWidthFac, mesh3D )

    !--
    call MOMX_tp_BL%Init( "MOMX_tp_BL", "kg.m-2.s-2", mesh3D)
    call MOMY_tp_BL%Init( "MOMY_tp_BL", "kg.m-2.s-2", mesh3D)
    call RHOT_tp_BL%Init( "RHOT_tp_BL", "kg.m-3.K/s", mesh3D)
    call RHOQ_tp_BL%Init( "RHOQ_tp_BL", "kg.m-3.K/s", mesh3D)
    call VViscCoef_BL%Init( "VViscCoef_BL", "m2/s", mesh3D )
    call VDiffCoef_BL%Init( "VDiffCoef_BL", "m2/s", mesh3D )
    call HFLX_BL%Init( "HFLX_BL", "m2/s", mesh3D )

    return
  end subroutine USER_sub_BLmixing_Init


!OCL SERIAL
  subroutine USER_sub_BLmixing_calc_tendency( MOMX_tp, MOMY_tp, RHOT_tp, RHOQv_tp, vars, vars2, &
    Dz, Lift, mesh3D )
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      RHOT_p  => PHYTEND_RHOT_ID,  & 
      DENS_VID => PRGVAR_DDENS_ID, &
      MOMX_VID => PRGVAR_MOMX_ID,  &
      MOMY_VID => PRGVAR_MOMY_ID,  &
      MOMZ_VID => PRGVAR_MOMZ_ID
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    implicit none
    class(MeshField3D), intent(inout) :: MOMX_tp
    class(MeshField3D), intent(inout) :: MOMY_tp
    class(MeshField3D), intent(inout) :: RHOT_tp
    class(MeshField3D), intent(inout) :: RHOQv_tp
    class(AtmosVarsContainer), intent(inout) :: vars
    class(AtmosVarsContainer), intent(inout) :: vars2
    class(SparseMat), intent(in) :: Dz, Lift    
    class(MeshBase3D), intent(in), target :: mesh3D

    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT, QV
    class(MeshField3D), pointer :: DDENS2, MOMX2, MOMY2, MOMZ2
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot

    integer :: ke
    !----------------------------------------------------------

    call vars2%PROGVARS_manager%Get3D( MOMZ_VID , MOMZ2  )
    call vars2%PROGVARS_manager%Get3D( MOMX_VID , MOMX2  )
    call vars2%PROGVARS_manager%Get3D( MOMY_VID , MOMY2  )
    call vars2%PROGVARS_manager%Get3D( DENS_VID , DDENS2 )

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      call AtmosVars_GetLocalMeshQTRC_Qv( n, mesh3D, &
        vars%QTRCVARS_manager, vars%PHYTENDS_manager, QV )

      call AtmosVars_GetLocalMeshPrgVars( n, mesh3D,  &
        vars%PROGVARS_manager, vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh   )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, mesh3D, &
        vars%AUXVARS_manager, PRES, PT                  )

      call BLmixing( RHOT_tp_BL%local(n)%val, RHOQ_tp_BL%local(n)%val,                    &
        MOMX_tp_BL%local(n)%val, MOMY_tp_BL%local(n)%val,                                 &
        VDiffCoef_BL%local(n)%val, VViscCoef_BL%local(n)%val, HFLX_BL%local(n)%val,       &
        DDENS%val, QV%val,                                                                &
        DDENS2%local(n)%val,  MOMX2%local(n)%val, MOMY2%local(n)%val, MOMZ2%local(n)%val, &
        DENS_hyd%val, PT%val, PRES%val, MOMX%val, MOMY%val,                               &
        Dz, Lift, lcmesh, lcmesh%refElem3D )    
    end do

    if (APPLY_NewFilter) then
      call newFilter%Apply( RHOT_tp_BL, mesh3D )
      call newFilter%Apply( RHOQ_tp_BL, mesh3D )
    else
      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh => mesh3D%lcmesh_list(n)

        !$omp parallel do
        do ke=lcmesh%NeS, lcmesh%NeE
          MOMX_tp%local(n)%val(:,ke) = MOMX_tp%local(n)%val(:,ke) + MOMX_tp_BL%local(n)%val(:,ke)
          MOMY_tp%local(n)%val(:,ke) = MOMY_tp%local(n)%val(:,ke) + MOMY_tp_BL%local(n)%val(:,ke)
          RHOT_tp%local(n)%val(:,ke) = RHOT_tp%local(n)%val(:,ke) + RHOT_tp_BL%local(n)%val(:,ke)
          RHOQv_tp%local(n)%val(:,ke) =  RHOQv_tp%local(n)%val(:,ke) + RHOQ_tp_BL%local(n)%val(:,ke)
        end do
      end do  
    end if

    return
  end subroutine USER_sub_BLmixing_calc_tendency

!--

!OCL SERIAL
  subroutine BLmixing( RHOT_dt, RHOQV_dt, MOMX_dt, MOMY_dt, &
    VDiffCoef_save, VViscCoef_save, HFLX_save,              &
    DDENS, QV,                                              &
    DDENS2, MOMX2, MOMY2, MOMZ2,                            &
    DENS_hyd, PT, PRES, MOMX, MOMY,                         &
    Dz, Lift, lcmesh, elem3D )
    use scale_sparsemat, only: &
      sparsemat, sparsemat_matmul
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: RHOT_dt(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOQV_dt(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX_dt(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY_dt(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: VDiffCoef_save(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: VViscCoef_save(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: HFLX_save(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: QV(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS2(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX2(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY2(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ2(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PT(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    class(SparseMat), intent(in) :: Dz, Lift

    real(RP) :: VDiffCoef(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: VViscCoef(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DiffFlux(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D,4)
    real(RP) :: U_(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: V_(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: PT_(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: QV_(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: za_tmp(elem3D%Nnode_h1D**2)
    real(RP) :: za(elem3D%Np,lcmesh%Ne2D)
    real(RP) :: Vabs_tmp(elem3D%Nnode_h1D**2)
    real(RP) :: Vabs(elem3D%Np,lcmesh%Ne2D)

    integer :: ke, ke_xy, ke_z

    integer :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer :: vmapP(elem3D%NfpTot,lcmesh%NeZ)

    real(RP) :: del_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,4)
    real(RP) :: Fz(elem3D%Np), LiftDelFlx(elem3D%Np)

    real(RP), parameter :: p_pbl = 850E2_RP
    real(RP), parameter :: p_strato = 100E2_RP

    real(RP) :: DENS(elem3D%Np)
    integer :: hslice(elem3D%Nnode_h1D**2)

    real(RP) :: Taper_Z_factor(elem3D%Np)

    ! real(RP) :: mass_check, mass_check2
    !-------------------------------------------------------------

    call lcmesh%GetVmapZ1D( vmapM, vmapP )

    !$omp parallel private(ke, hslice, za_tmp, Vabs_tmp, Taper_Z_factor)

    !$omp do
    do ke_xy=1, lcmesh%Ne2D
      ! hslice(:) = elem3D%Hslice(:,elem3D%Nnode_v)
      ! za_tmp(:) = lcmesh%zlev(hslice(:),ke_xy) / dble(elem3D%Nnode_v)
      hslice(:) = elem3D%Hslice(:,2)
      za_tmp(:) = 60.0_RP !lcmesh%zlev(hslice(:),ke_xy)
      za(:,ke_xy) = za_tmp(elem3D%IndexH2Dto3D(:))

      hslice(:) = elem3D%Hslice(:,1)
      ke = ke_xy
      Vabs_tmp(:) = sqrt( &
          MOMX2(hslice(:),ke) *  (lcmesh%G_ij(:,ke_xy,1,1) * MOMX2(hslice(:),ke) + lcmesh%G_ij(:,ke_xy,2,1) * MOMY2(hslice(:),ke)) &
        + MOMY2(hslice(:),ke) *  (lcmesh%G_ij(:,ke_xy,2,1) * MOMX2(hslice(:),ke) + lcmesh%G_ij(:,ke_xy,2,2) * MOMY2(hslice(:),ke)) &
        + MOMZ2(hslice(:),ke)**2 ) / ( DENS_hyd(hslice(:),ke) + DDENS(hslice(:),ke) )
      Vabs(:,ke_xy) = Vabs_tmp(elem3D%IndexH2Dto3D(:))
    end do

    !$omp do collapse(2)
    do ke_xy=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D

      where ( PRES(:,ke) > p_pbl )
        VDiffCoef(:,ke_z,ke_xy) = Vabs(:,ke_xy) * za(:,ke_xy)
      elsewhere
        VDiffCoef(:,ke_z,ke_xy) = Vabs(:,ke_xy) * za(:,ke_xy) * exp( - ((p_pbl-PRES(:,ke))/p_strato)**2 )
      endwhere
      where ( lcmesh%zlev(:,ke) < BL_TAPER_Z )
        Taper_Z_factor(:) = 0.5_RP * ( 1.0_RP - cos( PI * ( lcmesh%zlev(:,ke) ) / ( BL_TAPER_Z  ) ) )
        VDiffCoef(:,ke_z,ke_xy) = VDiffCoef(:,ke_z,ke_xy) * Taper_Z_factor(:) &
                                + BL_SFC_STAB_VDIFFCoef   * ( 1.0_RP - Taper_Z_factor(:) )
        VViscCoef(:,ke_z,ke_xy) = BL_SFC_STAB_VDIFFCoef   * ( 1.0_RP - Taper_Z_factor(:) )
      elsewhere
        VViscCoef(:,ke_z,ke_xy) = 0.0_RP
      end where

      PT_(:,ke_z,ke_xy) = PT(:,ke)
      QV_(:,ke_z,ke_xy) = QV(:,ke)
      U_(:,ke_z,ke_xy)  = MOMX(:,ke) / ( DENS_hyd(:,ke) + DDENS(:,ke) )
      V_(:,ke_z,ke_xy)  = MOMY(:,ke) / ( DENS_hyd(:,ke) + DDENS(:,ke) )
      nz (:,ke_z,ke_xy) = lcmesh%normal_fn(:,ke,3)

      VDiffCoef_save(:,ke) = VDiffCoef(:,ke_z,ke_xy)
      VViscCoef_save(:,ke) = VViscCoef(:,ke_z,ke_xy)
    end do
    end do

    !$omp end parallel

    call BLmixing_bnd_flux_grad( del_flux, &
      PT_, QV_, U_, V_, nz, vmapM, vmapP, lcmesh, elem3D )

    !$omp parallel do private(ke, Fz, LiftDelFlx, DENS) collapse(2)
    do ke_xy=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)

      call sparsemat_matmul(Dz, PT_(:,ke_z,ke_xy), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,1), LiftDelFlx)
      DiffFlux(:,ke_z,ke_xy,1) = &
        BL_CE_PT * DENS(:) * VDiffCoef(:,ke_z,ke_xy) * ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul(Dz, QV_(:,ke_z,ke_xy), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,2), LiftDelFlx)
      DiffFlux(:,ke_z,ke_xy,2) = &
        BL_CE_QV * DENS(:) * VDiffCoef(:,ke_z,ke_xy) * ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
        
        
      !-
      call sparsemat_matmul(Dz, U_(:,ke_z,ke_xy), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,3), LiftDelFlx)
      DiffFlux(:,ke_z,ke_xy,3) = &
        BL_CE_MOM * DENS(:) * VViscCoef(:,ke_z,ke_xy) * ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )        

      call sparsemat_matmul(Dz, V_(:,ke_z,ke_xy), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,4), LiftDelFlx)
      DiffFlux(:,ke_z,ke_xy,4) = &
        BL_CE_MOM * DENS(:) * VViscCoef(:,ke_z,ke_xy) * ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
    end do
    end do

    call BLmixing_bnd_flux( del_flux, &
      DiffFlux, nz, vmapM, vmapP, lcmesh, elem3D )

    !$omp parallel do private(ke, Fz, LiftDelFlx) collapse(2)
    do ke_xy=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke_xy + (ke_z-1)*lcmesh%Ne2D

      call sparsemat_matmul(Dz, DiffFlux(:,ke_z,ke_xy,1), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,1), LiftDelFlx)
      RHOT_dt(:,ke) = &
        + ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul(Dz, DiffFlux(:,ke_z,ke_xy,2), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,2), LiftDelFlx)
      RHOQV_dt(:,ke) = &
        + ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )

      call sparsemat_matmul(Dz, DiffFlux(:,ke_z,ke_xy,3), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,3), LiftDelFlx)
      MOMX_dt(:,ke) = &
        + ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
      
      call sparsemat_matmul(Dz, DiffFlux(:,ke_z,ke_xy,4), Fz)
      call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,4), LiftDelFlx)
      MOMY_dt(:,ke) = &
        + ( lcmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) )
    end do
    end do

    ! mass_check = 0.0_RP; mass_check2 = 0.0_RP
    ! !$omp parallel do private(ke, Fz, LiftDelFlx) reduction(+: mass_check, mass_check2)
    ! do ke_xy=1, lcmesh%Ne2D
    ! do ke_z=1, lcmesh%NeZ
    !   ke = ke_xy + (ke_z-1)*lcmesh%Ne2D
    !   call sparsemat_matmul(Dz, DiffFlux(:,ke_z,ke_xy,2), Fz)
    !   call sparsemat_matmul(Lift, lcmesh%Fscale(:,ke) * del_flux(:,ke_z,ke_xy,2), LiftDelFlx)
    !   mass_check = mass_check &
    !     + sum( lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) * elem3D%IntWeight_lgl(:) * ( lcmesh%Escale(:,ke,3,3) * ( Fz(:) ) ) )
    !   mass_check2 = mass_check2 &
    !     + sum( lcmesh%J(:,ke) * lcmesh%Gsqrt(:,ke) * elem3D%IntWeight_lgl(:) * ( LiftDelFlx(:) ) )
    ! end do
    ! end do
    ! LOG_INFO("BLM_mass_check: ++",*) mass_check + mass_check2      

    return
  end subroutine BLmixing

!OCL SERIAL
  subroutine BLmixing_bnd_flux_grad( bnd_flux, &
    PT, QV, U, V, nz, vmapM, vmapP, lcmesh, elem3D )
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: bnd_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,4)
    real(RP), intent(in) :: PT(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: QV(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: U(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: V(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem3D%NfpTot,lcmesh%NeZ)

    integer :: ke_z, ke2D, ke
    integer :: iM(elem3D%NfpTot), iP(elem3D%NfpTot)
    real(RP) :: nz_(elem3D%NfpTot)
    !--------------------------------------------------------

    !$omp parallel do private(iM, iP, ke2D, ke_z, ke, nz_) collapse(2)
    do ke2D=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      where ( abs(nz(:,ke_z,ke2D)) > 0.5_RP .and. iM(:) == iP(:) )
        nz_(:) = nz(:,ke_z,ke2D)
      elsewhere
        nz_(:) = ( 1.0_RP + sign(1.0_RP,nz(:,ke_z,ke2D)) ) * nz(:,ke_z,ke2D)
      endwhere
      !nz_(:) = nz(:,ke_z,ke2D)

      bnd_flux(:,ke_z,ke2D,1) = 0.5_RP * ( PT(iP,ke2D) - PT(iM,ke2D) ) * nz_(:)
      bnd_flux(:,ke_z,ke2D,2) = 0.5_RP * ( QV(iP,ke2D) - QV(iM,ke2D) ) * nz_(:)
      bnd_flux(:,ke_z,ke2D,3) = 0.5_RP * ( U(iP,ke2D) - U(iM,ke2D) ) * nz_(:)
      bnd_flux(:,ke_z,ke2D,4) = 0.5_RP * ( V(iP,ke2D) - V(iM,ke2D) ) * nz_(:)
    end do
    end do

    return
  end subroutine BLmixing_bnd_flux_grad

!OCL SERIAL
  subroutine BLmixing_bnd_flux( bnd_flux, &
    DiffFlux, nz, vmapM, vmapP, lcmesh, elem3D )
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: bnd_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,4)
    real(RP), intent(in) :: DiffFlux(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D,4)
    real(RP), intent(in) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem3D%NfpTot,lcmesh%NeZ)

    integer :: ke_z, ke2D, ke
    integer :: iM(elem3D%NfpTot), iP(elem3D%NfpTot)
    real(RP) :: DiffFlux_P(elem3D%NfpTot,4)
    real(RP) :: nz_(elem3D%NfpTot)
    !--------------------------------------------------------

    !$omp parallel do private(iM, iP, ke2D, ke_z, ke, DiffFlux_P, nz_) collapse(2)
    do ke2D=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      where ( abs(nz(:,ke_z,ke2D)) > 0.5_RP .and. iP(:) == iM(:) )
        DiffFlux_P(:,1) = - DiffFlux(iM,ke2D,1)
        DiffFlux_P(:,2) = - DiffFlux(iM,ke2D,2)
        DiffFlux_P(:,3) = - DiffFlux(iM,ke2D,3)
        DiffFlux_P(:,4) = - DiffFlux(iM,ke2D,4)
        nz_(:) = nz(:,ke_z,ke2D)
      elsewhere
        DiffFlux_P(:,1) = DiffFlux(iP,ke2D,1)
        DiffFlux_P(:,2) = DiffFlux(iP,ke2D,2)
        DiffFlux_P(:,3) = DiffFlux(iP,ke2D,3)
        DiffFlux_P(:,4) = DiffFlux(iP,ke2D,4)
        nz_(:) = ( 1.0_RP - sign(1.0_RP,nz(:,ke_z,ke2D)) ) * nz(:,ke_z,ke2D)
      endwhere

      bnd_flux(:,ke_z,ke2D,1) = 0.5_RP * ( DiffFlux_P(:,1) - DiffFlux(iM,ke2D,1) ) * nz_(:)
      bnd_flux(:,ke_z,ke2D,2) = 0.5_RP * ( DiffFlux_P(:,2) - DiffFlux(iM,ke2D,2) ) * nz_(:)
      bnd_flux(:,ke_z,ke2D,3) = 0.5_RP * ( DiffFlux_P(:,3) - DiffFlux(iM,ke2D,3) ) * nz_(:)
      bnd_flux(:,ke_z,ke2D,4) = 0.5_RP * ( DiffFlux_P(:,4) - DiffFlux(iM,ke2D,4) ) * nz_(:)
    end do
    end do

    return
  end subroutine BLmixing_bnd_flux

end module mod_user_sub_BLmixing
