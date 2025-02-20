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
  real(RP) :: BL_CE_PT = 0.0044_RP
  real(RP) :: BL_CE_QV = 0.0044_RP

  logical :: APPLY_NewFilter
  type(Filter) :: newFilter

  type(MeshField3D) :: RHOT_tp_BL
  type(MeshField3D) :: RHOQ_tp_BL

contains
  subroutine USER_sub_BLmixing_Init( mesh3D )
    implicit none
    class(MeshBase3D), intent(in) :: mesh3D

    character(len=H_SHORT) :: FilterShape
    real(RP) :: FilterWidthFac

    namelist / PARAM_USER_BLmixing / &
      BL_CE_PT, BL_CE_QV, &
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
    call RHOT_tp_BL%Init( "RHOT_tp_BL", "kg.m-3.K/s", mesh3D)
    call RHOQ_tp_BL%Init( "RHOQ_tp_BL", "kg.m-3.K/s", mesh3D)

    return
  end subroutine USER_sub_BLmixing_Init


!OCL SERIAL
  subroutine USER_sub_BLmixing_calc_tendency( vars, &
    Dz, Lift, mesh3D )
    use scale_atm_dyn_dgm_nonhydro3d_common, only: &
      RHOT_p  => PHYTEND_RHOT_ID      
    implicit none
    class(AtmosVars), intent(inout) :: vars
    class(SparseMat), intent(in) :: Dz, Lift    
    class(MeshBase3D), intent(in), target :: mesh3D

    class(LocalMesh3D), pointer :: lcmesh
    integer :: n
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT, QV
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    class(LocalMeshFieldBase), pointer :: RHOQv_tp

    integer :: ke
    !----------------------------------------------------------

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      call AtmosVars_GetLocalMeshQTRC_Qv( n, mesh3D, &
        vars%QTRCVARS_manager, vars%PHYTENDS_manager, QV, RHOQv_tp )

      call AtmosVars_GetLocalMeshPrgVars( n, mesh3D,  &
        vars%PROGVARS_manager, vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                  &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh   )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, mesh3D, &
        vars%AUXVARS_manager, PRES, PT                  )

      call BLmixing( RHOT_tp_BL%local(n)%val, RHOQ_tp_BL%local(n)%val, &
        DDENS%val, MOMX%val, MOMY%val, MOMZ%val, QV%val,                   &
        DENS_hyd%val, PT%val, PRES%val, Dz, Lift, lcmesh, lcmesh%refElem3D )
    end do

    if (APPLY_NewFilter) then
      call newFilter%Apply( vars%PHY_TEND(RHOT_p), RHOT_tp_BL, mesh3D )
      call newFilter%Apply( vars%PHY_TEND(RHOT_p+2), RHOQ_tp_BL, mesh3D )
    else
      do n=1, mesh3D%LOCAL_MESH_NUM
        lcmesh => mesh3D%lcmesh_list(n)
        call AtmosVars_GetLocalMeshQTRC_Qv( n, mesh3D, &
          vars%QTRCVARS_manager, vars%PHYTENDS_manager, QV, RHOQv_tp )

        !$omp parallel do
        do ke=lcmesh%NeS, lcmesh%NeE
          vars%PHY_TEND(RHOT_p)%local(n)%val(:,ke) = RHOT_tp_BL%local(n)%val(:,ke)
          RHOQv_tp%val(:,ke) =  RHOQ_tp_BL%local(n)%val(:,ke)
        end do
      end do  
    end if

    return
  end subroutine USER_sub_BLmixing_calc_tendency

!--

!OCL SERIAL
  subroutine BLmixing( RHOT_dt, RHOQV_dt, &
    DDENS, MOMX, MOMY, MOMZ, QV, &
    DENS_hyd, PT, PRES, &
    Dz, Lift, &
    lcmesh, elem3D )
    use scale_sparsemat, only: &
      sparsemat, sparsemat_matmul
    use scale_atm_dyn_dgm_nonhydro3d_rhot_hevi_common, only: &
      vi_gen_vmap => atm_dyn_dgm_nonhydro3d_rhot_hevi_common_gen_vmap      
    implicit none
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: RHOT_dt(elem3D%Np,lcmesh%NeA)
    real(RP), intent(inout) :: RHOQV_dt(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: QV(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PT(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES(elem3D%Np,lcmesh%NeA)
    class(SparseMat), intent(in) :: Dz, Lift

    real(RP) :: VDiffCoef(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D)
    real(RP) :: DiffFlux(elem3D%Np,lcmesh%NeZ,lcmesh%Ne2D,2)
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

    real(RP) :: del_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP) :: Fz(elem3D%Np), LiftDelFlx(elem3D%Np)

    real(RP), parameter :: p_pbl = 850E2_RP
    real(RP), parameter :: p_strato = 100E2_RP

    real(RP) :: G_11(elem3D%Np), G_12(elem3D%Np), G_22(elem3D%Np)
    real(RP) :: DENS(elem3D%Np)
    integer :: hslice(elem3D%Nnode_h1D**2)

    ! real(RP) :: mass_check, mass_check2
    !-------------------------------------------------------------

    call vi_gen_vmap( vmapM, vmapP, lcmesh, elem3D )

    !$omp parallel private(ke, hslice, za_tmp, Vabs_tmp )

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
          MOMX(hslice(:),ke) *  (lcmesh%G_ij(:,ke_xy,1,1) * MOMX(hslice(:),ke) + lcmesh%G_ij(:,ke_xy,2,1) * MOMY(hslice(:),ke)) &
        + MOMY(hslice(:),ke) *  (lcmesh%G_ij(:,ke_xy,2,1) * MOMX(hslice(:),ke) + lcmesh%G_ij(:,ke_xy,2,2) * MOMY(hslice(:),ke)) &
        + MOMZ(hslice(:),ke)**2 ) / ( DENS_hyd(hslice(:),ke) + DDENS(hslice(:),ke) )
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
      
      PT_(:,ke_z,ke_xy) = PT(:,ke)
      QV_(:,ke_z,ke_xy) = QV(:,ke)
      nz (:,ke_z,ke_xy) = lcmesh%normal_fn(:,ke,3)
    end do
    end do

    !$omp end parallel

    call BLmixing_bnd_flux_grad( del_flux, &
      PT_, QV_, nz, vmapM, vmapP, lcmesh, elem3D )

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

  subroutine BLmixing_bnd_flux_grad( bnd_flux, &
    PT, QV, nz, vmapM, vmapP, lcmesh, elem3D )
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: bnd_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP), intent(in) :: PT(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: QV(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D)
    real(RP), intent(in) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem3D%NfpTot,lcmesh%NeZ)

    integer :: ke_z, ke2D, ke
    integer :: iM(elem3D%NfpTot), iP(elem3D%NfpTot)
    !--------------------------------------------------------

    !$omp parallel do private(iM, iP, ke2D, ke_z, ke) collapse(2)
    do ke2D=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)
      bnd_flux(:,ke_z,ke2D,1) = 0.5_RP * ( PT(iP,ke2D) - PT(iM,ke2D) ) * nz(:,ke_z,ke2D)
      bnd_flux(:,ke_z,ke2D,2) = 0.5_RP * ( QV(iP,ke2D) - QV(iM,ke2D) ) * nz(:,ke_z,ke2D)
    end do
    end do

    return
  end subroutine BLmixing_bnd_flux_grad

  subroutine BLmixing_bnd_flux( bnd_flux, &
    DiffFlux, nz, vmapM, vmapP, lcmesh, elem3D )
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: bnd_flux(elem3D%NfpTot,lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP), intent(in) :: DiffFlux(elem3D%Np*lcmesh%NeZ,lcmesh%Ne2D,2)
    real(RP), intent(in) :: nz(elem3D%NfpTot,lcmesh%NeZ,lcmesh%NeX*lcmesh%NeY)
    integer, intent(in) :: vmapM(elem3D%NfpTot,lcmesh%NeZ)
    integer, intent(in) :: vmapP(elem3D%NfpTot,lcmesh%NeZ)

    integer :: ke_z, ke2D, ke
    integer :: iM(elem3D%NfpTot), iP(elem3D%NfpTot)
    real(RP) :: DiffFlux_P(elem3D%NfpTot,2)
    !--------------------------------------------------------

    !$omp parallel do private(iM, iP, ke2D, ke_z, ke, DiffFlux_P) collapse(2)
    do ke2D=1, lcmesh%Ne2D
    do ke_z=1, lcmesh%NeZ
      ke = ke2D + (ke_z-1)*lcmesh%Ne2D
      iM(:) = vmapM(:,ke_z); iP(:) = vmapP(:,ke_z)

      where ( abs(nz(:,ke_z,ke2D)) > 0.5_RP .and. iP(:) == iM(:) )
        DiffFlux_P(:,1) = - DiffFlux(iM,ke2D,1)
        DiffFlux_P(:,2) = - DiffFlux(iM,ke2D,2)
      elsewhere
        DiffFlux_P(:,1) = DiffFlux(iP,ke2D,1)
        DiffFlux_P(:,2) = DiffFlux(iP,ke2D,2)
      endwhere
      bnd_flux(:,ke_z,ke2D,1) = 0.5_RP * ( DiffFlux_P(:,1) - DiffFlux(iM,ke2D,1) ) * nz(:,ke_z,ke2D)
      bnd_flux(:,ke_z,ke2D,2) = 0.5_RP * ( DiffFlux_P(:,2) - DiffFlux(iM,ke2D,2) ) * nz(:,ke_z,ke2D)
    end do
    end do

    return
  end subroutine BLmixing_bnd_flux

end module mod_user_sub_BLmixing
