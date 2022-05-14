!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module for a benchmark test of dynamical core for Earth-like atmosphere proposed by Held and Suarez (1994).
!!          
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user

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
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  &      
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV  
  
  use mod_exp, only: experiment

  use mod_atmos_component, only: &
    AtmosComponent

  use scale_element_base, only: ElementBase3D
  use scale_element_line, only: LineElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: MeshField3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_mkinit
  public :: USER_setup
  public :: USER_calc_tendency
  public :: USER_update

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

  type, private, extends(experiment) :: Exp_Held_Suarez
  contains 
    procedure :: setInitCond_lc => exp_SetInitCond_Held_Suarez
    procedure :: geostrophic_balance_correction_lc => exp_geostrophic_balance_correction
  end type
  type(Exp_Held_Suarez), private :: exp_manager

  logical, private :: USER_do                   = .false. !< do user step?

  real(RP), private, parameter :: kf          = 1.0_RP / ( 86400.0_RP * 1.0_RP  )
  real(RP), private, parameter :: ka          = 1.0_RP / ( 86400.0_RP * 40.0_RP )
  real(RP), private, parameter :: ks          = 1.0_RP / ( 86400.0_RP * 4.0_RP  )
  real(RP), private, parameter :: TEMP_strato = 200.0_RP
  real(RP), private, parameter :: SFCTEMP_eq  = 315.0_RP
  real(RP), private, parameter :: DelT_y      = 60.0_RP
  real(RP), private, parameter :: DelPT_z     = 10.0_RP
  real(RP), private, parameter :: sigb        = 0.7_RP

  !-----------------------------------------------------------------------------

contains

!OCL SERIAL
  subroutine USER_mkinit ( atm )
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    call exp_manager%Init('Held_Suarez')

    call exp_manager%SetInitCond( atm%mesh,                &
      atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager, &
      atm%vars%QTRCVARS_manager                            )
    
    call exp_manager%Final()

    return
  end subroutine USER_mkinit

!OCL SERIAL  
  subroutine USER_setup( atm )
    implicit none
    
    class(AtmosComponent), intent(inout) :: atm

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr    
    !------------------------------------------


    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    return
  end subroutine USER_setup

!OCL SERIAL  
  subroutine USER_calc_tendency( atm )  
    implicit none

    class(AtmosComponent), intent(inout) :: atm
    !------------------------------------------

    return
  end subroutine USER_calc_tendency

!OCL SERIAL
  subroutine USER_update( atm )
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_in
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshPrgVars,    &
      AtmosVars_GetLocalMeshPhyAuxVars, &
      MOMX_p  => ATMOS_PHYTEND_MOMX_ID, &
      MOMY_p  => ATMOS_PHYTEND_MOMY_ID, &
      MOMZ_p  => ATMOS_PHYTEND_MOMZ_ID, &
      RHOH_p  => ATMOS_PHYTEND_RHOH_ID
    use scale_localmeshfield_base, only: LocalMeshFieldBase

    implicit none

    class(AtmosComponent), intent(inout) :: atm

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: PRES, PT
    class(LocalMeshFieldBase), pointer :: Rtot, CVtot, CPtot
    type(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D

    real(RP), allocatable :: DENS(:), T(:), Teq(:), sig(:), PRES_sfc(:)
    real(RP), allocatable :: rtauT(:), rtauV(:)
    real(RP), allocatable :: lat(:)

    real(DP) :: dt
    real(RP) :: Gamm
    !----------------------------------------------------------

    dt = atm%time_manager%dtsec
    gamm = CpDry / CvDry 

    do n=1, atm%mesh%ptr_mesh%LOCAL_MESH_NUM
      call AtmosVars_GetLocalMeshPrgVars( n, atm%mesh%ptr_mesh,  &
        atm%vars%PROGVARS_manager, atm%vars%AUXVARS_manager,     &
        DDENS, MOMX, MOMY, MOMZ, DRHOT,                          &
        DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot, lcmesh           )      

      call AtmosVars_GetLocalMeshPhyAuxVars( n, atm%mesh%ptr_mesh, &
        atm%vars%AUXVARS_manager, PRES, PT                         )
      
      elem3D => lcmesh%refElem3D

      allocate( DENS(elem3D%Np), T(elem3D%Np), Teq(elem3D%Np), sig(elem3D%Np) )
      allocate( PRES_sfc(elem3D%Nnode_h1D**2) )
      allocate( rtauT(elem3D%Np), rtauV(elem3D%Np) )
      allocate( lat(elem3D%Np) )

      !$omp parallel do private(DENS, T, Teq, PRES_sfc, sig, lat, rtauT, rtauV, ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)

        PRES_sfc(:) = PRES%val(elem3D%Hslice(:,1),ke2D)
        sig(:) = PRES%val(:,ke) / PRES_sfc(elem3D%IndexH2Dto3D)
        lat(:) = lcmesh%lat2D(elem3D%IndexH2Dto3D,ke2D)

        rtauT(:) = ka + (ks - ka) * max( 0.0_RP, (sig(:) - sigb)/(1.0_RP - sigb) ) * cos(lat(:))**4
        rtauV(:) = kf * max( 0.0_RP, (sig(:) - sigb)/(1.0_RP - sigb) )

        Teq(:) = max( TEMP_strato, &
          ( SFCTEMP_eq - DelT_y * sin(lat(:))**2 - DelPT_z * log(PRES%val(:,ke)/PRES00) * cos(lat(:))**2 ) &
          * (PRES%val(:,ke)/PRES00)**(Rdry/CPDry)                                                        )
        
        DENS(:) = DENS_hyd%val(:,ke) + DDENS%val(:,ke)
        T(:) =  PRES%val(:,ke) / ( Rdry * DENS(:) )

        MOMX%val(:,ke) = MOMX%val(:,ke) / ( 1.0_RP + dt * rtauV )
        MOMY%val(:,ke) = MOMY%val(:,ke) / ( 1.0_RP + dt * rtauV )

        !-- For the case of d DRHOT / dt = dens * Cv * ( Teq - T ) / tauT <- ( Ullrich and Jablonowski (2012, JCP) )
        ! DRHOT%val(:,ke) = DRHOT%val(:,ke) &
        !                 - dt * rtauT(:) / gamm * ( 1.0_RP - Teq(:) / T(:) ) * DENS(:) * PT%val(:,ke)     &
        !                 / ( 1.0_RP + dt * rtauT(:) / gamm * ( 1.0_RP + (gamm - 1.0_RP) * Teq(:) / T(:) ) )

        !- For the case of d DRHOT /dt = dens * Cp * ( Teq - T ) / tauT     
        !  <- It is based on the forcing form in Held and Surez in which the temperature evolution equation is assumed to be dT/dt = R/Cp * T/p * dp/dt + (Teq - T) / tauT
        DRHOT%val(:,ke) = DRHOT%val(:,ke) &
                        - dt * rtauT(:) * ( 1.0_RP - Teq(:) / T(:) ) * DENS(:) * PT%val(:,ke)     &
                        / ( 1.0_RP + dt * rtauT(:) * ( 1.0_RP + (gamm - 1.0_RP) * Teq(:) / T(:) ) )  
      end do

      deallocate( DENS, T, Teq, sig, PRES_sfc )
      deallocate( rtauT, rtauV )
      deallocate( lat )
    end do
    
    return
  end subroutine USER_update

  !------

!OCL SERIAL
  subroutine exp_SetInitCond_Held_Suarez( this,                            &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, tracer_field_list, &
    x, y, z, dom_xmin, dom_xmax, dom_ymin, dom_ymax, dom_zmin, dom_zmax,   &
    lcmesh, elem )
    
    use mod_exp, only: &
      TracerLocalMeshField_ptr

    use scale_const, only: &
      PRES00 => CONST_PRE00, &
      RPlanet => CONST_Radius, &
      CPdry => CONST_CPdry, &     
      CVdry => CONST_CVdry

    use scale_atm_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constT
  
    implicit none

    class(Exp_Held_Suarez), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
    type(TracerLocalMeshField_ptr), intent(inout) :: tracer_field_list(:)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: y(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_ymin, dom_ymax
    real(RP), intent(in) :: dom_zmin, dom_zmax
    
    real(RP) :: TEMP0           = 280.0_RP
    namelist /PARAM_EXP/ &
      TEMP0

    integer :: ke, ke2D
    integer :: p, p_h, p_v   
    integer :: ierr

    real(RP) :: DENS_ini(elem%Np)
    real(RP) :: PRES_ini(elem%Np)
    !-----------------------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("Held_Suarez_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("Held_Suarez_setup",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    call hydrostatic_calc_basicstate_constT( DENS_hyd, PRES_hyd, & ! (out)
      TEMP0, PRES00,                                             & ! (in)
      x, y, z, lcmesh, elem                                      ) ! (in)

    !$omp parallel do private( p_h, p_v, p, ke2D, &
    !$omp DENS_ini, PRES_ini                      )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)
      do p_v=1, elem%Nnode_v
      do p_h=1, elem%Nnode_h1D**2
        p = p_h + (p_v-1)*elem%Nnode_h1D**2

        call calc_hydrostatic_state_1pt( DENS_ini(p), PRES_ini(p),          & ! (out)
          lcmesh%zlev(p,ke), lcmesh%lon2D(p_h,ke2D), lcmesh%lat2D(p_h,ke2D) ) ! (in)
      end do
      end do

      DDENS(:,ke) = DENS_ini(:) - DENS_hyd(:,ke)
      DRHOT(:,ke) = PRES00 / Rdry * (  ( PRES_ini   (:) / PRES00 )**(CVdry/CPdry) &
                                     - ( PRES_hyd(:,ke) / PRES00 )**(CVdry/CPdry) )
    end do

    return
  end subroutine exp_SetInitCond_Held_Suarez

!OCL SERIAL
  subroutine calc_hydrostatic_state_1pt( dens, pres, z, lon, lat )
    implicit none 
    real(RP), intent(inout) :: dens    
    real(RP), intent(inout) :: pres
    real(RP), intent(in)  :: z
    real(RP), intent(in) :: lon    
    real(RP), intent(in) :: lat

    real(RP) :: geopot
    real(RP) :: temp
    real(RP) :: C1_temp
    real(RP) :: C2_temp
    real(RP) :: eta
    real(RP) :: ln_eta    
    real(RP) :: del_eta
    real(RP) :: exner
    real(RP) :: ln_exner
    real(RP) :: RdryOvCpDry
    real(RP) :: exner_Tstrato
    real(RP) :: ln_exner_Tstrato    
    real(RP) :: geopot_Tstrato

    integer :: itr
    integer,  parameter :: ITRMAX   = 1000
    real(RP), parameter :: CONV_EPS = 1E-15_RP 
    !------------------------------------------------

    C1_temp = SFCTEMP_eq - DelT_y * sin(lat)**2
    C2_temp = - CpDry / Rdry * DelPT_z * cos(lat)**2

    call calc_exner_at_strato( exner_Tstrato, & ! (out)
      lon, lat, C1_temp, C2_temp              ) ! (in)
    
    ln_exner_Tstrato = log(exner_Tstrato)

    geopot_Tstrato = - CpDry * (  &
                      C2_temp * exner_Tstrato * ln_exner_Tstrato         &
                    + ( exner_Tstrato - 1.0_RP ) * ( C1_temp - C2_temp ) &
                  )

    RdryOvCpDry = Rdry / CPDry

    !-- The loop for iteration
    itr     = 0
    eta     = 1.0E-8_RP ! Set first guess of eta
    del_eta = 1.0_RP

    do while( abs(del_eta) > CONV_EPS )
      ln_eta = log(eta)
      exner = eta**RdryOvCpDry
      ln_exner = log(exner)

      temp = ( C1_temp + C2_temp * ln_exner ) * exner  

      if ( temp <= TEMP_strato ) then
        geopot = geopot_Tstrato + CpDry * TEMP_strato * ( ln_exner_Tstrato - ln_exner )
        temp   = TEMP_strato
      else
      ! Pi = eta^(R/Cp)
      ! ln Pi = R/Cp ln eta
      ! d Phi / d Pi = - Cp PT
      ! = - Cp (C1_temp + C2_temp * ln_PI )
        geopot = - CpDry * (  &
                      C2_temp * exner * ln_exner                 &
                    + ( exner - 1.0_RP ) * ( C1_temp - C2_temp ) &
                   )
      end if

      del_eta = -  ( - Grav * z + geopot )        & ! <- F
                 * ( - eta / ( Rdry * temp ) )      ! <- (dF/deta)^-1

      eta = eta + del_eta
      itr = itr + 1

      if ( itr > ITRMAX ) then
        LOG_ERROR("Held_Suarez_setup_calc_hydrostatic_state_1pt",*) "Fail the convergence of iteration. Check!"
        LOG_ERROR_CONT(*) "* (lon,lat,z)=", lon, lat, z
        LOG_ERROR_CONT(*) "itr=", itr, "del_eta=", del_eta, "eta=", eta, "temp=", temp
        call PRC_abort
      end if                                   
    end do  !- End of loop for iteration ----------------------------

    pres = eta * PRES00
    dens = pres / ( Rdry * temp )

    return
  end subroutine calc_hydrostatic_state_1pt

!OCL SERIAL
  subroutine calc_exner_at_strato( exner, lon, lat, C1, C2 )
    implicit none 
    real(RP), intent(out) :: exner 
    real(RP), intent(in) :: lon    
    real(RP), intent(in) :: lat
    real(RP), intent(in) :: C1
    real(RP), intent(in) :: C2

    integer :: itr
    integer,  parameter :: ITRMAX   = 1000
    real(RP), parameter :: CONV_EPS = 1E-15_RP

    real(RP) :: del_exner
    real(RP) :: ln_exner
    !-----------------------------------------

    !-- The loop for iteration
    itr       = 0
    exner     = 1.0E-1_RP ! Set first guess of exner
    del_exner = 1.0_RP

    do while( abs(del_exner) > CONV_EPS )
      ln_exner = log(exner)
      del_exner = -  ( exner * ( C1 + C2 * ln_exner ) - TEMP_strato )  & ! <- F
                   / ( C2 * ln_exner + C1 + C2 )                         ! <- (dF/dexner)^-1

      exner = exner + del_exner
      itr = itr + 1

      if ( itr > ITRMAX ) then
        LOG_ERROR("Held_Suarez_setup_calc_eta_at_strato",*) "Fail the convergence of iteration. Check!"
        LOG_ERROR_CONT(*) "* (lon,lat,z)=", lon, lat
        LOG_ERROR_CONT(*) "itr=", itr, "del_exner=", del_exner, "eta=", exner
        call PRC_abort
      end if                                   
    end do  !- End of loop for iteration ----------------------------

    return
  end subroutine calc_exner_at_strato

!OCL SERIAL
  subroutine exp_geostrophic_balance_correction( this,  &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMY, MOMZ, DRHOT, &
    lcmesh, elem )
    
    implicit none

    class(Exp_Held_Suarez), intent(inout) :: this
    type(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: MOMY(elem%Np,lcmesh%NeA)    
    real(RP), intent(inout) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(inout) :: DRHOT(elem%Np,lcmesh%NeA)

    !---------------------------------------------------
    return
  end subroutine exp_geostrophic_balance_correction 

end module mod_user
