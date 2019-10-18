!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_dyn_driver
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof

  use scale_sparsemat  
  use scale_element_base
  use scale_element_quadrilateral
  use scale_localmesh_2d
  use scale_mesh_rectdom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D
  use mod_timeint_rk, only: timeint_rk

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: ATMOS_DYN_driver_setup
  public :: ATMOS_DYN_driver_finalize
  public :: ATMOS_DYN_driver

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public  :: ATMOS_DYN_CALC_DIFFVARS_FLAG = .false.
  character(len=H_SHORT), public :: ATMOS_DYN_TINTEG_TYPE = 'RK_TVD_3'

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  type(timeint_rk), private, allocatable :: tinteg_lc(:)
  
  real(RP), private :: ATMOS_DYN_VISCCOEF_H  = 0.0_RP
  real(RP), private :: ATMOS_DYN_VISCCOEF_V  = 0.0_RP
  real(RP), private :: ATMOS_DYN_DIFFCOEF_H  = 0.0_RP
  real(RP), private :: ATMOS_DYN_DIFFCOEF_V  = 0.0_RP

  logical, private :: ATMOS_DYN_EXPFILTER_FLAG   = .false.
  real(RP), private :: ATMOS_DYN_EXPFILTER_ETAC  = 1.0_RP
  real(RP), private :: ATMOS_DYN_EXPFILTER_ALPHA = 36.0_RP
  real(RP), private :: ATMOS_DYN_EXPFILTER_ORDER = 32

contains
  subroutine ATMOS_DYN_driver_setup()
    use scale_atm_dyn_nonhydro2d, only: &
      atm_dyn_nonhydro2d_Init,              &
      atm_dyn_nonhydro2d_prepair_expfilter
    use scale_time_manager, only: TIME_DTSEC

    use mod_atmos_mesh, only: mesh
    use mod_ATMOS_vars, only: PROG_VARS_NUM    
          
    implicit none

    namelist / PARAM_ATMOS_DYN / &
      ATMOS_DYN_TINTEG_TYPE,     &
      ATMOS_DYN_VISCCOEF_H,      &
      ATMOS_DYN_VISCCOEF_V,      &
      ATMOS_DYN_DIFFCOEF_H,      &
      ATMOS_DYN_DIFFCOEF_V,      &
      ATMOS_DYN_EXPFILTER_FLAG,  &
      ATMOS_DYN_EXPFILTER_ETAC,  &
      ATMOS_DYN_EXPFILTER_ALPHA, &
      ATMOS_DYN_EXPFILTER_ORDER
    
    integer :: ierr
    integer :: n

    type(LocalMesh2D), pointer :: lcmesh
    !-------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_DYN_driver_setup",*) 'Setup'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_DYN_driver_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_DYN_driver_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN)

    !
    call atm_dyn_nonhydro2d_Init( mesh )
    
    ! setup for time integrator
    allocate( tinteg_lc(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call tinteg_lc(n)%Init( ATMOS_DYN_TINTEG_TYPE, TIME_DTSEC, PROG_VARS_NUM,  &
                              2, (/ lcmesh%refElem%Np, lcmesh%NeA /)  )
    end do

    !
    if( ATMOS_DYN_VISCCOEF_H > 0.0_RP .or. ATMOS_DYN_VISCCOEF_V > 0.0_RP .or. &
        ATMOS_DYN_DIFFCOEF_H > 0.0_RP .or. ATMOS_DYN_DIFFCOEF_V > 0.0_RP        ) then
      ATMOS_DYN_CALC_DIFFVARS_FLAG = .true.
    else
      ATMOS_DYN_CALC_DIFFVARS_FLAG = .false.
    end if

    if (ATMOS_DYN_EXPFILTER_FLAG) then
      call atm_dyn_nonhydro2d_prepair_expfilter( mesh%lcmesh_list(1)%refElem2D,         &
        ATMOS_DYN_EXPFILTER_ETAC, ATMOS_DYN_EXPFILTER_ALPHA, ATMOS_DYN_EXPFILTER_ORDER )
    end if

    return
  end subroutine ATMOS_DYN_driver_setup

  subroutine ATMOS_DYN_driver_finalize()
    use scale_atm_dyn_nonhydro2d, only: &
      atm_dyn_nonhydro2d_Final
    
    use mod_atmos_mesh, only: mesh
    implicit none
    
    integer :: n
    !-------------------------------------------

    call atm_dyn_nonhydro2d_Final()

    do n=1, mesh%LOCAL_MESH_NUM
      call tinteg_lc(n)%Final()
    end do

    return
  end subroutine ATMOS_DYN_driver_finalize 

  !-------------------------------

  subroutine ATMOS_DYN_driver()   
    use mod_atmos_mesh, only: &
      mesh, Dx, Dz, Sx, Sz, Lift
    use mod_atmos_vars, only: &
      PROG_VARS_NUM, AUX_VARS_NUM, AUX_DIFFVARS_NUM,            &
      DDENS, MOMX, MOMZ, DRHOT,                                 &
      DPRES, DENS_hydro, PRES_hydro,                            &
      GxU, GzU, GxW, GzW, GxPT, GzPT,                           &
      VARS_DDENS_ID, VARS_MOMX_ID, VARS_MOMZ_ID, VARS_DRHOT_ID, &
      VARS_GxU_ID, VARS_GzU_ID, VARS_GxW_ID, VARS_GzW_ID,       &
      VARS_GxPT_ID, VARS_GzPT_ID,                               &
      PROG_VARS_comm, PROG_VARS_list,                           &
      AUX_DIFFVARS_comm, AUX_DIFFVARS_list
    use mod_atmos_bnd, only: &
      ATMOS_bnd_applyBC_prgvars, ATMOS_bnd_applyBC_auxvars
    use scale_atm_dyn_nonhydro2d, only: &
      atm_dyn_nonhydro2d_cal_tend,          &
      atm_dyn_nonhydro2d_cal_grad_diffVars, &
      atm_dyn_nonhydro2d_filter_prgvar
    
    implicit none

    integer :: rkstage
    integer :: tintbuf_ind

    type(LocalMesh2D), pointer :: lcmesh
    integer :: n
    !-------------------------------------------

    do rkstage=1, tinteg_lc(1)%nstage
      
      !* Exchange halo data
      call PROF_rapstart( 'ATM_DYN_exchange_prgv', 1)
      call PROG_VARS_comm%Put(PROG_VARS_list, 1)
      call PROG_VARS_comm%Exchange()
      call PROG_VARS_comm%Get(PROG_VARS_list, 1)
      call PROF_rapend( 'ATM_DYN_exchange_prgv', 1)

      !* Apply boundary conditions
      call PROF_rapstart( 'ATM_DYN_applyBC_prgv', 1)
      call ATMOS_bnd_applyBC_prgvars(    &
        DDENS, MOMX, MOMZ, DRHOT,        &
        DENS_hydro, PRES_hydro,          &
        mesh )
      call PROF_rapend( 'ATM_DYN_applyBC_prgv', 1)
      
      if (ATMOS_DYN_CALC_DIFFVARS_FLAG) then
        call PROF_rapstart( 'ATM_DYN_cal_grad_diffv', 1)
        do n=1, mesh%LOCAL_MESH_NUM
          lcmesh => mesh%lcmesh_list(n)
          call atm_dyn_nonhydro2d_cal_grad_diffVars( &
            GxU%local(n)%val, GzU%local(n)%val, GxW%local(n)%val, GzW%local(n)%val,        &
            GxPT%local(n)%val, GzPT%local(n)%val,                                          &
            DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,  &
            DENS_hydro%local(n)%val, PRES_hydro%local(n)%val,                              &
            Dx, Dz, Lift, lcmesh, lcmesh%refElem2D ) 
        end do
        call PROF_rapend( 'ATM_DYN_cal_grad_diffv', 1)
     

        !* Exchange halo data
        call PROF_rapstart( 'ATM_DYN_exchange_diffv', 1)
        call AUX_DIFFVARS_comm%Put(AUX_DIFFVARS_list, 1)
        call AUX_DIFFVARS_comm%Exchange()
        call AUX_DIFFVARS_comm%Get(AUX_DIFFVARS_list, 1)
        call PROF_rapend( 'ATM_DYN_exchange_diffv', 1)

        !* Apply boundary conditions
        call PROF_rapstart( 'ATM_DYN_applyBC_auxv', 1)
        call ATMOS_bnd_applyBC_auxvars(    &
          GxU, GzU, GxW, GzW, GxPT, GzPT,               &
          DENS_hydro, PRES_hydro,                       &
          ATMOS_DYN_VISCCOEF_H, ATMOS_DYN_VISCCOEF_V,   &
          ATMOS_DYN_DIFFCOEF_H, ATMOS_DYN_DIFFCOEF_V,   &
          mesh )
        call PROF_rapend( 'ATM_DYN_applyBC_auxv', 1)
      end if

      !* Update prognostic variables
      do n=1, mesh%LOCAL_MESH_NUM
        lcmesh => mesh%lcmesh_list(n)
        tintbuf_ind = tinteg_lc(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'ATM_DYN_cal_tend', 1)
        call atm_dyn_nonhydro2d_cal_tend( &
          tinteg_lc(n)%tend_buf2D(:,:,VARS_DDENS_ID,tintbuf_ind),                        &
          tinteg_lc(n)%tend_buf2D(:,:,VARS_MOMX_ID,tintbuf_ind),                         &
          tinteg_lc(n)%tend_buf2D(:,:,VARS_MOMZ_ID,tintbuf_ind),                         &
          tinteg_lc(n)%tend_buf2D(:,:,VARS_DRHOT_ID,tintbuf_ind),                        &
          DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,  &
          DENS_hydro%local(n)%val, PRES_hydro%local(n)%val,                              &
          GxU%local(n)%val, GzU%local(n)%val, GxW%local(n)%val, GzW%local(n)%val,        &
          GxPT%local(n)%val, GzPT%local(n)%val,                                          & 
          ATMOS_DYN_VISCCOEF_H, ATMOS_DYN_VISCCOEF_V,                                    &
          ATMOS_DYN_DIFFCOEF_H, ATMOS_DYN_DIFFCOEF_V,                                    &          
          Dx, Dz, Sx, Sz, Lift, lcmesh, lcmesh%refElem2D ) 
        call PROF_rapend( 'ATM_DYN_cal_tend', 1)
        
        call PROF_rapstart( 'ATM_DYN_tint_advance', 1)
        call tinteg_lc(n)%Advance( rkstage, DDENS%local(n)%val, VARS_DDENS_ID,     &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call tinteg_lc(n)%Advance( rkstage, MOMX%local(n)%val, VARS_MOMX_ID,       &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )  

        call tinteg_lc(n)%Advance( rkstage, MOMZ%local(n)%val, VARS_MOMZ_ID,       &
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )
        
        call tinteg_lc(n)%Advance( rkstage, DRHOT%local(n)%val, VARS_DRHOT_ID,     & 
                                   1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE )                                                                     
        call PROF_rapend('ATM_DYN_tint_advance', 1)
        
        call PROF_rapstart( 'ATM_DYN_tint_expfilter', 1)
        if (ATMOS_DYN_EXPFILTER_FLAG) then
          call atm_dyn_nonhydro2d_filter_prgvar( &
            DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val, &
            lcmesh, lcmesh%refElem2D )
        end if
      end do
    end do

    return
  end subroutine ATMOS_DYN_driver

end module mod_atmos_dyn_driver
