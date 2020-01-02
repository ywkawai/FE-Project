!-------------------------------------------------------------------------------
!> module ATMOSPHERE dynamical process
!!
!! @par Description
!!          Module for atmosphere dynamical process
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_dyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof
  use scale_const, only: &
    UNDEF => CONST_UNDEF8

  use scale_timeint_rk, only: TimeInt_RK
  
  use scale_mesh_base, only: MeshBase
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_localmesh_base, only: LocalMeshBase
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: &
    ElementBase, ElementBase2D, ElementBase3D
    
  use scale_meshfield_base, only: MeshFieldBase
  use scale_localmeshfield_base, only: LocalMeshFieldBase

  use scale_model_mesh_manager, only: ModelMeshBase
  use scale_model_var_manager, only: ModelVarManager
  use scale_model_component_proc, only:  ModelComponentProc

  use mod_atmos_dyn_bnd, only: AtmosDynBnd
  use mod_atmos_dyn_vars, only: &
    AtmosDynVars, AtmosDynVars_GetLocalMeshFields

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, extends(ModelComponentProc), public :: AtmosDyn
    type(TimeInt_RK), allocatable :: tint(:)
    type(AtmosDynBnd) :: boundary_cond
    type(AtmosDynVars) :: dyn_vars

    logical :: CALC_DIFFVARS_FLAG
    real(RP) :: VISCCOEF_H, VISCCOEF_V
    real(RP) :: DIFFCOEF_H, DIFFCOEF_V

    logical :: EXPFILTER_FLAG
    real(RP) :: EXPFILTER_ETAC
    real(RP) :: EXPFILTER_ALPHA
    integer :: EXPFILTER_ORDER
  contains
    procedure, public :: setup => AtmosDyn_setup 
    procedure, public :: calc_tendency => AtmosDyn_calc_tendency
    procedure, public :: update => AtmosDyn_update
    procedure, public :: finalize => AtmosDyn_finalize
  end type AtmosDyn

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
  !-----------------------------------------------------------------------------

contains

  subroutine AtmosDyn_setup( this, model_mesh )
    use mod_atmos_mesh, only: AtmosMesh
    use mod_atmos_vars, only: ATMOS_PROGVARS_NUM
    use scale_atm_dyn_nonhydro3d, only: &
      atm_dyn_nonhydro3d_prepair_expfilter
    implicit none

    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh

    character(len=H_SHORT) :: TINTEG_TYPE = 'RK_TVD_3'
    real(DP) :: DT_SEC     = 1.0_RP
    real(RP) :: VISCCOEF_H = 0.0_RP
    real(RP) :: VISCCOEF_V = 0.0_RP
    real(RP) :: DIFFCOEF_H = 0.0_RP
    real(RP) :: DIFFCOEF_V = 0.0_RP
    
    logical  :: EXPFILTER_FLAG = .false.
    real(RP) :: EXPFILTER_ETAC  = 2.0_RP/3.0_RP
    real(RP) :: EXPFILTER_ALPHA = 36.0_RP
    integer  :: EXPFILTER_ORDER = 16

    character(len=H_SHORT) :: coriolis_type = 'PLANE'   ! type of coriolis force: 'PLANE', 'SPHERE'
    real(RP) :: coriolis_f0         = 0.0_RP
    real(RP) :: coriolis_beta       = 0.0_RP
    real(RP) :: coriolis_y0         = UNDEF                  ! default is domain center    

    namelist / PARAM_ATMOS_DYN /       &
      TINTEG_TYPE,                            &
      DT_SEC,                                 &
      VISCCOEF_H, VISCCOEF_V,                 &
      DIFFCOEF_H, DIFFCOEF_V,                 &
      EXPFILTER_FLAG,                         &
      EXPFILTER_ETAC, EXPFILTER_ALPHA,        &
      EXPFILTER_ORDER,                        &
      CORIOLIS_TYPE,                          &
      CORIOLIS_f0, CORIOLIS_beta, CORIOLIS_y0
    
    integer :: ierr

    class(MeshBase), pointer :: ptr_mesh
    class(LocalMeshBase), pointer :: ptr_lcmesh
    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    !--------------------------------------------------

    if (.not. this%IsActivated()) return
    LOG_INFO('AtmosDyn_setup',*)

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN. Check!'
      call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN)
    
    !---------------------------------------------------

    call model_mesh%GetModelMesh( ptr_mesh )

    this%dtsec = DT_SEC

    if( VISCCOEF_H > 0.0_RP .or. VISCCOEF_V > 0.0_RP .or. &
        DIFFCOEF_H > 0.0_RP .or. DIFFCOEF_V > 0.0_RP      ) then
      this%CALC_DIFFVARS_FLAG = .true.
      this%VISCCOEF_H = VISCCOEF_H; this%VISCCOEF_V = VISCCOEF_V
      this%DIFFCOEF_H = DIFFCOEF_H; this%DIFFCOEF_V = DIFFCOEF_V
    else
      this%CALC_DIFFVARS_FLAG = .false.
    end if    

    this%EXPFILTER_FLAG = EXPFILTER_FLAG
    if ( this%EXPFILTER_FLAG ) then
      this%EXPFILTER_ALPHA = EXPFILTER_ALPHA
      this%EXPFILTER_ETAC = EXPFILTER_ETAC
      this%EXPFILTER_ORDER = EXPFILTER_ORDER

      call ptr_mesh%GetLocalMesh(1, ptr_lcmesh)
      select type( ptr_lcmesh )
      type is (LocalMesh3D)
        lcmesh3D => ptr_lcmesh
      end select
      call atm_dyn_nonhydro3d_prepair_expfilter( lcmesh3D%refElem3D,  &
        EXPFILTER_ETAC, EXPFILTER_ALPHA, EXPFILTER_ORDER )
    end if

    !-
    allocate( this%tint(ptr_mesh%LOCAL_MESH_NUM) )
    do n = 1, ptr_mesh%LOCAL_MESH_NUM
      call ptr_mesh%GetLocalMesh( n, ptr_lcmesh )
      call this%tint(n)%Init( TINTEG_TYPE, DT_SEC, ATMOS_PROGVARS_NUM, 2, &
        (/ ptr_mesh%refElem%Np, ptr_lcmesh%NeA /) )
    end do

    !- initialize an object to manage boundary conditions and variables for dynamical process
    call this%boundary_cond%Init()
    call this%boundary_cond%SetBCInfo( ptr_mesh )

    call this%dyn_vars%Init( model_mesh )
    call set_coriolis_parameter( this%dyn_vars, model_mesh, CORIOLIS_type, CORIOLIS_f0, CORIOLIS_beta, CORIOLIS_y0 )

    return  
  end subroutine AtmosDyn_setup

  subroutine AtmosDyn_calc_tendency( this, model_mesh, prgvars_list, auxvars_list )
    implicit none
    
    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list

    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    !LOG_INFO('AtmosDyn_tendency',*)

    return  
  end subroutine AtmosDyn_calc_tendency

  subroutine AtmosDyn_update( this, model_mesh, prgvars_list, auxvars_list )
    use scale_atm_dyn_nonhydro3d, only: &
      atm_dyn_nonhydro3d_cal_tend,          &
      atm_dyn_nonhydro3d_cal_grad_diffVars, &
      atm_dyn_nonhydro3d_filter_prgvar
    
    use mod_atmos_vars, only: &
      AtmosVars_GetLocalMeshFields, &
      ATMOS_PROGVARS_DDENS_ID, ATMOS_PROGVARS_DRHOT_ID,                      &
      ATMOS_PROGVARS_MOMX_ID, ATMOS_PROGVARS_MOMY_ID, ATMOS_PROGVARS_MOMZ_ID
     
    implicit none

    class(AtmosDyn), intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    class(ModelVarManager), intent(inout) :: prgvars_list
    class(ModelVarManager), intent(inout) :: auxvars_list

    integer :: rkstage
    integer :: tintbuf_ind

    class(MeshBase), pointer :: mesh
    class(MeshBase2D), pointer :: mesh2D    
    class(LocalMesh3D), pointer :: lcmesh
    integer :: n

    class(LocalMeshFieldBase), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    class(LocalMeshFieldBase), pointer :: &
      GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT
    class(LocalMeshFieldBase), pointer :: DENS_hyd, PRES_hyd
    class(LocalMeshFieldBase), pointer :: Coriolis


    integer :: v
    !--------------------------------------------------
    
    if (.not. this%IsActivated()) return
    LOG_INFO('AtmosDyn_update',*)
    
    call PROF_rapstart( 'ATM_DYN_update', 1)   

    call model_mesh%GetModelMesh( mesh )

    do rkstage=1, this%tint(1)%nstage

      !* Exchange halo data
      call PROF_rapstart( 'ATM_DYN_exchange_prgv', 2)
      call prgvars_list%MeshFieldComm_Exchange()
      call PROF_rapend( 'ATM_DYN_exchange_prgv', 2)

      do n=1, mesh%LOCAL_MESH_NUM
        call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshFields( n, &
          mesh, prgvars_list, auxvars_list,                               &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
          GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT,  &
          DENS_hyd, PRES_hyd, lcmesh                                      )
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)         
        
        !* Apply boundary conditions
        call PROF_rapstart( 'ATM_DYN_applyBC_prgv', 2)
        call this%boundary_cond%ApplyBC_PROGVARS_lc( n,                              & ! (in)
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                        & ! (inout)
          DENS_hyd%val, PRES_hyd%val,                                                & ! (in)
          lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3), & ! (in)
          lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D         ) ! (in)
        call PROF_rapend( 'ATM_DYN_applyBC_prgv', 2)
        
        if ( this%CALC_DIFFVARS_FLAG ) then
          call PROF_rapstart( 'ATM_DYN_cal_grad_diffv', 1)
          call atm_dyn_nonhydro3d_cal_grad_diffVars( &
            GxU%val, GyU%val, GzU%val, GxV%val, GyV%val, GzV%val, GxW%val, GyW%val, GzW%val, &
            GxPT%val, GyPT%val, GzPT%val,                                                    &
            DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                              &
            DENS_hyd%val, PRES_hyd%val,                                                      &
            model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3),          &
            model_mesh%LiftOptrMat,                                                          &
            lcmesh, lcmesh%refElem3D ) 
          call PROF_rapend( 'ATM_DYN_cal_grad_diffv', 2)
        end if
      end do

      if ( this%CALC_DIFFVARS_FLAG ) then
        !* Exchange halo data
        call PROF_rapstart( 'ATM_DYN_exchange_gradv', 2)
        call auxvars_list%MeshFieldComm_Exchange()
        call PROF_rapend( 'ATM_DYN_exchange_gradv', 2)      
        
        do n = 1, mesh%LOCAL_MESH_NUM
          call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
          call AtmosVars_GetLocalMeshFields( n, &
            mesh, prgvars_list, auxvars_list,                               &
            DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
            GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT,  &
            DENS_hyd, PRES_hyd, lcmesh                                      )
          call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)

          !* Apply boundary conditions
          call PROF_rapstart( 'ATM_DYN_applyBC_auxv', 1)
          call this%boundary_cond%ApplyBC_AUXVARS_lc( n,                                     & ! (in)
            GxU%val, GyU%val, GzU%val, GxV%val, GyV%val, GzV%val, GxW%val, GyW%val, GzW%val, & ! (inout)
            GxPT%val, GyPT%val, GzPT%val,                                                    & ! (inout)
            DENS_hyd%val, PRES_hyd%val,                                                      & ! (in)
            this%VISCCOEF_H, this%VISCCOEF_V, this%DIFFCOEF_H, this%DIFFCOEF_V,              & ! (in)
            lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), lcmesh%normal_fn(:,:,3),       & ! (in)
            lcmesh%vmapM, lcmesh%vmapP, lcmesh%vmapB, lcmesh, lcmesh%refElem3D               ) ! (in)
          call PROF_rapend( 'ATM_DYN_applyBC_auxv', 1)
        end do
      end if

      do n=1, mesh%LOCAL_MESH_NUM
        tintbuf_ind = this%tint(n)%tend_buf_indmap(rkstage)

        call PROF_rapstart( 'ATM_DYN_get_localmesh_ptr', 2)         
        call AtmosVars_GetLocalMeshFields( n, &
          mesh, prgvars_list, auxvars_list,                               &
          DDENS, MOMX, MOMY, MOMZ, DRHOT,                                 &
          GxU, GyU, GzU, GxV, GyV, GzV, GxW, GyW, GzW, GxPT, GyPT, GzPT,  &
          DENS_hyd, PRES_hyd, lcmesh                                      )
        
        call AtmosDynVars_GetLocalMeshFields( n,       &
          mesh, this%dyn_vars%AUXVARS_manager,         &
          Coriolis, lcmesh )
        call PROF_rapend( 'ATM_DYN_get_localmesh_ptr', 2)

        call PROF_rapstart( 'ATM_DYN_update_caltend', 2)
        call atm_dyn_nonhydro3d_cal_tend( &
          this%tint(n)%tend_buf2D(:,:,ATMOS_PROGVARS_DDENS_ID,tintbuf_ind),       &
          this%tint(n)%tend_buf2D(:,:,ATMOS_PROGVARS_MOMX_ID ,tintbuf_ind),       &
          this%tint(n)%tend_buf2D(:,:,ATMOS_PROGVARS_MOMY_ID ,tintbuf_ind),       &
          this%tint(n)%tend_buf2D(:,:,ATMOS_PROGVARS_MOMZ_ID ,tintbuf_ind),       &
          this%tint(n)%tend_buf2D(:,:,ATMOS_PROGVARS_DRHOT_ID,tintbuf_ind),       &
          DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val,                     &
          DENS_hyd%val, PRES_hyd%val,                                             &
          Coriolis%val,                                                           &
          GxU%val, GyU%val, GzU%val, GxV%val, GyV%val, GzV%val,                   &
          GxW%val, GyW%val, GzW%val, GxPT%val, GyPT%val, GzPT%val,                &
          this%VISCCOEF_H, this%VISCCOEF_V,                                       &
          this%DIFFCOEF_H, this%DIFFCOEF_V,                                       &          
          model_mesh%DOptrMat(1), model_mesh%DOptrMat(2), model_mesh%DOptrMat(3), &
          model_mesh%SOptrMat(1), model_mesh%SOptrMat(2), model_mesh%SOptrMat(3), &
          model_mesh%LiftOptrMat,                                                 &
          lcmesh, lcmesh%refElem3D, lcmesh%lcmesh2D, lcmesh%lcmesh2D%refElem2D ) 
        call PROF_rapend( 'ATM_DYN_update_caltend', 2)
 
        call PROF_rapstart( 'ATM_DYN_update_advance', 2)      
        call this%tint(n)%Advance( rkstage, DDENS%val, ATMOS_PROGVARS_DDENS_ID, &
                                  1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE  )
        
        call this%tint(n)%Advance( rkstage, MOMX%val, ATMOS_PROGVARS_MOMX_ID,   &
                                  1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE  )
        
        call this%tint(n)%Advance( rkstage, MOMY%val, ATMOS_PROGVARS_MOMY_ID,   &
                                  1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE  )

        call this%tint(n)%Advance( rkstage, MOMZ%val, ATMOS_PROGVARS_MOMZ_ID,   &
                                  1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE  )

        call this%tint(n)%Advance( rkstage, DRHOT%val, ATMOS_PROGVARS_DRHOT_ID, &
                                  1, lcmesh%refElem%Np, lcmesh%NeS, lcmesh%NeE  )
        call PROF_rapend( 'ATM_DYN_update_advance', 2)

        if (this%EXPFILTER_FLAG) then
          call PROF_rapstart( 'ATM_DYN_update_expfilter', 2)
          call atm_dyn_nonhydro3d_filter_prgvar(                & ! (inout)
            DDENS%val, MOMX%val, MOMY%val, MOMZ%val, DRHOT%val, & ! (in)
            lcmesh, lcmesh%refElem3D                            )
          call PROF_rapend( 'ATM_DYN_update_expfilter', 2)
        end if

      end do
    end do

    call PROF_rapend( 'ATM_DYN_update', 1)   


    return  
  end subroutine AtmosDyn_update

  subroutine AtmosDyn_finalize( this )
    implicit none
    class(AtmosDyn), intent(inout) :: this

    integer :: n
    !--------------------------------------------------
    if (.not. this%IsActivated()) return
    LOG_INFO('AtmosDyn_finalize',*)

    do n = 1, size(this%tint)
      call this%tint(n)%Final()
    end do
    deallocate( this%tint )

    call this%boundary_cond%Final()
    call this%dyn_vars%Final()

    return  
  end subroutine AtmosDyn_finalize  

  !--- private ---------------

  subroutine set_coriolis_parameter( this, model_mesh, &
    COLIORIS_type, f0, beta, y0_ )

    use mod_atmos_mesh, only: AtmosMesh
    use scale_const, only: &
      OHM     => CONST_OHM
    implicit none

    class(AtmosDynVars), target, intent(inout) :: this
    class(ModelMeshBase), intent(in) :: model_mesh
    character(*), intent(in) :: COLIORIS_type
    real(RP), intent(in) :: f0, beta, y0_

    class(AtmosMesh), pointer :: atm_mesh
    class(MeshBase2D), pointer :: mesh2D    

    class(LocalMeshFieldBase), pointer :: coriolis
    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: n
    real(RP) :: y0
    !-----------------------------------------------

    nullify( atm_mesh )
    select type(model_mesh)
    type is (AtmosMesh)
      atm_mesh => model_mesh
    end select
    call atm_mesh%mesh%GetMesh2D( mesh2D )

    if (y0_ == UNDEF) then
      y0 = 0.5_RP*(atm_mesh%mesh%ymax_gl -  atm_mesh%mesh%ymin_gl)
    else
      y0 = y0_
    end if

    do n = 1, mesh2D%LOCAL_MESH_NUM
      call AtmosDynVars_GetLocalMeshFields( n, atm_mesh%mesh, this%AUXVARS_manager, &
        coriolis, lcmesh3D )

      if ( trim(COLIORIS_type) == 'PLANE' ) then
        coriolis%val(:,:) = f0 + beta*(lcmesh3D%lcmesh2D%pos_en(:,:,2) - y0)
      else if ( trim(COLIORIS_type) == 'SPHERE' ) then
        coriolis%val(:,:) = 2.0_RP*OHM*sin(lcmesh3D%lat2D(:,:))
      else
        LOG_ERROR('AtmosDyn_set_colioris_parameter',*) 'Unexpected COLIORIS_type is specified. Check! COLIORIS_type=', COLIORIS_type
        call PRC_abort
      end if  
    end do

    return
  end subroutine set_coriolis_parameter

end module mod_atmos_dyn