!-------------------------------------------------------------------------------
!> module FElib / Fluid dyn solver / Atmosphere / Nonhydrostatic model / Common
!!
!! @par Description
!!      A coomon model for atmospheric nonhydrostatic dynamical core 
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module scale_atm_dyn_dgm_nonhydro3d_common
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    GRAV => CONST_GRAV,  &
    Rdry => CONST_Rdry,  &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00

  use scale_sparsemat  
  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_2d, only: LocalMesh2D  
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmeshfield_base, only: LocalMeshField3D
  use scale_meshfield_base, only: MeshField3D

  use scale_model_var_manager, only: &
    ModelVarManager, VariableInfo

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_dyn_dgm_nonhydro3d_common_Init
  public :: atm_dyn_dgm_nonhydro3d_common_Final
  public :: atm_dyn_dgm_nonhydro3d_common_setup_variables
  public :: atm_dyn_dgm_nonhydro3d_common_get_varinfo
  public :: atm_dyn_dgm_nonhydro3d_common_calc_pressure
  public :: atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES
  public :: atm_dyn_dgm_nonhydro3d_common_EnTot2PRES
  public :: atm_dyn_dgm_nonhydro3d_common_DRHOT2EnTot

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-
  integer, public, parameter :: PRGVAR_DDENS_ID   = 1
  integer, public, parameter :: PRGVAR_THERM_ID   = 2 ! Variable associated with energy equation (DRHOT or ETOT)
  integer, public, parameter :: PRGVAR_DRHOT_ID   = 2
  integer, public, parameter :: PRGVAR_ETOT_ID    = 2
  integer, public, parameter :: PRGVAR_MOMZ_ID    = 3
  integer, public, parameter :: PRGVAR_MOMX_ID    = 4
  integer, public, parameter :: PRGVAR_MOMY_ID    = 5
  integer, public, parameter :: PRGVAR_SCALAR_NUM = 3
  integer, public, parameter :: PRGVAR_HVEC_NUM   = 1
  integer, public, parameter :: PRGVAR_NUM        = 5

  integer, public, parameter :: PHYTEND_DENS_ID     = 1
  integer, public, parameter :: PHYTEND_MOMX_ID     = 2
  integer, public, parameter :: PHYTEND_MOMY_ID     = 3
  integer, public, parameter :: PHYTEND_MOMZ_ID     = 4
  integer, public, parameter :: PHYTEND_RHOT_ID     = 5
  integer, public, parameter :: PHYTEND_RHOH_ID     = 6
  integer, public, parameter :: PHYTEND_NUM         = 6

  !-
  integer, public, parameter :: AUXVAR_PRESHYDRO_ID     = 1
  integer, public, parameter :: AUXVAR_DENSHYDRO_ID     = 2
  integer, public, parameter :: AUXVAR_PRES_ID          = 3
  integer, public, parameter :: AUXVAR_PT_ID            = 4
  integer, public, parameter :: AUXVAR_Rtot_ID          = 5
  integer, public, parameter :: AUXVAR_CVtot_ID         = 6
  integer, public, parameter :: AUXVAR_CPtot_ID         = 7
  integer, public, parameter :: AUXVAR_Qdry_ID          = 8
  integer, public, parameter :: AUXVAR_PRESHYDRO_REF_ID = 9
  integer, public, parameter :: AUXVAR_NUM              = 9

  
  
  real(RP), public, allocatable :: IntrpMat_VPOrdM1(:,:)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
  !-------------------
  
contains
  subroutine atm_dyn_dgm_nonhydro3d_common_Init( mesh )

    implicit none
    class(MeshBase3D), intent(in) :: mesh

    integer :: p1, p2, p_
    type(ElementBase3D), pointer :: elem
    real(RP) :: invV_VPOrdM1(mesh%refElem3D%Np,mesh%refElem3D%Np)

    !--------------------------------------------

    elem => mesh%refElem3D
    allocate( IntrpMat_VPOrdM1(elem%Np,elem%Np) )
    
    InvV_VPOrdM1(:,:) = elem%invV
    do p2=1, elem%Nnode_h1D
    do p1=1, elem%Nnode_h1D
      p_ = p1 + (p2-1)*elem%Nnode_h1D + (elem%Nnode_v-1)*elem%Nnode_h1D**2
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem%V, invV_VPOrdM1)

    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_Init


  subroutine atm_dyn_dgm_nonhydro3d_common_Final()
    implicit none
    !--------------------------------------------
    
    deallocate( IntrpMat_VPOrdM1 )
    
    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_Final  

  subroutine atm_dyn_dgm_nonhydro3d_common_get_varinfo( &
    prgvar_info, auxvar_info, phytend_info              )

    implicit none
    type(VariableInfo), intent(out) :: prgvar_info(PRGVAR_NUM)
    type(VariableInfo), intent(out) :: auxvar_info(AUXVAR_NUM)
    type(VariableInfo), intent(out), optional :: phytend_info(PHYTEND_NUM)

    type(VariableInfo) :: PRGVAR_VARINFO(PRGVAR_NUM)
    DATA PRGVAR_VARINFO / &
      VariableInfo( PRGVAR_DDENS_ID, 'DDENS', 'deviation of density',         &
                    'kg/m3',  3, 'XYZ',  'air_density'                   ),   &
      VariableInfo( PRGVAR_THERM_ID, 'THERM', 'THERM',                        &
                    '-', 3, 'XYZ',  ''                                     ), &
      VariableInfo( PRGVAR_MOMZ_ID , 'MOMZ', 'momentum z',                    &
                    'kg/m2/s', 3, 'XYZ', 'northward_mass_flux_of_air'      ), &
      VariableInfo( PRGVAR_MOMX_ID , 'MOMX', 'momentum x',                    &
                    'kg/m2/s', 3, 'XYZ', 'upward_mass_flux_of_air'         ), &
      VariableInfo( PRGVAR_MOMY_ID , 'MOMY', 'momentum y',                    &
                    'kg/m2/s', 3, 'XYZ', 'eastward_mass_flux_of_air'       )  /

    type(VariableInfo) :: AUXVAR_VARINFO(AUXVAR_NUM)
    DATA AUXVAR_VARINFO / &
      VariableInfo( AUXVAR_PRESHYDRO_ID, 'PRES_hyd', 'hydrostatic part of pressure',             &
                      'Pa', 3, 'XYZ', ''                                                      ), &
      VariableInfo( AUXVAR_DENSHYDRO_ID, 'DENS_hyd', 'hydrostatic part of density',              &
                    'kg/m3', 3, 'XYZ', ''                                                     ), &
      VariableInfo( AUXVAR_PRES_ID     ,     'PRES', 'pressure',                                 &
                      'Pa', 3, 'XYZ', 'air_pressure'                                          ), &
      VariableInfo( AUXVAR_PT_ID       ,       'PT', 'potential temperature',                    &
                        'K', 3, 'XYZ', 'potential_temperature'                                ), &
      VariableInfo( AUXVAR_Rtot_ID     ,     'RTOT', 'Total gas constant',                       &
                        'J/kg/K', 3, 'XYZ', ''                                                ), &
      VariableInfo( AUXVAR_CVtot_ID    ,    'CVTOT', 'Total heat capacity',                      &
                        'J/kg/K', 3, 'XYZ', ''                                                ), &
      VariableInfo( AUXVAR_CPtot_ID    ,    'CPTOT', 'Total heat capacity',                      &
                        'J/kg/K', 3, 'XYZ', ''                                                ), &
      VariableInfo( AUXVAR_QDRY_ID     ,     'QDRY', 'dry air',                                  &
                        'kg/kg', 3, 'XYZ', ''                                                 ), &
      VariableInfo( AUXVAR_PRESHYDRO_REF_ID, 'PRES_hyd_REF', 'hydrostatic reference pressure',   &
                      'Pa', 3, 'XYZ', ''                                                        )/

    type(VariableInfo) :: PHYTEND_VARINFO(PHYTEND_NUM)
    DATA PHYTEND_VARINFO / &
      VariableInfo( PHYTEND_DENS_ID, 'DENS_tp', 'DENS_tp',                        &
                    'kg/m3/s',  3, 'XYZ',  'tendency of physical process for DENS' ),   &
      VariableInfo( PHYTEND_MOMX_ID, 'MOMX_tp', 'MOMX_tp',                        &
                    'kg/m2/s',  3, 'XYZ',  'tendency of physical process for MOMX' ),   &
      VariableInfo( PHYTEND_MOMY_ID, 'MOMY_tp', 'MOMY_tp',                        &
                    'kg/m2/s',  3, 'XYZ',  'tendency of physical process for MOMY' ),   &
      VariableInfo( PHYTEND_MOMZ_ID, 'MOMZ_tp', 'MOMZ_tp',                        &
                    'kg/m2/s',  3, 'XYZ',  'tendency of physical process for MOMZ' ),   &
      VariableInfo( PHYTEND_RHOT_ID, 'RHOT_tp', 'RHOT_tp',                        &
                    'kg/m3.K/s',  3, 'XYZ',  'tendency of physical process for RHOT' ), &
      VariableInfo( PHYTEND_RHOH_ID,  'RHOH_p',  'RHOH_p',                        &
                    'kg/m3.J/s',  3, 'XYZ',  'heating of physical process for THERM' )   /

    !----------------------------------------------------------

    prgvar_info(:) = PRGVAR_VARINFO
    auxvar_info(:) = AUXVAR_VARINFO
    if ( present(phytend_info) ) phytend_info(:) = PHYTEND_VARINFO

    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_get_varinfo

  subroutine atm_dyn_dgm_nonhydro3d_common_setup_variables( &
    prgvars, qtrcvars, auxvars, phytends,                             & ! (inout)
    prgvar_manager, qtrcvar_manager, auxvar_manager, phytend_manager, & ! (inout)
    PHYTEND_NUM_TOT, mesh3D,                                          & ! (in)
    PRGVAR_VARINFO )                                                    ! (out)

    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_dry
    use scale_tracer, only: &
      QA, TRACER_NAME, TRACER_DESC, TRACER_UNIT    
    use scale_mesh_base3d, only: MeshBase3D
    use scale_meshfield_base, only: MeshField3D
    implicit none
    integer, intent(in) :: PHYTEND_NUM_TOT
    type(MeshField3D), intent(inout) :: prgvars(PRGVAR_NUM)
    type(MeshField3D), intent(inout) :: qtrcvars(0:QA)    
    type(MeshField3D), intent(inout) :: auxvars(AUXVAR_NUM)
    type(MeshField3D), intent(inout) :: phytends(PHYTEND_NUM_TOT)
    type(ModelVarManager), intent(inout) :: prgvar_manager
    type(ModelVarManager), intent(inout) :: qtrcvar_manager
    type(ModelVarManager), intent(inout) :: auxvar_manager
    type(ModelVarManager), intent(inout) :: phytend_manager
    class(MeshBase3D), intent(in) :: mesh3D

    type(VariableInfo), intent(out) :: PRGVAR_VARINFO(PRGVAR_NUM)

    type(VariableInfo) :: AUXVAR_VARINFO(AUXVAR_NUM)
    type(VariableInfo) :: PHYTEND_VARINFO(PHYTEND_NUM)

    integer :: iv
    integer :: iq
    logical :: reg_file_hist

    type(VariableInfo) :: qtrc_dry_vinfo_tmp
    type(VariableInfo) :: qtrc_dry_tp_vinfo_tmp
    type(VariableInfo) :: qtrc_vinfo_tmp
    type(VariableInfo) :: qtrc_tp_vinfo_tmp
    !----------------------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_get_varinfo( PRGVAR_VARINFO, AUXVAR_VARINFO, PHYTEND_VARINFO ) ! (out)

    !- Initialize prognostic variables

    reg_file_hist = .true.    
    do iv = 1, PRGVAR_NUM
      call prgvar_manager%Regist(  &
        PRGVAR_VARINFO(iv), mesh3D,                           & ! (in) 
        prgvars(iv),                                          & ! (inout)
        reg_file_hist,  monitor_flag=.true., fill_zero=.true. ) ! (out)
    end do

    !- Initialize tracer variables

    reg_file_hist = .true.

    if ( ATMOS_HYDROMETEOR_dry ) then
      ! Dummy
      qtrc_dry_vinfo_tmp%ndims    = 3
      qtrc_dry_vinfo_tmp%dim_type = 'XYZ'
      qtrc_dry_vinfo_tmp%STDNAME  = ''

      qtrc_dry_vinfo_tmp%keyID = 0
      qtrc_dry_vinfo_tmp%NAME  = "QV"
      qtrc_dry_vinfo_tmp%DESC  = "Ratio of Water Vapor mass to total mass (Specific humidity)"
      qtrc_dry_vinfo_tmp%UNIT  = "kg/kg"
      call qtrcvar_manager%Regist( &
        qtrc_dry_vinfo_tmp, mesh3D,                      & ! (in) 
        qtrcvars(0),                                     & ! (inout)
        .false., monitor_flag=.false., fill_zero=.true.  ) ! (in)
    else
      qtrc_vinfo_tmp%ndims    = 3
      qtrc_vinfo_tmp%dim_type = 'XYZ'
      qtrc_vinfo_tmp%STDNAME  = ''
  
      do iq = 1, QA
        qtrc_vinfo_tmp%keyID = iq
        qtrc_vinfo_tmp%NAME  = TRACER_NAME(iq)
        qtrc_vinfo_tmp%DESC  = TRACER_DESC(iq)
        qtrc_vinfo_tmp%UNIT  = TRACER_UNIT(iq)
        call qtrcvar_manager%Regist( &
          qtrc_vinfo_tmp, mesh3D,                               & ! (in) 
          qtrcvars(iq),                                         & ! (inout)
          reg_file_hist, monitor_flag=.true., fill_zero=.true.  ) ! (in)
      end do
    end if

    !- Initialize auxiliary variables

    reg_file_hist = .true.
    do iv = 1, AUXVAR_NUM
      call auxvar_manager%Regist( &
        AUXVAR_VARINFO(iv), mesh3D,       & ! (in) 
        auxvars(iv),                      & ! (inout)
        reg_file_hist, fill_zero=.true.   ) ! (in)
    end do

    !- Initialize the tendency of physical processes

    reg_file_hist = .true.
    do iv = 1, PHYTEND_NUM
      call phytend_manager%Regist( &
        PHYTEND_VARINFO(iv), mesh3D,     & ! (in) 
        phytends(iv),                    & ! (inout)
        reg_file_hist, fill_zero=.true.  ) ! (in)
    end do

    if ( ATMOS_HYDROMETEOR_dry ) then
      ! Dummy
      qtrc_dry_tp_vinfo_tmp%ndims    = 3
      qtrc_dry_tp_vinfo_tmp%dim_type = 'XYZ'
      qtrc_dry_tp_vinfo_tmp%STDNAME  = ''

      iv = PHYTEND_NUM + 1
      qtrc_dry_tp_vinfo_tmp%keyID = 0
      qtrc_dry_tp_vinfo_tmp%NAME  = "QV_tp"
      qtrc_dry_tp_vinfo_tmp%DESC  = "tendency of physical process for QV"
      qtrc_dry_tp_vinfo_tmp%UNIT  = "kg/m3/s"
      call phytend_manager%Regist( &
        qtrc_dry_tp_vinfo_tmp, mesh3D,   & ! (in) 
        phytends(iv),                    & ! (inout)
        .false., fill_zero=.true.        ) ! (in)
    else
      qtrc_tp_vinfo_tmp%ndims    = 3
      qtrc_tp_vinfo_tmp%dim_type = 'XYZ'
      qtrc_tp_vinfo_tmp%STDNAME  = ''

      do iq = 1, QA
        iv = PHYTEND_NUM + iq 
        qtrc_tp_vinfo_tmp%keyID = iv
        qtrc_tp_vinfo_tmp%NAME  = trim(TRACER_NAME(iq))//'_tp'
        qtrc_tp_vinfo_tmp%DESC  = 'tendency of physical process for '//trim(TRACER_DESC(iq))
        qtrc_tp_vinfo_tmp%UNIT  = trim(TRACER_UNIT(iq))//'/s'

        call phytend_manager%Regist( &
          qtrc_tp_vinfo_tmp, mesh3D,       & ! (in) 
          phytends(iv),                    & ! (inout)
          reg_file_hist, fill_zero=.true.  ) ! (in)
      end do    
    end if

    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_setup_variables

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_common_calc_pressure( &
    PRES, DPRES,                               & ! (inout)
    DDENS, MOMX, MOMY, MOMZ, THERM,            & ! (in)
    PRES_hyd, DENS_hyd, Rtot, CVtot, CPtot,    & ! (in)
    mesh3D, ENTOT_CONSERVE_SCHEME_FLAG         ) ! (in)

    implicit none
    class(MeshField3D), intent(inout) :: PRES
    class(MeshField3D), intent(inout) :: DPRES
    class(MeshField3D), intent(in) :: DDENS
    class(MeshField3D), intent(in) :: MOMX
    class(MeshField3D), intent(in) :: MOMY
    class(MeshField3D), intent(in) :: MOMZ
    class(MeshField3D), intent(in) :: THERM
    class(MeshField3D), intent(in) :: PRES_hyd
    class(MeshField3D), intent(in) :: DENS_hyd
    class(MeshField3D), intent(in) :: Rtot
    class(MeshField3D), intent(in) :: CVtot
    class(MeshField3D), intent(in) :: CPtot
    class(MeshBase3D), intent(in), target :: mesh3D
    logical, intent(in) :: ENTOT_CONSERVE_SCHEME_FLAG
    
    integer :: n
    class(LocalMesh3D), pointer :: lcmesh3D
    !---------------------------

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D => mesh3D%lcmesh_list(n)

      if ( ENTOT_CONSERVE_SCHEME_FLAG ) then
        call atm_dyn_dgm_nonhydro3d_common_EnTot2PRES( PRES%local(n)%val, DPRES%local(n)%val,                  &
          DDENS%local(n)%val, MOMX%local(n)%val, MOMY%local(n)%val, MOMZ%local(n)%val, THERM%local(n)%val,     &
          PRES_hyd%local(n)%val, DENS_hyd%local(n)%val, Rtot%local(n)%val, CVtot%local(n)%val,                 &
          lcmesh3D, lcmesh3D%refElem3D )
      else
        call atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES( PRES%local(n)%val, DPRES%local(n)%val,                   &
          THERM%local(n)%val, PRES_hyd%local(n)%val, Rtot%local(n)%val, CVtot%local(n)%val, CPtot%local(n)%val, &
          lcmesh3D, lcmesh3D%refElem3D )
      end if
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_calc_pressure

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES( PRES, DPRES, &
    DRHOT, PRES_hyd, Rtot, CVtot, CPtot,                            &
    lcmesh, elem3D                                                  )

    use scale_const, only: &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
        
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: PRES(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(out) :: DPRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA) 

    integer :: ke
    real(RP) :: RHOT(elem3D%Np)
    real(RP) :: rP0
    !---------------------------------------------------------------

    rP0 = 1.0_RP / PRES00
    !$omp parallel do private( RHOT )
    do ke=lcmesh%NeS, lcmesh%NeE
      RHOT(:) = PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**(CvDry/CpDry) + DRHOT(:,ke)

      PRES(:,ke) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT(:) )**( CPtot(:,ke) / CVtot(:,ke) )
      DPRES(:,ke) = PRES(:,ke) - PRES_hyd(:,ke)
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_common_DRHOT2EnTot( EnTot,       &
    DDENS, MOMX, MOMY, MOMZ, DRHOT,                                  &
    DENS_hyd, PRES_hyd, Rtot, CVtot, CPtot,                          &
    lcmesh, elem3D                                                   )

    use scale_const, only: &
      GRAV => CONST_GRAV
        
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: EnTot(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA) 

    integer :: ke, ke2D

    real(RP) :: DENS(elem3D%Np)
    real(RP) :: mom_u1(elem3D%Np), mom_u2(elem3D%Np)

    real(RP) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP) :: DPRES(elem3D%Np,lcmesh%NeA)
    !---------------------------------------------------------------

    call atm_dyn_dgm_nonhydro3d_common_DRHOT2PRES( PRES, DPRES, &
      DRHOT, PRES_hyd, Rtot, CVtot, CPtot,                      &
      lcmesh, elem3D                                            )
    
    !$omp parallel do private( ke2D, DENS, mom_u1, mom_u2 )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      mom_u1(:) = lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,1,1) * MOMX(:,ke) + lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,1) * MOMY(:,ke)
      mom_u2(:) = lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,1) * MOMX(:,ke) + lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,2) * MOMY(:,ke)

      EnTot(:,ke) = PRES(:,ke) * CVtot(:,ke) / Rtot(:,ke)                                                  &
                  + 0.5_RP * ( MOMX(:,ke) * mom_u1(:) + MOMY(:,ke) * mom_u2(:) + MOMZ(:,ke)**2 ) / DENS(:) &
                  + Grav * DENS(:) * lcmesh%zlev(:,ke)
    end do

    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_DRHOT2EnTot

!OCL SERIAL
  subroutine atm_dyn_dgm_nonhydro3d_common_EnTot2PRES( PRES, DPRES, &
    DDENS, MOMX, MOMY, MOMZ, EnTot,         &
    PRES_hyd, DENS_hyd, Rtot, CVtot,        &
    lcmesh, elem3D )

    use scale_const, only: &
      GRAV => CONST_GRAV
    
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(elementbase3D), intent(in) :: elem3D
    real(RP), intent(out) :: PRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(out) :: DPRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: EnTot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)

    integer :: ke, ke2D

    real(RP) :: DENS(elem3D%Np)
    real(RP) :: mom_u1(elem3D%Np), mom_u2(elem3D%Np)
    !---------------------------------------------------------------

    !$omp parallel do private( ke2D, DENS, mom_u1, mom_u2 )
    do ke=lcmesh%NeS, lcmesh%NeE
      ke2D = lcmesh%EMap3Dto2D(ke)

      DENS(:) = DENS_hyd(:,ke) + DDENS(:,ke)
      mom_u1(:) = lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,1,1) * MOMX(:,ke) + lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,1) * MOMY(:,ke)
      mom_u2(:) = lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,1) * MOMX(:,ke) + lcmesh%G_ij(lcmesh%refElem3D%IndexH2Dto3D,ke2D,2,2) * MOMY(:,ke)

      PRES(:,ke) = (  EnTot(:,ke) - Grav * DENS(:) * lcmesh%zlev(:,ke)                                        &
                  - 0.5_RP * ( MOMX(:,ke) * mom_u1(:) + MOMY(:,ke) * mom_u2(:) + MOMZ(:,ke)**2 ) / DENS(:) &
                ) * Rtot(:,ke) / CVtot(:,ke)
      
      DPRES(:,ke) = PRES(:,ke) - PRES_hyd(:,ke)
    end do
    
    return
  end subroutine atm_dyn_dgm_nonhydro3d_common_EnTot2PRES

!-- private

end module scale_atm_dyn_dgm_nonhydro3d_common
