!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Preprocessing
!!
!! @par Description
!!          A module to provide preprocessing operations for variables before evaluating physics tendencies
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_atmos_phy_preproc
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_tracer, only: QA

  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_localmeshfield_base, only: &
    LocalMeshFieldBase, LocalMeshFieldBaseList
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D

  use scale_atm_dyn_dgm_nonhydro3d_common, only: &
    PRGVAR_NUM, AUXVAR_NUM, PHYTEND_NUM1 => PHYTEND_NUM,                                        &
    PRGVAR_DDENS_ID, PRGVAR_THERM_ID, PRGVAR_MOMZ_ID, PRGVAR_MOMX_ID, PRGVAR_MOMY_ID,           &
    AUXVAR_DENSHYDRO_ID, AUXVAR_PRESHYDRO_ID, AUXVAR_THERMHYDRO_ID, AUXVAR_PRESHYDRO_REF_ID,    &
    AUXVAR_Rtot_ID, AUXVAR_CPtot_ID, AUXVAR_CVtot_ID,                                           &
    AUXVAR_PRES_ID, AUXVAR_PT_ID, AUXVAR_Qdry_ID,                                               &
    PHYTEND_DENS_ID, PHYTEND_MOMX_ID, PHYTEND_MOMY_ID, PHYTEND_MOMZ_ID, PHYTEND_RHOT_ID,        &
    PHYTEND_RHOH_ID

  use scale_element_modalfilter, only: ModalFilter
  use scale_meshfield_filter_operation_3d, only: MeshFieldFilterOperation3D  

  use scale_element_operation_tensorprod3D, only: &
    ElementOperationTensorProd3D, &
    ElementOperationTensorprod3D_create

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Base type for preprocessing operations for variables before evaluating physics tendencies
  type, abstract, public :: PhysPreProcBase
    class(MeshBase3D), pointer :: mesh3D
  contains
    procedure(PhysPreProcModalBase_operate), public, deferred :: Operate
    procedure(PhysPreProcModalBase_Final), public, deferred :: Final
  end type PhysPreProcBase
  abstract interface 
    subroutine PhysPreProcModalBase_operate( this, &
      PRG_VARS, QTRC_VARS, AUX_VARS,    &
      PRG_VARS0, QTRC_VARS0, AUX_VARS0  )
      import PhysPreProcBase
      import MeshField3D
      import PRGVAR_NUM
      import AUXVAR_NUM
      import QA
      class(PhysPreProcBase), intent(inout) :: this
      type(MeshField3D), intent(inout) :: PRG_VARS(PRGVAR_NUM)
      type(MeshField3D), intent(inout) :: AUX_VARS(AUXVAR_NUM)
      type(MeshField3D), intent(inout) :: QTRC_VARS(0:QA)
      type(MeshField3D), intent(in) :: PRG_VARS0(PRGVAR_NUM)
      type(MeshField3D), intent(in) :: AUX_VARS0(AUXVAR_NUM)
      type(MeshField3D), intent(in) :: QTRC_VARS0(0:QA)
    end subroutine PhysPreProcModalBase_operate

    subroutine PhysPreProcModalBase_Final( this )
      import PhysPreProcBase
      class(PhysPreProcBase), intent(inout) :: this
    end subroutine PhysPreProcModalBase_Final
  end interface

  !> Derived type to represent no preprocessing operation for variables before evaluating physics tendencies  
  type, extends(PhysPreProcBase), public :: PhysPreProcNone
  contains
    procedure :: Init => PhysPreProcNone_Init
    procedure :: Final => PhysPreProcNone_Final
    procedure :: Operate => PhysPreProcNone_Operate
  end type PhysPreProcNone

  !> Derived type to represent modal filtering operation for variables before evaluating physics tendencies
  type, extends(PhysPreProcBase), public :: PhysPreProcModalFilter
    type(ModalFilter) :: mFilter
    class(ElementOperationTensorProd3D), allocatable :: elem_optr
  contains
    procedure :: Init => PhysPreProcModalFilter_Init
    procedure :: Final => PhysPreProcModalFilter_Final
    procedure :: Operate => PhysPreProcModalFilter_Operate
  end type PhysPreProcModalFilter

  !> Derived type to represent quasi-global filtering operation for variables before evaluating physics tendencies
  type, extends(PhysPreProcBase), public :: PhysPreProcGlobalFilter
    type(MeshFieldFilterOperation3D) :: gFilter_prgvar
    type(MeshFieldFilterOperation3D) :: gFilter_trcvar
    type(MeshField3D), allocatable :: RHOQ(:)
  contains
    procedure :: Init => PhysPreProcGlobalFilter_Init
    procedure :: Final => PhysPreProcGlobalFilter_Final
    procedure :: Operate => PhysPreProcGlobalFilter_Operate
  end type PhysPreProcGlobalFilter

  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures & variables
  !
contains
!--
  !> Initialize an object to represent no preprocessing operation before evaluating physics tendencies
!OCL SERIAL
  subroutine PhysPreProcNone_Init( this, mesh3D )
    implicit none
    class(PhysPreProcNone), intent(inout) :: this
    class(MeshBase3D), intent(in), target :: mesh3D
    !----------------------------------------------------
    this%mesh3D => mesh3D
    return
  end subroutine PhysPreProcNone_Init

  !> Finalize an object to represent no preprocessing operation before evaluating physics tendencies
!OCL SERIAL
  subroutine PhysPreProcNone_Final( this )
    implicit none
    class(PhysPreProcNone), intent(inout) :: this
    !----------------------------------------------------
    return
  end subroutine PhysPreProcNone_Final  

  !> Perform no preprocessing operation before evaluating physics tendencies
  !! This subroutine simply copies the input variables to the output variables without any modification.
!OCL SERIAL
  subroutine PhysPreProcNone_Operate( this, &
    PRG_VARS, QTRC_VARS, AUX_VARS,    &
    PRG_VARS0, QTRC_VARS0, AUX_VARS0  )
    use scale_tracer, only: QA
    implicit none
    class(PhysPreProcNone), intent(inout) :: this
    type(MeshField3D), intent(inout) :: PRG_VARS(PRGVAR_NUM)
    type(MeshField3D), intent(inout) :: AUX_VARS(AUXVAR_NUM)
    type(MeshField3D), intent(inout) :: QTRC_VARS(0:QA)
    type(MeshField3D), intent(in) :: PRG_VARS0(PRGVAR_NUM)
    type(MeshField3D), intent(in) :: AUX_VARS0(AUXVAR_NUM)
    type(MeshField3D), intent(in) :: QTRC_VARS0(0:QA)

    integer :: n
    integer :: iq, iqs
    !----------------------------------------------------

    do n=1, this%mesh3D%LOCAL_MESH_NUM
      call copy_data( &
        PRG_VARS(PRGVAR_DDENS_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMX_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMY_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMZ_ID)%local(n)%val, PRG_VARS(PRGVAR_THERM_ID)%local(n)%val,      & ! (inout)
        AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val, AUX_VARS(AUXVAR_PRESHYDRO_ID)%local(n)%val,                                                                                                                   & ! (inout)
        PRG_VARS0(PRGVAR_DDENS_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMX_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMY_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMZ_ID)%local(n)%val, PRG_VARS0(PRGVAR_THERM_ID)%local(n)%val, & ! (in)
        AUX_VARS0(AUXVAR_DENSHYDRO_ID)%local(n)%val, AUX_VARS0(AUXVAR_PRESHYDRO_ID)%local(n)%val,                                                                                                                 & ! (in)
        this%mesh3D%lcmesh_list(n), this%mesh3D%refElem3D ) ! (in)
    end do

    if ( QA == 0 ) then
      iqs = 0
    else
      iqs = 1
    end if
    do iq=iqs, QA
      do n=1, this%mesh3D%LOCAL_MESH_NUM
        call copy_data_qtrc( QTRC_VARS(iq)%local(n)%val,                                                                   & ! (inout)
          QTRC_VARS0(iq)%local(n)%val, AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val, PRG_VARS(PRGVAR_DDENS_ID)%local(n)%val, & ! (in)
          this%mesh3D%lcmesh_list(n), this%mesh3D%refElem3D ) ! (in)
      end do
    end do
    return
  end subroutine PhysPreProcNone_Operate

!--------------------------------------------
  !> Initialize an object to represent modal filtering operation before evaluating physics tendencies
!OCL SERIAL
  subroutine PhysPreProcModalFilter_Init( this,     &
    MF_ALPHA_h, MF_ORDER_h, MF_ALPHA_v, MF_ORDER_v, &
    mesh3D )
    implicit none
    class(PhysPreProcModalFilter), intent(inout) :: this
    real(RP), intent(in) :: MF_ALPHA_h
    integer, intent(in) :: MF_ORDER_h
    real(RP), intent(in) :: MF_ALPHA_v
    integer, intent(in) :: MF_ORDER_v
    class(MeshBase3D), intent(in), target :: mesh3D
    !----------------------------------------------------

    this%mesh3D => mesh3D

    call ElementOperationTensorprod3D_create( mesh3D%refElem3D, &
        this%elem_optr )
    
    call this%elem_optr%Setup_ModalFilter( 0.0_RP, MF_ALPHA_h, MF_ORDER_h, 0.0_RP, MF_ALPHA_v, MF_ORDER_v )
    call this%elem_optr%Setup_ModalFilter_tracer( 0.0_RP, MF_ALPHA_h, MF_ORDER_h, 0.0_RP, MF_ALPHA_v, MF_ORDER_v )    
    return
  end subroutine PhysPreProcModalFilter_Init

  !> Finalize an object to represent modal filtering operation before evaluating physics tendencies
!OCL SERIAL
  subroutine PhysPreProcModalFilter_Final( this )
    use scale_element_hexahedral, only: HexahedralElement
    implicit none
    class(PhysPreProcModalFilter), intent(inout) :: this
    !----------------------------------------------------

    call this%mfilter%Final()
    
    call this%elem_optr%Final()
    deallocate( this%elem_optr )
    return
  end subroutine PhysPreProcModalFilter_Final  

  !> Perform modal filtering operation before evaluating physics tendencies
!OCL SERIAL
  subroutine PhysPreProcModalFilter_Operate( this, &
    PRG_VARS, QTRC_VARS, AUX_VARS,    &
    PRG_VARS0, QTRC_VARS0, AUX_VARS0  )
    use scale_element_hexahedral, only: HexahedralElement
    use scale_tracer, only: QA
    implicit none
    class(PhysPreProcModalFilter), intent(inout) :: this
    type(MeshField3D), intent(inout) :: PRG_VARS(PRGVAR_NUM)
    type(MeshField3D), intent(inout) :: AUX_VARS(AUXVAR_NUM)
    type(MeshField3D), intent(inout) :: QTRC_VARS(0:QA)
    type(MeshField3D), intent(in) :: PRG_VARS0(PRGVAR_NUM)
    type(MeshField3D), intent(in) :: AUX_VARS0(AUXVAR_NUM)
    type(MeshField3D), intent(in) :: QTRC_VARS0(0:QA)

    integer :: n
    integer :: iq, iqs
    !----------------------------------------------------

    do n=1, this%mesh3D%LOCAL_MESH_NUM
      call apply_modal_filter_core( &
        PRG_VARS(PRGVAR_DDENS_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMX_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMY_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMZ_ID)%local(n)%val, PRG_VARS(PRGVAR_THERM_ID)%local(n)%val,      & ! (inout)
        AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val, AUX_VARS(AUXVAR_PRESHYDRO_ID)%local(n)%val,                                                                                                                   & ! (inout)
        PRG_VARS0(PRGVAR_DDENS_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMX_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMY_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMZ_ID)%local(n)%val, PRG_VARS0(PRGVAR_THERM_ID)%local(n)%val, & ! (in)
        AUX_VARS0(AUXVAR_DENSHYDRO_ID)%local(n)%val, AUX_VARS0(AUXVAR_PRESHYDRO_ID)%local(n)%val,                                                                                                                 & ! (in)
        this%elem_optr, this%mesh3D%lcmesh_list(n), this%mesh3D%refElem3D ) ! (in)
    end do

    if ( QA == 0 ) then
      iqs = 0
    else
      iqs = 1
    end if
    do iq=iqs, QA
      do n=1, this%mesh3D%LOCAL_MESH_NUM
        call apply_modal_filter_qtrc_core( QTRC_VARS(iq)%local(n)%val,                                                     & ! (inout)
          QTRC_VARS0(iq)%local(n)%val, AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val, PRG_VARS(PRGVAR_DDENS_ID)%local(n)%val, & ! (in)
          this%elem_optr, this%mesh3D%lcmesh_list(n), this%mesh3D%refElem3D ) ! (in)
      end do
    end do
    return
  end subroutine PhysPreProcModalFilter_Operate

!-- 
  !> Initialize an object to represent quasi-global filtering operation before evaluating physics tendencies
!OCL SERIAL
  subroutine PhysPreProcGlobalFilter_Init( this, &
    GLFilterOptrType, GLFilterShape, GLFilterWidthFac,  &
    Nnode_h1D_reconst, &
    mesh3D )
    implicit none
    class(PhysPreProcGlobalFilter), intent(inout) :: this
    character(len=*), intent(in) :: GLFilterOptrType
    character(len=*), intent(in) :: GLFilterShape
    real(RP), intent(in) :: GLFilterWidthFac
    integer, intent(in) :: Nnode_h1D_reconst
    class(MeshBase3D), intent(in), target :: mesh3D

    integer :: iq
    character(len=H_SHORT) :: varname
    !----------------------------------------------------
    
    this%mesh3D => mesh3D
    
    call this%gFilter_prgvar%Init( GLFilterOptrType, GLFilterShape, GLFilterWidthFac, Nnode_h1D_reconst, &
      PRGVAR_NUM, 0, 0, mesh3D )
    call this%gFilter_trcvar%Init( GLFilterOptrType, GLFilterShape, GLFilterWidthFac, Nnode_h1D_reconst, &
      max(1,QA), 0, 0, mesh3D )
    
    allocate( this%RHOQ(0:QA) )
    do iq=0, QA
      write(varname,'(a,I2.2)') 'RHOQ', iq
      call this%RHOQ(iq)%Init( trim(varname), "kg/m3", mesh3D )
    end do
    return
  end subroutine PhysPreProcGlobalFilter_Init

  !> Finalize an object to represent quasi-global filtering operation before evaluating physics tendencies
!OCL SERIAL
  subroutine PhysPreProcGlobalFilter_Final( this )
    use scale_element_hexahedral, only: HexahedralElement
    implicit none
    class(PhysPreProcGlobalFilter), intent(inout) :: this

    integer :: iq
    !----------------------------------------------------
    call this%gFilter_prgvar%Final()
    call this%gFilter_trcvar%Final()

    do iq=0, QA
      call this%RHOQ(iq)%Final()
    end do
    deallocate( this%RHOQ )
    return
  end subroutine PhysPreProcGlobalFilter_Final  

  !> Perform quasi-global filtering operation before evaluating physics tendencies
!OCL SERIAL
  subroutine PhysPreProcGlobalFilter_Operate( this, &
    PRG_VARS, QTRC_VARS, AUX_VARS,   &
    PRG_VARS0, QTRC_VARS0, AUX_VARS0 )
    use scale_element_hexahedral, only: HexahedralElement
    use scale_tracer, only: QA
    implicit none
    class(PhysPreProcGlobalFilter), intent(inout) :: this
    type(MeshField3D), intent(inout) :: PRG_VARS(PRGVAR_NUM)
    type(MeshField3D), intent(inout) :: AUX_VARS(AUXVAR_NUM)
    type(MeshField3D), intent(inout) :: QTRC_VARS(0:QA)
    type(MeshField3D), intent(in) :: PRG_VARS0(PRGVAR_NUM)
    type(MeshField3D), intent(in) :: AUX_VARS0(AUXVAR_NUM)
    type(MeshField3D), intent(in) :: QTRC_VARS0(0:QA)

    integer :: n
    integer :: iq, iqs
    !----------------------------------------------------

    do n=1, this%mesh3D%LOCAL_MESH_NUM
      call copy_data( &
        PRG_VARS(PRGVAR_DDENS_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMX_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMY_ID)%local(n)%val, PRG_VARS(PRGVAR_MOMZ_ID)%local(n)%val, PRG_VARS(PRGVAR_THERM_ID)%local(n)%val,      & ! (inout)
        AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val, AUX_VARS(AUXVAR_PRESHYDRO_ID)%local(n)%val,                                                                                                                   & ! (inout)
        PRG_VARS0(PRGVAR_DDENS_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMX_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMY_ID)%local(n)%val, PRG_VARS0(PRGVAR_MOMZ_ID)%local(n)%val, PRG_VARS0(PRGVAR_THERM_ID)%local(n)%val, & ! (in)
        AUX_VARS0(AUXVAR_DENSHYDRO_ID)%local(n)%val, AUX_VARS0(AUXVAR_PRESHYDRO_ID)%local(n)%val,                                                                                                                 & ! (in)
        this%mesh3D%lcmesh_list(n), this%mesh3D%refElem3D ) ! (in)
    end do

    call this%gFilter_prgvar%Apply( PRG_VARS, this%mesh3D )

    if ( QA == 0 ) then
      iqs = 0
    else
      iqs = 1
    end if
    do iq=iqs, QA
      do n=1, this%mesh3D%LOCAL_MESH_NUM
        call calc_data_rhoq( this%RHOQ(iq)%local(n)%val,                                                                    & ! (inout)
          QTRC_VARS0(iq)%local(n)%val, AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val, PRG_VARS0(PRGVAR_DDENS_ID)%local(n)%val, & ! (in)
          this%mesh3D%lcmesh_list(n), this%mesh3D%refElem3D ) ! (in)
      end do
    end do

    call this%gFilter_trcvar%Apply( this%RHOQ(iqs:QA), this%mesh3D )

    do iq=iqs, QA
      do n=1, this%mesh3D%LOCAL_MESH_NUM
        call calc_data_rhoq2qtrc( QTRC_VARS(iq)%local(n)%val,                                                             & ! (inout)
          this%RHOQ(iq)%local(n)%val, AUX_VARS(AUXVAR_DENSHYDRO_ID)%local(n)%val, PRG_VARS(PRGVAR_DDENS_ID)%local(n)%val, & ! (in)
          this%mesh3D%lcmesh_list(n), this%mesh3D%refElem3D ) ! (in)
      end do
    end do

    return
  end subroutine PhysPreProcGlobalFilter_Operate

!-- private --------------------------------
!OCL SERIAL  
  subroutine copy_data( &
    DDENS, MOMX, MOMY, MOMZ, THERM, DENS_hyd, PRES_hyd, &
    DDENS0, MOMX0, MOMY0, MOMZ0, THERM0, DENS_hyd0, PRES_hyd0, &
    lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DDENS(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: THERM(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMZ0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: THERM0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd0(elem%Np,lmesh%NeA)

    integer :: ke
    !--------------------------------------------------

    !$omp parallel do private(ke)
    do ke=lmesh%NeS, lmesh%NeE
      DDENS(:,ke) = DDENS0(:,ke)
      MOMX (:,ke) = MOMX0 (:,ke)
      MOMY (:,ke) = MOMY0 (:,ke)
      MOMZ (:,ke) = MOMZ0 (:,ke)
      THERM(:,ke) = THERM0(:,ke)

      DENS_hyd(:,ke) = DENS_hyd0(:,ke)
      PRES_hyd(:,ke) = PRES_hyd0(:,ke)
    end do    
    return
  end subroutine copy_data
!OCL SERIAL  
  subroutine copy_data_qtrc( QTRC, &
    QTRC0, DENS_hyd, DDENS, &
    lmesh, elem             )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout)  :: QTRC(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: QTRC0(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS(elem%Np,lmesh%NeA)

    integer :: ke
    !--------------------------------------------------

    !$omp parallel do private(ke)
    do ke=lmesh%NeS, lmesh%NeE
      QTRC(:,ke) = QTRC0(:,ke)
    end do    
    return
  end subroutine copy_data_qtrc

!OCL SERIAL  
  subroutine calc_data_rhoq( RHOQ, &
    QTRC0, DENS_hyd, DDENS0, &
    lmesh, elem              )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout)  :: RHOQ(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: QTRC0(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS0(elem%Np,lmesh%NeA)

    integer :: ke
    !--------------------------------------------------

    !$omp parallel do private(ke)
    do ke=lmesh%NeS, lmesh%NeE
      RHOQ(:,ke) = QTRC0(:,ke) * ( DDENS0(:,ke) + DENS_hyd(:,ke) )
    end do
    return
  end subroutine calc_data_rhoq

!OCL SERIAL  
  subroutine calc_data_rhoq2qtrc( QTRC, &
    RHOQ0, DENS_hyd, DDENS0, &
    lmesh, elem              )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout)  :: QTRC(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: RHOQ0(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS0(elem%Np,lmesh%NeA)

    integer :: ke
    !--------------------------------------------------

    !$omp parallel do private(ke)
    do ke=lmesh%NeS, lmesh%NeE
      QTRC(:,ke) = RHOQ0(:,ke) / ( DDENS0(:,ke) + DENS_hyd(:,ke) )
    end do
    return
  end subroutine calc_data_rhoq2qtrc

!-
!OCL SERIAL  
  subroutine apply_modal_filter_core( &
    DDENS, MOMX, MOMY, MOMZ, THERM, DENS_hyd, PRES_hyd, &
    DDENS0, MOMX0, MOMY0, MOMZ0, THERM0, DENS_hyd0, PRES_hyd0, &
    elem_optr, lmesh, elem )
    use scale_atm_dyn_dgm_modalfilter, only: &
        atm_dyn_dgm_modalfilter_apply
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: DDENS(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMY(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: MOMZ(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: THERM(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DDENS0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMZ0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: THERM0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd0(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd0(elem%Np,lmesh%NeA)
    class(ElementOperationTensorProd3D), intent(in) :: elem_optr

    integer :: ke
    !--------------------------------------------------

    !$omp parallel do private(ke)
    do ke=lmesh%NeS, lmesh%NeE
      DDENS(:,ke) = DDENS0(:,ke)
      MOMX (:,ke) = MOMX0 (:,ke)
      MOMY (:,ke) = MOMY0 (:,ke)
      MOMZ (:,ke) = MOMZ0 (:,ke)
      THERM(:,ke) = THERM0(:,ke)

      DENS_hyd(:,ke) = DENS_hyd0(:,ke)
      PRES_hyd(:,ke) = PRES_hyd0(:,ke)
    end do    
    call atm_dyn_dgm_modalfilter_apply( DDENS, MOMX, MOMY, MOMZ, THERM, &
        lmesh, elem, elem_optr, .true. )
    return
  end subroutine apply_modal_filter_core

!OCL SERIAL  
  subroutine apply_modal_filter_qtrc_core( QTRC, &
    QTRC0, DENS_hyd, DDENS, &
    elem_optr, lmesh, elem  )
    use scale_atm_dyn_dgm_modalfilter, only: &
        atm_dyn_dgm_tracer_modalfilter_apply
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout)  :: QTRC(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: QTRC0(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DDENS(elem%Np,lmesh%NeA)
    class(ElementOperationTensorProd3D), intent(in) :: elem_optr

    integer :: ke
    !--------------------------------------------------

    !$omp parallel do private(ke)
    do ke=lmesh%NeS, lmesh%NeE
      QTRC(:,ke) = QTRC0(:,ke)
    end do    
    call atm_dyn_dgm_tracer_modalfilter_apply( QTRC, & ! (inout)
      DENS_hyd, DDENS, lmesh, elem, elem_optr ) ! (in)
    return
  end subroutine apply_modal_filter_qtrc_core
end module mod_atmos_phy_preproc
