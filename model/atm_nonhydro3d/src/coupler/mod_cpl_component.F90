!-------------------------------------------------------------------------------
!> module Coupler component
!!
!! @par Description
!!          Coupler component module
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_cpl_component
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc
  use scale_model_component, only: ModelComponent

  use scale_mesh_base3d, only: MeshBase3D
  use scale_mesh_base2d, only: MeshBase2D

  use mod_cpl_vars, only: CouplerVars
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !

  !> Derived type to manage coupler component
  type, extends(ModelComponent), public :: CouplerComponent
    type(CouplerVars) :: vars             !< Object to mange variables with coupler component
  contains
    procedure, public :: setup => Coupler_setup
    procedure, public :: evaluate_activation => Coupler_evaluate_activation
    procedure, public :: setup_vars => Coupler_setup_vars
    procedure, public :: calc_tendency => Coupler_calc_tendency
    procedure, public :: update => Coupler_update
    procedure, public :: finalize => Coupler_finalize
  end type CouplerComponent

contains
  !> Setup coupler component 
  subroutine Coupler_setup( this )
    class(CouplerComponent), intent(inout), target :: this
    !------------------------------------------------------------------------------
    return
  end subroutine Coupler_setup

  !> Evaluate whether the coupler component is needed or not
  subroutine Coupler_evaluate_activation( this, &
    ocn_comp )
    implicit none
    class(CouplerComponent), intent(inout) :: this
    class(ModelComponent), intent(in) :: ocn_comp

    logical :: is_needed
    !------------------------------------------------------------------------------

    LOG_INFO( 'CouplerComponent_setup',*) 'Check whether the coupler component is needed'

    if ( ocn_comp%IsActivated()  ) then
      is_needed = .true.
    else
      is_needed = .false.
    end if
    call this%ModelComponent_Init('COUPLER', is_needed )
    return
  end subroutine Coupler_evaluate_activation

  subroutine Coupler_setup_vars( this, &
    atm_mesh, ocn_mesh )
    implicit none
    class(CouplerComponent), intent(inout) :: this
    class(MeshBase3D), intent(in) :: atm_mesh
    class(MeshBase3D), intent(in) :: ocn_mesh

    logical :: ACTIVATE_FLAG
    class(MeshBase2D), pointer :: atm_mesh2D
    class(MeshBase2D), pointer :: ocn_mesh2D
    !------------------------------------------------------------------------------
    
    call PROF_rapstart( 'Coupler_setup_vars', 1)
    call atm_mesh%GetMesh2D( atm_mesh2D )
    call ocn_mesh%GetMesh2D( ocn_mesh2D )

    call this%vars%Init( atm_mesh2D, ocn_mesh2D )
    call PROF_rapend( 'Coupler_setup_vars', 1)
    return
  end subroutine Coupler_setup_vars

  subroutine Coupler_finalize(this)
    class(CouplerComponent), intent(inout) :: this
    !------------------------------------------------------------------------------
    LOG_INFO('CouplerComponent_finalize',*)
    if ( .not. this%IsActivated() ) return
    return
  end subroutine Coupler_finalize


! Dummy subroutines for the coupler component -------------------------------

  subroutine Coupler_calc_tendency( this, force )
    class(CouplerComponent), intent(inout) :: this
    logical, intent(in) :: force
    !------------------------------------------------------------------------------
    return
  end subroutine Coupler_calc_tendency

  subroutine Coupler_update( this )
    class(CouplerComponent), intent(inout) :: this
    !------------------------------------------------------------------------------
    return
  end subroutine Coupler_update

end module mod_cpl_component