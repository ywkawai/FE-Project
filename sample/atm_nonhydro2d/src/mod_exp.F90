!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_exp
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc

  use scale_const, only: &
    PI => CONST_PI,       &
    GRAV => CONST_GRAV,   &
    Rdry => CONST_Rdry,   &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00
    
  use scale_meshfield_base, only: MeshField2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D

  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  type, abstract, public :: experiment
    character(len=H_SHORT) :: label
  contains
    procedure, public :: Init => exp_Init
    procedure, public :: Final => exp_Final
    procedure, public :: SetInitCond => exp_SetInitCond
    procedure(exp_SetInitCond_lc), private, deferred :: setInitCond_lc
  end type

  abstract interface
    subroutine exp_SetInitCond_lc( &
      this, DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
      x, z, dom_xmin, dom_xmax, dom_zmin, dom_zmax, lcmesh, elem )

      import experiment
      import LocalMesh2D 
      import QuadrilateralElement
      import RP

      class(experiment), intent(inout) :: this
      type(LocalMesh2D), intent(in) :: lcmesh
      class(QuadrilateralElement), intent(in) :: elem
      real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
      real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: x(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: z(elem%Np,lcmesh%NeA)
      real(RP), intent(in) :: dom_xmin, dom_xmax
      real(RP), intent(in) :: dom_zmin, dom_zmax
    end subroutine exp_SetInitCond_lc
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  subroutine exp_Init( this, exp_name )
    implicit none
    class(experiment), intent(inout) :: this
    character(len=*), intent(in) :: exp_name
    !----------------------------------------------------------------------

    this%label = exp_name
    return
  end subroutine exp_Init

  subroutine exp_Final( this )
    implicit none
    class(experiment), intent(inout) :: this
    !----------------------------------------------------------------------

    return
  end subroutine exp_Final
  
  subroutine exp_SetInitCond( this )
    use mod_atmos_vars, only: &
      DENS_hydro, PRES_hydro,   &
      DDENS, MOMX, MOMZ, DRHOT
    use mod_atmos_mesh, only: &
      mesh, refElem
    implicit none

    class(experiment), intent(inout) :: this

    integer :: n
    type(LocalMesh2D), pointer :: lcmesh
    !----------------------------------------------------------------------
    
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      call this%setInitCond_lc( &
        DENS_hydro%local(n)%val, PRES_hydro%local(n)%val,                                   & ! (out)
        DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,       & ! (out)
        lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),                                         & ! (in)
        mesh%xmin_gl, mesh%xmax_gl, mesh%ymin_gl, mesh%ymax_gl, lcmesh, refElem )             ! (in) 
    end do

    return
  end subroutine exp_SetInitCond
  
  subroutine exp_SetInitCond_gravwave( &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
    x, z, dom_xmin, dom_xmax, dom_zmin, dom_zmax, lcmesh, elem )
    
    implicit none

    type(LocalMesh2D), intent(in) :: lcmesh
    class(QuadrilateralElement), intent(in) :: elem
    real(RP), intent(out) :: DENS_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: PRES_hyd(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DDENS(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMX(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: MOMZ(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: DRHOT(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: x(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: z(elem%Np,lcmesh%Ne)
    real(RP), intent(in) :: dom_xmin, dom_xmax
    real(RP), intent(in) :: dom_zmin, dom_zmax
    
    real(RP) :: T0
    real(RP) :: DTemp
    real(RP) :: x_c, r_d
    real(RP) :: U0

    namelist /PARAM_EXP/ &
      T0, DTemp,        &
      x_c, r_d,         &
      U0

    integer :: k
    real(RP) :: DENS00
    real(RP) :: H0, Htop
    real(RP) :: Theta_hyd(elem%Np)
    real(RP) :: DENSb(elem%Np), Tb(elem%Np), DT(elem%Np), DENS(elem%Np)
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.5_RP * (dom_xmin + dom_xmax)
    r_d = (dom_xmax - dom_xmin)/60.0_RP
    T0    = 250.0_RP
    DTemp = 0.01_RP
    U0 = 0.0_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("exp_SetInitCond_gravwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("exp_SetInitCond_gravwave",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---
    H0 = Rdry*T0/Grav
    Htop = dom_zmax - dom_zmin
    DENS00 = PRES00/(H0*Grav)
    do k=1, lcmesh%Ne
      
      DENS_hyd(:,k) = DENS00*exp(-z(:,k)/H0)
      PRES_hyd(:,k) = PRES00*exp(-z(:,k)/H0)
      Theta_hyd(:) = T0*exp(z(:,k)/H0)**(Rdry/CPdry)

      Tb(:) = DTemp * exp(- ((x(:,k) - x_c)/r_d)**2) * sin(PI*z(:,k)/Htop)
      DT(:) = Tb(:) * exp(0.5_RP*z(:,k)/H0)

      DENSb(:) = - PRES00*Tb(:)/(Rdry*T0**2)
      DDENS(:,k) = DENSb(:) * exp(- 0.5_RP*z(:,k)/H0)
      DENS(:) = DENS_hyd(:,k) + DDENS(:,k)

      MOMX(:,k) = DENS(:)*U0
      MOMZ(:,k) = 0.0_RP
      DRHOT(:,k) =   DENS(:) * (T0 + DT(:)) * (PRES00/PRES_hyd(:,k))**(Rdry/CPdry) &
                   - DENS_hyd(:,k)*Theta_hyd(:)
    end do

    return
  end subroutine exp_SetInitCond_gravwave
  
end module mod_exp