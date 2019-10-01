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
  use scale_localmesh_2d, only: LocalMesh2D

  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: exp_Init
  public :: exp_Final
  public :: exp_SetInitCond

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

  character(len=H_SHORT) :: expName
  
  character(*), parameter :: EXP_GRAVITYWAVE_LABEL = 'GravityWave'
  character(*), parameter :: EXP_DENSCURRENT_LABEL = 'DensityCurrent'

contains
  subroutine exp_Init( exp_name )
    implicit none

    character(len=*), intent(in) :: exp_name
    !----------------------------------------------------------------------

    expName = exp_name
  end subroutine exp_Init

  subroutine exp_Final()
  end subroutine exp_Final
  
  subroutine exp_SetInitCond( &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
    mesh )
    implicit none

    type(MeshField2D), intent(inout) :: DENS_hyd
    type(MeshField2D), intent(inout) :: PRES_hyd
    type(MeshField2D), intent(inout) :: DDENS
    type(MeshField2D), intent(inout) :: MOMX
    type(MeshField2D), intent(inout) :: MOMZ
    type(MeshField2D), intent(inout) :: DRHOT
    type(MeshRectDom2D), intent(in), target :: mesh

    integer :: n
    type(LocalMesh2D), pointer :: lcmesh

    abstract interface
      subroutine SetInitCond( &
        DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
        x, z, dom_xmin, dom_xmax, dom_zmin, dom_zmax, lcmesh, elem )
    
        import LocalMesh2D 
        import ElementBase2D
        import RP

        type(LocalMesh2D), intent(in) :: lcmesh
        class(ElementBase2D), intent(in) :: elem
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
      end subroutine SetInitCond
    end interface
    procedure(SetInitCond), pointer :: SetInitCond_p => NULL()

    !----------------------------------------------------------------------

    select case(expName)
    case (EXP_GRAVITYWAVE_LABEL)
      SetInitCond_p => exp_SetInitCond_gravwave
    case (EXP_DENSCURRENT_LABEL)
      SetInitCond_p => exp_SetInitCond_densitycurrent     
    case default
      LOG_ERROR("exp_SetInitCond",*) trim(expName), ' is not supported. Check!'
      call PRC_abort
    end select
    
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)

      call SetInitCond_p( &
        DENS_hyd%local(n)%val, PRES_hyd%local(n)%val,                                       & ! (out)
        DDENS%local(n)%val, MOMX%local(n)%val, MOMZ%local(n)%val, DRHOT%local(n)%val,       & ! (out)
        lcmesh%pos_en(:,:,1), lcmesh%pos_en(:,:,2),                                         & ! (in)
        mesh%xmin_gl, mesh%xmax_gl, mesh%ymin_gl, mesh%ymax_gl, lcmesh, lcmesh%refElem2D )    ! (in) 
    end do

    return
  end subroutine exp_SetInitCond
  
  subroutine exp_SetInitCond_gravwave( &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
    x, z, dom_xmin, dom_xmax, dom_zmin, dom_zmax, lcmesh, elem )
    
    implicit none

    type(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
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
  
  subroutine exp_SetInitCond_densitycurrent( &
    DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
    x, z, dom_xmin, dom_xmax, dom_zmin, dom_zmax, lcmesh, elem )
    
    implicit none

    type(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
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
    
    real(RP) :: THETA0
    real(RP) :: DTHETA
    real(RP) :: x_c, z_c, r_x, r_z

    namelist /PARAM_EXP/ &
      THETA0, DTHETA,        &
      x_c, z_c, r_x, r_z

    integer :: k
    real(RP) :: THETA(elem%Np), DENS(elem%Np), dens_zfunc(elem%Np), RHOT(elem%Np)
    real(RP) :: r(elem%Np)
    integer :: ierr
    !-----------------------------------------------------------------------------

    x_c = 0.0_RP; z_c = 3.0E3_RP
    r_x = 4.0E3_RP; r_z = 2.0E3_RP
    THETA0    = 300.0_RP
    DTHETA    = -15.0_RP

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_EXP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("exp_SetInitCond_densitycurrent",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("exp_SetInitCond_densitycurrent",*) 'Not appropriate names in namelist PARAM_EXP. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_EXP)

    !---

    do k=1, lcmesh%Ne
      dens_zfunc(:) = (1.0_RP - Grav*z(:,k)/(CpDry*THETA0))**(CVdry/Rdry)
      DENS_hyd(:,k) = PRES00/(THETA0*Rdry) * dens_zfunc(:)
      PRES_hyd(:,k) = PRES00 * (Rdry*DENS_hyd(:,k)*THETA0/PRES00)**(CPdry/Cvdry)

      r(:) = min(1.0_RP, sqrt(((x(:,k) - x_c)/r_x)**2 + ((z(:,k) - z_c)/r_z)**2))
      THETA(:) = THETA0 + DTHETA*0.5_RP*(1.0_RP + cos(PI*r(:)))

      DENS(:) = PRES00/(THETA(:)*Rdry) * dens_zfunc(:)
      DDENS(:,k) = DENS(:) - DENS_hyd(:,k)

      DRHOT(:,k) = DENS(:)*THETA(:) - DENS_hyd(:,k)*THETA0


      MOMX(:,k) = 0.0_RP
      MOMZ(:,k) = 0.0_RP
    end do

    return
  end subroutine exp_SetInitCond_densitycurrent

end module mod_exp