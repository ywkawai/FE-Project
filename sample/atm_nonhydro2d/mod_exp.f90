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
  !++ Public procedures
  !
  public :: exp_setup
  public :: exp_finalize
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

  logical :: InitCond_GalerkinProjFlag 

contains
  subroutine exp_setup( exp_name )
    implicit none

    character(len=*), intent(in) :: exp_name
    !----------------------------------------------------------------------

    expName = exp_name
    InitCond_GalerkinProjFlag = .false.

    return
  end subroutine exp_setup

  subroutine exp_finalize()
    implicit none

    return
  end subroutine exp_finalize
  
  subroutine exp_SetInitCond()
    use mod_atmos_vars, only: &
      DENS_hydro, PRES_hydro,   &
      DDENS, MOMX, MOMZ, DRHOT
    use mod_atmos_mesh, only: &
      mesh, refElem
    implicit none

    integer :: n
    type(LocalMesh2D), pointer :: lcmesh

    abstract interface
      subroutine SetInitCond( &
        DENS_hyd, PRES_hyd, DDENS, MOMX, MOMZ, DRHOT, &
        x, z, dom_xmin, dom_xmax, dom_zmin, dom_zmax, lcmesh, elem )
    
        import LocalMesh2D 
        import QuadrilateralElement
        import RP

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
  
  subroutine exp_SetInitCond_densitycurrent( &
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
    
    real(RP) :: THETA0
    real(RP) :: DTHETA
    real(RP) :: x_c, z_c, r_x, r_z

    namelist /PARAM_EXP/ &
      THETA0, DTHETA,            &
      x_c, z_c, r_x, r_z,        &
      InitCond_GalerkinProjFlag


    integer :: k
    real(RP) :: THETA(elem%Np), DENS(elem%Np), dens_zfunc(elem%Np), RHOT(elem%Np)
    real(RP) :: r(elem%Np)

    integer, parameter :: IntrpPolyOrder = 12
    type(QuadrilateralElement) :: elem_intrp
    real(RP), allocatable :: x_intrp(:), z_intrp(:)
    real(RP) :: vx(elem%Nv), vz(elem%Nv)
    real(RP), allocatable :: IntrpMat(:,:), InvV_intrp(:,:)
    integer :: p1, p2, p_, p_intrp

    real(RP), allocatable :: r_intrp(:)
    real(RP), allocatable :: THETA_intrp(:)
  
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
    call elem_intrp%Init( IntrpPolyOrder, .false. )
    allocate( IntrpMat(elem%Np,elem_intrp%Np) )
    allocate( InvV_intrp(elem%Np,elem_intrp%Np) )
    allocate( x_intrp(elem_intrp%Np), z_intrp(elem_intrp%Np) )
  
    allocate( r_intrp(elem_intrp%Np) )
    allocate( THETA_intrp(elem_intrp%Np) )


    InvV_intrp(:,:) = 0.0_RP
    do p2=1, elem%PolyOrder+1
    do p1=1, elem%PolyOrder+1
      p_ = p1 + (p2-1)*(elem%PolyOrder + 1)
      p_intrp = p1 + (p2-1)*(elem_intrp%PolyOrder + 1)
      InvV_intrp(p_,:) = elem_intrp%invV(p_intrp,:)
    end do
    end do    
    IntrpMat(:,:) = matmul(elem%V, InvV_intrp)

    !----
    do k=1, lcmesh%Ne
      vx(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),1)
      vz(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),2)
      x_intrp(:) = vx(1) + 0.5_RP*(elem_intrp%x1(:) + 1.0_RP)*(vx(2) - vx(1))
      z_intrp(:) = vz(1) + 0.5_RP*(elem_intrp%x2(:) + 1.0_RP)*(vz(3) - vz(1))

      dens_zfunc(:) = (1.0_RP - Grav*z(:,k)/(CpDry*THETA0))**(CVdry/Rdry)
      DENS_hyd(:,k) = PRES00/(THETA0*Rdry) * dens_zfunc(:)
      PRES_hyd(:,k) = PRES00 * (Rdry*DENS_hyd(:,k)*THETA0/PRES00)**(CPdry/Cvdry)

      r(:) = min(1.0_RP, sqrt(((x(:,k) - x_c)/r_x)**2 + ((z(:,k) - z_c)/r_z)**2))
      r_intrp(:) = min(1.0_RP, sqrt(((x_intrp(:) - x_c)/r_x)**2 + ((z_intrp(:) - z_c)/r_z)**2))
      
      THETA_intrp(:) = THETA0                          &
         + DTHETA*0.5_RP*(1.0_RP + cos(PI*r_intrp(:))) &
           / (1.0_RP - Grav*z_intrp(:)/(CpDry*THETA0)) 
      !THETA(:) = THETA0 + DTHETA*0.5_RP*(1.0_RP + cos(PI*r(:)))
      THETA(:) = matmul(IntrpMat, THETA_intrp)

      DENS(:) = PRES00/(THETA(:)*Rdry) * dens_zfunc(:)
      DDENS(:,k) = DENS(:) - DENS_hyd(:,k)

      DRHOT(:,k) = DENS(:)*THETA(:) - DENS_hyd(:,k)*THETA0


      MOMX(:,k) = 0.0_RP
      MOMZ(:,k) = 0.0_RP
    end do
    
    call elem_intrp%Final()

    return
  end subroutine exp_SetInitCond_densitycurrent

end module mod_exp