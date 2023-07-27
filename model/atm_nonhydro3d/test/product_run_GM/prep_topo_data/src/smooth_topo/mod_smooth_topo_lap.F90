#include "scalelib.h"
module mod_smooth_topo_lap
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
    PRC_abort
  use scale_const, only: &
    RPlanet => CONST_RADIUS, &
    PI => CONST_PI

  use scale_element_base, only: ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D

  use scale_sparsemat, only: &
    SparseMat, sparsemat_matmul  
  use scale_meshfield_base, only: &
    MeshField2D

  use scale_meshfieldcomm_base, only: &
    MeshFieldContainer
  use scale_meshfieldcomm_cubedspheredom2d, only: &
    MeshFieldCommCubedSphereDom2D

  use scale_timeint_rk, only: timeint_rk  

  implicit none
  private

  public :: smooth_topo_lap_Init
  public :: smooth_topo_lap_apply
  public :: smooth_topo_lap_Final

  type(MeshFieldCommCubedSphereDom2D) :: vars_comm
  type(MeshFieldContainer) :: comm_vars_ptr(1)

  type(MeshFieldCommCubedSphereDom2D) :: grad_vars_comm
  type(MeshFieldContainer) :: comm_grad_vars_ptr(2)

  type(MeshField2D), target :: tmp_topography
  type(MeshField2D), target :: grad_topo_alph
  type(MeshField2D), target :: grad_topo_beta

  type(SparseMat) :: Dx, Dy
  type(SparseMat) :: Lift

  type(timeint_rk), allocatable :: tinteg_lc(:)

  real(RP), parameter :: DELT = 0.001_RP
  integer :: SMOOTH_LAPLACIAN_ITER_NUM
  real(RP) :: visc_coef

  real(RP) :: STAB_TAU

contains
!OCL SERIAL
  subroutine smooth_topo_lap_Init( mesh2D )
    implicit none
    class(MeshCubedSphereDom2D), intent(in), target :: mesh2D

    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase2D), pointer :: elem
    integer :: domid

    real(RP) :: dx_eff
    
    integer :: EFOLD_GRID_NUM

    namelist / PARAM_SMOOTH_TOPO_LAP / &
        STAB_TAU, &
        EFOLD_GRID_NUM

    integer :: ierr
    !------------------------------------------------

    STAB_TAU = 1.E-11_RP
    EFOLD_GRID_NUM = 2

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SMOOTH_TOPO_LAP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("SMOOTH_TOPO_LAP_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("SMOOTH_TOPO_LAP_setup",*) 'Not appropriate names in namelist PARAM_SMOOTH_TOPO_LAP. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_SMOOTH_TOPO_LAP)

    !--

    call tmp_topography%Init( "tmp_topo", "m", mesh2D )
    call grad_topo_alph%Init( "dalpha_topo", "m", mesh2D )
    call grad_topo_beta%Init( "dbeta_topo", "m", mesh2D )

    !-- Initialize data commnication
    call vars_comm%Init( 1, 0, mesh2D )
    comm_vars_ptr(1)%field2d => tmp_topography

    call grad_vars_comm%Init( 0, 1, mesh2D )
    comm_grad_vars_ptr(1)%field2d => grad_topo_alph
    comm_grad_vars_ptr(2)%field2d => grad_topo_beta

    elem => mesh2D%refElem2D
    call Dx%Init( elem%Dx1, storage_format='ELL' )
    call Dy%Init( elem%Dx2, storage_format='ELL' )
    call Lift%Init( elem%Lift, storage_format='ELL' )

    allocate( tinteg_lc(mesh2D%LOCAL_MESH_NUM) )

    do domid = 1, mesh2D%LOCAL_MESH_NUM
      lcmesh2D => mesh2D%lcmesh_list(domid)

      call tinteg_lc(domid)%Init( 'ERK_SSP_4s3o', DELT, 1,  &
        2, (/ lcmesh2D%refElem%Np, lcmesh2D%NeA /)  )
    end do

    dx_eff = sqrt( 4.0_RP * PI * RPlanet**2 / ( 6.0_RP * mesh2D%NeGX * mesh2D%NeGY * elem%Np ) )
    visc_coef =  ( dx_eff / ( PI ) )**2 / 1.0_RP

    SMOOTH_LAPLACIAN_ITER_NUM = int( dble(0.5_RP * EFOLD_GRID_NUM)**2 / DELT )

    return
  end subroutine smooth_topo_lap_Init

!OCL SERIAL
  subroutine smooth_topo_lap_apply( sm_topo, ori_topo )
    implicit none

    type(MeshField2D), intent(inout), target :: sm_topo
    type(MeshField2D), intent(in), target :: ori_topo

    class(MeshBase2D), pointer :: mesh
    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: domid, ke

    integer :: itr
    integer :: rkstage
    integer :: tintbuf_ind
    !------------------------------------------------

    mesh => sm_topo%mesh

    do domid = 1, mesh%LOCAL_MESH_NUM
      lcmesh2D => mesh%lcmesh_list(domid)

      !$omp parallel do
      do ke = lcmesh2D%NeS, lcmesh2D%NeE
        tmp_topography%local(domid)%val(:,ke) = ori_topo%local(domid)%val(:,ke)
      end do
    end do    

    !----

    do itr = 1, SMOOTH_LAPLACIAN_ITER_NUM
      LOG_INFO("smooth_topo_apply_laplacian",*) "itr=", itr, "/", SMOOTH_LAPLACIAN_ITER_NUM

      do rkstage = 1, tinteg_lc(1)%nstage 
        call vars_comm%Put( comm_vars_ptr, 1 )
        call vars_comm%Exchange()
        call vars_comm%Get( comm_vars_ptr, 1 )

        do domid = 1, mesh%LOCAL_MESH_NUM
          lcmesh2D => mesh%lcmesh_list(domid)

          call calculate_grad( grad_topo_alph%local(domid)%val, grad_topo_beta%local(domid)%val, & ! (out)
            tmp_topography%local(domid)%val, lcmesh2D, lcmesh2D%refElem2D )
        end do

        call grad_vars_comm%Put( comm_grad_vars_ptr, 1 )
        call grad_vars_comm%Exchange()
        call grad_vars_comm%Get( comm_grad_vars_ptr, 1 )

        do domid = 1, mesh%LOCAL_MESH_NUM
          lcmesh2D => mesh%lcmesh_list(domid)

          tintbuf_ind = tinteg_lc(domid)%tend_buf_indmap(rkstage)

          call calculate_div( tinteg_lc(domid)%tend_buf2D_ex(:,:,1,tintbuf_ind),     & ! (out)
            grad_topo_alph%local(domid)%val, grad_topo_beta%local(domid)%val,        & ! (in)
            tmp_topography%local(domid)%val, visc_coef, lcmesh2D, lcmesh2D%refElem2D ) ! (in)

          call tinteg_lc(domid)%Advance( rkstage, tmp_topography%local(domid)%val, 1, &
            1, lcmesh2D%refElem2D%Np, lcmesh2D%NeS, lcmesh2D%NeE )
        end do
      end do
    end do

    do domid = 1, mesh%LOCAL_MESH_NUM
      lcmesh2D => mesh%lcmesh_list(domid)

      !$omp parallel do
      do ke = lcmesh2D%NeS, lcmesh2D%NeE
        sm_topo%local(domid)%val(:,ke) = tmp_topography%local(domid)%val(:,ke)
      end do
    end do        
    
    return
  end subroutine smooth_topo_lap_apply

!OCL SERIAL
  subroutine smooth_topo_lap_Final()
    implicit none

    integer :: domid
    !------------------------------------------------

    do domid = 1, size(tinteg_lc)
      call tinteg_lc(domid)%Final()
    end do
    deallocate( tinteg_lc )

    call vars_comm%Final()
    call grad_vars_comm%Final()

    call tmp_topography%Final()
    call grad_topo_alph%Final()
    call grad_topo_beta%Final()    

    call Dx%Final(); call Dy%Final()
    call Lift%Final()

    return
  end subroutine smooth_topo_lap_Final

!-----
  subroutine calculate_div( div, &
    grad_alph, grad_beta, q, vis_coef, lcmesh, elem )

    implicit none
    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: div(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: grad_alph(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: grad_beta(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: q(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: vis_coef

    integer :: ke, p

    real(RP) :: Fx(elem%Np)
    real(RP) :: Fy(elem%Np)
    real(RP) :: LiftDelFlx(elem%Np)
    real(RP) :: dqdx1(elem%Np)
    real(RP) :: dqdx2(elem%Np)

    real(RP) :: bnd_flux(elem%NfpTot,lcmesh%Ne)
    !--------------------------------------------

    call calc_bnd_flux_div( bnd_flux, &
      grad_alph, grad_beta, q, vis_coef,                              &
      lcmesh%Gsqrt, lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), &
      lcmesh%VMapM, lcmesh%VMapP, lcmesh, elem          )

    !$omp parallel do private( ke, Fx, Fy, LiftDelFlx, dqdx1, dqdx2 )
    do ke=lcmesh%NeS, lcmesh%NeE
      call sparsemat_matmul( Dx, grad_alph(:,ke), Fx )
      call sparsemat_matmul( Dy, grad_beta(:,ke), Fy )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * bnd_flux(:,ke), LiftDelFlx)

      div(:,ke) = vis_coef * &
                  ( lcmesh%Escale(:,ke,1,1) * Fx(:) &
                  + lcmesh%Escale(:,ke,2,2) * Fy(:) &
                  + LiftDelFlx(:) ) / lcmesh%Gsqrt(:,ke)
    end do

    return
  end subroutine calculate_div

!OCL SERIAL
  subroutine calc_bnd_flux_div( bnd_flux, &
    grad_alph, grad_beta, q, vis_coef,    &
    Gsqrt, nx, ny, vmapM, vmapP,          &
    lcmesh, elem )
    implicit none

    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: bnd_flux(elem%NfpTot,lcmesh%Ne)
    real(RP), intent(in) :: grad_alph(elem%Np*lcmesh%NeA)
    real(RP), intent(in) :: grad_beta(elem%Np*lcmesh%NeA)
    real(RP), intent(in) :: q(elem%Np*lcmesh%NeA)
    real(RP), intent(in) :: vis_coef
    real(RP), intent(in) :: Gsqrt(elem%Np*lcmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lcmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lcmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lcmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lcmesh%Ne)

    integer :: ke
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    !--------------------------------------------

    !$omp parallel do private( iM, iP )
    do ke=lcmesh%NeS, lcmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      bnd_flux(:,ke) = 0.5_RP * ( &
      !   ( grad_alph(iP) - grad_alph(iM) ) * ( 1.0_RP + nx(:,ke) ) &
      ! + ( grad_beta(iP) - grad_beta(iM) ) * ( 1.0_RP + ny(:,ke) ) &
        ( grad_alph(iP) - grad_alph(iM) ) * nx(:,ke)  &
      + ( grad_beta(iP) - grad_beta(iM) ) * ny(:,ke)  &
      + STAB_TAU * ( Gsqrt(iM) * q(iP) - Gsqrt(iM) * q(iM) ) )
    end do
    
    return
  end subroutine calc_bnd_flux_div

  subroutine calculate_grad( grad_alph, grad_beta, q, lcmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: grad_alph(elem%Np,lcmesh%NeA)
    real(RP), intent(out) :: grad_beta(elem%Np,lcmesh%NeA)
    real(RP), intent(in) :: q(elem%Np,lcmesh%NeA)

    integer :: ke, p

    real(RP) :: Fx(elem%Np)
    real(RP) :: Fy(elem%Np)
    real(RP) :: LiftDelFlx(elem%Np)
    real(RP) :: dqdx1(elem%Np)
    real(RP) :: dqdx2(elem%Np)

    real(RP) :: bnd_flux(elem%NfpTot,lcmesh%Ne,2)
    !--------------------------------------------

    call calc_bnd_flux_grad( bnd_flux, &
      q, lcmesh%normal_fn(:,:,1), lcmesh%normal_fn(:,:,2), &
      lcmesh%VMapM, lcmesh%VMapP, lcmesh, elem )

    !$omp parallel do private( ke, Fx, Fy, LiftDelFlx, dqdx1, dqdx2 )
    do ke=lcmesh%NeS, lcmesh%NeE
      call sparsemat_matmul( Dx, q(:,ke), Fx )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * bnd_flux(:,ke,1), LiftDelFlx)
      dqdx1(:) = lcmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul( Dy, q(:,ke), Fy )
      call sparsemat_matmul( Lift, lcmesh%Fscale(:,ke) * bnd_flux(:,ke,2), LiftDelFlx)
      dqdx2(:) = lcmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:)

      grad_alph(:,ke) = lcmesh%Gsqrt(:,ke) * ( lcmesh%GIJ(:,ke,1,1) * dqdx1(:) + lcmesh%GIJ(:,ke,2,1) * dqdx2(:) )
      grad_beta(:,ke) = lcmesh%Gsqrt(:,ke) * ( lcmesh%GIJ(:,ke,2,1) * dqdx1(:) + lcmesh%GIJ(:,ke,2,2) * dqdx2(:) )
    end do

    return
  end subroutine calculate_grad

!OCL SERIAL
  subroutine calc_bnd_flux_grad( bnd_flux, &
    q, nx, ny, vmapM, vmapP, lcmesh, elem )
    implicit none

    class(LocalMesh2D), intent(in) :: lcmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: bnd_flux(elem%NfpTot,lcmesh%Ne,2)
    real(RP), intent(in) :: q(elem%Np*lcmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lcmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lcmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lcmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lcmesh%Ne)

    integer :: ke
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: del(elem%NfpTot)
    !--------------------------------------------

    !$omp parallel do private( iM, iP, del )
    do ke=lcmesh%NeS, lcmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)

      del(:) = 0.5_RP * ( q(iP) - q(iM) )

      ! bnd_flux(:,ke,1) = del(:) * ( 1.0_RP - nx(:,ke) )
      ! bnd_flux(:,ke,2) = del(:) * ( 1.0_RP - ny(:,ke) )
      bnd_flux(:,ke,1) = del(:) * nx(:,ke)
      bnd_flux(:,ke,2) = del(:) * ny(:,ke)
    end do
    
    return
  end subroutine calc_bnd_flux_grad

end module mod_smooth_topo_lap
