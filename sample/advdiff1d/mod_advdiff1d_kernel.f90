#include "scaleFElib.h"
module mod_advdiff1d_kernel
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof

  use scale_localmesh_1d, only: LocalMesh1D
  use scale_element_base, only: ElementBase1D
  use scale_sparsemat, only: SparseMat
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public types & variables
  !
  public :: advdiff1d_kernel_cal_aux
  public :: advdiff1d_kernel_cal_tend

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: cal_dqdt
  private :: cal_elembnd_flux_aux
  private :: cal_elembnd_flux_aux_gpu
  private :: cal_elembnd_flux
  private :: cal_elembnd_flux_gpu  
  private :: cal_dqdt_gpu

contains

  !> Calculate the auxiliary variable for the tendency
  !!
  subroutine advdiff1d_kernel_cal_aux( dqdx, & ! (out)
    q_, Dx, Lift, lmesh, elem      ) ! (in)
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdx(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Lift

    real(RP) :: ebnd_flux(elem%NfpTot,lmesh%Ne)
    !------------------------------------------------------------------------

    !$acc data create(ebnd_flux)
    call PROF_rapstart( 'cal_tend_bndflux_aux', 2)
#ifndef _OPENACC
    call cal_elembnd_flux_aux( ebnd_flux,    & ! (out)
      q_, lmesh%normal_fn(:,:,1),            & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem  ) ! (in)
#else
    call cal_elembnd_flux_aux_gpu( ebnd_flux, & ! (out)
      q_, lmesh%normal_fn(:,:,1),             & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem   ) ! (in)
#endif
    call PROF_rapend( 'cal_tend_bndflux_aux', 2)

    call PROF_rapstart( 'cal_grad_interior', 2)
#ifndef _OPENACC
    call cal_grad( dqdx,    & ! (out)
      q_, ebnd_flux,        & ! (in)
      Dx, Lift, lmesh, elem ) ! (in)
#else
    call cal_grad_gpu( dqdx, & ! (out)
      q_, ebnd_flux,         & ! (in)
      Dx, Lift, lmesh, elem  ) ! (in)
#endif
    call PROF_rapend( 'cal_grad_interior', 2)
    !$acc end data
    return
  end subroutine advdiff1d_kernel_cal_aux

  !> Calculate the tendency
  !! dqdt = - Dx ( uq ) + L ( <u q>_numflx - uq )
  !!
  subroutine advdiff1d_kernel_cal_tend( dqdt,       & ! (out)
    q_, dqdx_, u_, DIFF_COEF, Dx, Lift, lmesh, elem ) ! (in)
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: dqdx_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: DIFF_COEF
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Lift

    real(RP) :: ebnd_flux(elem%NfpTot,lmesh%Ne)
    !------------------------------------------------------------------------

    !$acc data create(ebnd_flux)
    call PROF_rapstart( 'cal_tend_bndflux', 2)
#ifndef _OPENACC
    call cal_elembnd_flux( ebnd_flux,                   & ! (out)
      q_, dqdx_, u_, DIFF_COEF, lmesh%normal_fn(:,:,1), & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem             ) ! (in)
#else
    call cal_elembnd_flux_gpu( ebnd_flux,               & ! (out)
      q_, dqdx_, u_, DIFF_COEF, lmesh%normal_fn(:,:,1), & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem             ) ! (in)
#endif
    call PROF_rapend( 'cal_tend_bndflux', 2)

    call PROF_rapstart( 'cal_tend_interior', 2)
#ifndef _OPENACC
    call cal_dqdt( dqdt,                   & ! (out)
      q_, dqdx_, u_, DIFF_COEF, ebnd_flux, & ! (in)
      Dx, Lift, lmesh, elem                ) ! (in)
#else
    call cal_dqdt_gpu( dqdt,               & ! (out)
      q_, dqdx_, u_, DIFF_COEF, ebnd_flux, & ! (in)
      Dx, Lift, lmesh, elem                ) ! (in)
#endif
    call PROF_rapend( 'cal_tend_interior', 2)
    !$acc end data

    ! !$acc update host( dqdt )
    ! write(*,*) "dqdt=", dqdt(1:elem%Np,lmesh%NeS:lmesh%NeE)
    return
  end subroutine advdiff1d_kernel_cal_tend

!-- private --------------------

  subroutine cal_grad( dqdx, q_, ebnd_flux, Dx, Lift, lmesh, elem )
    use scale_sparsemat, only: sparsemat_matmul
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ebnd_flux(elem%NfpTot,lmesh%Ne)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Lift

    integer :: ke
    real(RP) :: Fx(elem%Np)
    real(RP) :: LiftBndFlux(elem%Np)
    !------------------------------------------------------------------------

    !$omp parallel do private(ke, Fx,LiftBndFlux)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul( Dx, q_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * ebnd_flux(:,ke), LiftBndFlux )

      dqdx(:,ke) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                 + LiftBndFlux(:)
    end do
    return
  end subroutine cal_grad

  subroutine cal_dqdt( dqdt, q_, dqdx_, u_, DIFF_COEF, ebnd_flux, Dx, Lift, lmesh, elem )
    use scale_sparsemat, only: sparsemat_matmul
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dqdx_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DIFF_COEF
    real(RP), intent(in) :: ebnd_flux(elem%NfpTot,lmesh%Ne)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Lift

    integer :: ke
    real(RP) :: Fx(elem%Np)
    real(RP) :: LiftBndFlux(elem%Np)
    !------------------------------------------------------------------------

    !$omp parallel do private(ke, Fx,LiftBndFlux)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul( Dx, q_(:,ke) * u_(:,ke) - DIFF_COEF * dqdx_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * ebnd_flux(:,ke), LiftBndFlux )

      dqdt(:,ke) = - (  lmesh%Escale(:,ke,1,1) * Fx(:) &
                      + LiftBndFlux(:) )
    end do
    return
  end subroutine cal_dqdt

  subroutine cal_elembnd_flux_aux( ebnd_flux, q_, nx, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  ebnd_flux(elem%NfpTot*lmesh%Ne) !< Flux at element boundaries
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)          !< Normal vector at element boundaries
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)        !< Mapping array to convert the node id for boundary data (outside own element) into that of all nodes
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)        !< Mapping array to convert the node id for boundary data (inside) own element) into that of all nodes
     
    integer :: i, iP, iM
    real(RP) :: delVar
    !------------------------------------------------------------------------

    !$omp parallel do private(i, iM, iP, delVar)
    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)
      delVar = 0.5_RP * ( q_(iP) - q_(iM) )
      ebnd_flux(i) = ( 1.0_RP + nx(i) ) * delVar * nx(i)
    end do
    return
  end subroutine cal_elembnd_flux_aux  

  !> Calculate the contribution at element boundaries: 
  !! 0.5 * [ ( [qu]^+ [qu]^- ) - ( [qu]^+ [qu]^- ) ] - [qu]^-
  subroutine cal_elembnd_flux( ebnd_flux,                     & ! (out)
      q_, dqdx_, u_, DIFF_COEF, nx, vmapM, vmapP, lmesh, elem ) ! (in)
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  ebnd_flux(elem%NfpTot,lmesh%Ne) !< Flux at element boundaries
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  dqdx_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DIFF_COEF
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)          !< Normal vector at element boundaries
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (outside own element) into that of all nodes
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (inside) own element) into that of all nodes

    integer :: ke
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: uP(elem%NfpTot), uM(elem%NfpTot)
    real(RP) :: qP(elem%NfpTot), qM(elem%NfpTot)
    real(RP) :: dqdxP(elem%NfpTot), dqdxM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
    !------------------------------------------------------------------------

    !$omp parallel do private(iP,iM,uP,uM,qP,qM,dqdxP,dqdxM,alpha)
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      uM(:) = u_(iM(:)); uP(:) = u_(iP(:))
      qM(:) = q_(iM(:)); qP(:) = q_(iP(:))
      dqdxM(:) = dqdx_(iM(:)); dqdxP(:) = dqdx_(iP(:))

      alpha = 0.5_RP * abs( uP(:) + uM(:) )
      ebnd_flux(:,ke) = 0.5_RP * ( &  
          ( qP(:) * uP(:) - qM(:) * uM(:) ) * nx(:,ke) &
           - alpha(:) * ( qP(:) - qM(:) )              &
           - DIFF_COEF * ( 1.0_RP - nx(:,ke) )         &
              * ( dqdxP(:) - dqdxM(:) ) * nx(:,ke)     )  
    end do
    return
  end subroutine cal_elembnd_flux

!* -- For GPU ---------------------------------------------------------------
  
  subroutine cal_grad_gpu( dqdx, q_, ebnd_flux, Dx, Lift, lmesh, elem )
    use scale_sparsemat, only: sparsemat_matmul
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdx(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ebnd_flux(elem%NfpTot,lmesh%Ne)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Lift

    integer :: ke, p
    real(RP) :: Fx(elem%Np)
    real(RP) :: LiftBndFlux(elem%Np)
    !------------------------------------------------------------------------

    !$acc parallel loop gang private(Fx,LiftBndFlux) present(dqdx, q_, ebnd_flux, lmesh, elem)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul( Dx, q_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke), ebnd_flux(:,ke), LiftBndFlux )
      !$acc loop vector independent
      do p=1, elem%Np
        dqdx(p,ke) = lmesh%Escale(p,ke,1,1) * Fx(p) &
                   + LiftBndFlux(p)
      end do
    end do
    return
  end subroutine cal_grad_gpu

  subroutine cal_dqdt_gpu( dqdt, q_, dqdx_, u_, DIFF_COEF, ebnd_flux, Dx, Lift, lmesh, elem )
    use scale_sparsemat, only: sparsemat_matmul
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: dqdx_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DIFF_COEF
    real(RP), intent(in) :: ebnd_flux(elem%NfpTot,lmesh%Ne)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Lift

    integer :: ke, p
    real(RP) :: Flx(elem%Np)
    real(RP) :: Fx(elem%Np)
    real(RP) :: LiftBndFlux(elem%Np)
    !------------------------------------------------------------------------

    !$acc parallel loop gang private(Flx,Fx,LiftBndFlux) present(dqdt, q_,dqdx_, u_, ebnd_flux, lmesh, elem)
    do ke=lmesh%NeS, lmesh%NeE
      !$acc loop vector independent
      do p=1, elem%Np
        Flx(p) = q_(p,ke) * u_(p,ke) - DIFF_COEF * dqdx_(p,ke)
      end do
      call sparsemat_matmul( Dx, Flx, Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke), ebnd_flux(:,ke), LiftBndFlux )
      !$acc loop vector independent
      do p=1, elem%Np
        dqdt(p,ke) = - (  lmesh%Escale(p,ke,1,1) * Fx(p) &
                        + LiftBndFlux(p) )
      end do
    end do
    return
  end subroutine cal_dqdt_gpu

  subroutine cal_elembnd_flux_aux_gpu( ebnd_flux, q_, nx, vmapM, vmapP, lmesh, elem )
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  ebnd_flux(elem%NfpTot,lmesh%Ne) !< Flux at element boundaries
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)          !< Normal vector at element boundaries
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (outside own element) into that of all nodes
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (inside) own element) into that of all nodes
     
    integer :: ke, p
    integer :: iP, iM
    real(RP) :: qP, qM

    real(RP) :: delVar
    !------------------------------------------------------------------------

    !$acc parallel loop collapse(2) present(ebnd_flux, q_, nx, vmapM, vmapP, lmesh, elem)
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%NfpTot
      iM = vmapM(p,ke); iP = vmapP(p,ke)
      qM = q_(iM); qP = q_(iP)

      delVar = 0.5_RP * ( qP - qM )
      ebnd_flux(p,ke) = ( 1.0_RP + nx(p,ke) ) * delVar * nx(p,ke)
    end do
    end do
    return
  end subroutine cal_elembnd_flux_aux_gpu

  subroutine cal_elembnd_flux_gpu( ebnd_flux,               & ! (out)
    q_, dqdx_, u_, DIFF_COEF, nx, vmapM, vmapP, lmesh, elem ) ! (in)
    implicit none
    class(LocalMesh1D), intent(in) :: lmesh
    class(ElementBase1D), intent(in) :: elem  
    real(RP), intent(out) ::  ebnd_flux(elem%NfpTot,lmesh%Ne) !< Flux at element boundaries
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  dqdx_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DIFF_COEF
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)          !< Normal vector at element boundaries
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (outside own element) into that of all nodes
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (inside) own element) into that of all nodes

    integer :: ke, p
    integer :: iP, iM
    real(RP) :: uP, uM
    real(RP) :: qP, qM
    real(RP) :: dqdxP, dqdxM
    real(RP) :: alpha
    !------------------------------------------------------------------------

    !$acc parallel loop collapse(2) present(ebnd_flux, q_, dqdx_, u_, nx, vmapM, vmapP, lmesh, elem)
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%NfpTot
      iM = vmapM(p,ke); iP = vmapP(p,ke)
      uM = u_(iM); uP = u_(iP)
      qM = q_(iM); qP = q_(iP)
      dqdxM = dqdx_(iM); dqdxP = dqdx_(iP)

      alpha = 0.5_RP * abs( uP + uM )
      ebnd_flux(p,ke) = 0.5_RP * ( &  
          ( qP * uP - qM * uM ) * nx(p,ke) &
           - alpha * ( qP - qM )           &
           - DIFF_COEF * ( 1.0_RP - nx(p,ke) ) &
              * ( dqdxP - dqdxM ) * nx(p,ke)   )
    end do
    end do
    return
  end subroutine cal_elembnd_flux_gpu
end module mod_advdiff1d_kernel