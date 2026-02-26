#include "scaleFElib.h"
module mod_advect2d_kernel
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof

  use scale_localmesh_2d, only: LocalMesh2D
  use scale_element_base, only: ElementBase2D
  use scale_sparsemat, only: SparseMat
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public types & variables
  !
  public :: advect2d_kernel_cal_tend

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: cal_dqdt
  private :: cal_elembnd_flux
  private :: cal_elembnd_flux_gpu  
  private :: cal_dqdt_gpu

contains

  !> Calculate the tendency
  !! dqdt = - Dx ( uq ) - Dy ( vq ) + L ( < vec q>_numflx - vec q ).n
  !!
  subroutine advect2d_kernel_cal_tend( dqdt, & ! (out)
    q_, u_, v_, Dx, Dy, Lift, lmesh, elem    ) ! (in)
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in)  :: v_(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift

    real(RP) :: Fx(elem%Np), Fy(elem%Np), LiftBndFlux(elem%Np)
    real(RP) :: ebnd_flux(elem%NfpTot,lmesh%Ne)

    integer :: ke
    !------------------------------------------------------------------------

    !$acc data create(ebnd_flux)
    call PROF_rapstart( 'cal_tend_bndflux', 2)
#ifndef _OPENACC
    call cal_elembnd_flux( ebnd_flux,                             & ! (out)
      q_, u_, v_, lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem  ) ! (in)
#else
    call cal_elembnd_flux_gpu( ebnd_flux,                         & ! (out)
      q_, u_, v_, lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), & ! (in)
      lmesh%vmapM, lmesh%vmapP, lmesh, elem  ) ! (in)
#endif
    call PROF_rapend( 'cal_tend_bndflux', 2)

    call PROF_rapstart( 'cal_tend_interior', 2)
#ifndef _OPENACC
    call cal_dqdt( dqdt,        & ! (out)
      q_, u_, v_, ebnd_flux,    & ! (in)
      Dx, Dy, Lift, lmesh, elem ) ! (in)
#else
    call cal_dqdt_gpu( dqdt,       & ! (out)
        q_, u_, v_, ebnd_flux,     & ! (in)
        Dx, Dy, Lift, lmesh, elem  ) ! (in)
#endif
    call PROF_rapend( 'cal_tend_interior', 2)
    !$acc end data
    return
  end subroutine advect2d_kernel_cal_tend

!-- private --------------------

  subroutine cal_dqdt( dqdt, q_, u_, v_, ebnd_flux, Dx, Dy, Lift, lmesh, elem )
    use scale_sparsemat, only: sparsemat_matmul
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: v_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ebnd_flux(elem%NfpTot,lmesh%Ne)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift

    integer :: ke
    real(RP) :: Fx(elem%Np)
    real(RP) :: Fy(elem%Np)
    real(RP) :: LiftBndFlux(elem%Np)
    !------------------------------------------------------------------------

    !$omp parallel do private(ke)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul( Dx, q_(:,ke) * u_(:,ke), Fx )
      call sparsemat_matmul( Dy, q_(:,ke) * v_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * ebnd_flux(:,ke), LiftBndFlux )

      dqdt(:,ke) = - (  lmesh%Escale(:,ke,1,1) * Fx(:) &
                      + lmesh%Escale(:,ke,2,2) * Fy(:) &
                      + LiftBndFlux(:) )
    end do
    return
  end subroutine cal_dqdt

  !> Calculate the contribution at element boundaries: 
  !! 0.5 * [ ( [qu]^+ [qu]^- ) - ( [qu]^+ [qu]^- ) ] - [qu]^-
  subroutine cal_elembnd_flux( ebnd_flux,           & ! (out)
      q_, u_, v_, nx, ny, vmapM, vmapP, lmesh, elem ) ! (in)
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem  
    real(RP), intent(out) ::  ebnd_flux(elem%NfpTot,lmesh%Ne) !< Flux at element boundaries
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  v_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)          !< Normal vector at element boundaries
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)          !< Normal vector at element boundaries
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (outside own element) into that of all nodes
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (inside) own element) into that of all nodes

    integer :: ke
    integer :: iP(elem%NfpTot), iM(elem%NfpTot)
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot)
    real(RP) :: qP(elem%NfpTot), qM(elem%NfpTot)
    real(RP) :: alpha(elem%NfpTot)
    !------------------------------------------------------------------------

    !$omp parallel do private(iP,iM,VelP,VelM,qP,qM,alpha)
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      VelM(:) = u_(iM(:)) * nx(:,ke) + v_(iM(:)) * ny(:,ke)
      VelP(:) = u_(iP(:)) * nx(:,ke) + v_(iP(:)) * ny(:,ke)
      qM(:) = q_(iM(:)); qP(:) = q_(iP(:))

      alpha = 0.5_RP * abs( VelP(:) + VelM(:) )
      ebnd_flux(:,ke) = 0.5_RP * ( &  
          (    qP(:) * VelP(:) - qM(:) * VelM(:) ) &
           - alpha(:) * ( qP(:) - qM(:) )          )  
    end do
    return
  end subroutine cal_elembnd_flux

!* -- For GPU ---------------------------------------------------------------
  
  subroutine cal_dqdt_gpu( dqdt, q_, u_, v_, ebnd_flux, Dx, Dy, Lift, lmesh, elem )
    use scale_sparsemat, only: sparsemat_matmul
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: dqdt(elem%Np,lmesh%Ne)
    real(RP), intent(in) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: v_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: ebnd_flux(elem%NfpTot,lmesh%Ne)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Lift

    integer :: ke, p
    real(RP) :: Fx(elem%Np)
    real(RP) :: Fy(elem%Np)
    real(RP) :: LiftBndFlux(elem%Np)
    !------------------------------------------------------------------------

    !$acc parallel loop gang private(Fx,Fy,LiftBndFlux) present(dqdt, q_, u_, v_, ebnd_flux, lmesh, elem)
    do ke=lmesh%NeS, lmesh%NeE
      call sparsemat_matmul( Dx, q_(:,ke), u_(:,ke), Fx )
      call sparsemat_matmul( Dy, q_(:,ke), v_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke), ebnd_flux(:,ke), LiftBndFlux )
      !$acc loop vector independent
      do p=1, elem%Np
        dqdt(p,ke) = - (  lmesh%Escale(p,ke,1,1) * Fx(p) &
                        + lmesh%Escale(p,ke,2,2) * Fy(p) &
                        + LiftBndFlux(p) )
      end do
    end do
    return
  end subroutine cal_dqdt_gpu

  subroutine cal_elembnd_flux_gpu( ebnd_flux, & ! (out)
    q_, u_, v_, nx, ny, vmapM, vmapP, lmesh, elem     ) ! (in)
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem  
    real(RP), intent(out) ::  ebnd_flux(elem%NfpTot,lmesh%Ne) !< Flux at element boundaries
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  v_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)          !< Normal vector at element boundaries
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)          !< Normal vector at element boundaries
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (outside own element) into that of all nodes
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)        !< Mapping array to convert the node id for boundary data (inside) own element) into that of all nodes

    integer :: ke, p
    integer :: iP, iM
    real(RP) :: velP, velM
    real(RP) :: qP, qM
    real(RP) :: alpha
    !------------------------------------------------------------------------

    !$acc parallel loop collapse(2) present(ebnd_flux, q_, u_, v_, nx, ny, vmapM, vmapP, lmesh, elem)
    do ke=lmesh%NeS, lmesh%NeE
    do p=1, elem%NfpTot
      iM = vmapM(p,ke); iP = vmapP(p,ke)
      velM = u_(iM) * nx(p,ke) + v_(iM) * ny(p,ke)
      velP = u_(iP) * nx(p,ke) + v_(iP) * ny(p,ke)
      qM = q_(iM); qP = q_(iP)

      alpha = 0.5_RP * abs( velP + velM )
      ebnd_flux(p,ke) = 0.5_RP * ( &  
            ( qP * velP - qM * velM ) &
           - alpha * ( qP - qM )      )  
    end do
    end do
    return
  end subroutine cal_elembnd_flux_gpu
end module mod_advect2d_kernel