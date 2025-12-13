#include "scaleFElib.h"
module mod_poisson2d_smoother
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale_prc
  use scale_io

  use scale_sparsemat, only: &
    SparseMat, SparseMat_matmul
  
  use scale_element_base, only: ElementBase2D
  use scale_element_quadrilateral, only: QuadrilateralElement
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_mesh_base2d, only: MeshBase2D

  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  public :: Poisson2d_smoother_Init
  public :: Poisson2d_smoother_Final

  public :: Poisson2d_smoother_advance_itr_1step
  ! public :: poisson_smoother_evaluate_error_norm

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !

contains
  !> Initialization
  subroutine Poisson2d_smoother_Init()
    implicit none
    !---------------------------------------------------------------------------
    return
  end subroutine Poisson2d_smoother_Init

  !> Finalization
  subroutine Poisson2d_smoother_Final()
    implicit none
    !---------------------------------------------------------------------------
    return
  end subroutine Poisson2d_smoother_Final

  !> Smoother
  subroutine Poisson2d_smoother_advance_itr_1step( &
    q, res, & 
    f, qx, qy,      &
    itr, var_comm, aux_comm, Dx, Dy, mesh2D )
    implicit none

    class(MeshField2D), intent(inout), target :: q
    class(MeshField2D), intent(inout) :: res
    class(MeshField2D), intent(in) :: f
    class(MeshField2D), intent(inout), target :: qx
    class(MeshField2D), intent(inout), target :: qy
    integer, intent(in) :: itr
    class(MeshFieldCommRectDom2D), intent(inout) :: var_comm
    class(MeshFieldCommRectDom2D), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    class(MeshBase2D), intent(in), target :: mesh2D

    integer :: ldomID
    class(LocalMesh2D), pointer :: lmesh2D

    type(MeshFieldContainer) :: var_comm_list(1)
    type(MeshFieldContainer) :: aux_comm_list(2) 
    !---------------------------------------------------------------------------

    var_comm_list(1)%field2d => q
    aux_comm_list(1)%field2d => qx
    aux_comm_list(2)%field2d => qy

    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    !-
    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldomID)
      call cal_grad_lc( qx%local(ldomID)%val, qy%local(ldomID)%val,   & ! (out)
        1, q%local(ldomID)%val, Dx, Dy, lmesh2D, lmesh2D%refElem2D    ) ! (in)
    end do

    call aux_comm%Put( aux_comm_list, 1 )
    call aux_comm%Exchange()
    call aux_comm%Get( aux_comm_list, 1 )

    !-
    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldomID)
      call cal_q_lc( q%local(ldomID)%val, res%local(ldomID)%val,                     &
        1, .false., f%local(ldomID)%val, qx%local(ldomID)%val, qy%local(ldomID)%val, &
        lmesh2D%VMapM, lmesh2D%VMapP, lmesh2D, lmesh2D%refElem2D )
    end do

    
    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    return
  end subroutine Poisson2d_smoother_advance_itr_1step

!- private
  subroutine cal_grad_lc( qx, qy, &
    color_id, q, Dx, Dy, lmesh, elem )
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qy(elem%Np,lmesh%NeA)
    integer, intent(in) :: color_id
    real(RP), intent(in) :: q(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy

    integer :: ke
    real(RP) :: q1(elem%Np)
    real(RP) :: q2(elem%Np)
    !----------------------------
    !$omp parallel do private(qx, qy)
    do ke=lmesh%NeS, lmesh%NeE
      call Sparsemat_matmul(Dx, q(:,ke), q1)
      call Sparsemat_matmul(Dy, q(:,ke), q2)

      qx(:,ke) = lmesh%G_ij(:,ke,1,1) * q1(:) +  lmesh%G_ij(:,ke,2,1) * q2(:)
      qy(:,ke) = lmesh%G_ij(:,ke,1,2) * q1(:) +  lmesh%G_ij(:,ke,2,2) * q2(:)
    end do
    return
  end subroutine cal_grad_lc

  subroutine cal_q_lc( q, res,                &
    color_id, cal_res, rhs, qx, qy,           &
    vmapM, vmapP, lmesh, elem                 )
    use scale_linalgebra, only: Linalgebra_SolveLinEq
    implicit none
    class(LocalMesh2D), intent(in) :: lmesh
    class(ElementBase2D), intent(in) :: elem
    real(RP), intent(inout) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(out) :: res(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rhs(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np*lmesh%NeA)    
    integer, intent(in) :: color_id
    logical, intent(in) :: cal_res
    integer, intent(in) :: vmapM(elem%Nfp,elem%Nfaces,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%Nfp,elem%Nfaces,lmesh%Ne)

    integer :: ke, ke2
    integer :: i, j

    integer :: p

    integer :: is
    integer :: i_int

    integer :: f1, f2

    real(RP) :: E11, E22

    real(RP) :: Dx(elem%Np,elem%Np)
    real(RP) :: Dy(elem%Np,elem%Np)
    real(RP) :: Sx(elem%Np,elem%Np)
    real(RP) :: Sy(elem%Np,elem%Np)
    real(RP) :: Dx2(elem%Np,elem%Np)
    real(RP) :: Dy2(elem%Np,elem%Np)

    real(RP) :: GsqrtDx(elem%Np,elem%Np)
    real(RP) :: GsqrtDy(elem%Np,elem%Np)
    real(RP) :: GsqrtI(elem%Np,elem%Np)

    real(RP) :: Dn1(elem%Np,elem%Np)    
    real(RP) :: GsqrtDn1(elem%Np,elem%Np)

    real(RP) :: OPT11(elem%Np,elem%Np)
    real(RP) :: OPT12(elem%Np,elem%Np)
    real(RP) :: mmE(elem%Np,elem%Np)   
    real(RP) :: Lift_(elem%Np,elem%Np)
    real(RP) :: Lift_x(elem%Np,elem%Np)       
    real(RP) :: Lift_y(elem%Np,elem%Np)
    real(RP) :: Lift_x_q(elem%Np)       
    real(RP) :: Lift_y_q(elem%Np)
    real(RP) :: Msrc(elem%Np)

    real(RP) :: massEdge(elem%Np,elem%Np,elem%Nfaces)
    real(RP) :: Emat(elem%Np,elem%Nfp*elem%Nfaces)
    integer :: Fm(elem%Nfp), Fm1(elem%Nfp)
    integer :: vidP(elem%Nfp), vidM(elem%Nfp)
    integer :: id
    real(RP) :: lnx, lny, lsJ    

    real(RP) :: gtau, hinv

    real(RP) :: q_(elem%Np)
    !---------------------------

    Emat(:,:) = matmul(elem%M, elem%Lift)
    massEdge(:,:,:)  = 0d0
    do f1=1, elem%Nfaces
      Fm = elem%Fmask(:,f1)
      massEdge(Fm,Fm,f1) = Emat(Fm, (f1-1)*elem%Nfp+1:f1*elem%Nfp)
    end do

    do j=1, lmesh%NeY
      if (cal_res) then
        is = 1; i_int = 1
      else
        is = 1+mod(j+color_id,2); i_int = 2
      end if      
      do i=is, lmesh%NeX, i_int
        ke = i + (j-1)*lmesh%NeX
        
        E11 = lmesh%Escale(1,ke,1,1)
        E22 = lmesh%Escale(1,ke,2,2)

        Dx(:,:) = E11 * elem%Dx1
        Dy(:,:) = E22 * elem%Dx2
        Sx(:,:) = matmul(transpose(Dx),elem%M)
        Sy(:,:) = matmul(transpose(Dy),elem%M)

        do p=1, elem%Np
          GsqrtI(:,p) = lmesh%Gsqrt(:,ke)
          GsqrtDx(:,p) = GsqrtI(:,p) * ( lmesh%GIJ(:,ke,1,1) * Dx(:,p) + lmesh%GIJ(:,ke,1,2) * Dy(:,p) )
          GsqrtDy(:,p) = GsqrtI(:,p) * ( lmesh%GIJ(:,ke,2,1) * Dx(:,p) + lmesh%GIJ(:,ke,2,2) * Dy(:,p) )
        end do

        Msrc(:) = matmul(elem%M, rhs(:,ke) * lmesh%Gsqrt(:,ke) * lmesh%J(1,ke) )
        OPT11(:,:) = lmesh%J(1,ke) * ( &
          - ( matmul(Sx, GsqrtDx) + matmul(Sy, GsqrtDy) ) &
          )
        

        do f1=1, 4        
          ke2 = lmesh%EToE(ke,f1); f2 = lmesh%EToF(ke,f1)
          
          id = 1 + (f1-1)*elem%Nfp
          lnx = lmesh%normal_fn(id,ke,1)
          lny = lmesh%normal_fn(id,ke,2)
          lsJ = lmesh%sJ(id,ke) 
          Fm1(:) = elem%Fmask(:,f1)
          
          Dn1(:,:) = lnx * Dx(:,:) + lny * Dy(:,:)
          GsqrtDn1(:,:) = lnx * GsqrtDx(:,:) + lny * GsqrtDy(:,:)
          
          mmE(:,:) = lsJ * massEdge(:,:,f1)  
          hinv = lmesh%Fscale(id,ke) 
          gtau = 1.0_RP * dble(elem%Np) * hinv
                  
          Lift_(:,:) = matmul(elem%invM, mmE)
          do p=1, elem%Np
            Lift_x(:,p) = lmesh%Gsqrt(:,ke) * ( lnx * lmesh%GIJ(:,ke,1,1) + lny * lmesh%GIJ(:,ke,1,2) ) * Lift_(:,p)
            Lift_y(:,p) = lmesh%Gsqrt(:,ke) * ( lnx * lmesh%GIJ(:,ke,2,1) + lny * lmesh%GIJ(:,ke,2,2) ) * Lift_(:,p)
          end do
          Lift_x(:,:) = matmul(elem%M, Lift_x)
          Lift_y(:,:) = matmul(elem%M, Lift_y)
          
          OPT11(:,:) = OPT11(:,:) + 0.5_RP * ( &
            - gtau * mmE(:,:) + matmul(mmE, GsqrtDn1) &
            + ( matmul(transpose(Dx),Lift_x) + matmul(transpose(Dy),Lift_y) ) &
            )
          
          Lift_x_q(:) = matmul( Lift_x(:,Fm1), q(VMapP(:,f1,ke)) )
          Lift_y_q(:) = matmul( Lift_y(:,Fm1), q(VMapP(:,f1,ke)) )
          Msrc(:)   = Msrc(:) + 0.5_RP * ( &
            & - gtau*matmul(mmE(:,Fm1), q(VMapP(:,f1,ke)))                       &
            & + matmul(transpose(Dx),Lift_x_q) + matmul(transpose(Dy),Lift_y_q)  &
            & )
          
          Msrc(Fm1) = Msrc(Fm1) &
            & - 0.5_RP * matmul(mmE(Fm1,Fm1), lmesh%Gsqrt(Fm1,ke)*(lnx*qx(VMapP(:,f1,ke)) + lny*qy(VMapP(:,f1,ke))))
        end do 
        
        if (cal_res) then
          res(:,ke) = matmul(elem%invM, Msrc - matmul(OPT11,q(1+(ke-1)*elem%Np:ke*elem%Np)))
        else
          call LinAlgebra_SolveLinEq( OPT11, Msrc, q_ )
          !pold(:) = p(1+(k1-1)*elem%Np:k1*elem%Np)
          q(1 + (ke-1)*elem%Np:ke*elem%Np) = q_(:)
        end if        
      end do
    end do
    return
  end subroutine cal_q_lc  
end module mod_poisson2d_smoother
