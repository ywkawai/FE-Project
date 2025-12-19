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

  !--

contains
  !> Initialization
  subroutine Poisson2d_smoother_Init( mesh )
    implicit none
    class(MeshBase2D), intent(in), target :: mesh
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
!OCL SERIAL
  subroutine Poisson2d_smoother_advance_itr_1step( &
    q, res,                                  & 
    f, qx, qy,                               &
    itr, var_comm, aux_comm, Dx, Dy, mesh2D, &
    cal_res_flag, zero_initial_guess, is_mg_top_level )
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
    logical, intent(in) :: cal_res_flag
    logical, intent(in) :: zero_initial_guess
    logical, intent(in) :: is_mg_top_level

    integer :: ldomID
    class(LocalMesh2D), pointer :: lmesh2D
    integer :: ke

    real(RP) :: res_lcdom_int0
    real(RP) :: res_lcdom_int

    type(MeshFieldContainer) :: var_comm_list(1)
    !---------------------------------------------------------------------------

    if ( zero_initial_guess ) then
      do ldomID=1, mesh2D%LOCAL_MESH_NUM
        lmesh2D => mesh2D%lcmesh_list(ldomID)
        !$omp parallel do
        do ke=lmesh2D%NeS, lmesh2D%NeE
          q%local(ldomID)%val(:,ke) = 0.0_RP
        end do
      end do
    end if

    if ( cal_res_flag ) then
      call Poisson2d_smoother_advance_itr_1step_color( &
        q, res,                                  & 
        f, qx, qy,                               &
        itr, var_comm, aux_comm, Dx, Dy, mesh2D, &
        cal_res_flag, 0, is_mg_top_level )

      res_lcdom_int0 = 0.0_RP
      do ldomID=1, mesh2D%LOCAL_MESH_NUM
        lmesh2D => mesh2D%lcmesh_list(ldomID)
        !$omp parallel do reduction(+:res_lcdom_int0)
        do ke=lmesh2D%NeS, lmesh2D%NeE
          res_lcdom_int0 = res_lcdom_int0 + sum(res%local(ldomID)%val(:,ke)**2)
        end do
      end do
    end if

    ! Red
    call Poisson2d_smoother_advance_itr_1step_color( &
      q, res,                                  & 
      f, qx, qy,                               &
      itr, var_comm, aux_comm, Dx, Dy, mesh2D, &
      .false., 0, is_mg_top_level )

    ! Red
    call Poisson2d_smoother_advance_itr_1step_color( &
      q, res,                                  & 
      f, qx, qy,                               &
      itr, var_comm, aux_comm, Dx, Dy, mesh2D, &
      .false., 1, is_mg_top_level )
      
    if ( cal_res_flag ) then
      call Poisson2d_smoother_advance_itr_1step_color( &
        q, res,                                  & 
        f, qx, qy,                               &
        itr, var_comm, aux_comm, Dx, Dy, mesh2D, &
        cal_res_flag, 0, is_mg_top_level )

      res_lcdom_int = 0.0_RP
      do ldomID=1, mesh2D%LOCAL_MESH_NUM
        lmesh2D => mesh2D%lcmesh_list(ldomID)
        !$omp parallel do reduction(+:res_lcdom_int)
        do ke=lmesh2D%NeS, lmesh2D%NeE
          res_lcdom_int = res_lcdom_int + sum(res%local(ldomID)%val(:,ke)**2)
        end do
      end do
      LOG_INFO("Poisson2d_smoother_advance_itr_1step",*) "itr=", itr, ": lc_res0=", sqrt(res_lcdom_int0), ", lc_res1=", sqrt(res_lcdom_int)
    end if
    
    var_comm_list(1)%field2d => q
    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    return
  end subroutine Poisson2d_smoother_advance_itr_1step

!OCL SERIAL
  subroutine Poisson2d_smoother_advance_itr_1step_color( &
    q, res,                                  & 
    f, qx, qy,                               &
    itr, var_comm, aux_comm, Dx, Dy, mesh2D, &
    cal_res_flag, color_id, is_mg_top_level )
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
    logical, intent(in) :: cal_res_flag
    integer, intent(in) :: color_id
    logical, intent(in) :: is_mg_top_level

    integer :: ldomID
    class(LocalMesh2D), pointer :: lmesh2D

    type(MeshFieldContainer) :: var_comm_list(1)
    type(MeshFieldContainer) :: aux_comm_list(2)
    !---------------------------------------------------------------------------

    !-
    var_comm_list(1)%field2d => q

    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldomID)
      call cal_grad_lc( qx%local(ldomID)%val, qy%local(ldomID)%val,       & ! (out)
        color_id, q%local(ldomID)%val, Dx, Dy, lmesh2D, lmesh2D%refElem2D ) ! (in)
    end do

    !- 
    aux_comm_list(1)%field2d => qx
    aux_comm_list(2)%field2d => qy

    call aux_comm%Put( aux_comm_list, 1 )
    call aux_comm%Exchange()
    call aux_comm%Get( aux_comm_list, 1 )

    do ldomID=1, mesh2D%LOCAL_MESH_NUM
      lmesh2D => mesh2D%lcmesh_list(ldomID)
      call cal_q_lc( q%local(ldomID)%val, res%local(ldomID)%val,                                 &
        color_id, cal_res_flag, is_mg_top_level,f%local(ldomID)%val, qx%local(ldomID)%val, qy%local(ldomID)%val, &
        lmesh2D%VMapM, lmesh2D%VMapP, lmesh2D, lmesh2D%refElem2D )
    end do
    
    return
  end subroutine Poisson2d_smoother_advance_itr_1step_color

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
    !$omp parallel do private(q1, q2)
    do ke=lmesh%NeS, lmesh%NeE
      call Sparsemat_matmul(Dx, q(:,ke), q1)
      call Sparsemat_matmul(Dy, q(:,ke), q2)

      q1(:) = lmesh%Escale(:,ke,1,1) * q1(:)
      q2(:) = lmesh%Escale(:,ke,2,2) * q2(:)

      qx(:,ke) = lmesh%G_ij(:,ke,1,1) * q1(:) &
              +  lmesh%G_ij(:,ke,2,1) * q2(:)
      qy(:,ke) = lmesh%G_ij(:,ke,1,2) * q1(:) &
              +  lmesh%G_ij(:,ke,2,2) * q2(:)
    end do
    return
  end subroutine cal_grad_lc

!OCL SERIAL
  subroutine cal_q_lc0( q, res,                       &
    color_id, cal_res, is_mg_top_level, rhs, qx, qy, &
    vmapM, vmapP, lmesh, elem                        )
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
    logical, intent(in) :: is_mg_top_level
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
    real(RP) :: Sx_tr(elem%Np,elem%Np)
    real(RP) :: Sy_tr(elem%Np,elem%Np)
    real(RP) :: Dx2(elem%Np,elem%Np)
    real(RP) :: Dy2(elem%Np,elem%Np)

    real(RP) :: GsqrtDx(elem%Np,elem%Np)
    real(RP) :: GsqrtDy(elem%Np,elem%Np)
    real(RP) :: GsqrtI(elem%Np,elem%Np)

    real(RP) :: Dn1(elem%Np,elem%Np)    
    real(RP) :: GsqrtDn1(elem%Np,elem%Np)
    real(RP) :: GsqrtDn2(elem%Nfp)

    real(RP) :: OPT11(elem%Np,elem%Np)
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
    real(RP) :: lnx, lny  

    real(RP) :: q_(elem%Np)
    real(RP) :: q_P(elem%Nfp)
    real(RP) :: qx_P(elem%Nfp)
    real(RP) :: qy_P(elem%Nfp)

    real(RP) :: gtau_
    !---------------------------

    Emat(:,:) = matmul(elem%M, elem%Lift)
    massEdge(:,:,:)  = 0.0_RP
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
        Sx_tr(:,:) = E11 * matmul(transpose(Dx),elem%M)
        Sy_tr(:,:) = E22 * matmul(transpose(Dy),elem%M)

        do p=1, elem%Np
          GsqrtI(:,p) = lmesh%Gsqrt(:,ke)
          GsqrtDx(:,p) = GsqrtI(:,p) * ( lmesh%GIJ(:,ke,1,1) * Dx(:,p) + lmesh%GIJ(:,ke,1,2) * Dy(:,p) )
          GsqrtDy(:,p) = GsqrtI(:,p) * ( lmesh%GIJ(:,ke,2,1) * Dx(:,p) + lmesh%GIJ(:,ke,2,2) * Dy(:,p) )
        end do

        if ( is_mg_top_level ) then
          Msrc(:) = matmul(elem%M, rhs(:,ke) * lmesh%Gsqrt(:,ke) )
        else
          Msrc(:) = rhs(:,ke)
        end if

        OPT11(:,:) = - ( matmul(Sx_tr, GsqrtDx) + matmul(Sy_tr, GsqrtDy) )
        

        do f1=1, 4        
          ke2 = lmesh%EToE(ke,f1); f2 = lmesh%EToF(ke,f1)
          
          id = 1 + (f1-1)*elem%Nfp
          lnx = lmesh%normal_fn(id,ke,1)
          lny = lmesh%normal_fn(id,ke,2)
          Fm1(:) = elem%Fmask(:,f1)
          
          Dn1(:,:) = lnx * Dx(:,:) + lny * Dy(:,:)
          GsqrtDn1(:,:) = lnx * GsqrtDx(:,:) + lny * GsqrtDy(:,:)
          
          mmE(:,:) = lmesh%Fscale(id,ke) * massEdge(:,:,f1)  
                           
          Lift_(:,:) = matmul(elem%invM, mmE)
          do p=1, elem%Np
            Lift_x(:,p) = lmesh%Gsqrt(:,ke) * ( lnx * lmesh%GIJ(:,ke,1,1) + lny * lmesh%GIJ(:,ke,1,2) ) * Lift_(:,p)
            Lift_y(:,p) = lmesh%Gsqrt(:,ke) * ( lnx * lmesh%GIJ(:,ke,2,1) + lny * lmesh%GIJ(:,ke,2,2) ) * Lift_(:,p)
          end do
          
          gtau_ = elem%PolyOrder * (elem%PolyOrder + 1 ) * lmesh%Fscale(id,ke) * 0.5_RP * 7.0_RP
          
          OPT11(:,:) = OPT11(:,:) + 0.5_RP * ( &
            - gtau_ * mmE(:,:) + matmul(mmE, GsqrtDn1)         &
            + ( matmul(Sx_tr,Lift_x) + matmul(Sy_tr,Lift_y) ) )
          
          q_P(:) = q(VMapP(:,f1,ke))
          Lift_x_q(:) = matmul(Lift_x(:,Fm1), q_P(:))
          Lift_y_q(:) = matmul(Lift_y(:,Fm1), q_P(:))
          Msrc(:) = Msrc(:) + 0.5_RP * ( &
            + matmul(Sx_tr,Lift_x_q) + matmul(Sy_tr,Lift_y_q) )

          qx_P(:) = qx(VMapP(:,f1,ke))
          qy_P(:) = qy(VMapP(:,f1,ke))
          GsqrtDn2(:) = lmesh%Gsqrt(Fm1(:),ke) * ( lnx * qx_P(:) + lny * qy_P(:) )
          Msrc(Fm1) = Msrc(Fm1) - 0.5_RP * (  &
            matmul(mmE(Fm1,Fm1), gtau_ * q_P(:) + GsqrtDn2(:)) )
        end do 
        
        if ( cal_res ) then
          do p=1, elem%Np
            q_(p) = q(p+(ke-1)*elem%Np)
          end do
          q_(:) = Msrc(:) - matmul(OPT11(:,:), q_(:))
          ! res(:,ke) = matmul(elem%invM, q_)
          res(:,ke) = q_
        else
          call LinAlgebra_SolveLinEq( OPT11, Msrc, q_ )
          !pold(:) = p(1+(k1-1)*elem%Np:k1*elem%Np)
          do p=1, elem%Np
            q(p+(ke-1)*elem%Np) = q_(p)
          end do
        end if        
      end do
    end do
    return
  end subroutine cal_q_lc0  

!OCL SERIAL
  subroutine cal_q_lc( q, res,                       &
    color_id, cal_res, is_mg_top_level, rhs, qx, qy, &
    vmapM, vmapP, lmesh, elem                        )
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
    logical, intent(in) :: is_mg_top_level
    integer, intent(in) :: vmapM(elem%Nfp,elem%Nfaces,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%Nfp,elem%Nfaces,lmesh%Ne)

    integer :: ke, ke2
    integer :: i, j

    integer :: p

    integer :: j_int

    integer :: f1, f2

    real(RP) :: E11, E22

    real(RP) :: Dx(elem%Np,elem%Np)
    real(RP) :: Dy(elem%Np,elem%Np)
    real(RP) :: Sx_tr0(elem%Np,elem%Np)
    real(RP) :: Sy_tr0(elem%Np,elem%Np)
    real(RP) :: Sx_tr(elem%Np,elem%Np)
    real(RP) :: Sy_tr(elem%Np,elem%Np)

    real(RP) :: GsqrtDx(elem%Np,elem%Np)
    real(RP) :: GsqrtDy(elem%Np,elem%Np)
    real(RP) :: GsqrtI(elem%Np,elem%Np)

    real(RP) :: Dn1(elem%Np,elem%Np)    
    real(RP) :: GsqrtDn1(elem%Np,elem%Np)

    real(RP) :: OPT11(elem%Np,elem%Np)
    real(RP) :: mmE(elem%Np,elem%Np)   
    real(RP) :: Lift_(elem%Np,elem%Np)
    real(RP) :: Lift_x(elem%Np,elem%Np)       
    real(RP) :: Lift_y(elem%Np,elem%Np)
    real(RP) :: Lift_x_q(elem%Np)       
    real(RP) :: Lift_y_q(elem%Np)
    real(RP) :: Msrc(elem%Np)

    real(RP) :: massEdge(elem%Np,elem%Np,elem%Nfaces)
    real(RP) :: Emat(elem%Np,elem%Nfp*elem%Nfaces)
    integer :: Fm1(elem%Nfp)
    integer :: id
    real(RP) :: lnx, lny  

    real(RP) :: q_(elem%Np)
    real(RP) :: q_P(elem%Nfp)
    real(RP) :: qx_P(elem%Nfp)
    real(RP) :: qy_P(elem%Nfp)
    real(RP) :: GsqrtDn2(elem%Nfp)

    real(RP) :: gtau_

    integer :: parity, parity_max
    integer :: j0
    integer :: ii, ncol
    !---------------------------

    Emat(:,:) = matmul(elem%M, elem%Lift)
    massEdge(:,:,:)  = 0.0_RP
    do f1=1, elem%Nfaces
      Fm1 = elem%Fmask(:,f1)
      massEdge(Fm1,Fm1,f1) = Emat(Fm1, (f1-1)*elem%Nfp+1:f1*elem%Nfp)
    end do

    parity_max = merge(0, 1, cal_res)

    Sx_tr0(:,:) = matmul(transpose(elem%Dx1),elem%M)
    Sy_tr0(:,:) = matmul(transpose(elem%Dx2),elem%M)

   !$omp parallel private( j0, j_int, ncol, i, ke, &
   !$omp E11, E22, Dx, Dy, Sx_tr, Sy_tr, GsqrtI, GsqrtDx, GsqrtDy,           &
   !$omp Dn1, GsqrtDn1, GsqrtDn2, OPT11, mmE, Lift_, Lift_x, Lift_y, Lift_x_q, Lift_y_q,   &
   !$omp Msrc, f1, ke2, f2, id, lnx, lny, Fm1, q_, q_P, qx_P, qy_P, gtau_, p     )

    do parity=0, parity_max

      if (cal_res) then
        j0 = 1; j_int = 1; ncol = lmesh%NeX
      else
        j0 = 1 + modulo(parity-color_id, 2)
        j_int = 2
        ncol = (lmesh%NeX+1-parity)/2
      end if

      !$omp do collapse(2)
      do j=j0, lmesh%NeY, j_int
      do ii=1, ncol
        
        if (cal_res) then
          i = ii
        else
          i = 2*ii - 1 + parity
        end if
        ke = i + (j-1)*lmesh%NeX
        
        ! Assume E11 and E22 are constant in an element
        E11 = lmesh%Escale(1,ke,1,1)
        E22 = lmesh%Escale(1,ke,2,2)

        Dx(:,:) = E11 * elem%Dx1
        Dy(:,:) = E22 * elem%Dx2
        Sx_tr(:,:) = E11 * ( E11 * Sx_tr0(:,:) )
        Sy_tr(:,:) = E22 * ( E22 * Sy_tr0(:,:) )

        do p=1, elem%Np
          GsqrtI(:,p) = lmesh%Gsqrt(:,ke)
          GsqrtDx(:,p) = GsqrtI(:,p) * ( lmesh%GIJ(:,ke,1,1) * Dx(:,p) + lmesh%GIJ(:,ke,1,2) * Dy(:,p) )
          GsqrtDy(:,p) = GsqrtI(:,p) * ( lmesh%GIJ(:,ke,2,1) * Dx(:,p) + lmesh%GIJ(:,ke,2,2) * Dy(:,p) )
        end do

        if ( is_mg_top_level ) then
          Msrc(:) = matmul(elem%M, rhs(:,ke) * lmesh%Gsqrt(:,ke) )
        else
          Msrc(:) = rhs(:,ke)
        end if

        OPT11(:,:) = - ( matmul(Sx_tr, GsqrtDx) + matmul(Sy_tr, GsqrtDy) )

        do f1=1, 4        
          ke2 = lmesh%EToE(ke,f1); f2 = lmesh%EToF(ke,f1)
          
          id = 1 + (f1-1)*elem%Nfp
          lnx = lmesh%normal_fn(id,ke,1)
          lny = lmesh%normal_fn(id,ke,2)
          Fm1(:) = elem%Fmask(:,f1)
          
          Dn1(:,:) = lnx * Dx(:,:) + lny * Dy(:,:)
          GsqrtDn1(:,:) = lnx * GsqrtDx(:,:) + lny * GsqrtDy(:,:)
          
          mmE(:,:) = lmesh%Fscale(id,ke) * massEdge(:,:,f1)  
                          
          Lift_(:,:) = matmul(elem%invM, mmE)
          do p=1, elem%Np
            Lift_x(:,p) = lmesh%Gsqrt(:,ke) * ( lnx * lmesh%GIJ(:,ke,1,1) + lny * lmesh%GIJ(:,ke,1,2) ) * Lift_(:,p)
            Lift_y(:,p) = lmesh%Gsqrt(:,ke) * ( lnx * lmesh%GIJ(:,ke,2,1) + lny * lmesh%GIJ(:,ke,2,2) ) * Lift_(:,p)
          end do
          
          gtau_ = elem%PolyOrder * (elem%PolyOrder + 1 ) * lmesh%Fscale(id,ke) * 0.5_RP * 7.0_RP
          
          OPT11(:,:) = OPT11(:,:) + 0.5_RP * ( &
            - gtau_ * mmE(:,:) + matmul(mmE, GsqrtDn1)         &
            + ( matmul(Sx_tr,Lift_x) + matmul(Sy_tr,Lift_y) ) )
          
          q_P(:) = q(VMapP(:,f1,ke))
          Lift_x_q(:) = matmul(Lift_x(:,Fm1), q_P(:))
          Lift_y_q(:) = matmul(Lift_y(:,Fm1), q_P(:))
          Msrc(:) = Msrc(:) + 0.5_RP * ( &
            + matmul(Sx_tr,Lift_x_q) + matmul(Sy_tr,Lift_y_q) )

          qx_P(:) = qx(VMapP(:,f1,ke))
          qy_P(:) = qy(VMapP(:,f1,ke))
          GsqrtDn2(:) = lmesh%Gsqrt(Fm1(:),ke) * ( lnx * qx_P(:) + lny * qy_P(:) )
          Msrc(Fm1) = Msrc(Fm1) - 0.5_RP * (  &
            matmul(mmE(Fm1,Fm1), gtau_ * q_P(:) + GsqrtDn2(:)) )
        end do 
        
        if ( cal_res ) then
          do p=1, elem%Np
            q_(p) = q(p+(ke-1)*elem%Np)
          end do
          res(:,ke) = Msrc(:) - matmul(OPT11(:,:), q_(:))
        else
          call LinAlgebra_SolveLinEq( OPT11, Msrc, q_ )
          do p=1, elem%Np
            q(p+(ke-1)*elem%Np) = q_(p)
          end do
        end if        
      end do
      end do
    end do
    !$omp end parallel
    return
  end subroutine cal_q_lc  
end module mod_poisson2d_smoother
