#include "scaleFElib.h"
module mod_poisson3d_smoother
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_precision
  use scale_prc
  use scale_io

  use scale_sparsemat, only: &
    SparseMat, SparseMat_matmul
  
  use scale_element_base, only: ElementBase3D
  use scale_element_hexahedral, only: HexahedralElement
  use scale_localmesh_3d, only: LocalMesh3D
  use scale_mesh_base3d, only: MeshBase3D

  use scale_meshfield_base, only: MeshField3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldCommBase, MeshFieldContainer
  use scale_meshfieldcomm_cubedom3d, only: MeshFieldCommCubeDom3D

  use scale_multigrid_smoother_base, only: &
    MGSmoother3DBase, MGSmoother_PRE_ID
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type, public, extends(MGSmoother3DBase) :: MGSmoother_Poisson3D
  contains
    procedure :: Init => MGSmoother_Poisson3D_Init
    procedure :: Final => MGSmoother_Poisson3D_Final
    procedure :: Advance_itr_1step => MGSmoother_Poisson3D_advance_itr_1step
  end type MGSmoother_Poisson3D

  ! public :: poisson_smoother_evaluate_error_norm
  integer, parameter, public :: MGSmoother_Possion3D_AUX_SCALAR_NUM  = 1
  integer, parameter, public :: MGSmoother_Possion3D_AUX_HVEC_NUM    = 1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: Poisson3d_smoother_advance_itr_1step_color
  private :: cal_grad_lc
  private :: cal_q_lc

  integer, parameter :: AUXVAR_QZ_ID = 1
  integer, parameter :: AUXVAR_QX_ID = 2
  integer, parameter :: AUXVAR_QY_ID = 3
  
contains
  !> Initialization
  subroutine MGSmoother_Poisson3D_Init( this, mesh )
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    class(MeshBase3D), intent(in), target :: mesh
    !---------------------------------------------------------------------------
    return
  end subroutine MGSmoother_Poisson3D_Init

  !> Finalization
  subroutine MGSmoother_Poisson3D_Final(this)
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    !---------------------------------------------------------------------------
    return
  end subroutine MGSmoother_Poisson3D_Final

  !> Smoother
!OCL SERIAL
  subroutine MGSmoother_Poisson3D_advance_itr_1step( this, &
    q, res,                                      & 
    f, aux_var,                                  &
    itr, var_comm, aux_comm, Dx, Dy, Dz, mesh3D, &
    cal_res_flag, zero_initial_guess,            &
    mg_p_level, mg_h_level, pre_or_post_smooth )

    use scale_mesh_hierarchy_base, only: &
      pMG_FINEST_LEVEL => MESH_HIERARCHY_pMG_FINEST_LEVEL
    implicit none
    class(MGSmoother_Poisson3D), intent(inout) :: this
    type(MeshField3D), intent(inout), target :: q
    type(MeshField3D), intent(inout) :: res
    type(MeshField3D), intent(in) :: f
    type(MeshField3D), intent(inout), target :: aux_var(:)
    integer, intent(in) :: itr
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: DZ
    class(MeshBase3D), intent(in), target :: mesh3D
    logical, intent(in) :: cal_res_flag
    logical, intent(in) :: zero_initial_guess
    integer, intent(in) :: mg_p_level
    integer, intent(in) :: mg_h_level
    integer, intent(in) :: pre_or_post_smooth

    integer :: ldomID
    class(LocalMesh3D), pointer :: lmesh3D
    integer :: ke

    real(RP) :: res_lcdom_int0
    real(RP) :: res_lcdom_int

    type(MeshFieldContainer) :: var_comm_list(1)
    logical :: is_mg_top_level
    integer :: color_id_offset
    !---------------------------------------------------------------------------

    is_mg_top_level = ( mg_p_level == pMG_FINEST_LEVEL )

    if ( pre_or_post_smooth == MGSmoother_PRE_ID ) then
      color_id_offset = 2
    else
      color_id_offset = 1
    end if

    if ( zero_initial_guess ) then
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)
        !$omp parallel do
        do ke=lmesh3D%NeS, lmesh3D%NeE
          q%local(ldomID)%val(:,ke) = 0.0_RP
        end do
      end do
    end if

    ! if ( cal_res_flag ) then
    !   call Poisson3d_smoother_advance_itr_1step_color( &
    !     q, res,                                                                  & 
    !     f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID),  &
    !     itr, var_comm, aux_comm, Dx, Dy, Dz, mesh3D,                             &
    !     cal_res_flag, 0, is_mg_top_level )

    !   res_lcdom_int0 = 0.0_RP
    !   do ldomID=1, mesh3D%LOCAL_MESH_NUM
    !     lmesh3D => mesh3D%lcmesh_list(ldomID)
    !     !$omp parallel do reduction(+:res_lcdom_int0)
    !     do ke=lmesh3D%NeS, lmesh3D%NeE
    !       res_lcdom_int0 = res_lcdom_int0 + sum(res%local(ldomID)%val(:,ke)**2)
    !     end do
    !   end do
    !   LOG_INFO("Poisson3d_smoother_advance_itr_1step",*) "itr=", itr, ": lc_res0=", sqrt(res_lcdom_int0)
    ! end if

    ! Red
    call Poisson3d_smoother_advance_itr_1step_color( &
      q, res,                                                                  & 
      f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID),  &
      itr, var_comm, aux_comm, Dx, Dy, Dz, mesh3D,                             &
      .false., mod(color_id_offset,2), is_mg_top_level )
    
    ! Red
    call Poisson3d_smoother_advance_itr_1step_color( &
      q, res,                                                                  & 
      f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID),  &
      itr, var_comm, aux_comm, Dx, Dy, Dz, mesh3D,                             &
      .false., mod(color_id_offset+1,2), is_mg_top_level )
     
    if ( cal_res_flag ) then
      call Poisson3d_smoother_advance_itr_1step_color( &
        q, res,                                                                  & 
        f, aux_var(AUXVAR_QX_ID), aux_var(AUXVAR_QY_ID), aux_var(AUXVAR_QZ_ID),  &
        itr, var_comm, aux_comm, Dx, Dy, Dz, mesh3D,                             &
        cal_res_flag, 0, is_mg_top_level )

      res_lcdom_int = 0.0_RP
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        lmesh3D => mesh3D%lcmesh_list(ldomID)
        !$omp parallel do reduction(+:res_lcdom_int)
        do ke=lmesh3D%NeS, lmesh3D%NeE
          res_lcdom_int = res_lcdom_int + sum(res%local(ldomID)%val(:,ke)**2)
        end do
      end do
      ! LOG_INFO("Poisson3d_smoother_advance_itr_1step",*) "itr=", itr, ": lc_res0=", sqrt(res_lcdom_int0), ", lc_res1=", sqrt(res_lcdom_int)
      LOG_INFO("Poisson3d_smoother_advance_itr_1step",*) "itr=", itr, ": lc_res1=", sqrt(res_lcdom_int)
    end if
    
    var_comm_list(1)%field3d => q
    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    return
  end subroutine MGSmoother_Poisson3D_advance_itr_1step

!-- private -----------------------

!OCL SERIAL
  subroutine Poisson3d_smoother_advance_itr_1step_color( &
    q, res,                                      & 
    f, qx, qy, qz,                               &
    itr, var_comm, aux_comm, Dx, Dy, Dz, mesh3D, &
    cal_res_flag, color_id, is_mg_top_level )
    implicit none

    type(MeshField3D), intent(inout), target :: q
    type(MeshField3D), intent(inout) :: res
    type(MeshField3D), intent(in) :: f
    type(MeshField3D), intent(inout), target :: qx
    type(MeshField3D), intent(inout), target :: qy
    type(MeshField3D), intent(inout), target :: qz
    integer, intent(in) :: itr
    class(MeshFieldCommBase), intent(inout) :: var_comm
    class(MeshFieldCommBase), intent(inout) :: aux_comm
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz
    class(MeshBase3D), intent(in), target :: mesh3D
    logical, intent(in) :: cal_res_flag
    integer, intent(in) :: color_id
    logical, intent(in) :: is_mg_top_level

    integer :: ldomID
    class(LocalMesh3D), pointer :: lmesh3D

    type(MeshFieldContainer) :: var_comm_list(1)
    type(MeshFieldContainer) :: aux_comm_list(3)

    integer :: color_id_offset
    !---------------------------------------------------------------------------

    !-
    var_comm_list(1)%field3d => q

    call var_comm%Put( var_comm_list, 1 )
    call var_comm%Exchange()
    call var_comm%Get( var_comm_list, 1 )

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(ldomID)
      call cal_grad_lc( qx%local(ldomID)%val, qy%local(ldomID)%val,  qz%local(ldomID)%val, & ! (out)
        color_id, q%local(ldomID)%val, Dx, Dy, Dz, lmesh3D, lmesh3D%refElem3D ) ! (in)
    end do

    !- 
    aux_comm_list(AUXVAR_QZ_ID)%field3d => qz
    aux_comm_list(AUXVAR_QX_ID)%field3d => qx
    aux_comm_list(AUXVAR_QY_ID)%field3d => qy

    call aux_comm%Put( aux_comm_list, 1 )
    call aux_comm%Exchange()
    call aux_comm%Get( aux_comm_list, 1 )

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(ldomID)
      call cal_q_lc( q%local(ldomID)%val, res%local(ldomID)%val,          &
        color_id, cal_res_flag, is_mg_top_level,f%local(ldomID)%val,      &
        qx%local(ldomID)%val, qy%local(ldomID)%val, qy%local(ldomID)%val, &
        lmesh3D%VMapM, lmesh3D%VMapP, lmesh3D, lmesh3D%refElem3D )
    end do
    
    return
  end subroutine Poisson3d_smoother_advance_itr_1step_color

!- private
  subroutine cal_grad_lc( qx, qy, qz, &
    color_id, q, Dx, Dy, Dz, lmesh, elem )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: qx(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qy(elem%Np,lmesh%NeA)
    real(RP), intent(out) :: qz(elem%Np,lmesh%NeA)
    integer, intent(in) :: color_id
    real(RP), intent(in) :: q(elem%Np,lmesh%NeA)
    type(SparseMat), intent(in) :: Dx
    type(SparseMat), intent(in) :: Dy
    type(SparseMat), intent(in) :: Dz

    integer :: ke
    real(RP) :: q1(elem%Np)
    real(RP) :: q2(elem%Np)
    real(RP) :: q3(elem%Np)
    !----------------------------
    
    !$omp parallel do private(q1, q2, q3)
    do ke=lmesh%NeS, lmesh%NeE
      call Sparsemat_matmul(Dx, q(:,ke), q1)
      call Sparsemat_matmul(Dy, q(:,ke), q2)
      call Sparsemat_matmul(Dz, q(:,ke), q3)

      qx(:,ke) = lmesh%Escale(:,ke,1,1) * q1(:)
      qy(:,ke) = lmesh%Escale(:,ke,2,2) * q2(:)
      qz(:,ke) = lmesh%Escale(:,ke,3,3) * q3(:)
    end do
    return
  end subroutine cal_grad_lc

!OCL SERIAL
  subroutine cal_q_lc( q, res,                          &
    color_id, cal_res, is_mg_top_level, rhs, qx, qy,qz, &
    vmapM, vmapP, lmesh, elem                           )
    use scale_prc
    use scale_linalgebra, only: Linalgebra_SolveLinEq
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: q(elem%Np*lmesh%NeA)
    real(RP), intent(out) :: res(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: rhs(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: qx(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: qy(elem%Np*lmesh%NeA)    
    real(RP), intent(in) :: qz(elem%Np*lmesh%NeA)    
    integer, intent(in) :: color_id
    logical, intent(in) :: cal_res
    logical, intent(in) :: is_mg_top_level
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)

    integer :: ke, ke2
    integer :: i, j, k

    integer :: p

    integer :: f1, f2

    real(RP) :: E11, E22, E33

    real(RP) :: Dx(elem%Np,elem%Np)
    real(RP) :: Dy(elem%Np,elem%Np)
    real(RP) :: Dz(elem%Np,elem%Np)
    real(RP) :: Sx_tr0(elem%Np,elem%Np)
    real(RP) :: Sy_tr0(elem%Np,elem%Np)
    real(RP) :: Sz_tr0(elem%Np,elem%Np)
    real(RP) :: Sx_tr(elem%Np,elem%Np)
    real(RP) :: Sy_tr(elem%Np,elem%Np)
    real(RP) :: Sz_tr(elem%Np,elem%Np)

    real(RP) :: GsqrtDx(elem%Np,elem%Np)
    real(RP) :: GsqrtDy(elem%Np,elem%Np)
    real(RP) :: GsqrtDz(elem%Np,elem%Np)
    real(RP) :: GsqrtI(elem%Np,elem%Np)
  
    real(RP) :: GsqrtDn1(elem%Np,elem%Np)

    real(RP) :: OPT11(elem%Np,elem%Np)
    real(RP) :: mmE(elem%Np,elem%Np)   
    real(RP) :: Lift_(elem%Np,elem%Np)
    real(RP) :: Lift_x(elem%Np,elem%Np)       
    real(RP) :: Lift_y(elem%Np,elem%Np)
    real(RP) :: Lift_z(elem%Np,elem%Np)
    real(RP) :: Lift_x_q(elem%Np)       
    real(RP) :: Lift_y_q(elem%Np)
    real(RP) :: Lift_z_q(elem%Np)
    real(RP) :: Msrc(elem%Np)

    real(RP) :: Emat(elem%Np,elem%NfpTot)
    real(RP) :: massEdge_h(elem%Np,elem%Np,elem%Nfaces_h)
    real(RP) :: massEdge_v(elem%Np,elem%Np,elem%Nfaces_v)
    integer :: Fm1_h(elem%Nfp_h)
    integer :: Fm1_v(elem%Nfp_v)
    integer :: vid_s, vid_e
    real(RP) :: lnx, lny, lnz

    real(RP) :: q_(elem%Np)
    real(RP) :: q_P_h(elem%Nfp_h)
    real(RP) :: q_P_v(elem%Nfp_v)
    real(RP) :: gradQ_P_h(elem%Nfp_h,2)
    real(RP) :: gradQ_P_v(elem%Nfp_v)
    real(RP) :: GsqrtDn2_h(elem%Nfp_h)
    real(RP) :: GsqrtDn2_v(elem%Nfp_v)

    real(RP) :: gtau_

    integer :: parity, parity_max
    integer :: i0, i_int
    integer :: ii, ncol

    integer :: nfp_os

    real(RP) :: q_old
    real(RP), parameter :: omg = 1.0_RP
    !---------------------------


    Emat(:,:) = matmul(elem%M, elem%Lift)

    massEdge_h(:,:,:)  = 0.0_RP
    do f1=1, elem%Nfaces_h
      Fm1_h(:) = elem%Fmask_h(:,f1)
      massEdge_h(Fm1_h,Fm1_h,f1) = Emat(Fm1_h, (f1-1)*elem%Nfp_h+1:f1*elem%Nfp_h)
    end do

    massEdge_v(:,:,:)  = 0.0_RP
    nfp_os = elem%Nfaces_h * elem%Nfp_h
    do f1=1, elem%Nfaces_v
      Fm1_v(:) = elem%Fmask_v(:,f1)
      massEdge_v(Fm1_v,Fm1_v,f1) = Emat(Fm1_v, nfp_os+(f1-1)*elem%Nfp_v+1:nfp_os+f1*elem%Nfp_v)
    end do

    parity_max = merge(0, 1, cal_res)

    Sx_tr0(:,:) = matmul(transpose(elem%Dx1),elem%M)
    Sy_tr0(:,:) = matmul(transpose(elem%Dx2),elem%M)
    Sz_tr0(:,:) = matmul(transpose(elem%Dx3),elem%M)

!    !$omp parallel private( i0, i_int, ncol, i, ke, &
!    !$omp E11, E22, Dx, Dy, Dz, Sx_tr, Sy_tr, Sz_tr, GsqrtI, GsqrtDx, GsqrtDy, GsqrtDz,   &
!    !$omp Dn1, GsqrtDn1, GsqrtDn2_h, GsqrtDn2_v,                                          &
!    !$omp OPT11, mmE, Lift_, Lift_x, Lift_y, Lift_z, Lift_x_q, Lift_y_q,  Lift_z_q,       &
!    !$omp Msrc, f1, ke2, f2, vid_s, vid_e, lnx, lny, lnz, Fm1_h, Fm1_v,                   &
!    !$omp q_, q_P_h, q_P_v, gradQ_P_h, gradQ_P_v, gtau_, p     )
    do parity=0, parity_max

!      !$omp do collapse(2)
      do k=1, lmesh%NeZ
      do j=1, lmesh%NeY

        if ( cal_res ) then
          i0 = 1
          ncol = lmesh%NeX
        else
          i0 = 1 + modulo(parity - modulo(j + k + color_id, 2), 2)  ! 1 or 2
          ncol = (lmesh%NeX - i0)/2 + 1
          if ( ncol < 1 ) cycle
        end if

        do ii=1, ncol
          if (cal_res) then
            i = ii
          else
            i = i0 + (ii-1)*2
          end if
          ke = i + (j-1)*lmesh%NeX + (k-1)*lmesh%NeX*lmesh%NeY
          
          ! Assume E11, E22, and E33 are constant in an element
          E11 = lmesh%Escale(1,ke,1,1)
          E22 = lmesh%Escale(1,ke,2,2)
          E33 = lmesh%Escale(1,ke,3,3)

          Dx(:,:) = E11 * elem%Dx1
          Dy(:,:) = E22 * elem%Dx2
          Dz(:,:) = E33 * elem%Dx3
          Sx_tr(:,:) = E11 * Sx_tr0(:,:)
          Sy_tr(:,:) = E22 * Sy_tr0(:,:)
          Sz_tr(:,:) = E33 * Sz_tr0(:,:)

          do p=1, elem%Np
            GsqrtI(:,p) = lmesh%Gsqrt(:,ke)
            GsqrtDx(:,p) = GsqrtI(:,p) * Dx(:,p)
            GsqrtDy(:,p) = GsqrtI(:,p) * Dy(:,p)
            GsqrtDz(:,p) = GsqrtI(:,p) * Dz(:,p)
          end do

          if ( is_mg_top_level ) then
            Msrc(:) = matmul(elem%M, rhs(:,ke) * lmesh%Gsqrt(:,ke) )
          else
            Msrc(:) = rhs(:,ke)
          end if

          OPT11(:,:) = - ( matmul(Sx_tr, GsqrtDx) + matmul(Sy_tr, GsqrtDy) + matmul(Sz_tr, GsqrtDz) )

          do f1=1, elem%Nfaces_h       
            ke2 = lmesh%EToE(ke,f1); f2 = lmesh%EToF(ke,f1)
            
            vid_s = 1 + (f1-1)*elem%Nfp_h
            vid_e = vid_s + elem%Nfp_h - 1

            lnx = lmesh%normal_fn(vid_s,ke,1)
            lny = lmesh%normal_fn(vid_s,ke,2)
            Fm1_h(:) = elem%Fmask_h(:,f1)
                        
            mmE(:,:) = lmesh%Fscale(vid_s,ke) * massEdge_h(:,:,f1)  
                            
            Lift_(:,:) = matmul(elem%invM, mmE)
            do p=1, elem%Np
              Lift_x(:,p) = lmesh%Gsqrt(:,ke) * lnx * Lift_(:,p)
              Lift_y(:,p) = lmesh%Gsqrt(:,ke) * lny * Lift_(:,p)
            end do
            
            gtau_ = elem%PolyOrder_h * ( elem%PolyOrder_h + 1 ) * lmesh%Fscale(vid_s,ke) * 0.5_RP * 1.2_RP
            
            GsqrtDn1(:,:) = lnx * GsqrtDx(:,:) + lny * GsqrtDy(:,:)

            OPT11(:,:) = OPT11(:,:) + 0.5_RP * ( &
              - gtau_ * mmE(:,:) + matmul(mmE, GsqrtDn1)         &
              + ( matmul(Sx_tr,Lift_x) + matmul(Sy_tr,Lift_y) ) )
            
            q_P_h(:) = q(VMapP(vid_s:vid_e,ke))
            Lift_x_q(:) = matmul(Lift_x(:,Fm1_h), q_P_h(:))
            Lift_y_q(:) = matmul(Lift_y(:,Fm1_h), q_P_h(:))
            Msrc(:) = Msrc(:) + 0.5_RP * ( &
              + matmul(Sx_tr,Lift_x_q) + matmul(Sy_tr,Lift_y_q) )

            gradQ_P_h(:,1) = qx(VMapP(vid_s:vid_e,ke))
            gradQ_P_h(:,2) = qy(VMapP(vid_s:vid_e,ke))
            GsqrtDn2_h(:) = lmesh%Gsqrt(Fm1_h(:),ke) * ( lnx * gradQ_P_h(:,1) + lny * gradQ_P_h(:,2) )
            Msrc(Fm1_h) = Msrc(Fm1_h) - 0.5_RP * (  &
              matmul(mmE(Fm1_h,Fm1_h), gtau_ * q_P_h(:) + GsqrtDn2_h(:)) )
          end do 

          do f1=1, elem%Nfaces_v       
            ke2 = lmesh%EToE(ke,elem%Nfaces_h+f1); f2 = lmesh%EToF(ke,elem%Nfaces_h+f1)
            
            vid_s = 1 + elem%Nfaces_h*elem%Nfp_h + (f1-1)*elem%Nfp_v
            vid_e = vid_s + elem%Nfp_v - 1

            lnz = lmesh%normal_fn(vid_s,ke,3)
            Fm1_v(:) = elem%Fmask_v(:,f1)
            
            GsqrtDn1(:,:) = lnz * GsqrtDz(:,:)
            
            mmE(:,:) = lmesh%Fscale(vid_s,ke) * massEdge_v(:,:,f1)  
                            
            Lift_(:,:) = matmul(elem%invM, mmE)
            do p=1, elem%Np
              Lift_z(:,p) = lmesh%Gsqrt(:,ke) * Lift_(:,p)
            end do
            
            gtau_ = elem%PolyOrder_v * (elem%PolyOrder_v + 1 ) * lmesh%Fscale(vid_s,ke) * 0.5_RP * 1.2_RP
            
            OPT11(:,:) = OPT11(:,:) + 0.5_RP * ( &
              - gtau_ * mmE(:,:) + matmul(mmE, GsqrtDn1) &
              + matmul(Sz_tr,Lift_z)                     )
            
            q_P_v(:) = q(VMapP(vid_s:vid_e,ke))
            Lift_z_q(:) = matmul(Lift_z(:,Fm1_v), q_P_v(:))
            Msrc(:) = Msrc(:) + 0.5_RP * ( &
              + matmul(Sz_tr,Lift_z_q)     )

            gradQ_P_v(:) = qz(VMapP(vid_s:vid_e,ke))
            GsqrtDn2_v(:) = lmesh%Gsqrt(Fm1_v(:),ke) * ( lnz * gradQ_P_v(:) )
            Msrc(Fm1_v) = Msrc(Fm1_v) - 0.5_RP * (  &
              matmul(mmE(Fm1_v,Fm1_v), gtau_ * q_P_v(:) + GsqrtDn2_v(:)) )
          end do 
          
          if ( cal_res ) then
            do p=1, elem%Np
              q_(p) = q(p+(ke-1)*elem%Np)
            end do
            res(:,ke) = Msrc(:) - matmul(OPT11(:,:), q_(:))
          else
            call LinAlgebra_SolveLinEq( OPT11, Msrc, q_ )
            do p=1, elem%Np
              q_old = q(p+(ke-1)*elem%Np)
              q(p+(ke-1)*elem%Np) = q_old + omg * ( q_(p) - q_old )
            end do
          end if        
        end do
      end do
      end do
    end do
!    !$omp end parallel
    return
  end subroutine cal_q_lc  
end module mod_poisson3d_smoother

