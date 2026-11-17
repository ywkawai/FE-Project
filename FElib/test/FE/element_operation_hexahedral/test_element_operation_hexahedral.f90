#include "scaleFElib.h"
program test_element_operation_hexahedral
  use scale_precision
  use scale_prc
  use scale_io  
  use scale

  use scale_sparsemat
  use scale_element_base
  use scale_element_hexahedral

  use scale_element_operation_base
  use scale_element_operation_general
  use scale_element_operation_tensorprod3D
  use scale_unittest_util, only: UNITTEST_init, UNITTEST_final
  use scale_unittest_assert, only: UNITTEST_assert_true

  implicit none

  integer :: p
  integer, parameter :: Pmax=7
  !-------------------------------------------------

  call init()
  do p=1, Pmax
    call check_operation( p, .false. )
  end do
  call final()

contains
  subroutine check_operation( porder, LumpedMassMatFlag )
    use scale_element_operation_tensorprod3D, only: &
      ElementOperationTensorprod3D_create
    implicit none
    integer, intent(in) :: porder
    logical, intent(in) :: LumpedMassMatFlag

    type(HexahedralElement) :: refElem
    type(ElementOperationGeneral) :: elem_oper_general
    class(ElementOperationTensorProd3D), allocatable :: elem_oper_tensorprod

    type(SparseMat) :: Dx, Dy, Dz, Lift
    character(len=H_SHORT) :: SpMV_StorageFormat = 'ELL'
    !-----------------------------------------------

    LOG_INFO("check_operation",*)  'Initialize HexahedralElement .. porder=', porder

    call refElem%Init( porder, porder, LumpedMassMatFlag )
    
    !---------------------

    LOG_INFO("check_operation",*)  'Elementwise operations using SpMV .. porder=', porder

    call Dx%Init( refElem%Dx1, storage_format=SpMV_StorageFormat )
    call Dy%Init( refElem%Dx2, storage_format=SpMV_StorageFormat )
    call Dz%Init( refElem%Dx3, storage_format=SpMV_StorageFormat )
    call Lift%Init( refElem%Lift, storage_format=SpMV_StorageFormat )

    call elem_oper_general%Init( refElem, Dx, Dy, Dz, Lift )
    call check_operation_p( elem_oper_general, refElem, 'general' )

    LOG_INFO("check_operation",*)  'Elementwise operations using low computational cost kernels assuming tensor product elements .. porder=', porder

    call ElementOperationTensorprod3D_create( refElem, elem_oper_tensorprod )
    call check_operation_p( elem_oper_tensorprod, refElem, 'tensorprod3D' )

    !---------------------
    call elem_oper_general%Final()
    call elem_oper_tensorprod%Final()

    call refElem%Final()
    return
  end subroutine check_operation

  subroutine check_operation_p( elem_oper, elem, oper_name )
    implicit none
    class(ElementOperationBase3D), intent(in) :: elem_oper
    class(ElementBase3D), intent(in) :: elem
    character(len=*), intent(in) :: oper_name

    character(len=H_MID) :: lbl

    real(RP) :: dat_in(elem%Np)
    real(RP) :: dat_in_vec(elem%Np,3)

    real(RP) :: dat_out_grad(elem%Np,4)
    real(RP) :: dat_out_grad_ans(elem%Np,3)

    real(RP) :: dat_in_f(elem%NfpTot)
    real(RP) :: dat_out_lift(elem%Np)
    real(RP) :: dat_out_lift_ans(elem%Np)

    real(RP) :: dat_out_div(elem%Np)
    real(RP) :: dat_out_div_ans(elem%Np)

    real(RP) :: Gsqrt(elem%Np)
    real(RP) :: Escale(3,elem%Np)

    integer :: f, fp, fps
    !----------------------------------------------------------

    call gen_dat( elem, dat_in, 1.0_RP )

    call gen_grad( elem, dat_out_grad_ans )

    do f=1, 4
      fps = (f-1)*elem%Nfp_h
      do fp=1, elem%Nfp_h
        dat_in_f(fps+fp) = dat_in(elem%Fmask_h(fp,f))
      end do
    end do
    do f=1, 2
      fps = 4*elem%Nfp_h + (f-1)*elem%Nfp_v
      do fp=1, elem%Nfp_v
        dat_in_f(fps+fp) = dat_in(elem%Fmask_v(fp,f))
      end do
    end do
    dat_out_lift_ans(:) = matmul(elem%Lift, dat_in_f(:))

    Escale(1,:) = 1.0_RP; Escale(2,:) = 2.0_RP; Escale(3,:) = 0.2_RP; 
    Gsqrt(:) = 100.0_RP - elem%x1(:)**2

    call gen_dat( elem, dat_in_vec(:,1), 1.0_RP )
    call gen_dat( elem, dat_in_vec(:,2), 2.0_RP )
    call gen_dat( elem, dat_in_vec(:,3), 3.0_RP )
    call gen_div( elem, dat_out_lift_ans, Escale, Gsqrt, dat_out_div_ans )

    write(lbl,'(a,a,i0,a)') trim(oper_name), ' (porder=', elem%PolyOrder_h, '): '

    !--
    call elem_oper%Dx( dat_in, dat_out_grad(:,1) )
    call assert( dat_out_grad(:,1), dat_out_grad_ans(:,1), trim(lbl)//' Dx' )
    
    call elem_oper%Dy( dat_in, dat_out_grad(:,2) )
    call assert( dat_out_grad(:,2), dat_out_grad_ans(:,2), trim(lbl)//' Dy' )

    call elem_oper%Dz( dat_in, dat_out_grad(:,3) )
    call assert( dat_out_grad(:,3), dat_out_grad_ans(:,3), trim(lbl)//' Dz' )

    !--
    call elem_oper%Lift( dat_in_f, dat_out_lift )
    call assert( dat_out_lift, dat_out_lift_ans, trim(lbl)//' Lift' )

    !--
    call elem_oper%Div( dat_in_vec, dat_in_f, &
      dat_out_grad )
    dat_out_div(:) = ( Escale(1,:) * dat_out_grad(:,1) &
                     + Escale(2,:) * dat_out_grad(:,2) &
                     + Escale(3,:) * dat_out_grad(:,3) &
                     + dat_out_grad(:,4) ) / Gsqrt(:)
    call assert( dat_out_div, dat_out_div_ans(:), trim(lbl)//' Div' )

    return
  end subroutine check_operation_p

  subroutine gen_dat( elem, dat, fac )
    implicit none
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: dat(elem%Np)
    real(RP), intent(in), optional :: fac
    integer :: pord
    !----------------------------------------------------------

    pord = elem%PolyOrder_h
    dat(:) = ( 4.0_RP * elem%x1(:)**pord &
             + 3.0_RP * elem%x2(:)**pord &
             + 2.0_RP * elem%x3(:)**pord ) * fac
    return
  end subroutine gen_dat

  subroutine gen_grad( elem, dat )
    implicit none
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(out) :: dat(elem%Np,3)

    integer :: pord
    !----------------------------------------------------------
    pord = elem%PolyOrder_h

    dat(:,1) = 4.0_RP * elem%x1(:)**(pord-1) * real(pord,kind=RP)
    dat(:,2) = 3.0_RP * elem%x2(:)**(pord-1) * real(pord,kind=RP)
    dat(:,3) = 2.0_RP * elem%x3(:)**(pord-1) * real(pord,kind=RP)
    return
  end subroutine gen_grad

  subroutine gen_div( elem, dat_lift, Escale, Gsqrt, dat )
    implicit none
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(in) :: dat_lift(elem%Np)
    real(RP), intent(in) :: Escale(3,elem%Np)
    real(RP), intent(in) :: Gsqrt(elem%Np)
    real(RP), intent(out) :: dat(elem%Np)

    integer :: pord
    !----------------------------------------------------------
    pord = elem%PolyOrder_h

    dat(:) = ( &
        Escale(1,:) * 4.0_RP * elem%x1(:)**(pord-1) * real(pord,kind=RP) * 1.0_RP &
      + Escale(2,:) * 3.0_RP * elem%x2(:)**(pord-1) * real(pord,kind=RP) * 2.0_RP &
      + Escale(3,:) * 2.0_RP * elem%x3(:)**(pord-1) * real(pord,kind=RP) * 3.0_RP &
      + dat_lift(:) ) / Gsqrt(:)
    
    return
  end subroutine gen_div

  subroutine init()
    implicit none
    !----------------------------------------------

    call UNITTEST_init( "test_element_operation_hexahedral" )

    return
  end subroutine init

  subroutine final()
    implicit none
    !----------------------------------------------

    call UNITTEST_final()
    return
  end subroutine final

  !--------------------------

  !> Check that the squared L2 norm of (vals - ans) is less than or equal to EPS
  subroutine assert( vals, ans, label )
    implicit none
    real(RP), intent(in) :: vals(:)
    real(RP), intent(in) :: ans(size(vals))
    character(*), intent(in) :: label

    real(RP), parameter :: EPS = 1.0E-15_RP
    real(RP) :: l2_error
    character(len=H_MID) :: msg
    !--------------------------------------

    l2_error = sum( (vals(:) - ans(:))**2 )
    write(msg,'(a,es12.4)') 'The error is too large: l2error=', sqrt(l2_error)
    call UNITTEST_assert_true( l2_error <= EPS, label, msg )

    return
  end subroutine assert

end program test_element_operation_hexahedral
