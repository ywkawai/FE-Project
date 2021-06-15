program test_eigen_analysis  
  use scale_precision
  use scale_const, only: &
    PI => CONST_PI
  use scale_io
  implicit none

  integer, parameter :: pmax = 7
  integer :: p

  real(RP), parameter :: ADV_VEL      = 1.0_RP
  real(RP), parameter :: BETA_UPWIND  = 1.0_RP
  real(RP), parameter :: BETA_CENTRAL = 0.0_RP
  real(RP), parameter :: Helem        = 1.0_RP

  integer, parameter :: nout = 10
  !----------------------------------

  do p=1, pmax
    call eigen_analysis( p, BETA_UPWIND, 'modal', 'weak' )
    call eigen_analysis( p, BETA_UPWIND, 'nodal', 'weak' )
  end do

contains

  subroutine eigen_analysis( porder, beta, basis_type, form_type )
    use scale_element_line, only: LineElement
    use scale_polynominal, only: &
      polynominal_genLegendrePoly, &
      polynominal_genDLegendrePoly
    implicit none
    integer, intent(in) :: porder
    real(RP), intent(in) :: beta
    character(len=*), intent(in) :: basis_type
    character(len=*), intent(in) :: form_type

    type(LineElement) :: elem
    real(RP) :: KC(porder+1,porder+1)
    real(RP) :: KM(porder+1,porder+1)
    real(RP) :: KP(porder+1,porder+1)
    real(RP) :: Stiff(porder+1,porder+1)
    real(RP) :: Minv(porder+1,porder+1)
    real(RP) :: phiM1(porder+1)
    real(RP) :: phiP1(porder+1)
    real(RP) :: P1D_ (porder+1,porder+1)
    real(RP) :: DP1D_(porder+1,porder+1)

    integer :: p1, p2

    integer :: k
    integer, parameter :: kmax = 400
    real(RP) :: dwnum
    real(RP) :: wnum
    complex(RP) :: Mat(porder+1,porder+1)
    real(RP) :: cs_kh, si_kh
    complex(RP), parameter :: ei = (0.0_RP, 1.0_RP)

    complex(RP) :: dummy(1,1)
    complex(RP) :: vr(porder+1,porder+1), w(porder+1), rwork(2*(porder+1))
    real(RP) :: w_r(porder+1), w_i(porder+1), dw_r(porder+1), dw_i(porder+1)
    complex(RP), allocatable :: work(:)
    integer :: lwork
    integer :: info

    real(RP) :: wnum_nondim(kmax)
    real(RP) :: wnum_nondim_num(2,kmax,porder+1)
    character(len=H_MID) :: filename
    !------------------------------------

    write(*,'(a,i2,a,a)') "Eigen analysis: porder=", porder, ', basis_type: ', trim(basis_type)

    call elem%Init( porder, .false. )

    if (basis_type == 'nodal' ) then
      Minv(:,:) = elem%invM(:,:)
      Stiff(:,:) = matmul(elem%M, elem%Sx1(:,:))
      phiM1(:) = 0.0_RP; phiP1(:) = 0.0_RP
      phiM1(1) = 1.0_RP; phiP1(porder+1) = 1.0_RP
    else if (basis_type == 'modal' ) then
      P1D_ (:,:) = Polynominal_GenLegendrePoly( elem%PolyOrder, elem%x1(:) )
      DP1D_(:,:) = polynominal_genDLegendrePoly( elem%PolyOrder, elem%x1(:), P1D_(:,:) )
      ! normalization
      do p1=1, porder+1
        P1D_ (:,p1) =  P1D_(:,p1) * sqrt(dble(p1-1) + 0.5_RP)
        DP1D_(:,p1) = DP1D_(:,p1) * sqrt(dble(p1-1) + 0.5_RP)
      end do
      Minv(:,:) = 0.0_RP
      do p2=1, porder+1
        Minv(p2,p2) = 1.0_RP
        phiM1(p2) = P1D_(1,p2)        
        phiP1(p2) = P1D_(porder+1,p2)
        do p1=1, porder+1
          Stiff(p1,p2) = sum( elem%IntWeight_lgl(:) * P1D_(:,p2) * DP1D_(:,p1) )
        end do
      end do
    else
      write(*,*) 'Check! Invalid basis_type:', basis_type
      stop
    end if

    !--
    call cosntruct_mat( KC, KP, KM,           & ! (out)
      Stiff, Minv, phiM1, phiP1, porder, beta ) ! (in)

    !-- Calculate eigenvalues

    ! Use routine workspace query to get optimal workspace.    
    lwork = -1
    call zgeev( 'No left vectors', 'Vectors (right)', porder+1, Mat, porder+1, w, dummy, 1, &
      vr, porder+1, dummy, lwork, rwork, info )
    lwork = max( (64+1)*(porder+1), nint(real(dummy(1,1))) )
    allocate( work(lwork) )
  
    dwnum = PI * dble(porder + 1) / Helem / dble( kmax )
    dw_r(:) = 0.0_RP; dw_i(:) = 0.0_RP
    do k=0, kmax-1
      wnum = dble(k) * dwnum
      cs_kh = cos(wnum * Helem)
      si_kh = sin(wnum * Helem)
      do p2=1, porder+1
      do p1=1, porder+1
        Mat(p1,p2) = 2.0_RP * ( ( ( KM(p1,p2) + KP(p1,p2) ) * cs_kh + KC(p1,p2) ) * ei - ( - KM(p1,p2) + KP(p1,p2) ) * si_kh )
      end do
      end do

      call zgeev( 'No left vectors', 'Vectors (right)', porder+1, Mat, porder+1, w, dummy, 1, &
      vr, porder+1, work, lwork, rwork, info )
    
      if ( info /= 0 ) then
        write(*,*) 'Failure in ZGEEV.  INFO = ', info
        stop
      end if
      
      w_r(:) =  real(w(:)) / dble(porder + 1)
      w_i(:) = aimag(w(:)) / dble(porder + 1)
      if ( k > 0 ) then
        call sort_eigenval( w_r, w_i, dw_r, dw_i,                  &
          wnum_nondim_num(1,k,:), wnum_nondim_num(2,k,:), porder+1 )
      end if
      wnum_nondim(k+1) = wnum * Helem / dble(porder + 1) 
      wnum_nondim_num(1,k+1,:) = w_r(:)
      wnum_nondim_num(2,k+1,:) = w_i(:)

      ! write(*,*) 'K = ', wnum_nondim(k+1)
      ! write(*,*) 'Eigenvalues: Re', wnum_nondim_num(1,k+1,:)
      ! write(*,*) 'Eigenvalues: Im', wnum_nondim_num(2,k+1,:)
    end do
    
    !--
    do p1=1, porder+1
      write(filename,'(a,a,a,a,a,i2.2,a,i2.2,a)') './data/eigenval_', trim(basis_type), '_', trim(form_type), '_P', porder, '_Mode', p1, '.dat'
      open( nout, file=trim(filename), status='replace' )
      write( nout, '(a)') "# K Re(Km) Im(Km)"
      do k=1, kmax
        write( nout, '(3f15.8)') wnum_nondim(k), wnum_nondim_num(1,k,p1), wnum_nondim_num(2,k,p1)
      end do
      close( nout )
    end do

    call elem%Final()  

    return
  end subroutine eigen_analysis

  subroutine cosntruct_mat( K, Kp, Km,        &
      Stiff, Minv, phiM1, phiP1, porder, beta )
    implicit none
    integer, intent(in) :: porder
    real(RP), intent(out) :: K (porder+1,porder+1)
    real(RP), intent(out) :: Kp(porder+1,porder+1)
    real(RP), intent(out) :: Km(porder+1,porder+1)
    real(RP), intent(in) :: Stiff(porder+1,porder+1)
    real(RP), intent(in) :: Minv(porder+1,porder+1)
    real(RP), intent(in) :: phiM1(porder+1)
    real(RP), intent(in) :: phiP1(porder+1)
    real(RP), intent(in) :: beta

    real(RP) :: bM, bP
    real(RP) :: tmpMat (porder+1,porder+1)
    real(RP) :: tmpMatM(porder+1,porder+1)
    real(RP) :: tmpMatP(porder+1,porder+1)

    integer :: l, m
    !------------------------------------------

    bM = 0.5_RP * ( 1.0_RP - sign(beta, ADV_VEL) )
    bP = 0.5_RP * ( 1.0_RP + sign(beta, ADV_VEL) )

    do m=1, porder+1
    do l=1, porder+1
      tmpMat (l,m) = Stiff(l,m) + bM * phiM1(m) * phiM1(l) - bP * phiP1(m) * phiP1(l)
      tmpMatM(l,m) =   bP * phiP1(m) * phiM1(l) 
      tmpMatP(l,m) = - bM * phiM1(m) * phiP1(l) 
    end do
    end do

    K (:,:) = matmul(Minv, tmpMat )
    Km(:,:) = matmul(Minv, tmpMatM)
    Kp(:,:) = matmul(Minv, tmpMatP)

    return
  end subroutine cosntruct_mat

  subroutine sort_eigenval( w_r, w_i, dw_r, dw_i, &
      wr_prev, wi_prev, n )
    use scale_quicksort, only: QUICKSORT_exec_with_idx
    implicit none
    integer, intent(in) :: n
    real(RP), intent(inout) :: w_r(n)
    real(RP), intent(inout) :: w_i(n)
    real(RP), intent(inout) :: dw_r(n)
    real(RP), intent(inout) :: dw_i(n)
    real(RP), intent(in) :: wr_prev(n)
    real(RP), intent(in) :: wi_prev(n)

    integer :: p1, p2, nn
    real(RP) :: wr_ori(n)    
    real(RP) :: wi_ori(n)
    real(RP) :: dw_ori(2,n)
    integer :: loc_min
    real(RP) :: diff(n)
    logical :: flag(n)
    integer :: sorted_ind(n)
    !----------------------
    
    wr_ori(:) = w_r(:)
    wi_ori(:) = w_i(:)
    dw_ori(1,:) = dw_r(:)
    dw_ori(2,:) = dw_i(:)
    flag(:) = .false.

    do p1=1, n
      forall(nn=1:n) sorted_ind(nn) = nn
      diff(:) = ( (wr_ori(:) - wr_prev(p1)) - dw_ori(1,p1) )**2 &
              + ( (wi_ori(:) - wi_prev(p1)) - dw_ori(2,p1) )**2
      call QUICKSORT_exec_with_idx(n, diff(:), sorted_ind(:))
      do p2=1, n
        loc_min = sorted_ind(p2)
        if ( .not. flag(loc_min) ) then
          flag(loc_min) = .true.
          w_r(p1) = wr_ori(loc_min)
          w_i(p1) = wi_ori(loc_min)
          dw_r(p1) = w_r(p1) - wr_prev(p1)
          dw_i(p1) = w_i(p1) - wi_prev(p1)
          exit
        end if
      end do
    end do

    return
  end subroutine sort_eigenval

end program test_eigen_analysis