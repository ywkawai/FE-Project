#include "scaleFElib.h"
module mod_dg_eigen_analysis
  use scale_precision
  use scale_const, only: &
    PI => CONST_PI
  use scale_io

  use scale_polynominal, only: &
    Polynominal_GenLagrangePoly,  &
    polynominal_genDLegendrePoly, &
    Polynominal_GenLegendrePoly
  use scale_element_line, only: LineElement

  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: eigen_analysis

  !-----------------------------------------------------------------------------
  !++ Public parameters & variables
  !
  complex(RP), parameter :: ei = (0.0_RP, 1.0_RP)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: cosntruct_mat
  private :: construct_full_discrete_mat
  private :: calc_muhat
  private :: sort_eigenval

contains

  subroutine eigen_analysis( &
    output_dir, porder, Helem, beta, ADV_VEL, basis_type, form_type, &
    MF_alph, MF_order, MF_alph_lbl,                                  &
    tlev_num, tlev_slot,                                             &
    tscheme_name, courant_num                                        )

    use scale_element_modalfilter, only: &
      ModalFilter
    implicit none
    character(len=*), intent(in) :: output_dir
    integer, intent(in) :: porder
    real(RP), intent(in) :: Helem
    real(RP), intent(in) :: beta
    real(RP), intent(in) :: ADV_VEL
    character(len=*), intent(in) :: basis_type
    character(len=*), intent(in) :: form_type  ! weak or strong
    real(RP), intent(in) :: MF_alph
    integer, intent(in) :: MF_order
    character(len=*), intent(in) :: MF_alph_lbl
    integer, intent(in) :: tlev_num
    integer, intent(in) :: tlev_slot(tlev_num)
    character(len=*), intent(in), optional :: tscheme_name
    real(RP), intent(in), optional :: courant_num


    type(LineElement) :: elem
    real(RP) :: KC(porder+1,porder+1)
    real(RP) :: KM(porder+1,porder+1)
    real(RP) :: KP(porder+1,porder+1)
    real(RP) :: Stiff(porder+1,porder+1)
    real(RP) :: Mass(porder+1,porder+1)
    real(RP) :: Minv(porder+1,porder+1)
    real(RP) :: phiM1(porder+1)
    real(RP) :: phiP1(porder+1)
    real(RP) :: P1D_ (porder+1,porder+1)
    real(RP) :: DP1D_(porder+1,porder+1)
    type(ModalFilter) :: MFilter

    integer :: p1, p2

    integer :: k
    integer, parameter :: kmax = 2000
    real(RP) :: dwnum
    real(RP) :: wnum
    complex(RP) :: Mat(porder+1,porder+1)    
    real(RP) :: cs_kh, si_kh

    complex(RP) :: dummy(1,1)
    complex(RP) :: vr(porder+1,porder+1)
    complex(RP) :: lamb(porder+1)
    complex(RP) :: vl(porder+1,porder+1)
    real(DP) :: rwork(2*(porder+1))
    real(DP) :: rconde(porder+1), rcondv(porder+1), scale(porder+1)
    complex(RP), allocatable :: work(:)
    real(RP) :: abnrm
    integer :: lwork
    integer :: ihi, ilo, info    
    integer :: ipiv(porder+1)

    real(RP) :: w_r(porder+1), w_i(porder+1), dw_r(porder+1), dw_i(porder+1)
    integer :: ipiv_eigenval(porder+1)

    complex(RP) :: mu_hat(porder+1)
    complex(RP) :: theta(porder+1)
    complex(RP) :: EvecMat(porder+1,porder+1)

    real(RP) :: wnum_nondim(kmax)
    real(RP) :: wnum_nondim_num(2,kmax,porder+1)
    real(RP) :: eigensol_weight(kmax,porder+1)    
    real(RP) :: AmplifiFac_comb(tlev_num,kmax)
    real(RP) :: PhaseError_comb(tlev_num,kmax)
    character(len=H_MID) :: filename

    complex(RP) :: Ue    (porder+1,tlev_num)
    complex(RP) :: Ue_num(porder+1,tlev_num)
    complex(RP) :: tmpvec1(porder+1)
    complex(RP) :: tmpvec2(porder+1)
    real(RP) :: gex, gnum
    real(RP) :: fac
    real(RP) :: phase(porder+1)
    integer :: n
    complex(RP) :: phase_error_tmp

    integer :: nout
    !------------------------------------

    LOG_INFO("eigen_analysis",*) "Eigen analysis: porder=", porder, ', basis_type: ', trim(basis_type)

    if ( present(tscheme_name) ) then
      LOG_INFO("eigen_analysis",*) " full discrete: tscheme:" // trim(tscheme_name) // ', courant_num=', courant_num
    end if

    call elem%Init( porder, .false. )

    if ( basis_type == 'nodal' ) then
      Mass(:,:) = elem%M(:,:)
      Minv(:,:) = elem%invM(:,:)
      Stiff(:,:) = matmul(elem%M, elem%Sx1(:,:))
      phiM1(:) = 0.0_RP; phiP1(:) = 0.0_RP
      phiM1(1) = 1.0_RP; phiP1(porder+1) = 1.0_RP
    else if ( basis_type == 'modal' ) then
      P1D_ (:,:) = Polynominal_GenLegendrePoly( elem%PolyOrder, elem%x1(:) )
      DP1D_(:,:) = polynominal_genDLegendrePoly( elem%PolyOrder, elem%x1(:), P1D_(:,:) )
      ! normalization
      do p1=1, porder+1
        P1D_ (:,p1) =  P1D_(:,p1) * sqrt(dble(p1-1) + 0.5_RP)
        DP1D_(:,p1) = DP1D_(:,p1) * sqrt(dble(p1-1) + 0.5_RP)
      end do
      Mass(:,:) = 0.0_RP
      do p2=1, porder+1
        Mass(p2,p2) = 1.0_RP
        phiM1(p2) = P1D_(1,p2)        
        phiP1(p2) = P1D_(porder+1,p2)
        do p1=1, porder+1
          Stiff(p1,p2) = sum( elem%IntWeight_lgl(:) * P1D_(:,p2) * DP1D_(:,p1) )
        end do
      end do
      Minv(:,:) = Mass(:,:)
    else
      LOG_ERROR("eigen_analysis",*) 'Check! Invalid basis_type:', basis_type
      stop
    end if

    !--
    call cosntruct_mat( KC, KP, KM,                     & ! (out)
      Stiff, Minv, phiM1, phiP1, porder, beta, ADV_VEL, & ! (in) 
      form_type ) ! (in)

    if ( basis_type == 'nodal' ) then
      call MFilter%Init( elem, 0.0_RP, MF_alph, MF_order )
    else if ( basis_type == 'modal' ) then
      call MFilter%Init( elem, 0.0_RP, MF_alph, MF_order )

      MFilter%FilterMat(:,:) = 0.0_RP
      do p1=1, porder+1
        MFilter%FilterMat(p1,p1) = exp( - MF_alph * ( ((dble(p1-1)/dble(porder) - 0.0_RP)/(1.0_RP - 0.0_RP))**MF_order ) )
        !write(*,*) MFilter%FilterMat(p1,:)
      end do
    end if

    !-- Calculate eigenvalues

    ! Use routine workspace query to get optimal workspace.    
    lwork = -1
    call zgeevx( 'Balance', 'Vectors (left)', 'Vectors (right)', 'Both reciprocal condition numbers', &
      porder+1, Mat, porder+1, lamb, vl, porder+1, &
      vr, porder+1, ilo, ihi, scale, abnrm, rconde, rcondv, dummy, lwork, rwork, info )

    lwork = max( (128+1)*(porder+1), nint(real(vl(1,1))) )
    allocate( work(lwork) )

    dwnum = PI * dble(porder + 1) / Helem / dble( kmax-1 )
    dw_r(:) = 0.0_RP; dw_i(:) = 0.0_RP
    do k=0, kmax-1
      wnum = dble(k) * dwnum
      cs_kh = cos( wnum * Helem )
      si_kh = sin( wnum * Helem )

      call calc_muhat( mu_hat(:),                 & ! (out)
        wnum, Helem, elem, basis_type, dble(Minv) ) ! (in)

      do n=1, tlev_num
        do p1=1, porder+1
          phase(1) = mod( dble(tlev_slot(n)) * courant_num * wnum * Helem, 2.0_RP * PI )
          Ue(p1,n) = mu_hat(p1) * ( cos(phase(1)) - ei * sin(phase(1)) ) ! * exp(ei * k * 0.0_RP) .. xe=0
        end do
      end do

      do p2=1, porder+1
      do p1=1, porder+1
        Mat(p1,p2) = 2.0_RP * ( ( ( KM(p1,p2) + KP(p1,p2) ) * cs_kh + KC(p1,p2) ) + ( - KM(p1,p2) + KP(p1,p2) ) * si_kh * ei )
      end do
      end do
      if ( present(tscheme_name) ) then
        call construct_full_discrete_mat( Mat, & ! (inout)
          tscheme_name, courant_num, porder    )  ! (in)        
      end if

      if ( MF_alph > 0.0_RP ) Mat(:,:) = matmul(MFilter%FilterMat(:,:), Mat(:,:))

      call zgeevx( 'Balance', 'Vectors (left)', 'Vectors (right)', 'Both reciprocal condition numbers', &
        porder+1, Mat, porder+1, lamb, vl, porder+1, vr, porder+1, &
        ilo, ihi, scale, abnrm, rconde, rcondv, work, lwork, rwork, info )

      if ( info /= 0 ) then
        LOG_ERROR("eigen_analysis",*) 'Failure in ZGEEVX.  INFO = ', info
        stop
      end if

      theta(:) = mu_hat(:)
      EvecMat(:,:) = vr(:,:)
      call zgesv( porder+1, 1, vr, porder+1, ipiv, theta, porder+1, info )
      
      if ( present(tscheme_name) ) then
        fac =  dble(porder + 1) * courant_num
        w_r(:) = - aimag(log(lamb(:))) / fac
        w_i(:) =    real(log(lamb(:))) / fac
      else
        w_r(:) = - aimag(lamb(:)) / dble(porder + 1)
        w_i(:) =    real(lamb(:)) / dble(porder + 1)
      end if
      
      do n=1, tlev_num
        fac = dble(porder + 1) * courant_num * dble(tlev_slot(n))
        phase(:) = mod( fac * w_r(:), 2.0_RP * PI )
        do p1=1, porder+1
          Ue_num(p1,n) = sum(  EvecMat(p1,:) * theta(:) * exp( fac * w_i(:) )  &
            * ( cos( phase(:) ) - ei * sin( phase(:) ) )                       ) ! * exp(ei * k * 0.0_RP) .. xe=0
        end do
      end do   

      if ( k > 0 ) then
        call sort_eigenval( w_r, w_i, dw_r, dw_i, ipiv_eigenval,   & ! (out)
          wnum_nondim_num(1,k,:), wnum_nondim_num(2,k,:), porder+1 ) ! (in)
      else
        forall(n=1:porder+1) ipiv_eigenval(n) = n
      end if
      wnum_nondim(k+1) = wnum * Helem / dble(porder + 1) 
      wnum_nondim_num(1,k+1,:) = w_r(:)
      wnum_nondim_num(2,k+1,:) = w_i(:)
      eigensol_weight(k+1,:) = abs(theta(ipiv_eigenval(:)))

      do n=1, tlev_num        
        tmpvec1(:) = matmul(Mass, Ue    (:,n))
        tmpvec2(:) = matmul(Mass, Ue_num(:,n))
        gex  = sqrt( real(sum( conjg(Ue    (:,n)) * tmpvec1(:) )) )
        gnum = sqrt( real(sum( conjg(Ue_num(:,n)) * tmpvec2(:) )) )
        phase_error_tmp = sum( conjg(Ue_num(:,n)) * tmpvec1(:) )

        AmplifiFac_comb(n,k+1)  = gnum / gex
        PhaseError_comb(n,k+1) = aimag( log( phase_error_tmp ) ) 
        ! Note:  We do not multiply this phase error by a factor of 1.0_RP / dble( porder+1 )
        ! This reason is becasue we would like to evaluate the phase error at a timing for the various p cases, 
        ! when remaining a value of the Courant number relatated to the effective grid spacing (he/(p+1)). 

        ! if (k==200) then
        !   fac = dble(porder + 1) * courant_num * dble(tlev_slot(n))
        !   write(*,*) "tlev_num=", tlev_slot(n), ":", abs(log(gnum/gex))/dble(tlev_slot(n)), &
        !     ":", w_i(:)*courant_num*(porder+1)
        ! end if
      end do
    end do

    !-- output the results for each mode

    do p1=1, porder+1
      if ( present( tscheme_name ) ) then 
        write(filename,'(a,a,a,a,a,i2.2,a,i2.2,a)') trim(output_dir)//'eigenval_', trim(basis_type), '_', trim(form_type), '_P', porder, '_Mode', p1, '_'//trim(tscheme_name)
      else
        write(filename,'(a,a,a,a,a,i2.2,a,i2.2)') trim(output_dir)//'eigenval_', trim(basis_type), '_', trim(form_type), '_P', porder, '_Mode', p1
      end if
      if ( MF_alph > 0.0_RP ) then
        write(filename,'(3a)') trim(filename), '_MF'//trim(MF_alph_lbl), '.dat' 
      else
        filename = trim(filename)//".dat"
      end if        

      nout = IO_get_available_fid()
      open( nout, file=trim(filename), status='replace' )
      write( nout, '(a)') "# K Re(Km) Im(Km) Weight"
      do k=1, kmax
        write( nout, '(4f20.14)') wnum_nondim(k), wnum_nondim_num(1,k,p1), wnum_nondim_num(2,k,p1), eigensol_weight(k,p1)
      end do
      close( nout )
    end do

    !-- output the results for combined mode

    if ( present( tscheme_name ) ) then 
      write(filename,'(a,a,a,a,a,i2.2,a,a)') trim(output_dir)//'eigenval_', trim(basis_type), '_', trim(form_type), '_P', porder, '_combinedmode_', trim(tscheme_name)
    else
      write(filename,'(a,a,a,a,a,i2.2,a)') trim(output_dir)//'eigenval_', trim(basis_type), '_', trim(form_type), '_P', porder, '_combinedmode'
    end if
    if ( MF_alph > 0.0_RP ) then
      write(filename,'(3a)') trim(filename), '_MF'//trim(MF_alph_lbl), '.dat' 
    else
      filename = trim(filename) // ".dat"
    end if

    nout = IO_get_available_fid()
    open( nout, file=trim(filename), status='replace' )
    write( nout, '(a)') "# K G(n=1) phi(n=1) G(n=10) phi(n=10) G(n=50) phi(n=50) G(n=100) phi(n=100) G(n=200) phi(n=200) G(n=500) phi(n=500) G(n=1000) phi(n=1000) G(n=3000) phi(n=3000)"
    do k=1, kmax
      write( nout, '(17e30.15E3)') wnum_nondim(k), &
        AmplifiFac_comb(1,k), PhaseError_comb(1,k), AmplifiFac_comb(2,k), PhaseError_comb(2,k), &
        AmplifiFac_comb(3,k), PhaseError_comb(3,k), AmplifiFac_comb(4,k), PhaseError_comb(4,k), &
        AmplifiFac_comb(5,k), PhaseError_comb(5,k), AmplifiFac_comb(6,k), PhaseError_comb(6,k), &
        AmplifiFac_comb(7,k), PhaseError_comb(7,k), AmplifiFac_comb(8,k), PhaseError_comb(8,k)
    end do
    close( nout )
    !-

    call MFilter%Final()
    call elem%Final()  

    return
  end subroutine eigen_analysis

  subroutine cosntruct_mat( K, Kp, Km,                &
    Stiff, Minv, phiM1, phiP1, porder, beta, ADV_VEL, &
    form_type )
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
    real(RP), intent(in) :: ADV_VEL
    character(len=*), intent(in) :: form_type

    real(RP) :: bM, bP
    real(RP) :: tmpMat (porder+1,porder+1)
    real(RP) :: tmpMatM(porder+1,porder+1)
    real(RP) :: tmpMatP(porder+1,porder+1)

    integer :: l, m
    !------------------------------------------

    bM = 0.5_RP * ( 1.0_RP - sign(beta, ADV_VEL) )
    bP = 0.5_RP * ( 1.0_RP + sign(beta, ADV_VEL) )

    if ( form_type == 'weak' ) then
      do m=1, porder+1
      do l=1, porder+1
        tmpMat (l,m) = Stiff(l,m) + bM * phiM1(m) * phiM1(l) - bP * phiP1(m) * phiP1(l)
        tmpMatM(l,m) =   bP * phiP1(m) * phiM1(l) 
        tmpMatP(l,m) = - bM * phiM1(m) * phiP1(l) 
      end do
      end do
    else if ( form_type == 'strong' ) then
      do m=1, porder+1
      do l=1, porder+1
        tmpMat (l,m) = - Stiff(m,l) - ( 1.0_RP - bM ) * phiM1(m) * phiM1(l) + ( 1.0_RP - bP ) * phiP1(m) * phiP1(l)
        tmpMatM(l,m) =   bP * phiP1(m) * phiM1(l) 
        tmpMatP(l,m) = - bM * phiM1(m) * phiP1(l) 
      end do
      end do      
    end if

    K (:,:) = matmul(Minv, tmpMat )
    Km(:,:) = matmul(Minv, tmpMatM)
    Kp(:,:) = matmul(Minv, tmpMatP)

    return
  end subroutine cosntruct_mat

  subroutine construct_full_discrete_mat( G, & ! (inout)
    tscheme_name, courant_num, porder )
    implicit none

    integer, intent(in) :: porder
    complex(RP), intent(inout) :: G(porder+1,porder+1)
    character(len=*), intent(in) :: tscheme_name
    real(RP), intent(in) :: courant_num

    complex(RP) :: A(porder+1,porder+1)
    complex(RP) :: Id(porder+1,porder+1)
    complex(RP) :: tmpMat(porder+1,porder+1)
    complex(RP) :: tmpMat2(porder+1,porder+1)

    integer :: rk_nstage
    integer :: rk_s
    integer :: m
    !---------------------------------

    select case( tscheme_name )
    case ( 'RK2' )
      rk_nstage = 2
    case ( 'RK3' )
      rk_nstage = 3
    case ( 'RK4' )
      rk_nstage = 4
    case ( 'RKo4s10' )
      rk_nstage = 10
    case default
      write(*,*) 'Check! Invalid tscheme_name:', tscheme_name
      stop
    end select

    Id(:,:) = 0.0_RP
    do m=1, porder+1
      Id(m,m) = 1.0_RP
    end do

    A(:,:) = G(:,:)
    G(:,:) = Id(:,:)
    tmpMat(:,:) = Id(:,:)
    tmpMat2(:,:) = Id(:,:)
    do rk_s=1, rk_nstage
      tmpMat(:,:) = courant_num / real(rk_s, kind=RP) * matmul(tmpMat, A)
      tmpMat2(:,:) = courant_num * matmul(tmpMat2, A)
      G(:,:) = G(:,:) + tmpMat(:,:)
      if ( tscheme_name == 'RKo4s10' .and. rk_s == 4 ) exit
    end do
    if ( tscheme_name == 'RKo4s10' ) then
      tmpMat2(:,:) = courant_num * matmul(tmpMat2, A)
      G(:,:) = G(:,:) + tmpMat2(:,:) * ( 17.0_RP / 2160.0_RP )

      tmpMat2(:,:) = matmul(courant_num * tmpMat2, A)
      G(:,:) = G(:,:) + tmpMat2(:,:) * ( 7.0_RP / 6480.0_RP )

      tmpMat2(:,:) = matmul(courant_num * tmpMat2, A)
      G(:,:) = G(:,:) + tmpMat2(:,:) / 9720.0_RP

      tmpMat2(:,:) = matmul(courant_num * tmpMat2, A)
      G(:,:) = G(:,:) + tmpMat2(:,:) / 155520.0_RP

      tmpMat2(:,:) = matmul(courant_num * tmpMat2, A)
      G(:,:) = G(:,:) + tmpMat2(:,:) / 4199040.0_RP

      tmpMat2(:,:) = matmul(courant_num * tmpMat2, A)
      G(:,:) = G(:,:) + tmpMat2(:,:) / 251942400.0_RP
    end if

    return
  end subroutine construct_full_discrete_mat


  subroutine calc_muhat( muhat, wnum, Helem, elem, basis_type, Minv )
    implicit none
    type(LineElement), intent(in) :: elem
    complex(RP), intent(out) :: muhat(elem%Np)
    real(RP), intent(in) :: wnum
    real(RP), intent(in) :: Helem
    character(len=*), intent(in) :: basis_type
    real(RP), intent(in) :: Minv(elem%Np,elem%Np)

    complex(RP) :: tmp(elem%Np)
    integer, parameter :: Porder_int = 12
    real(RP) :: basis(Porder_int+1,elem%Np)
    type(LineElement) :: elem_int

    integer :: pp
    real(RP) :: phase(Porder_int+1)
    !--------------------------

    call elem_int%Init(Porder_int, .false.)

    if (basis_type == 'nodal' ) then
      basis(:,:) = Polynominal_GenLagrangePoly(elem%PolyOrder, elem%x1(:), elem_int%x1(:))
    else if (basis_type == 'modal' ) then
      basis(:,:) = Polynominal_GenLegendrePoly(elem%PolyOrder, elem_int%x1(:))
      ! normalization
      do pp=1, elem%PolyOrder+1
        basis(:,pp) = basis(:,pp) * sqrt(dble(pp-1) + 0.5_RP)
      end do      
    end if 

    do pp=1, elem%PolyOrder+1
      phase(:) = 0.5_RP * wnum * Helem * elem_int%x1(:)
      tmp(pp) = sum( elem_int%IntWeight_lgl(:) * ( cos(phase) + ei * sin(phase) ) * basis(:,pp) )
    end do
    muhat(:) = matmul( Minv, tmp )

    call elem_int%Final()
    return
  end subroutine calc_muhat

  subroutine sort_eigenval( w_r, w_i, dw_r, dw_i, ipiv, &
    wr_prev, wi_prev, n )
    use scale_quicksort, only: QUICKSORT_exec_with_idx
    implicit none
    integer, intent(in) :: n
    real(RP), intent(inout) :: w_r(n)
    real(RP), intent(inout) :: w_i(n)
    real(RP), intent(inout) :: dw_r(n)
    real(RP), intent(inout) :: dw_i(n)
    integer, intent(out) :: ipiv(n)
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
    ipiv(:) = -9999

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
          ipiv(p1) = loc_min
          exit
        end if
      end do
    end do

    return
  end subroutine sort_eigenval

end module mod_dg_eigen_analysis