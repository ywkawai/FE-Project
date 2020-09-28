#include "scalelib.h"
program test_polynominal
  use scale_precision
  use scale_io
  use scale_polynominal
  implicit none

  !----------------------------------------------
  
  !* n=1
  real(RP), parameter :: lgl_pts_n1(2) = (/ -1.0_RP, 1.0_RP /)
  real(RP), parameter :: lgl_intw_n1(2) = (/ 1.0_RP, 1.0_RP /)
  
  real(RP), parameter :: gl_pts_n1(1) = (/ 0.0_RP /)
  real(RP), parameter :: gl_intw_n1(1) = (/ 2.0_RP /)

  !* n=2
  real(RP), parameter :: lgl_pts_n2(3) = (/ -1.0_RP, 0.0_RP, 1.0_RP /)
  real(RP), parameter :: lgl_intw_n2(3) = (/ 1.0_RP/3.0_RP, 4.0_RP/3.0_RP, 1.0_RP/3.0_RP /)
  
  real(RP), parameter :: gl_pts_n2(2) = (/ -1.0_RP/sqrt(3.0_RP), 1.0_RP/sqrt(3.0_RP) /)
  real(RP), parameter :: gl_intw_n2(2) = (/ 1.0_RP, 1.0_RP /)

  !* n=3
  real(RP), parameter :: lgl_pts_n3(4) = (/ -1.0_RP, - sqrt(1.0_RP/5.0_RP), + sqrt(1.0_RP/5.0_RP), 1.0_RP /)
  real(RP), parameter :: lgl_intw_n3(4) = (/ 1.0_RP/6.0_RP, 5.0_RP/6.0_RP, 5.0_RP/6.0_RP, 1.0_RP/6.0_RP /)

  real(RP), parameter :: gl_pts_n3(3) = (/ -sqrt(3.0_RP)/sqrt(5.0_RP), 0.0_RP, sqrt(3.0_RP)/sqrt(5.0_RP) /)
  real(RP), parameter :: gl_intw_n3(3) = (/ 5.0_RP/9.0_RP, 8.0_RP/9.0_RP, 5.0_RP/9.0_RP /)

  !* n=4
  real(RP), parameter :: lgl_pts_n4(5) = (/ -1.0_RP,  -sqrt(21.0_RP)/7.0_RP, 0.0_RP, +sqrt(21.0_RP)/7.0_RP, 1.0_RP /)
  real(RP), parameter :: lgl_intw_n4(5) = (/ 1.0_RP/10.0_RP, 49.0_RP/90.0_RP, 32.0_RP/45.0_RP, 49.0_RP/90.0_RP, 1.0_RP/10.0_RP /)

  real(DP), parameter :: gl_pts_n4(4) = (/ &
    -sqrt(3.0_RP + 2.0_RP*sqrt(6.0_RP/5.0_RP))/sqrt(7.0_RP), -sqrt(3.0_RP - 2.0_RP*sqrt(6.0_RP/5.0_RP))/sqrt(7.0_RP), &
    +sqrt(3.0_RP - 2.0_RP*sqrt(6.0_RP/5.0_RP))/sqrt(7.0_RP), +sqrt(3.0_RP + 2.0_RP*sqrt(6.0_RP/5.0_RP))/sqrt(7.0_RP) /)
  real(DP), parameter :: gl_intw_n4(4) = (/ &
    (18.0_RP - sqrt(30.0_RP))/36.0_RP, (18.0_RP + sqrt(30.0_RP))/36.0_RP, &
    (18.0_RP + sqrt(30.0_RP))/36.0_RP, (18.0_RP - sqrt(30.0_RP))/36.0_RP /)  

 !----------------------------------------------------------------------------------------

  call chechk_value(1, lgl_pts_n1, lgl_intw_n1, gl_pts_n1, gl_intw_n1)
  call chechk_value(2, lgl_pts_n2, lgl_intw_n2, gl_pts_n2, gl_intw_n2)
  call chechk_value(3, lgl_pts_n3, lgl_intw_n3, gl_pts_n3, gl_intw_n3)
  call chechk_value(4, lgl_pts_n4, lgl_intw_n4, gl_pts_n4, gl_intw_n4)

  call output_expandfunc(4, 200)
  write(*,*) "test_polynominal has been succeeded!"
contains
  subroutine chechk_value(Norder, lgl_pts_ans, lgl_intw_ans, gl_pts_ans, gl_intw_ans)
    integer, intent(in) :: Norder
    real(RP), intent(in) :: lgl_pts_ans(Norder+1)
    real(RP), intent(in) :: lgl_intw_ans(Norder+1)
    real(RP), intent(in) :: gl_pts_ans(Norder)
    real(RP), intent(in) :: gl_intw_ans(Norder)

    real(RP) :: LGLpts(Norder+1)
    real(RP) :: Legendre_LGLPts(size(LGLpts),Norder+1)
    real(RP) :: LGL_intw(Norder+1)
    real(RP) :: GLpts(Norder)
    real(RP) :: GL_intw(Norder)    
    integer :: l
    
    real(RP), parameter :: CHECK_EPS = 5.0E-15_RP

   !--------------------------------------------------------------------

    LGLpts(:) = Polynominal_GenGaussLobattoPt(Norder)
    Legendre_LGLPts(:,:) = Polynominal_GenLegendrePoly(Norder, LGLpts)
    LGL_intw(:) = Polynominal_GenGaussLobattoPtIntWeight(Norder)

    GLpts(:) = Polynominal_GenGaussLegendrePt(Norder)
    GL_intw(:) = Polynominal_GenGaussLegendrePtIntWeight(Norder)

    write(*,*) "******** Norder=", Norder
    
    write(*,*) "- LGLpts:", LGLpts(:)
    if (maxval(abs(LGLpts - lgl_pts_ans)) > CHECK_EPS) then
      LOG_ERROR('chechk_value',*) 'Error of LGL points are large. Check!', &
        "answer=", lgl_pts_ans
      stop
    end if

    write(*,*) "--"
    write(*,*) "- Legendre (at LGLpts): "
    do l=1, Norder+1
      write(*,*) l-1, ":", Legendre_LGLPts(:,l)
    end do

    write(*,*) "--"
    write(*,*) "LGL weight:",  LGL_intw(:)
    if (maxval(abs(LGL_intw - lgl_intw_ans)) > CHECK_EPS) then
      LOG_ERROR('chechk_value',*) 'Error of LGL weights are large. Check!', &
        "answer=", lgl_intw_ans
      stop
    end if

    write(*,*) "--"    
    write(*,*) "- GLpts:", GLpts(:)
    if (maxval(abs(GLpts - gl_pts_ans)) > CHECK_EPS) then
      LOG_ERROR('chechk_value',*) 'Error of GL points are large. Check!', &
        "answer=", gl_pts_ans
      stop
    end if

    write(*,*) "--"
    write(*,*) "GL weight:",  GL_intw(:)
    if (maxval(abs(GL_intw - gl_intw_ans)) > CHECK_EPS) then
      LOG_ERROR('chechk_value',*) 'Error of GL weights are large. Check!', &
        "answer=", gl_intw_ans
      stop
    end if    
    write(*,*) "**************************************"

  end subroutine chechk_value

  subroutine output_expandfunc(Norder, Nnode_intrp)
    integer, intent(in) :: Norder
    integer, intent(in) :: Nnode_intrp

    real(RP) :: pts(Nnode_intrp)
    real(RP) :: lgl_pts(Norder+1)
    real(RP) :: val_Legendre(Nnode_intrp,Norder+1)
    real(RP) :: val_Lagrange(Nnode_intrp,Norder+1)
    
    integer :: i, n
    character(len=32) :: file_name
    !--------------------------------------------------------------------

    pts(1) = -1.0_RP
    do i=2,Nnode_intrp
      pts(i) = pts(i-1) + 2.0_RP/dble(Nnode_intrp-1)
    end do 

    val_Legendre(:,:) = Polynominal_GenLegendrePoly(Norder, pts)
    do n=0, Norder
      write(*,*) 'output legendre polynominal: N=', n
      write(file_name,'(A,I2.2,A)') 'legendreN',n,'.dat'
      open(10, file=trim(file_name))
      write(10,*) '#-- Legende polynominal: Norder=', n
      write(10,*) '# x   Pn(x)'
      ! output
      do i=1, Nnode_intrp
        write(10,'(E18.8,E18.8)') pts(i), val_Legendre(i,n+1)
      end do
      close(10)
    end do

    lgl_pts(:) = Polynominal_GenGaussLobattoPt(Norder)
    val_Lagrange(:,:) = Polynominal_GenLagrangePoly(Norder, lgl_pts, pts)
    do n=0, Norder
      write(*,*) 'output lagrange polynominal: N=', n
      write(file_name,'(A,I2.2,A)') 'lagrangeN',n,'.dat'
      open(10, file=trim(file_name))
      write(10,*) '#-- Lagrange polynominal: n=', n
      write(10,*) '# x   phi_n(x)'
      ! output
      do i=1, Nnode_intrp
        write(10,'(E18.8,E18.8)') pts(i), val_Lagrange(i,n+1)
      end do
      close(10)
    end do    
  end subroutine output_expandfunc

end program test_polynominal
