#include "scalelib.h"
program test_modalfilter
    use scale_precision
    use scale_element_line
    use scale_element_hexahedral
    use scale_element_modalfilter
    implicit none

    type(LineElement) :: elem1D
    type(HexahedralElement) :: elem
    type(ModalFilter) :: filter
    type(ModalFilter) :: filter1D

    integer, parameter :: porder=7
    integer, parameter :: Np1D=porder+1
    integer, parameter :: Np=(porder+1)**3

    integer :: x, y, z, kp
    real(RP) :: q_in(Np)
    real(RP) :: q_out_ori(Np)
    real(RP) :: q_out_new(Np)
    real(RP) :: error(Np)
    real(RP) :: rmse

    real(RP) :: filter1D_tr(Np1D, Np1D)

    call elem%Init(porder, porder, .false.)
    call elem1D%Init(porder, .false.)

    call filter%Init(elem, 0D0, 0D0, 16, 0D0, 0D0, 16)
    call filter1D%Init(elem1D, 0D0, 0D0, 16)

    filter1D_tr(:,:) = transpose(filter1D%FilterMat(:,:))

    do z = 1, 8
    do y = 1, 8
    do x = 1, 8
        kp = x + (y - 1) * 8 + (z - 1) * 8 * 8
        q_in(kp) = sin(dble(x)) * sin(dble(y)) * sin(dble(z)) 
    end do
    end do
    end do

    q_out_ori(:) = matmul(filter%FilterMat, q_in)

    q_out_new(:) = 0.0_RP

    call apply_filter_xyz(filter1D%FilterMat, filter1D_tr, q_in, q_out_new)

    !-- calculate error: err = q_out_new - q_out_ori 
    do kp=1, Np
        error(kp) = q_out_new(kp) - q_out_ori(kp)
    end do

    !-- save error
    open(10, file="modalfilter_err.dat", access='stream', form='unformatted', status='new')
    write(10) error(:)
    close(10)

    !-- calculate rmse
    do kp=1, Np
        rmse = rmse + (q_out_new(kp) - q_out_ori(kp)) ** 2
    end do

    rmse = sqrt(rmse / Np)

    print *, "RMSE: ", rmse
        
   
contains
    subroutine apply_filter_xyz(filterMat, filterMat_tr, q_in, q_out_new)
        implicit none

        real(RP), intent(in) :: filterMat(8, 8)
        real(RP), intent(in) :: filterMat_tr(8, 8)
        real(RP), intent(inout) :: q_in(8,8,8)
        real(RP), intent(inout) :: q_out_new(8,8,8)

        integer :: i, j, k

        !-- x direction
        do k=1, 8
        do j=1, 8
        do i=1, 8
            q_out_new(i,j,k) = filterMat(i,1) * q_in(1,j,k) + &
                               filterMat(i,2) * q_in(2,j,k) + & 
                               filterMat(i,3) * q_in(3,j,k) + & 
                               filterMat(i,4) * q_in(4,j,k) + & 
                               filterMat(i,5) * q_in(5,j,k) + & 
                               filterMat(i,6) * q_in(6,j,k) + & 
                               filterMat(i,7) * q_in(7,j,k) + & 
                               filterMat(i,8) * q_in(8,j,k) 
        end do
        end do
        end do

        !-- y direction
        do k=1, 8
        do j=1, 8
        do i=1, 8
            q_in(i,j,k) = q_out_new(i,1,k) * filterMat_tr(1,j) + &
                          q_out_new(i,2,k) * filterMat_tr(2,j) + &
                          q_out_new(i,3,k) * filterMat_tr(3,j) + &
                          q_out_new(i,4,k) * filterMat_tr(4,j) + &
                          q_out_new(i,5,k) * filterMat_tr(5,j) + &
                          q_out_new(i,6,k) * filterMat_tr(6,j) + &
                          q_out_new(i,7,k) * filterMat_tr(7,j) + &
                          q_out_new(i,8,k) * filterMat_tr(8,j)
        end do
        end do
        end do

        !-- z direction
        do k=1, 8
        do j=1, 8
        do i=1, 8
            q_out_new(i,j,k) = q_in(i,j,1) * filterMat_tr(1,k) + &
                               q_in(i,j,2) * filterMat_tr(2,k) + & 
                               q_in(i,j,3) * filterMat_tr(3,k) + & 
                               q_in(i,j,4) * filterMat_tr(4,k) + & 
                               q_in(i,j,5) * filterMat_tr(5,k) + & 
                               q_in(i,j,6) * filterMat_tr(6,k) + & 
                               q_in(i,j,7) * filterMat_tr(7,k) + & 
                               q_in(i,j,8) * filterMat_tr(8,k)
        end do
        end do
        end do

    end subroutine apply_filter_xyz
end program test_modalfilter