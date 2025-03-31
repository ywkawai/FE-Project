!> module FElib / Data / Utility
!!
!! @par Description
!!           A module for providing utility routines associated with DG field data
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scaleFElib.h"
module scale_meshfield_fvm_util
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io

  use scale_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MeshFieldFVMUtil_interp_FVtoDG_3D
  public :: MeshFieldFVMUtil_interp_FVtoDG_2D
  interface MeshFieldFVMUtil_interp_FVtoDG
    module procedure :: MeshFieldFVMUtil_interp_FVtoDG_2D
    module procedure :: MeshFieldFVMUtil_interp_FVtoDG_3D
  end interface
  public :: MeshFieldFVMUtil_interp_FVtoDG

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  
contains

!> Interpolate 3D FV data into DG nodes
!! @pre Before calling this subroutine, we need to setup scale_comm_cartesC module.
!!
!OCL SERIAL
  subroutine MeshFieldFVMUtil_interp_FVtoDG_3D( out_dg, &
      in_fv, lcmesh3D, elem3D_, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      interp_ord, out_dg_flush_flag )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait    
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D_
    integer, intent(in) :: IS, IE, IA, IHALO
    integer, intent(in) :: JS, JE, JA, JHALO
    integer, intent(in) :: KS, KE, KA, KHALO
    real(RP), intent(inout) :: out_dg(elem3D_%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: in_fv(KA,IA,JA)
    integer, intent(in) :: interp_ord
    logical, intent(in), optional :: out_dg_flush_flag !< Default: true

    integer :: ke_x, ke_y, ke_z
    integer :: i, j, k
    integer :: kelem 

    real(RP) :: tend_fv_cp(KA,IA,JA)
    real(RP) :: dqdx(KA,IA,JA), dqdy(KA,IA,JA), dqdz(KA,IA,JA)
    real(RP) :: dqdxx(KA,IA,JA), dqdyy(KA,IA,JA), dqdzz(KA,IA,JA)
    real(RP) :: dqdxy(KA,IA,JA), dqdxz(KA,IA,JA), dqdyz(KA,IA,JA)

    logical :: out_dg_flush_flag_
    !----------------------------------

    if ( present(out_dg_flush_flag) ) then
      out_dg_flush_flag_ = out_dg_flush_flag
    else
      out_dg_flush_flag_ = .true.
    end if
    
    tend_fv_cp(:,:,:) = in_fv(:,:,:)
    if ( interp_ord > 0 ) then
      call COMM_vars8( tend_fv_cp(:,:,:), 1 )
      call COMM_wait( tend_fv_cp(:,:,:), 1, .true. )
    end if

    !$omp parallel private(ke_x,ke_y,ke_z,kelem, i,j,k) 

    if ( out_dg_flush_flag_ ) then
      !$omp workshare
      out_dg(:,:) = 0.0_RP
      !$omp end workshare
    end if

    if ( interp_ord > 0 ) then
      !$omp do collapse(2)
      do j=JS-1, JE+1
      do i=IS, IE
      do k=KS, KE
        dqdx(k,i,j) = 0.25_RP * ( tend_fv_cp(k,i+1,j) - tend_fv_cp(k,i-1,j) )
      enddo
      enddo
      enddo
      !$omp do collapse(2)
      do j=JS, JE
      do i=IS-1, IE+1
      do k=KS, KE
        dqdy(k,i,j) = 0.25_RP * ( tend_fv_cp(k,i,j+1) - tend_fv_cp(k,i,j-1) )
      enddo
      enddo
      enddo
      !$omp do collapse(2)
      do j=JS, JE
      do i=IS, IE
        do k=KS+1, KE-1
          dqdz(k,i,j) = 0.25_RP * ( tend_fv_cp(k+1,i,j) - tend_fv_cp(k-1,i,j) )
        enddo
        dqdz(KS,i,j) = 0.25_RP * ( - tend_fv_cp(KS+2,i,j) + 4.0_RP * tend_fv_cp(KS+1,i,j) - 3.0_RP * tend_fv_cp(KS,i,j) )
        dqdz(KE,i,j) = 0.25_RP * ( 3.0_RP * tend_fv_cp(KE,i,j) - 4.0_RP * tend_fv_cp(KE-1,i,j) + tend_fv_cp(KE-2,i,j) )
      enddo
      enddo       
    end if
    if ( interp_ord > 1 ) then
      !$omp do collapse(2)
      do j=JS, JE
      do i=IS, IE
      do k=KS, KE
        dqdxx(k,i,j) = 0.25_RP * ( tend_fv_cp(k,i+1,j) - 2.0_RP * tend_fv_cp(k,i,j) + tend_fv_cp(k,i-1,j) )
      enddo
      enddo
      enddo
      !$omp do collapse(2)
      do j=JS, JE
      do i=IS, IE
      do k=KS, KE
        dqdyy(k,i,j) = 0.25_RP * ( tend_fv_cp(k,i,j+1) - 2.0_RP * tend_fv_cp(k,i,j) + tend_fv_cp(k,i,j-1) )
      enddo
      enddo
      enddo
      !$omp do collapse(2)
      do j=JS, JE
      do i=IS, IE
      do k=KS, KE
        dqdxy(k,i,j) = 0.25_RP * ( dqdx(k,i,j+1) - dqdx(k,i,j-1) )
      enddo
      enddo
      enddo 
      !$omp do collapse(2)
      do j=JS, JE
      do i=IS, IE
        do k=KS+1, KE-1
          dqdzz(k,i,j) = 0.25_RP * ( tend_fv_cp(k+1,i,j) - 2.0_RP * tend_fv_cp(k,i,j) + tend_fv_cp(k-1,i,j) )
          dqdxz(k,i,j) = 0.25_RP * ( dqdx(k+1,i,j) - dqdx(k-1,i,j) )
          dqdyz(k,i,j) = 0.25_RP * ( dqdy(k+1,i,j) - dqdy(k-1,i,j) )
        enddo
        dqdzz(KS,i,j) = dqdzz(KS+1,i,j); dqdzz(KE,i,j) = dqdzz(KE-1,i,j); 
        dqdxz(KS,i,j) = 0.5_RP * ( dqdx(KS+1,i,j) - dqdx(KS,i,j) ); dqdxz(KE,i,j) = 0.5_RP * ( dqdx(KE,i,j) - dqdx(KE-1,i,j) ); 
        dqdyz(KS,i,j) = 0.5_RP * ( dqdy(KS+1,i,j) - dqdy(KS,i,j) ); dqdyz(KE,i,j) = 0.5_RP * ( dqdy(KE,i,j) - dqdy(KE-1,i,j) ); 
      enddo
      enddo
    end if

    select case(interp_ord)
    case(0) !- Piecewise constant
      !$omp do collapse(2)
      do ke_z=1, lcmesh3D%NeZ
      do ke_y=1, lcmesh3D%NeY
      do ke_x=1, lcmesh3D%NeX
        kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
        i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z
        out_dg(:,kelem) = out_dg(:,kelem) &
                         + in_fv(k,i,j)
      enddo
      enddo
      enddo    
    case(1) !- Piecewise linear reconstruction
      !$omp do collapse(2)
      do ke_z=1, lcmesh3D%NeZ
      do ke_y=1, lcmesh3D%NeY
      do ke_x=1, lcmesh3D%NeX
        kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
        i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z
        out_dg(:,kelem) = out_dg(:,kelem) &
                          + in_fv(k,i,j)             &
                          + dqdx(k,i,j) * elem3D_%x1(:) &
                          + dqdy(k,i,j) * elem3D_%x2(:) &
                          + dqdz(k,i,j) * elem3D_%x3(:)
      enddo
      enddo
      enddo          
    case(2) !- Piecewise quadratic reconstruction
      !$omp do collapse(2)
      do ke_z=1, lcmesh3D%NeZ
      do ke_y=1, lcmesh3D%NeY
      do ke_x=1, lcmesh3D%NeX
        kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
        i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z
        out_dg(:,kelem) = out_dg(:,kelem) &
                          + in_fv(k,i,j)      &
                          + dqdx(k,i,j) * elem3D_%x1(:) &
                          + dqdy(k,i,j) * elem3D_%x2(:) &
                          + dqdz(k,i,j) * elem3D_%x3(:) &
                          + dqdxy(k,i,j) * elem3D_%x1(:) * elem3D_%x2(:) &
                          + dqdxz(k,i,j) * elem3D_%x1(:) * elem3D_%x3(:) &
                          + dqdyz(k,i,j) * elem3D_%x2(:) * elem3D_%x3(:) &
                          + 0.5_RP * dqdxx(k,i,j) * ( elem3D_%x1(:)**2 - 1.0_RP / 3.0_RP ) &
                          + 0.5_RP * dqdyy(k,i,j) * ( elem3D_%x2(:)**2 - 1.0_RP / 3.0_RP ) &
                          + 0.5_RP * dqdzz(k,i,j) * ( elem3D_%x3(:)**2 - 1.0_RP / 3.0_RP )
      enddo
      enddo
      enddo 
    end select
    !$omp end parallel
    return
  end subroutine MeshFieldFVMUtil_interp_FVtoDG_3D

!> Interpolate 2D FV data into DG nodes
!! @pre Before calling this subroutine, we need to setup scale_comm_cartesC module.
!!
!OCL SERIAL
  subroutine MeshFieldFVMUtil_interp_FVtoDG_2D( out_dg, &
    in_fv, lcmesh2D, elem2D_, IS, IE, IA, IHALO, JS, JE, JA, JHALO, &
    interp_ord, out_dg_flush_flag )
  use scale_comm_cartesC, only: &
     COMM_vars8, &
     COMM_wait    
  implicit none
  class(LocalMesh2D), intent(in) :: lcmesh2D
  class(ElementBase2D), intent(in) :: elem2D_
  integer, intent(in) :: IS, IE, IA, IHALO
  integer, intent(in) :: JS, JE, JA, JHALO
  real(RP), intent(inout) :: out_dg(elem2D_%Np,lcmesh2D%NeA)
  real(RP), intent(in) :: in_fv(IA,JA)
  integer, intent(in) :: interp_ord
  logical, intent(in), optional :: out_dg_flush_flag !< Default: true

  integer :: ke_x, ke_y
  integer :: i, j
  integer :: kelem 

  real(RP) :: tend_fv_cp(IA,JA)
  real(RP) :: dqdx(IA,JA), dqdy(IA,JA)
  real(RP) :: dqdxx(IA,JA), dqdyy(IA,JA)
  real(RP) :: dqdxy(IA,JA)
  
  logical :: out_dg_flush_flag_
  !----------------------------------

  if ( present(out_dg_flush_flag) ) then
    out_dg_flush_flag_ = out_dg_flush_flag
  else
    out_dg_flush_flag_ = .true.
  end if

  tend_fv_cp(:,:) = in_fv(:,:)
  if ( interp_ord > 0 ) then
    call COMM_vars8( tend_fv_cp(:,:), 1 )
    call COMM_wait( tend_fv_cp(:,:), 1, .false. )
  end if

  !$omp parallel private(ke_x,ke_y,kelem, i,j) 

  if ( out_dg_flush_flag_ ) then
    !$omp workshare
    out_dg(:,:) = 0.0_RP
    !$omp end workshare
  end if

  if ( interp_ord > 0 ) then
    !$omp do
    do j=JS-1, JE+1
    do i=IS, IE
      dqdx(i,j) = 0.25_RP * ( tend_fv_cp(i+1,j) - tend_fv_cp(i-1,j) )
    enddo
    enddo
    !$omp do
    do j=JS, JE
    do i=IS-1, IE+1
      dqdy(i,j) = 0.25_RP * ( tend_fv_cp(i,j+1) - tend_fv_cp(i,j-1) )
    enddo
    enddo
  end if
  if ( interp_ord > 1 ) then
    !$omp do
    do j=JS, JE
    do i=IS, IE
      dqdxx(i,j) = 0.25_RP * ( tend_fv_cp(i+1,j) - 2.0_RP * tend_fv_cp(i,j) + tend_fv_cp(i-1,j) )
    enddo
    enddo
    !$omp do
    do j=JS, JE
    do i=IS, IE
      dqdyy(i,j) = 0.25_RP * ( tend_fv_cp(i,j+1) - 2.0_RP * tend_fv_cp(i,j) + tend_fv_cp(i,j-1) )
    enddo
    enddo
    !$omp do
    do j=JS, JE
    do i=IS, IE
      dqdxy(i,j) = 0.25_RP * ( dqdx(i,j+1) - dqdx(i,j-1) )
    enddo
    enddo 
  end if

  select case(interp_ord)
  case(0) !- Piecewise constant
    !$omp do
    do ke_y=1, lcmesh2D%NeY
    do ke_x=1, lcmesh2D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh2D%NeX
      i = IHALO + ke_x; j = JHALO + ke_y;
      out_dg(:,kelem) = out_dg(:,kelem) &
                       + in_fv(i,j)
    enddo
    enddo    
  case(1) !- Piecewise linear reconstruction
    !$omp do
    do ke_y=1, lcmesh2D%NeY
    do ke_x=1, lcmesh2D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh2D%NeX 
      i = IHALO + ke_x; j = JHALO + ke_y
      out_dg(:,kelem) = out_dg(:,kelem) &
                        + in_fv(i,j)    &
                        + dqdx(i,j) * elem2D_%x1(:) &
                        + dqdy(i,j) * elem2D_%x2(:)
    enddo
    enddo
  case(2) !- Piecewise quadratic reconstruction
    !$omp do
    do ke_y=1, lcmesh2D%NeY
    do ke_x=1, lcmesh2D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh2D%NeX
      i = IHALO + ke_x; j = JHALO + ke_y
      out_dg(:,kelem) = out_dg(:,kelem) &
                        + in_fv(i,j)    &
                        + dqdx(i,j) * elem2D_%x1(:) &
                        + dqdy(i,j) * elem2D_%x2(:) &
                        + dqdxy(i,j) * elem2D_%x1(:) * elem2D_%x2(:) &
                        + 0.5_RP * dqdxx(i,j) * ( elem2D_%x1(:)**2 - 1.0_RP / 3.0_RP ) &
                        + 0.5_RP * dqdyy(i,j) * ( elem2D_%x2(:)**2 - 1.0_RP / 3.0_RP )
    enddo
    enddo
  end select
  !$omp end parallel
  return
end subroutine MeshFieldFVMUtil_interp_FVtoDG_2D

!- private ----------------------


end module scale_meshfield_fvm_util