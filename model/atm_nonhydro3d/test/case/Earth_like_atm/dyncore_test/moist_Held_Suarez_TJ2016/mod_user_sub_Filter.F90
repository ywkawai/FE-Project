#include "scalelib.h"
module mod_user_sub_Filter

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: PRC_abort 
  use scale_const, only: &
    Rdry => CONST_Rdry,    &
    Rvap => CONST_Rvap,    &
    CPdry => CONST_CPdry,  & 
    CVdry => CONST_CVdry,  & 
    LHV0 => CONST_LHV0,    &
    PRES00 => CONST_PRE00, &
    Grav => CONST_GRAV, &
    OHM => CONST_OHM,   &
    RPlanet => CONST_RADIUS, &
    PI => CONST_PI
  use scale_tracer, only: &
    TRACER_inq_id

  use scale_element_base, only: ElementBase1D, ElementBase3D
  use scale_element_line, only: LineElement
  use scale_element_hexahedral, only: HexahedralElement
  use scale_mesh_base2d, only: MeshBase2D
  use scale_mesh_base3d, only: MeshBase3D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_localmesh_3d, only: LocalMesh3D  
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_localmeshfield_base, only: LocalMeshFieldBase
  use scale_mesh_cubedspheredom3d, only: &
    MeshCubedSphereDom3D
  use scale_meshfieldcomm_cubedspheredom3d, only: &
    MeshFieldCommCubedSphereDom3D
  use scale_meshfieldcomm_base, only: &
    MeshFieldContainer

  use mod_user_sub_data, only: &
    MeshFieldCommCubedSphereDom3D_2
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  !
  type, public :: Filter
    real(RP), allocatable :: FilterMat_h1D(:,:)
  
    type(MeshFieldCommCubedSphereDom3D) :: vars_comm
    type(MeshFieldCommCubedSphereDom3D_2) :: vars_comm_2
  contains
    procedure :: Init => USER_sub_Filter_Init
    procedure :: Apply => USER_sub_apply_filter
    procedure, private :: Prepair_FilterMat => prepair_filter_matrix
  end type

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
contains
  subroutine USER_sub_Filter_Init( this, FilterShape, FilterWidthFac, mesh3D )
    implicit none
    class(Filter), intent(inout) :: this
    class(MeshBase3D), intent(in) :: mesh3D
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: FilterWidthFac
    !------------------------------------------

    LOG_NEWLINE
    LOG_INFO("USER_Filter_setup",*) 'Setup'

    !--
    call this%Prepair_FilterMat( mesh3D%refElem3D, FilterShape, FilterWidthFac )

    select type(mesh3D)
    type is (MeshCubedSphereDom3D)
      call this%vars_comm%Init( 1, 0, 0, mesh3D )
      call this%vars_comm_2%Init( 1, 0, 0, mesh3D )
    end select
    
    return
  end subroutine USER_sub_Filter_Init

!OCL serial
  subroutine USER_sub_apply_filter( this, q, &
    mesh3D )
    implicit none
    class(Filter), intent(inout) :: this
    type(MeshField3D), intent(inout), target :: q
    class(MeshBase3D), intent(in), target :: mesh3D

    class(LocalMesh3D), pointer :: lmesh3D
    class(ElementBase3D), pointer :: elem3D
    real(RP), allocatable :: tmp3D(:,:,:,:)

    integer :: n
    type(MeshFieldContainer) :: comm_vars_list(1)
    !-----------------------------------------------------

    comm_vars_list(1)%field3d => q
    ! call this%vars_comm%Put(comm_vars_list, 1)
    ! call this%vars_comm%Exchange()
    ! call this%vars_comm%Get(comm_vars_list, 1)
    call this%vars_comm_2%Put(comm_vars_list, 1)
    call this%vars_comm_2%Exchange()
    call this%vars_comm_2%Get(comm_vars_list, 1)

    do n=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(n)
      elem3D => lmesh3D%refElem3D
!      allocate( tmp3D(0:elem3D%Nnode_h1D+1,0:elem3D%Nnode_h1D+1,elem3D%Nnode_v,lmesh3D%Ne) )
      allocate( tmp3D(-elem3D%Nnode_h1D+1:2*elem3D%Nnode_h1D,-elem3D%Nnode_h1D+1:2*elem3D%Nnode_h1D,elem3D%Nnode_v,lmesh3D%Ne) )

      call extract_tmp3D( tmp3D, &
        q%local(n)%val, q%local(n)%val, lmesh3D, elem3D, lmesh3D%VMapP, 0 )
      call apply_filter1d_x( q%local(n)%val, &
        tmp3D, this%FilterMat_h1D, elem3D%Nnode_h1D, elem3D%Nnode_h1D, elem3D%Nnode_v, lmesh3D%Ne, lmesh3D, elem3D )

      deallocate( tmp3D )
    end do

    comm_vars_list(1)%field3d => q
    ! call this%vars_comm%Put(comm_vars_list, 1)
    ! call this%vars_comm%Exchange()
    ! call this%vars_comm%Get(comm_vars_list, 1)
    call this%vars_comm_2%Put(comm_vars_list, 1)
    call this%vars_comm_2%Exchange()
    call this%vars_comm_2%Get(comm_vars_list, 1)

    do n=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(n)
      elem3D => lmesh3D%refElem3D

      allocate( tmp3D(-elem3D%Nnode_h1D+1:2*elem3D%Nnode_h1D,-elem3D%Nnode_h1D+1:2*elem3D%Nnode_h1D,elem3D%Nnode_v,lmesh3D%Ne) )

      call extract_tmp3D( tmp3D, &
        q%local(n)%val, q%local(n)%val, lmesh3D, elem3D, lmesh3D%VMapP, 1 )
      call apply_filter1d_y( q%local(n)%val, &
        tmp3D, this%FilterMat_h1D, elem3D%Nnode_h1D, elem3D%Nnode_h1D, elem3D%Nnode_v, lmesh3D%Ne, lmesh3D, elem3D )

      deallocate( tmp3D )
    end do
    
    return
  end subroutine USER_sub_apply_filter

!- Private subroutines -------------------
  subroutine prepair_filter_matrix( this, elem3D, FilterShape, FilterWidthFac )
    use scale_polynominal, only: &
      Polynominal_GenLagrangePoly
    implicit none
    class(Filter), intent(inout) :: this
    class(ElementBase3D), intent(in) :: elem3D
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: FilterWidthFac

    integer :: p1, p2

    real(RP) :: filterW
    real(RP) :: filterW2
    type(LineElement) :: elem1D
    type(LineElement) :: elem1D_intrp

    real(RP), allocatable :: lag(:,:)
    real(RP), allocatable :: filter_func(:)
    !--------------------------

    call elem1D%Init( elem3D%PolyOrder_h, .false. )
    call elem1D_intrp%Init( min(2*elem3D%PolyOrder_h, 11), .false. )
    LOG_INFO("prepair_filter_matrix Pos:",*) elem1D%x1(:)
    LOG_INFO("prepair_filter_matrix IntWeight:",*) elem1D%IntWeight_lgl(:)

!    allocate( this%FilterMat_h1D(elem1D%Np,0:elem1D%Np+1) )
    allocate( this%FilterMat_h1D(elem1D%Np,-elem1D%Np+1:elem1D%Np+elem1D%Np) )
    allocate( lag(elem1D_intrp%Np,elem1D%PolyOrder+1) )
    allocate( filter_func(elem1D_intrp%Np) )
    lag(:,:) = Polynominal_GenLagrangePoly( elem1D%PolyOrder, elem1D%x1, elem1D_intrp%x1 )

    this%FilterMat_h1D(:,:) = 0.0_RP

    filterW = FilterWidthFac * 2.0_RP / real(elem1D%Np,kind=RP)
    LOG_INFO("prepair_filter_matrix FilterW:",*) filterW

    do p1=1, elem1D%Np

      !--
      call calc_filter_kenrnel( filter_func, &
        FilterShape, elem1D_intrp%x1(:) - elem1D%x1(p1), filterW, elem1D_intrp%Np )      
      do p2=1, elem1D%Np
        this%FilterMat_h1D(p1,p2) = sum( elem1D_intrp%IntWeight_lgl(:) * lag(:,p2) * filter_func(:) )
      end do
      filterW2 = sum(elem1D_intrp%IntWeight_lgl(:) * filter_func(:))

      !--
      call calc_filter_kenrnel( filter_func, &
        FilterShape, elem1D_intrp%x1(:) - 2.0_RP - elem1D%x1(p1), filterW, elem1D_intrp%Np )
      
      ! this%FilterMat_h1D(p1,0) = sum( elem1D_intrp%IntWeight_lgl(:) * filter_func(:) )
      do p2=1, elem1D%Np
        this%FilterMat_h1D(p1,-elem1D%Np+p2) = sum( elem1D_intrp%IntWeight_lgl(:) * lag(:,p2) * filter_func(:) )
      end do
      filterW2 = filterW2 + sum(elem1D_intrp%IntWeight_lgl(:) * filter_func(:))

      !--
      call calc_filter_kenrnel( filter_func, &
        FilterShape, elem1D_intrp%x1(:) + 2.0_RP - elem1D%x1(p1), filterW, elem1D_intrp%Np )
      
      ! this%FilterMat_h1D(p1,elem1D%Np+1) = sum( elem1D_intrp%IntWeight_lgl(:) * filter_func(:) )
      do p2=1, elem1D%Np
        this%FilterMat_h1D(p1,elem1D%Np+p2) = sum( elem1D_intrp%IntWeight_lgl(:) * lag(:,p2) * filter_func(:) )
      end do
      filterW2 = filterW2 + sum(elem1D_intrp%IntWeight_lgl(:) * filter_func(:))
      
      LOG_INFO("prepair_filter_matrix",*) p1, ":", filterW2, ":", this%FilterMat_h1D(p1,:)
      this%FilterMat_h1D(p1,:) = this%FilterMat_h1D(p1,:) / filterW2 ! normalization
    end do

    call elem1D%Final()
    call elem1D_intrp%Final()

    return
  end subroutine prepair_filter_matrix
!OCL serial
  subroutine calc_filter_kenrnel( filter_kernel, &
      FilterShape, x, filter_width, Np )
    implicit none
    integer, intent(in) :: Np
    real(RP), intent(out) :: filter_kernel(Np)
    character(*), intent(in) :: FilterShape
    real(RP), intent(in) :: x(Np)
    real(RP), intent(in) :: filter_width
    !------------------------------------------------

    ! 
    select case(FilterShape)
    case( "GAUSSIAN")
      filter_kernel(:) = exp( - (x(:)/filter_width)**2 )
    case( "TOPHAT")
      filter_kernel(:) = 0.5_RP * ( sign(1.0_RP, x(:) + 0.5_RP * filter_width) - sign(1.0_RP, x(:) - 0.5_RP * filter_width) )
    case default
      LOG_ERROR("calc_filter_kenrnel",*) "The specified FilterShape is not supported. Check!", FilterShape
      call PRC_abort
    end select

    return
  end subroutine calc_filter_kenrnel

!OCL SERIAL
  subroutine apply_filter1d_x( q, q0, Filter1D, Npx, Npy, Npz, Ne, lmesh3D, elem3D )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    integer, intent(in) :: Npx, Npy, Npz, Ne
    real(RP), intent(out) :: q(Npx,Npy,Npz,lmesh3D%NeA)
!    real(RP), intent(in) :: q0(0:Npx+1,0:Npy+1,Npz,Ne)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,-Npy+1:2*Npy,Npz,Ne)
!    real(RP), intent(in) :: Filter1D(elem3D%Nnode_h1D,0:Npx+1)
    real(RP), intent(in) :: Filter1D(elem3D%Nnode_h1D,-Npx+1:2*Npx)

    integer :: ke, pz, py, px, p
    real(RP) :: tmp
    !-------------------------------------

    !$omp parallel do private(ke,px,py,pz,p,tmp) collapse(2)
    do ke=1, Ne
    do pz=1, Npz
      do py=1, Npy
      do px=1, Npx
        tmp = 0.0_RP
        do p=-Npx+1, Npx+Npx
          tmp = tmp + Filter1D(px,p) * q0(p,py,pz,ke)
        end do
        q(px,py,pz,ke) = tmp
      end do
      end do
    end do
    end do

    return
  end subroutine apply_filter1d_x
!OCL SERIAL
  subroutine apply_filter1d_y( q, q0, Filter1D, Npx, Npy, Npz, Ne, lmesh3D, elem3D )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    integer, intent(in) :: Npx, Npy, Npz, Ne
    real(RP), intent(out) :: q(Npx,Npy,Npz,lmesh3D%NeA)
!    real(RP), intent(in) :: q0(0:Npx+1,0:Npy+1,Npz,Ne)
    real(RP), intent(in) :: q0(-Npx+1:2*Npx,-Npy+1:2*Npy,Npz,Ne)
!    real(RP), intent(in) :: Filter1D(elem3D%Nnode_h1D,0:Npx+1)
    real(RP), intent(in) :: Filter1D(elem3D%Nnode_h1D,-Npy+1:2*Npy)

    integer :: ke, px, py, pz, p
    real(RP) :: tmp
    !-------------------------------------

    !$omp parallel do private(ke,px,py,pz,p,tmp) collapse(2)
    do ke=1, Ne
    do pz=1, Npz
      do py=1, Npy
      do px=1, Npx
        tmp = 0.0_RP
        do p=-Npy+1, Npy+Npy
          tmp = tmp + Filter1D(py,p) * q0(px,p,pz,ke)
        end do
        q(px,py,pz,ke) = tmp
      end do
      end do
    end do
    end do

    return
  end subroutine apply_filter1d_y

!OCL SERIAL
  subroutine extract_tmp3D( tmp3D, q0, q0_, lmesh3D, elem3D, vmapP, x_or_y )
    implicit none
    class(LocalMesh3D), intent(in) :: lmesh3D
    class(ElementBase3D), intent(in) :: elem3D
!    real(RP), intent(out) :: tmp3D(0:elem3D%Nnode_h1D+1,0:elem3D%Nnode_h1D+1,elem3D%Nnode_v,lmesh3D%Ne)
    real(RP), intent(out) :: tmp3D(-elem3D%Nnode_h1D+1:2*elem3D%Nnode_h1D,-elem3D%Nnode_h1D+1:2*elem3D%Nnode_h1D,elem3D%Nnode_v,lmesh3D%Ne)
    real(RP), intent(in) :: q0(elem3D%Np*lmesh3D%NeA) 
    real(RP), intent(in) :: q0_(elem3D%Nnode_h1D,elem3D%Nnode_h1D,elem3D%Nnode_v,lmesh3D%NeA) 
    integer, intent(in) :: vmapP(elem3D%NfpTot,lmesh3D%Ne) 
    integer, intent(in) :: x_or_y !< x: 0, y: 1

    integer :: ke, ke_h, ke_x, ke_y, ke_z, px, py, ph, pz, i, f
    integer :: fso, fs, fe
    integer :: iP(elem3D%Nnode_h1D)
    real(RP) :: halo_h(elem3D%Nnode_h1D,elem3D%Nnode_v,elem3D%Nnode_h1D,lmesh3D%NeX,lmesh3D%NeZ,4) 
    !------------------------------

    !$omp parallel private(ke, ke_x, ke_y, ke_h, ke_z, i, px, py, ph, pz, fso, fs, fe, iP, f)

    !$omp do collapse(3)
    do f=1, 4
      do ke_z=1, lmesh3D%NeZ
      do ke_h=1, lmesh3D%NeX
        fso = lmesh3D%Ne * elem3D%Np + (f-1)*elem3D%Nfp_h*elem3D%Nnode_h1D*lmesh3D%NeX*lmesh3D%NeZ &
            + (ke_h-1)*elem3D%Nfp_h*elem3D%Nnode_h1D + (ke_z-1)*elem3D%Nfp_h*elem3D%Nnode_h1D*lmesh3D%NeX
        do ph=1, elem3D%Nnode_h1D
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D + (ph-1)*elem3D%Nfp_h
          fe = fso + elem3D%Nnode_h1D - 1
          halo_h(:,pz,ph,ke_h,ke_z,f) = q0(fs:fe)
        end do
        end do
      end do
      end do
    end do

    !$omp do collapse(3)
    do ke_z=1, lmesh3D%NeZ
    do ke_y=1, lmesh3D%NeY
    do ke_x=1, lmesh3D%NeX
      ke = ke_x + (ke_y-1) * lmesh3D%NeX + (ke_z-1) * lmesh3D%NeX * lmesh3D%NeY

      do pz=1, elem3D%Nnode_v
      do py=1, elem3D%Nnode_h1D
      do px=1, elem3D%Nnode_h1D
        tmp3D(px,py,pz,ke) = q0_(px,py,pz,ke)
      end do
      end do
      end do

      if ( x_or_y == 1 ) then ! y-direction
        ! Face 1
        fso = 0
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D; fe = fs + elem3D%Nnode_h1D - 1
          iP(:) = vmapP(fs:fe,ke)
          if ( iP(1) <= lmesh3D%Ne * elem3D%Np ) then
            do ph=1, elem3D%Nnode_h1D
              tmp3D(1:elem3D%Nnode_h1D,-ph+1,pz,ke) = q0(iP(:) - (ph-1)*elem3D%Nnode_h1D)
            end do
          else
            do ph=1, elem3D%Nnode_h1D
              tmp3D(1:elem3D%Nnode_h1D,-ph+1,pz,ke) = halo_h(:,pz,ph,ke_x,ke_z,1)
            end do
          end if
        enddo

        ! Face 3
        fso = 2 * elem3D%Nfp_h
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D; fe = fs + elem3D%Nnode_h1D - 1
          iP(:) = vmapP(fs:fe,ke)  
          if ( iP(1) <= lmesh3D%Ne * elem3D%Np ) then
            do ph=1, elem3D%Nnode_h1D
              tmp3D(1:elem3D%Nnode_h1D,ph+elem3D%Nnode_h1D,pz,ke) = q0(iP(:) + (ph-1)*elem3D%Nnode_h1D)
            end do
          else
            do ph=1, elem3D%Nnode_h1D
              tmp3D(1:elem3D%Nnode_h1D,ph+elem3D%Nnode_h1D,pz,ke) = halo_h(:,pz,ph,ke_x,ke_z,3)
            end do
          end if          
        enddo
      else if ( x_or_y == 0 ) then ! x-direction
        ! Face 2
        fso = elem3D%Nfp_h
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D; fe = fs + elem3D%Nnode_h1D - 1
          iP(:) = vmapP(fs:fe,ke)  
          if ( iP(1) <= lmesh3D%Ne * elem3D%Np ) then
            do ph=1, elem3D%Nnode_h1D
              tmp3D(ph+elem3D%Nnode_h1D,1:elem3D%Nnode_h1D,pz,ke) = q0(iP(:) + (ph-1))
            end do
          else
            do ph=1, elem3D%Nnode_h1D
              tmp3D(ph+elem3D%Nnode_h1D,1:elem3D%Nnode_h1D,pz,ke) = halo_h(:,pz,ph,ke_y,ke_z,2)
            end do
          end if           
        enddo

        ! Face 4
        fso = 3 * elem3D%Nfp_h
        do pz=1, elem3D%Nnode_v
          fs = fso + 1 + (pz-1)*elem3D%Nnode_h1D; fe = fs + elem3D%Nnode_h1D - 1
          iP(:) = vmapP(fs:fe,ke)  
          if ( iP(1) <= lmesh3D%Ne * elem3D%Np ) then
            do ph=1, elem3D%Nnode_h1D
              tmp3D(-ph+1,1:elem3D%Nnode_h1D,pz,ke) = q0(iP(:) - (ph-1))
            end do
          else
            do ph=1, elem3D%Nnode_h1D
              tmp3D(-ph+1,1:elem3D%Nnode_h1D,pz,ke) = halo_h(:,pz,ph,ke_y,ke_z,4)
            end do
          end if           
        enddo
      end if
    
    end do
    end do
    end do

    !$omp end parallel

    return
  end subroutine extract_tmp3D

end module mod_user_sub_Filter
