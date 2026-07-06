!-------------------------------------------------------------------------------
!> module Vertical interpolation for spectral analysis
!!
!! @author Yuta Kawai, Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_spectral_analysis_vintrp
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: PRC_abort

  use scale_localmesh_3d, only: LocalMesh3D
  use scale_element_base, only: ElementBase3D
  use scale_meshfield_base, only: MeshField2D, MeshField3D
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !

  !> Derived type to manage vertical interpolation for spectral analysis
  type, public :: SpectralAnalysisVIntrp
    integer               :: ZintrpType_ID      !< ID of vertical interpolation type
    real(RP), allocatable :: IntrpMat(:,:,:)     !< Vertical interpolation matrix for each target level
    real(RP), allocatable :: coef1(:)            !< Interpolation coefficient for the lower node, for each target level
    integer,  allocatable :: keZ1(:), keZ2(:)    !< Vertical element indices of the lower/upper nodes used for interpolation, for each target level
  contains
    procedure :: Init          => VIntrp_Init
    procedure :: Final         => VIntrp_Final
    procedure :: GetLevelValue => VIntrp_GetLevelValue
  end type SpectralAnalysisVIntrp

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: ZintrpType_Nearest_ID                 = 1  !< ID of vertical interpolation type: Nearest
  integer, parameter :: ZintrpType_SamplingUniPt_LinInterp_ID = 2  !< ID of vertical interpolation type: Sampling at uniform points with linear interpolation

contains
  !> Setup vertical interpolation strategy and precompute per-level interpolation tables
!OCL SERIAL
  subroutine VIntrp_Init( this, zintrp_type_name, levels, lmesh3D )
    implicit none
    class(SpectralAnalysisVIntrp), intent(inout) :: this
    character(len=*), intent(in) :: zintrp_type_name  !< Name of vertical interpolation type
    real(RP), intent(in) :: levels(:)                  !< Array of vertical levels
    class(LocalMesh3D), intent(in) :: lmesh3D

    integer :: k, LevelNum
    !----------------------------------------------

    select case( trim(zintrp_type_name) )
    case( "Nearest" )
      this%ZintrpType_ID = ZintrpType_Nearest_ID
    case( "SamplingUniPt_LinearIntrp" )
      this%ZintrpType_ID = ZintrpType_SamplingUniPt_LinInterp_ID
    case default
      LOG_INFO("SpectralAnalysisVIntrp_Init",*) "The specified zintrp_type_name is invalid. Check!", trim(zintrp_type_name)
      call PRC_abort
    end select

    LevelNum = size(levels)
    allocate( this%IntrpMat(lmesh3D%refElem3D%Nnode_v,2,LevelNum) )
    allocate( this%keZ1(LevelNum), this%keZ2(LevelNum) )
    allocate( this%coef1(LevelNum) )

    do k=1, LevelNum
      call prep_VIntrp( this%IntrpMat(:,:,k), this%keZ1(k), this%keZ2(k), this%coef1(k), &
        levels(k), lmesh3D, lmesh3D%refElem3D )
    end do

    return
  end subroutine VIntrp_Init

  !> Finalize vertical interpolation strategy
!OCL SERIAL
  subroutine VIntrp_Final( this )
    implicit none
    class(SpectralAnalysisVIntrp), intent(inout) :: this
    !----------------------------

    if ( allocated(this%IntrpMat) ) then
      deallocate( this%IntrpMat, this%keZ1, this%keZ2, this%coef1 )
    end if
    return
  end subroutine VIntrp_Final

  !> Interpolate a 3D field onto the k-th target level, writing the result into a 2D slice
!OCL SERIAL
  subroutine VIntrp_GetLevelValue( this, k, target_lev, lcmesh3D, elem3D, tmp_field3D, gvar_out )
    implicit none
    class(SpectralAnalysisVIntrp), intent(in) :: this
    integer, intent(in) :: k                              !< Target level index
    real(RP), intent(in) :: target_lev                    !< Target level
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    type(MeshField3D), intent(in) :: tmp_field3D          !< Source 3D field
    type(MeshField2D), intent(inout) :: gvar_out          !< Destination 2D slice at the target level

    integer :: ke2D, ke, ke2
    integer :: ke_z1, ke_z2
    integer :: p, pp, ph, pz
    real(RP), allocatable :: tmp3D_zxy(:,:,:)
    real(RP) :: coef1, coef2
    !---------------------------------------------------

    if ( this%ZintrpType_ID == ZintrpType_Nearest_ID ) then
      do ke = lcmesh3D%NeS, lcmesh3D%NeE
          ke2D = lcmesh3D%EMap3Dto2D(ke)
          do pz=1, elem3D%Nnode_v
          do ph=1, elem3D%Nnode_h1D**2
            p = ph + (pz-1)*elem3D%Nnode_h1D**2
            pp = ph + min(pz,elem3D%Nnode_v-1)*elem3D%Nnode_h1D**2
            if ( lcmesh3D%pos_en(p,ke,3) <= target_lev   &
              .and. lcmesh3D%pos_en(pp,ke,3) >= target_lev ) then
              gvar_out%local(1)%val(ph,ke2D) = tmp_field3D%local(1)%val(p,ke)
            end if
          end do
          end do
        end do
    else if ( this%ZintrpType_ID == ZintrpType_SamplingUniPt_LinInterp_ID ) then
      allocate( tmp3D_zxy(elem3D%Nnode_v,elem3D%Nnode_h1D**2,2) )

      coef1 = this%coef1(k); coef2 = 1.0_RP - coef1
      !$omp parallel do private(ke2D, ke_z1, ke_z2, ke, ke2, ph, tmp3D_zxy)
      do ke2D=1, lcmesh3D%Ne2D
        ke_z1 = this%keZ1(k); ke_z2 = this%keZ2(k)
        ke = ke2D + (ke_z1-1)*lcmesh3D%Ne2D
        ke2 = ke2D + (ke_z2-1)*lcmesh3D%Ne2D
        do ph=1, elem3D%Nnode_h1D**2
          tmp3D_zxy(:,ph,1) = tmp_field3D%local(1)%val(elem3D%Colmask(:,ph),ke)
          tmp3D_zxy(:,ph,2) = tmp_field3D%local(1)%val(elem3D%Colmask(:,ph),ke2)
        end do
        do ph=1, elem3D%Nnode_h1D**2
          gvar_out%local(1)%val(ph,ke2D) = &
              coef1 * sum( this%IntrpMat(:,1,k) * tmp3D_zxy(:,ph,1) ) &
            + coef2 * sum( this%IntrpMat(:,2,k) * tmp3D_zxy(:,ph,2) )
        end do
      end do

      deallocate( tmp3D_zxy )
    end if

    return
  end subroutine VIntrp_GetLevelValue

!-- private subroutines ---------------------------------------------------

!OCL SERIAL
  subroutine prep_VIntrp( IntrpMat, VIntrp_k1, VIntrp_k2, VIntrp_coef1, target_lev, lmesh3D, elem3D )
    use scale_polynomial, only: &
      Polynomial_GenLagrangePoly
    implicit none
    class(ElementBase3D), intent(in) :: elem3D
    class(LocalMesh3D), intent(in) :: lmesh3D
    real(RP), intent(out) :: IntrpMat(elem3D%Nnode_v,2)
    integer, intent(out) :: VIntrp_k1, VIntrp_k2
    real(RP), intent(out) :: VIntrp_coef1
    real(RP), intent(in) :: target_lev

    integer :: ke, ke_z, ke_z2
    integer :: pz, pz1, pz2
    integer :: k, kk
    real(RP) :: vz(8)
    real(RP) :: delz
    real(RP) :: lev_uniform(elem3D%Nnode_v*lmesh3D%NeZ)
    real(RP) :: lev_uniform_xi(elem3D%Nnode_v*lmesh3D%NeZ)

    real(RP) :: lag_poly(2,elem3D%Nnode_v)
    !---------------------------------

    !$omp parallel do private(ke,pz,k,vz,delz)
    do ke_z=1, lmesh3D%NeZ
      ke = 1 + (ke_z-1)*lmesh3D%Ne2D
      vz(:) = lmesh3D%pos_ev(lmesh3D%EToV(ke,:),3)
      delz = ( vz(5) - vz(1) ) / real(elem3D%Nnode_v, kind=RP)
      do pz=1, elem3D%Nnode_v
        k = pz + (ke_z-1)*elem3D%Nnode_v
        lev_uniform(k) = vz(1) + (pz - 0.5_RP)*delz
        lev_uniform_xi(k) = - 1.0_RP + 2.0_RP * ( lev_uniform(k) - vz(1) ) / ( vz(5) - vz(1) )
      end do
    end do
    VIntrp_k1 = -1
    VIntrp_k2 = -1
    do ke_z=1, lmesh3D%NeZ
    do pz=1, elem3D%Nnode_v
      k = pz+(ke_z-1)*elem3D%Nnode_v
      kk = min(k+1, size(lev_uniform))
      if ( lev_uniform(k) <= target_lev .and. target_lev <= lev_uniform(kk) ) then
        if ( pz+1 > elem3D%Nnode_v ) then
          ke_z2 = ke_z + 1; pz2 = 1
        else
          ke_z2 = ke_z; pz2 = pz + 1
        end if
        VIntrp_k1 = ke_z; VIntrp_k2 = ke_z2; pz1 = pz
        exit
      end if
    end do
    end do

    if ( VIntrp_k1 < 0 .or. VIntrp_k2 < 0 ) then
      LOG_INFO("prep_VIntrp",*) "The specified level is invalid. Check!", target_lev
      call PRC_abort
    end if

    !-
    k = pz1 + (VIntrp_k1-1)*elem3D%Nnode_v
    kk = pz2 + (VIntrp_k2-1)*elem3D%Nnode_v

    lag_poly(:,:) = Polynomial_GenLagrangePoly( elem3D%PolyOrder_v, elem3D%x3(elem3D%Colmask(:,1)), (/ lev_uniform_xi(k), lev_uniform_xi(kk) /) )
    IntrpMat(:,:) = transpose(lag_poly)
    VIntrp_coef1 = ( lev_uniform(kk) - target_lev ) / ( lev_uniform(kk) - lev_uniform(k) )

    LOG_INFO("prep_VIntrp",*) "k1,k2=", VIntrp_k1, VIntrp_k2, ": pz1, pz2=", pz1, pz2
    LOG_INFO("prep_VIntrp",*) "lev=", lev_uniform(k), lev_uniform(kk)
    LOG_INFO("prep_VIntrp",*) "lev_xi=", lev_uniform_xi(k), lev_uniform_xi(kk)
    LOG_INFO("prep_VIntrp",*) "VIntrp_coef1=", VIntrp_coef1

    return
  end subroutine prep_VIntrp

end module mod_spectral_analysis_vintrp
