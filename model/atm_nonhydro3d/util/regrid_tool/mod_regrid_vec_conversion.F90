!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_regrid_vec_conversion
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_file_h
  use scale_prof
  use scale_prc, only: &
    PRC_myrank, PRC_abort

  use scale_mesh_base2d, only: &
    MeshBase2D
  use scale_mesh_base3d, only: &
    MeshBase3D
  use scale_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_mesh_cubedspheredom2d, only: MeshCubedSphereDom2D
  use scale_mesh_rectdom2d, only: MeshRectDom2D
  use scale_mesh_cubedspheredom3d, only: MeshCubedSphereDom3D
  use scale_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_file_base_meshfield, only: FILE_base_meshfield

  use mod_regrid_mesh_base, only: &
    regrid_mesh_base
  use mod_regrid_nodemap, only: &
    regrid_nodemap
  
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedures
  !
  public :: regrid_vec_conversion_Init
  public :: regrid_vec_conversion_Final
  public :: regrid_vec_conversion_Do
  
  type(MeshField3D), public :: out_veclon
  type(MeshField3D), public :: out_veclat

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  character(len=H_LONG) :: out_basename      = ''       ! Basename of the input  file 

  type(MeshField3D) :: out_vec_comp1
  type(MeshField3D) :: out_vec_comp2


contains
!OCL SERIAL
  subroutine regrid_vec_conversion_Init( out_mesh )
    implicit none
    class(MeshBase3D), intent(in) :: out_mesh
        
    integer :: ierr
    !----------------------------------------------

    call out_vec_comp1%Init( "vec_comp1", "", out_mesh )
    call out_vec_comp2%Init( "vec_comp2", "", out_mesh )
    call out_veclon%Init( "veclon_comp1", "", out_mesh )
    call out_veclat%Init( "veclat_comp2", "", out_mesh )

    return
  end subroutine regrid_vec_conversion_Init

!OCL SERIAL
  subroutine regrid_vec_conversion_Do( &
    istep, out_mesh, nodeMap_list,     &
    GP_flag, out_mesh_GP, GPMat        )

    use scale_meshutil_cubedsphere2d, only: &
      MeshUtilCubedSphere2D_getPanelID
    use mod_regrid_interp_field, only: &
      regrid_interp_field_Interpolate
    
    use scale_localmesh_2d, only: LocalMesh2D
    use scale_localmesh_3d, only: LocalMesh3D
    use scale_element_base, only: ElementBase3D
    use scale_cubedsphere_coord_cnv, only: &
      CubedSphereCoordCnv_LonLat2CSPos
    
    use scale_const, only: &
      PI => CONST_PI,         &
      RPlanet => CONST_radius
    implicit none
    integer, intent(in) :: istep
    class(regrid_mesh_base), intent(in), target :: out_mesh
    type(regrid_nodemap), intent(in) :: nodeMap_list(:)
    logical, intent(in) :: GP_flag
    class(regrid_mesh_base), intent(in), target :: out_mesh_GP
    real(RP), intent(in) :: GPMat(:,:)

    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: n
    integer, allocatable :: inPanelID(:,:)
    integer, allocatable :: inPanelID3D(:,:)
    integer :: Np1D, Ne2D

    real(RP) :: X, Y, del2
    integer :: ke, ke2D
    class(elementbase3D), pointer :: elem3D

    real(RP), allocatable :: out_x(:,:)
    real(RP), allocatable :: out_y(:,:)
    real(RP), allocatable :: out_x3D(:,:)
    real(RP), allocatable :: out_y3D(:,:)

    integer :: p_h, ke_h
    real(RP) :: out_lon(1)
    real(RP) :: out_lat(1)    
    !----------------------------------------------


    call regrid_interp_field_Interpolate( out_vec_comp1, &
      istep, "U", out_mesh, nodeMap_list,                &
      GP_flag, out_mesh_GP, GPMat                        )
    call regrid_interp_field_Interpolate( out_vec_comp2, &
      istep, "V", out_mesh, nodeMap_list,                &
      GP_flag, out_mesh_GP, GPMat                        )

    do n=1, out_mesh%mesh3D%LOCAL_MESH_NUM
      lcmesh => out_mesh%mesh3D%lcmesh_list(n)
      elem3D => lcmesh%refElem3D

      lcmesh2D => lcmesh%lcmesh2D
      Np1D = lcmesh2D%refElem2D%Nfp
      Ne2D = lcmesh%NeX * lcmesh%NeY

      allocate( inPanelID(Np1D**2,Ne2D) )
      allocate( inPanelID3D(elem3D%Np,lcmesh%Ne) )
      allocate( out_x(Np1D**2,Ne2D), out_y(Np1D**2,Ne2D) )
      allocate( out_x3D(elem3D%Np,lcmesh%Ne) )
      allocate( out_y3D(elem3D%Np,lcmesh%Ne) )

      call MeshUtilCubedSphere2D_getPanelID ( &
        inPanelID(:,:),                             & ! (out)
        lcmesh2D%pos_en(:,:,1) * PI / 180.0_RP,       & ! (in)
        lcmesh2D%pos_en(:,:,2) * PI / 180.0_RP,       & ! (in)
        Np1D**2 * Ne2D                              ) ! (in)      

      do ke_h=1, Ne2D
      do p_h=1, Np1D**2
        out_lon(1) = lcmesh2D%pos_en(p_h,ke_h,1) * PI / 180.0_RP
        out_lat(1) = lcmesh2D%pos_en(p_h,ke_h,2) * PI / 180.0_RP

        call CubedSphereCoordCnv_LonLat2CSPos( &
          inPanelID(p_h,ke_h),                       & ! (in)
          out_lon(1), out_lat(1), 1,                 & ! (in)
          out_x(p_h,ke_h), out_y(p_h,ke_h)           ) ! (out)
      end do
      end do
    
      !$omp parallel do private(ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)
        inPanelID3D(:,ke) = inPanelID(elem3D%IndexH2Dto3D(:),ke2D)
        out_x3D(:,ke) = out_x(elem3D%IndexH2Dto3D(:),ke2D)
        out_y3D(:,ke) = out_y(elem3D%IndexH2Dto3D(:),ke2D)
      end do

      call CS2LonLatVec( inPanelID3D, out_x3D, out_y3D, lcmesh%gam, &
        elem3D%Np, lcmesh%Ne, lcmesh%NeA,               &
        out_vec_comp1%local(n)%val(:,:),                &
        out_vec_comp2%local(n)%val(:,:),                &
        out_veclon%local(n)%val(:,:),             &
        out_veclat%local(n)%val(:,:)              )
      
      !$omp parallel do private(ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)
        out_veclon%local(n)%val(:,ke) = &
          out_veclon%local(n)%val(:,ke) * cos(lcmesh%pos_en(:,ke,2) * PI / 180.0_RP )
      end do
        
      deallocate( inPanelID, inPanelID3D )
      deallocate( out_x, out_y )
      deallocate( out_x3D, out_y3D )
    end do

    return
  end subroutine regrid_vec_conversion_Do

!OCL SERIAL
  subroutine regrid_vec_conversion_Final()
    implicit none
    !----------------------------------------------

    call out_vec_comp1%Final()
    call out_vec_comp2%Final()

    return
  end subroutine regrid_vec_conversion_Final

!-----------------------------

!OCL SERIAL
  subroutine CS2LonLatVec( &
    panelID, alpha, beta, gam, Np, Ne, NeA,    & ! (in)
    VecAlpha, VecBeta,                         & ! (in)
    VecLon, VecLat                             ) ! (out)

    use scale_const, only: &
      EPS => CONST_EPS,       &
      RPlanet => CONST_RADIUS
    implicit none

    integer, intent(in) :: Np
    integer, intent(in) :: Ne
    integer, intent(in) :: NeA
    integer, intent(in) :: panelID(Np,Ne)
    real(RP), intent(in) :: alpha(Np,Ne)
    real(RP), intent(in) :: beta (Np,Ne)
    real(RP), intent(in) :: gam(Np,NeA)
    real(DP), intent(in) :: VecAlpha(Np,NeA)
    real(DP), intent(in) :: VecBeta (Np,NeA)
    real(RP), intent(out) :: VecLon(Np,NeA)
    real(RP), intent(out) :: VecLat(Np,NeA)

    integer :: p, ke
    real(RP) :: X ,Y, del2
    real(RP) :: s
    real(RP) :: radius
    !-----------------------------------------------------------------------------

    !$omp parallel do collapse(2) private( p, X, Y, del2, s, radius )
    do ke=1, Ne
    do p=1, Np
      X = tan( alpha(p,ke) )
      Y = tan( beta (p,ke) )
      del2 = 1.0_RP + X**2 + Y**2
      radius = RPlanet * gam(p,ke)

      select case( panelID(p,ke) )
      case( 1, 2, 3, 4 )
        VecLon(p,ke) = VecAlpha(p,ke) * radius
        VecLat(p,ke) = ( - X * Y * VecAlpha(p,ke) + ( 1.0_RP + Y**2 ) * VecBeta(p,ke) ) &
                  * radius * sqrt( 1.0_RP + X**2 ) / del2 
      case ( 5 )
        s = radius
      case ( 6 )
        s = - radius
      end select
      select case( panelID(p,ke) )
      case( 5, 6 )
        VecLon(p,ke) = (- Y * ( 1.0 + X**2 ) * VecAlpha(p,ke) + X * ( 1.0_RP + Y**2 ) * VecBeta(p,ke) ) &
                     * s / ( X**2 + Y**2 + EPS )
        VecLat(p,ke) = (- X * ( 1.0 + X**2 ) * VecAlpha(p,ke) - Y * ( 1.0_RP + Y**2 ) * VecBeta(p,ke) ) &
                     * s / ( del2 * sqrt( X**2 + Y**2 ) + EPS )
      end select
    end do
    end do

    return
  end subroutine CS2LonLatVec

end module mod_regrid_vec_conversion