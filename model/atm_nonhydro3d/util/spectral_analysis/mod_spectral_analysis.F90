#include "scalelib.h"
module mod_spectral_analysis
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
    PRC_abort
  use scale_const, only: &
    PI => CONST_PI, &
    RPlanet => CONST_RADIUS
  
  use scale_mesh_cubedom3d, only: &
    MeshCubeDom3D
  use scale_element_base, only: &
    ElementBase2D
  use scale_localmesh_2d, only: LocalMesh2D
  use scale_meshfield_base, only: &
    MeshField2D, &
    MeshField2DList

  use scale_meshfield_spectral_transform, only: &
    MeshField_SpetralTransform2D, &
    ST_EVALTYPE_SAMPLE_UNIFORM_PTS, &
    ST_EVALTYPE_L2PROJECTION_1
  implicit none
  private

  public :: spectral_analysis_Init
  public :: spectral_analysis_do

  type(MeshField_SpetralTransform2D) :: ST_tool
  integer :: varNumTot
  
contains
!OCL SERIAL
  subroutine spectral_analysis_Init( eval_type, ks, ke, ls, le, var_num, LayerNum, var2D_num, mesh3D, &
    GLQuadOrd, NsamplePerElem1D )
    use scale_element_quadrilateral, only: QuadrilateralElement
    implicit none
    character(len=*), intent(in) :: eval_type
    integer, intent(in) :: ks, ke
    integer, intent(in) :: ls, le
    integer, intent(in) :: var_num
    integer, intent(in) :: LayerNum
    integer, intent(in) :: var2D_num
    class(MeshCubeDom3D), intent(in) :: mesh3D
    integer, intent(in) :: GLQuadOrd
    integer, intent(in) :: NsamplePerElem1D
    !----------------------------------------------

    varNumTot = var_num * LayerNum + var2D_num
    select case(trim(eval_type))
    case( "Uniform_Sampling" )
      call ST_tool%Init( ST_EVALTYPE_SAMPLE_UNIFORM_PTS, ks, ke, ls, ke, mesh3D%mesh2D, varNumTot, NsamplePtPerElem1D=NsamplePerElem1D )
    case( "L2_Projection")
      call ST_tool%Init( ST_EVALTYPE_L2PROJECTION_1, ks, ke, ls, ke, mesh3D%mesh2D, varNumTot, GLQuadOrd=GLQuadOrd )
    case default
      LOG_INFO("spectral_analysis_Init",*) "Invalid eval_type is specified. Check!"
      call PRC_abort
    end select
    return
  end subroutine spectral_analysis_Init

!OCL SERIAL
  subroutine spectral_analysis_do( g_var_list, g_var2D_list, &
    var_num, var2D_num, mesh_num_x, mesh_num_y, LayerNum,    &
    s_ks, s_ke, s_ls, s_le, &
    s_var, s_var2D )
    implicit none
    integer, intent(in) :: var_num
    integer, intent(in) :: var2D_num
    integer, intent(in) :: mesh_num_x
    integer, intent(in) :: mesh_num_y
    integer, intent(in) :: LayerNum
    integer, intent(in) :: s_ks, s_ke
    integer, intent(in) :: s_ls, s_le
    type(MeshField2D), intent(in), target :: g_var_list(LayerNum,var_num,mesh_num_x,mesh_num_y)
    type(MeshField2D), intent(in), target :: g_var2D_list(var2D_num,mesh_num_x,mesh_num_y)
    real(RP), intent(out) :: s_var(s_ks:s_ke,s_ls:s_le,2,LayerNum,var_num)
    real(RP), intent(out) :: s_var2D(s_ks:s_ke,s_ls:s_le,2,var2D_num)

    type(MeshField2DList) :: var_list(varNumTot,mesh_num_x,mesh_num_y)
    integer :: v, k
    integer :: vv
    integer :: mx, my

    class(LocalMesh2D), pointer :: lmesh2D
    !----------------------------------------------

    do my=1, mesh_num_y
    do mx=1, mesh_num_x
      do v=1, var_num
      do k=1, LayerNum
        vv = k + (v-1)*LayerNum
        var_list(vv,mx,my)%ptr => g_var_list(k,v,mx,my)
      end do
      end do
      do v=1, var2D_num
        vv = v + LayerNum * var_num
        var_list(vv,mx,my)%ptr => g_var2D_list(v,mx,my)
      end do
    end do
    end do

    ! do my=1, mesh_num_y
    ! do mx=1, mesh_num_x
    !   lmesh2D => g_var_list(1,1,mx,1)%mesh%lcmesh_list(1)
    !   LOG_INFO("CHECK",*) "mx=", mx, lmesh2D%xmin, lmesh2D%xmax, lmesh2D%ymin, lmesh2D%ymax 
    !   LOG_INFO("CHECK",*) "mx,my=", mx,my, "U=", g_var_list(1,1,mx,my)%local(1)%val(:,lmesh2D%NeS:lmesh2D%NeE)
    ! end do
    ! end do
    call ST_tool%Transform( var_list, mesh_num_x, mesh_num_y )

    ! do v=s_ks,s_ke
    !   LOG_INFO("CHECK",*) "spectra k=", v, ST_tool%spectral_coef(v,0,:,1)
    ! end do

    do v=1, var_num
    do k=1, LayerNum
      vv = k + (v-1)*LayerNum
      s_var(:,:,1,k,v) = ST_tool%spectral_coef(:,:,1,vv)
      s_var(:,:,2,k,v) = ST_tool%spectral_coef(:,:,2,vv)
    end do
    end do
    do v=1, var2D_num
      vv = v + LayerNum * var_num
      s_var2D(:,:,1,v) = ST_tool%spectral_coef(:,:,1,vv)
      s_var2D(:,:,2,v) = ST_tool%spectral_coef(:,:,2,vv)
    end do
    return
  end subroutine spectral_analysis_do
end module mod_spectral_analysis