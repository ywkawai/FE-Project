#include "scalelib.h"
program test_advect2d
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !  
  use scale_const, only: &
    PI  => CONST_PI
  use scale_precision
  use scale_prc
  use scale_io
  use scale_prof

  use scale_sparsemat  
  use scale_element_base
  use scale_element_quadrial
  use scale_localmesh_2d
  use scale_mesh_rectdom2d

  use scale_localmeshfield_base, only: LocalMeshField2D
  use scale_meshfield_base, only: MeshField2D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_rectdom2d, only: MeshFieldCommRectDom2D

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_put,   &
    FILE_HISTORY_meshfield_write

  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

  use scale_time_manager, only: &
    TIME_manager_advance,                              &
    TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP,            &
    TIME_DTSEC, TIME_NSTEP 

  !-----------------------------------------------------------------------------
  implicit none

  integer :: NeGX
  integer :: NeGY
  integer, parameter :: NLocalMeshPerPrc = 1

  ! sin, cosbell, hat
  character(len=H_SHORT) :: InitShapeName = 'cosbell'

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: dom_ymin = -1.0_RP
  real(RP), parameter :: dom_ymax = +1.0_RP
  
  real(RP), parameter :: ADV_VELX = sqrt(1.0_RP)
  real(RP), parameter :: ADV_VELY = sqrt(1.0_RP)

  type(QuadrialElement) :: refElem
  integer :: PolyOrder
  logical, parameter :: DumpedMassMatFlag = .false.
  integer, parameter :: PolyOrderErrorCheck = 6
  type(sparsemat) :: Dx, Sx, Dy, Sy, Lift
  
  type(MeshRectDom2D), target :: mesh
  type(MeshField2D), target :: q, q0, qexact  
  type(MeshField2D), target :: u, v
  type(MeshFieldCommRectDom2D) :: fields_comm
  type(MeshFieldContainer) :: field_list(3)  
  integer :: HST_ID(2)

  integer :: n, k, p
  type(LocalMesh2D), pointer :: lcmesh
  
  integer :: nowstep
  integer :: rkstage
  integer, parameter :: nrkstage = 3
  real(RP) :: tsec_
  real(RP) :: rkcoef1(nrkstage) = (/ 0.0_RP, 3.0_RP/4.0_RP, 1.0_RP/3.0_RP /)
  real(RP) :: rkcoef2(nrkstage) = (/ 1.0_RP, 1.0_RP/4.0_RP, 2.0_RP/3.0_RP /)

  real(RP), allocatable :: IntrpMat(:,:)
  real(RP) :: intw_intrp(PolyOrderErrorCheck**2)
  real(RP) :: x_intrp(PolyOrderErrorCheck**2)
  real(RP) :: y_intrp(PolyOrderErrorCheck**2)

  integer :: nstep_eval_error
  !-------------------------------------------------------

  call init()
  call set_initcond()

  field_list(1)%field2d => q
  field_list(2)%field2d => u
  field_list(3)%field2d => v

  do nowstep=1, TIME_NSTEP
  
    do n=1, mesh%LOCAL_MESH_NUM
      q0%local(n)%val(:,:) = q%local(n)%val(:,:)
    end do

    do rkstage=1, nrkstage
      !* Exchange halo data
      call PROF_rapstart( 'exchange_halo', 1)

      call fields_comm%Put(field_list, 1)
      call fields_comm%Exchange()
      call fields_comm%Get(field_list, 1)

      call PROF_rapend( 'exchange_halo', 1)

      !* Update prognostic variables
      call PROF_rapstart( 'update_dyn', 1)
      do n=1, mesh%LOCAL_MESH_NUM
        call update_dyn( &
          TIME_DTSEC, rkcoef1(rkstage), rkcoef2(rkstage), &
          mesh%lcmesh_list(n), refElem,                   &
          q%local(n)%val, q0%local(n)%val,                &
          u%local(n)%val, v%local(n)%val )        
      end do
      call PROF_rapend( 'update_dyn', 1)

    end do
    
    !* Advance time
    call TIME_manager_advance()

    tsec_ = TIME_NOWDATE(6) + TIME_NOWMS
    if (mod(nowstep,nstep_eval_error) == 0) then 
      LOG_PROGRESS('(A,F13.5,A)') "t=", real(tsec_), "[s]"
      call evaluate_error(tsec_)
    end if
    call FILE_HISTORY_set_nowdate( TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP )

    !* Output
    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_write()
  end do

  call final()

contains
  subroutine update_dyn( dt, rkcoef_1, rkcoef_2, lmesh, elem, q_, q0_, u_, v_)
    implicit none

    real(RP), intent(in) :: dt
    real(RP), intent(in) :: rkcoef_1, rkcoef_2
    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem
    real(RP), intent(inout) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)    :: q0_(elem%Np,lmesh%NeA)
    real(RP), intent(in)    :: u_(elem%Np,lmesh%NeA)
    real(RP), intent(in)    :: v_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), Fy(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: dqdt(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)
    !------------------------------------------------------------------------

    call PROF_rapstart( 'update_dyn_cal_bndflux', 2)
    call cal_del_flux_dyn( del_flux,                              & ! (out)
      q_, u_, v_, lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                   & ! (in)
      lmesh, elem )                                                 ! (in)
    
    call PROF_rapend( 'update_dyn_cal_bndflux', 2)

    !-----
    call PROF_rapstart( 'update_dyn_cal_interior', 2)
    
    do k = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,k)*u_(:,k), Fx)
      call sparsemat_matmul(Dy, q_(:,k)*v_(:,k), Fy)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k), LiftDelFlx)
      
      dqdt(:) = - (  lmesh%Escale(:,k,1,1) * Fx(:) &
                   + lmesh%Escale(:,k,2,2) * Fy(:) &
                   + LiftDelFlx )
      q_(:,k) = rkcoef_1*q0_(:,k) + rkcoef_2*(q_(:,k) + dt * dqdt)
    end do

    call PROF_rapend( 'update_dyn_cal_interior', 2)

    return
  end subroutine update_dyn

  subroutine cal_del_flux_dyn( del_flux, q_, u_, v_, nx, ny, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh2D), intent(in) :: lmesh
    class(elementbase2D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  v_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: VelP, VelM, alpha
    !------------------------------------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      VelM = u_(iM)*nx(i) + v_(iM)*ny(i)
      VelP = u_(iP)*nx(i) + v_(iP)*ny(i)

      alpha = 0.5_RP*abs(VelM + VelP)
      del_flux(i) = 0.5_RP*(               &
          ( q_(iP)*VelP - q_(iM)*VelM )    &
        - alpha*(q_(iP) - q_(iM))        )
    end do

    return
  end subroutine cal_del_flux_dyn

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine evaluate_error(tsec)

    implicit none

    real(DP), intent(in) :: tsec

    real(RP) :: l2error
    real(RP) :: linferror
    real(RP) :: q_intrp(PolyOrderErrorCheck**2)
    real(RP) :: qexact_intrp(PolyOrderErrorCheck**2)
    real(RP) :: x_uwind(refElem%Np), y_vwind(refElem%Np)
    real(RP) :: x_uwind_intrp(PolyOrderErrorCheck**2), y_vwind_intrp(PolyOrderErrorCheck**2)
    real(RP) :: pos_intrp(PolyOrderErrorCheck**2,2)
    real(RP) vx(4), vy(4)

    !------------------------------------------------------------------------

    linferror = 0.0_RP
    l2error = 0.0_RP
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE

        x_uwind(:) = get_upwind_pos(lcmesh%pos_en(:,k,1), ADV_VELX, tsec, dom_xmin, dom_xmax)
        y_vwind(:) = get_upwind_pos(lcmesh%pos_en(:,k,2), ADV_VELY, tsec, dom_ymin, dom_ymax)

        vx(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),1)
        vy(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),2)
        pos_intrp(:,1) = vx(1) + 0.5_RP*(x_intrp(:) + 1.0_RP)*(vx(2) - vx(1))
        pos_intrp(:,2) = vy(1) + 0.5_RP*(y_intrp(:) + 1.0_RP)*(vy(3) - vy(1))
        x_uwind_intrp(:) = get_upwind_pos(pos_intrp(:,1), ADV_VELX, tsec, dom_xmin, dom_xmax)
        y_vwind_intrp(:) = get_upwind_pos(pos_intrp(:,2), ADV_VELY, tsec, dom_ymin, dom_ymax)

        qexact%local(n)%val(:,k) = get_profile(InitShapeName, x_uwind, y_vwind)
        qexact_intrp(:) = get_profile(InitShapeName, x_uwind_intrp, y_vwind_intrp)
        q_intrp(:) = matmul(IntrpMat, q%local(n)%val(:,k))

        l2error = l2error &
          + sum(   lcmesh%J(1,k) * intw_intrp(:)          &
                * ( q_intrp(:) - qexact_intrp(:) )**2 )
        
        linferror = max(linferror, maxval(abs(q%local(n)%val(:,k) - qexact%local(n)%val(:,k))))
      end do
    end do

    LOG_INFO("evaluate_error_l2",*), sqrt(l2error)/( (dom_xmax - dom_xmin) * (dom_ymax - dom_ymin) )
    LOG_INFO("evaluate_error_linf",*) linferror

  end subroutine evaluate_error

  function get_upwind_pos(pos, ADV_VEL, nowtime, dom_min, dom_max) result(upos)
    real(RP), intent(in) :: pos(:)
    real(RP), intent(in) :: ADV_VEL
    real(RP), intent(in) :: nowtime
    real(RP), intent(in) :: dom_min, dom_max
    real(RP) :: upos(size(pos))

    integer :: period
    !-------

    period = ADV_VEL*nowtime/(dom_max - dom_min)
    
    upos(:) = pos(:) - (ADV_VEL*nowtime - dble(period)*(dom_max - dom_min))
    where (upos < dom_min)
      upos = dom_max + (upos - dom_min)
    end where
  end function get_upwind_pos


  subroutine set_initcond()
    implicit none
    !------------------------------------------------------------------------

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE
        q%local(n)%val(:,k) = get_profile(InitShapeName, lcmesh%pos_en(:,k,1), lcmesh%pos_en(:,k,2))
        qexact%local(n)%val(:,k) = q%local(n)%val(:,k)
        u%local(n)%val(:,k) = ADV_VELX
        v%local(n)%val(:,k) = ADV_VELY
      end do
    end do

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_write()   
  
    LOG_PROGRESS('(A,F13.5,A)') "t=", real(0.0_RP), "[s]"
    call evaluate_error(0.0_RP)
    return
  end subroutine set_initcond

  function get_profile(profile_name, x, y) result(profile)
    implicit none

    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(:)
    real(RP), intent(in) :: y(:)
    real(RP) :: profile(size(x))

    real(RP) :: half_width = 0.15_RP
    real(RP) :: dist(size(x))
    !------------------------------------------------------------------------

    profile(:) = 0.0_RP
    dist(:) = sqrt(x(:)**2 + y(:)**2)

    select case(InitShapeName)
    case ('sin')
      profile(:) = sin( PI*x(:) )
    case ('cosbell')
      where( dist <= half_width )
        profile(:) = (1.0_RP + cos(PI*dist(:)/half_width))*0.5_RP
      end where
    case ('hat')
      where( dist <= half_width )
        profile(:) = 1.0_RP
      end where
    end select

    return
  end function get_profile

  subroutine init()

    use scale_calendar, only: CALENDAR_setup
    use scale_time_manager, only: TIME_manager_Init 
    use scale_file_history_meshfield, only: FILE_HISTORY_meshfield_setup  
    use scale_file_history, only: FILE_HISTORY_reg 
    
    use scale_polynominal, only: &
      Polynominal_genLegendrePoly,    &
      Polynominal_GenGaussLegendrePt, &
      Polynominal_GenGaussLegendrePtIntWeight
        
    implicit none

    namelist /PARAM_TEST/ &
      NeGX, NeGY, PolyOrder, InitShapeName, &
      nstep_eval_error
    
    integer :: comm, myrank, nprocs
    logical :: ismaster

    real(RP) :: r_int1D_i(PolyOrderErrorCheck)
    real(RP) :: r_int1Dw_i(PolyOrderErrorCheck)
    real(RP), allocatable :: P_int1D_ori(:,:)
    real(RP), allocatable :: Vint(:,:)
    
    integer :: p1, p2, p1_, p2_
    integer :: n_, l_

    integer :: ierr
    !------------------------------------------------------------------------

    call PRC_MPIstart( comm )
     
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call IO_setup( "test", allow_noconf = .true. )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
  
    !--- read namelist

    NeGX = 2; NeGY = 2; PolyOrder = 1 
    InitShapeName    = 'sin'
    nstep_eval_error = 5

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TEST,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("init",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("init",*) 'Not appropriate names in namelist PARAM_TEST. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_TEST)

    ! setup profiler
    call PROF_setup
    call PROF_rapstart( "total", 0 )
    call PROF_rapstart( "init", 1 )

    ! setup calendar & initial time
    call CALENDAR_setup
    call TIME_manager_Init

    !------   
    
    call refElem%Init(PolyOrder, DumpedMassMatFlag)
    call Dx%Init(refElem%Dx1)
    call Sx%Init(refElem%Sx1)
    call Dy%Init(refElem%Dx2)
    call Sy%Init(refElem%Sx2)
    call Lift%Init(refElem%Lift)
  
    call mesh%Init( &
      NeGX, NeGY,                             &
      dom_xmin, dom_xmax, dom_ymin, dom_ymax, &
      .true., .true.,                         &
      refElem, NLocalMeshPerPrc )

    call mesh%Generate()

    !---
    call q%Init( "q", "1", mesh )
    call q0%Init( "q0", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call u%Init( "u", "m/s", mesh )
    call v%Init( "v", "m/s", mesh )
    call fields_comm%Init(3, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( mesh2d_=mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='XY')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='XY')
    !---

    allocate( P_int1D_ori(PolyOrderErrorCheck,Polyorder+1) )
    allocate( Vint(PolyOrderErrorCheck**2,(PolyOrder+1)**2) )

    r_int1D_i(:) = Polynominal_GenGaussLegendrePt( PolyOrderErrorCheck )
    r_int1Dw_i(:) = Polynominal_GenGaussLegendrePtIntWeight( PolyOrderErrorCheck )
    P_int1D_ori(:,:) = Polynominal_GenLegendrePoly(refElem%PolyOrder, r_int1D_i)

    do p2_=1, PolyOrderErrorCheck
    do p1_=1, PolyOrderErrorCheck
      n_= p1_ + (p2_-1)*PolyOrderErrorCheck
      x_intrp(n_) = r_int1D_i(p1_)
      y_intrp(n_) = r_int1D_i(p2_)
      intw_intrp(n_) = r_int1Dw_i(p1_) * r_int1Dw_i(p2_)

      do p2=1, refElem%Nfp
      do p1=1, refElem%Nfp
        l_ = p1 + (p2-1)*refElem%Nfp
        Vint(n_,l_) =  P_int1D_ori(p1_,p1) * sqrt(real(p1-1,kind=RP) + 0.5_RP) &
                     * P_int1D_ori(p2_,p2) * sqrt(real(p2-1,kind=RP) + 0.5_RP)
      end do
      end do
    end do
    end do

    allocate( IntrpMat(PolyOrderErrorCheck**2,(PolyOrder+1)**2) )
    IntrpMat(:,:) = matmul(Vint, refElem%invV)

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()

    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none
    !------------------------------------------------------------------------

    call PROF_rapstart( "final", 1 )

    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
    call q0%Final()
    call qexact%Final()
    call u%Final()
    call V%Final()

    call fields_comm%Final()
    call mesh%Final()
    
    call Dx%Final()
    call Sx%Final()
    call Dy%Final()
    call Sy%Final()
    call Lift%Final()
    call refElem%Final()
    
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()

    return
  end subroutine final

end program test_advect2d
