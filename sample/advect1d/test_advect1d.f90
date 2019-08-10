#include "scalelib.h"
program test_advect1d
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
  use scale_element_line
  use scale_localmesh_1d
  use scale_mesh_linedom1d

  use scale_localmeshfield_base, only: LocalMeshField1D
  use scale_meshfield_base, only: MeshField1D
  use scale_meshfieldcomm_base, only: MeshFieldContainer
  use scale_meshfieldcomm_1d, only: MeshFieldComm1D

  use scale_file_history_meshfield, only: &
    FILE_HISTORY_meshfield_put,           &
    FILE_HISTORY_meshfield_write

  use scale_file_history, only: &
    FILE_HISTORY_set_nowdate

  use scale_time_manager, only: &
    TIME_manager_advance,                              &
    TIME_NOWDATE, TIME_NOWMS, TIME_NOWSTEP,            &
    TIME_DTSEC, TIME_NSTEP 

  !-----------------------------------------------------------------------------
  implicit none

  integer, parameter :: NeGX = 32
  integer, parameter :: NLocalMeshPerPrc = 1

  ! sin, cosbell, hat
  character(len=H_SHORT) :: InitShapeName = 'sin'

  real(RP), parameter :: dom_xmin = -1.0_RP
  real(RP), parameter :: dom_xmax = +1.0_RP
  real(RP), parameter :: ADV_VEL  = 1.0_RP

  type(LineElement) :: refElem
  integer, parameter :: PolyOrder = 3
  logical, parameter :: DumpedMassMatFlag = .false.
  type(sparsemat) :: Dx, Sx, Lift
  integer, parameter :: PolyOrderErrorCheck = 6

  type(MeshLineDom1D), target :: mesh
  type(MeshField1D), target :: q, q0, qexact  
  type(MeshField1D), target :: u
  type(MeshFieldComm1D) :: fields_comm
  type(MeshFieldContainer) :: field_list(2)  
  integer :: HST_ID(2)

  integer :: n, k, p
  type(LocalMesh1D), pointer :: lcmesh
  
  integer :: nowstep
  integer :: rkstage
  integer, parameter :: nrkstage = 3
  real(RP) :: tsec_
  real(RP) :: rkcoef1(nrkstage) = (/ 0.0_RP, 3.0_RP/4.0_RP, 1.0_RP/3.0_RP /)
  real(RP) :: rkcoef2(nrkstage) = (/ 1.0_RP, 1.0_RP/4.0_RP, 2.0_RP/3.0_RP /)

  real(RP) :: IntrpMat(PolyOrderErrorCheck,PolyOrder+1)
  real(RP) :: intw_intrp(PolyOrderErrorCheck)
  real(RP) :: x_intrp(PolyOrderErrorCheck)

  !-------------------------------------------------------

  call init()
  call set_initcond()

  field_list(1)%field1d => q
  field_list(2)%field1d => u

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
          q%local(n)%val, q0%local(n)%val, u%local(n)%val )        
      end do
      call PROF_rapend( 'update_dyn', 1)

    end do
    
    !* Advance time
    call TIME_manager_advance()

    tsec_ = TIME_NOWDATE(6) + TIME_NOWMS
    if (mod(tsec_,0.25_RP) == 0.0_RP) then 
      write(*,*) "t=", real(tsec_), "[s]"
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
  subroutine update_dyn( dt, rkcoef_1, rkcoef_2, lmesh, elem, q_, q0_, u_ )
    implicit none

    real(RP), intent(in) :: dt
    real(RP), intent(in) :: rkcoef_1, rkcoef_2
    class(LocalMesh1D), intent(in) :: lmesh
    class(elementbase1D), intent(in) :: elem
    real(RP), intent(inout) :: q_(elem%Np,lmesh%NeA)
    real(RP), intent(in)    :: q0_(elem%Np,lmesh%NeA)
    real(RP), intent(in)    :: u_(elem%Np,lmesh%NeA)

    real(RP) :: Fx(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: dqdt(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne)

    !------------------------------------------------------------------------

    call PROF_rapstart( 'update_dyn_cal_bndflux', 2)
    call cal_del_flux_dyn( del_flux,                            & ! (out)
      q_, u_, lmesh%normal_fn(:,:,1), lmesh%vmapM, lmesh%vmapP, & ! (in)
      lmesh, elem )                                               ! (in)
    call PROF_rapend( 'update_dyn_cal_bndflux', 2)

    !-----
    call PROF_rapstart( 'update_dyn_cal_interior', 2)

    do k = lmesh%NeS, lmesh%NeE
      call sparsemat_matmul(Dx, q_(:,k)*u_(:,k), Fx)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,k)*del_flux(:,k), LiftDelFlx)
      
      dqdt(:) = -(lmesh%Escale(:,k,1,1) * Fx(:) + LiftDelFlx)
      q_(:,k) = rkcoef_1*q0_(:,k) + rkcoef_2*(q_(:,k) + dt * dqdt)
    end do

    call PROF_rapend( 'update_dyn_cal_interior', 2)

    return
  end subroutine update_dyn

  subroutine cal_del_flux_dyn( del_flux, q_, u_, nx, vmapM, vmapP, lmesh, elem )
    implicit none

    class(LocalMesh1D), intent(in) :: lmesh
    class(elementbase1D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) ::  q_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  u_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
     
    integer :: i, iP, iM
    real(RP) :: alpha
    !------------------------------------------------------------------------

    do i=1, elem%NfpTot*lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      alpha = 0.5_RP*abs(u_(iP) + u_(iM))
      del_flux(i) = 0.5_RP*(                       &
          (q_(iP)*u_(iP) - q_(iM)*u_(iM))*nx(i)    &
        - alpha*(q_(iP) - q_(iM))              )
    end do

    return
  end subroutine cal_del_flux_dyn


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine evaluate_error(tsec)
    use scale_polynominal, only: Polynominal_genLegendrePoly
    implicit none
    real(DP), intent(in) :: tsec


    real(RP) :: l2error
    real(RP) :: q_intrp(PolyOrderErrorCheck)
    real(RP) :: qexact_intrp(PolyOrderErrorCheck)
    real(RP) :: x_uwind(refElem%Np)
    real(RP) :: x_uwind_intrp(PolyOrderErrorCheck)
    real(RP) :: pos_intrp(PolyOrderErrorCheck)
    real(RP) vx(2)
    !------------------------------------------------------------------------

    l2error = 0.0_RP
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE

        x_uwind(:) = get_upwind_pos(lcmesh%pos_en(:,k,1), ADV_VEL, tsec, dom_xmin, dom_xmax)

        vx(:) = lcmesh%pos_ev(lcmesh%EToV(k,:),1)
        pos_intrp(:) = vx(1) + 0.5_RP*(x_intrp(:) + 1.0_RP)*(vx(2) - vx(1))
        x_uwind_intrp(:) = get_upwind_pos(pos_intrp(:), ADV_VEL, tsec, dom_xmin, dom_xmax)

        qexact%local(n)%val(:,k) = get_profile(InitShapeName, x_uwind)
        qexact_intrp(:) = get_profile(InitShapeName, x_uwind_intrp)
        q_intrp(:) = matmul(IntrpMat, q%local(n)%val(:,k))

        l2error = l2error &
          + sum(   lcmesh%J(1,k) * intw_intrp(:)          &
                * ( q_intrp(:) - qexact_intrp(:) )**2 ) 
      end do
    end do

    write(*,*) "L2 error:", sqrt(l2error)/(dom_xmax - dom_xmin)
  end subroutine evaluate_error

  function get_upwind_pos(pos, ADV_VEL, nowtime, dom_min, dom_max) result(upos)
    real(RP), intent(in) :: pos(:)
    real(RP), intent(in) :: ADV_VEL
    real(RP), intent(in) :: nowtime
    real(RP), intent(in) :: dom_min, dom_max
    real(RP) :: upos(size(pos))

    real(RP) :: period
    !-------

    period = ADV_VEL*nowtime/(dom_max - dom_min)
    
    upos(:) = pos(:) - (ADV_VEL*nowtime - dble(period)*(dom_max - dom_min))
    where (upos < dom_min)
      upos = dom_max + (upos - dom_min)
    end where
  end function get_upwind_pos

  subroutine set_initcond()
    implicit none

    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh => mesh%lcmesh_list(n)
      do k=lcmesh%NeS, lcmesh%NeE
        q%local(n)%val(:,k) = get_profile(InitShapeName, lcmesh%pos_en(:,k,1))
        qexact%local(n)%val(:,k) = q%local(n)%val(:,k)
        u%local(n)%val(:,k) = ADV_VEL
      end do
    end do

    call FILE_HISTORY_meshfield_put(HST_ID(1), q)
    call FILE_HISTORY_meshfield_put(HST_ID(2), qexact)
    call FILE_HISTORY_meshfield_write()   
  
    return
  end subroutine set_initcond

  function get_profile(profile_name, x) result(profile)
    implicit none

    character(*), intent(in) :: profile_name
    real(RP), intent(in) :: x(:)
    real(RP) :: profile(size(x))

    real(RP) :: half_width = 0.3_RP

    profile(:) = 0.0_RP

    select case(InitShapeName)
    case ('sin')
      profile(:) = sin( PI*x(:) )
    case ('cosbell')
      where( abs(x) <= half_width )
        profile(:) = (1.0_RP + cos(PI*x/half_width))*0.5_RP
      end where
    case ('hat')
      where( abs(x) <= half_width )
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

    integer :: comm, myrank, nprocs
    logical :: ismaster
    
    real(RP) :: Vint(PolyOrderErrorCheck,PolyOrder+1)
    real(RP) :: P_int1D_ori(PolyOrderErrorCheck,PolyOrder+1)
    integer :: n_, l_

    !----------------------------------------------
    
    call PRC_MPIstart( comm )
     
    call PRC_SINGLECOM_setup( comm,    & ! [IN]
      nprocs, myrank, ismaster )         ! [OUT]
    
    call PRC_ERRHANDLER_setup( .false., ismaster ) ! [IN]
    
    ! setup scale_io
    call IO_setup( "test", allow_noconf = .true. )
    
    ! setup log
    call IO_LOG_setup( myrank, ismaster )   
    
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
    call Lift%Init(refElem%Lift)

    call mesh%Init( &
      NeGX,                              &
      dom_xmin, dom_xmax,                &
      refElem, NLocalMeshPerPrc )
    
    call mesh%Generate()

    !---
    call q%Init( "q", "1", mesh )
    call q0%Init( "q0", "1", mesh )
    call qexact%Init( "qexact", "1", mesh )
    call u%Init( "u", "m/s", mesh )
    call fields_comm%Init(2, 0, mesh)
    
    call FILE_HISTORY_meshfield_setup( mesh )
    call FILE_HISTORY_reg( q%varname, "q", q%unit, HST_ID(1), dim_type='X')
    call FILE_HISTORY_reg( qexact%varname, "qexact", q%unit, HST_ID(2), dim_type='X')

    !------------------
    x_intrp(:) = Polynominal_GenGaussLegendrePt( PolyOrderErrorCheck )
    intw_intrp(:) = Polynominal_GenGaussLegendrePtIntWeight( PolyOrderErrorCheck )
    P_int1D_ori(:,:) = Polynominal_GenLegendrePoly(refElem%PolyOrder, x_intrp)

    do n_=1, PolyOrderErrorCheck  
      do l_=1, refElem%Np
        Vint(n_,l_) =  P_int1D_ori(n_,l_) * sqrt(real(l_-1,kind=RP) + 0.5_RP)
      end do
    end do
    IntrpMat(:,:) = matmul(Vint, refElem%invV)   

    call PROF_rapend( "init", 1 )
    return
  end subroutine init

  subroutine final()
    use scale_file_history_meshfield, only: &
      FILE_HISTORY_meshfield_finalize
    use scale_time_manager, only: TIME_manager_Final    
    implicit none

    call PROF_rapstart( "final", 1 )
    call FILE_HISTORY_meshfield_finalize()

    call q%Final()
    call q0%Final()
    call qexact%Final()
    call u%Final()

    call fields_comm%Final()
    call mesh%Final()
    
    call Dx%Final()
    call Sx%Final()
    call Lift%Final()
    call refElem%Final()
    
    call TIME_manager_Final

    call PROF_rapend( "final", 1 )
    call PROF_rapend( "total", 0 )
    call PROF_rapreport
    call PRC_MPIfinish()

    return
  end subroutine final

end program test_advect1d
