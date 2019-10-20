!-------------------------------------------------------------------------------
#include "scaleFElib.h"
module mod_operator_fvm
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc


  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type, procedures
  !
  type, public :: operator_fvm
    integer :: flux_scheme_id

    integer :: KS, KE, KA
    integer :: IS, IE, IA
    integer :: JS, JE, JA
  contains
    procedure :: Init => operator_fvm_Init
    procedure :: Final => operator_fvm_Final
    procedure :: C_flux_XYW => operator_fvm_C_flux_XYW
    procedure :: C_flux_UYZ => operator_fvm_C_flux_UYZ
    procedure :: C_flux_XVZ => operator_fvm_C_flux_XVZ
  end type

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-------------------

  integer, parameter :: FLUX_SCHEME_ID_UD1 = 1
  integer, parameter :: FLUX_SCHEME_ID_CD2 = 2
  integer, parameter :: FLUX_SCHEME_ID_UD3 = 3
  integer, parameter :: FLUX_SCHEME_ID_CD4 = 4
  integer, parameter :: FLUX_SCHEME_ID_UD5 = 5


  real(RP), parameter :: F2  = 0.5_RP
  
  real(RP), parameter :: F31  = -1.0_RP/12.0_RP
  real(RP), parameter :: F32  =  7.0_RP/12.0_RP
  real(RP), parameter :: F33  =  3.0_RP/12.0_RP

  real(RP), parameter :: F41 =   7.0_RP/12.0_RP
  real(RP), parameter :: F42 = - 1.0_RP/12.0_RP

  real(RP), parameter :: F51  =  1.0_RP/60.0_RP
  real(RP), parameter :: F52  = -8.0_RP/60.0_RP
  real(RP), parameter :: F53  = 37.0_RP/60.0_RP
  real(RP), parameter :: F54  = -5.0_RP/60.0_RP
  real(RP), parameter :: F55  = 10.0_RP/60.0_RP

  integer, parameter :: DIM_X_ID  = 1
  integer, parameter :: DIM_Y_ID  = 2
  integer, parameter :: DIM_Z_ID  = 3
contains

!----------------

  subroutine operator_fvm_Init( this,    &
    flux_scheme_name,                    &
    KS, KE, KA, IS, IE, IA, JS, JE, JA )
    
    class(operator_fvm), intent(inout) :: this
    character(*), intent(in) :: flux_scheme_name
    integer, intent(in) :: KS, KE, KA
    integer, intent(in) :: IS, IE, IA
    integer, intent(in) :: JS, JE, JA
    !----------------------------------------

    select case( flux_scheme_name )
    case ('UD1')
      this%flux_scheme_id = FLUX_SCHEME_ID_UD1
    case ('CD2')
      this%flux_scheme_id = FLUX_SCHEME_ID_CD2
    case ('UD3')
      this%flux_scheme_id = FLUX_SCHEME_ID_UD3
    case ('CD4')
      this%flux_scheme_id = FLUX_SCHEME_ID_CD4
    case ('UD5')
      this%flux_scheme_id = FLUX_SCHEME_ID_UD5
    case default
      LOG_ERROR("operator_fvm_Init",*) trim(flux_scheme_name)//' is not supported. Check!'
      call PRC_abort
    end select

    this%KS = KS; this%KE = KE; this%KA = KA
    this%IS = IS; this%IE = IE; this%IA = IA
    this%JS = JS; this%JE = JE; this%JA = JA

    return
  end subroutine operator_fvm_Init

  subroutine operator_fvm_Final( this )
    class(operator_fvm), intent(inout) :: this
    !----------------------------------------
    
    return
  end subroutine operator_fvm_Final

  !----------------

  subroutine operator_fvm_C_flux_XYW( this, flux, vel, q )

    class(operator_fvm), intent(inout) :: this
    real(RP), intent(out) :: flux(this%KA, this%IA, this%JA)
    real(RP), intent(in) :: vel(this%KA, this%IA, this%JA)
    real(RP), intent(in) :: q(this%KA, this%IA, this%JA)

    integer :: i, j, k
    !----------------------------------------

    select case( this%flux_scheme_id )
    case (FLUX_SCHEME_ID_UD1)

    case (FLUX_SCHEME_ID_CD2)
      do j=this%JS, this%JE
      do i=this%IS, this%IE
      do k=this%KS-1, this%KE
        flux(k,i,j) = vel(k,i,j)*( F2*(q(k,i,j) + q(k+1,i,j)) )
      end do
      end do
      end do
    case (FLUX_SCHEME_ID_UD3)
      do j=this%JS, this%JE
      do i=this%IS, this%IE
      do k=this%KS-1, this%KE
        flux(k,i,j) = vel(k,i,j) &
         * ( ( F31 * ( q(k+2,i,j)+q(k-1,i,j) ) + F32 * ( q(k+1,i,j)+q(k,i,j) ) ) &
           - ( F31 * ( q(k+2,i,j)-q(k-1,i,j) ) + F33 * ( q(k+1,i,j)-q(k,i,j) ) ) * sign(1.0_RP,vel(k,i,j)) )
      end do
      end do
      end do 
    case (FLUX_SCHEME_ID_CD4)
      do j=this%JS, this%JE
      do i=this%IS, this%IE
      do k=this%KS-1, this%KE
        flux(k,i,j) = vel(k,i,j)*( F41*(q(k,i,j) + q(k+1,i,j)) + F42*(q(k+2,i,j) + q(k-1,i,j)) )
      end do
      end do
      end do      
    case (FLUX_SCHEME_ID_UD5)
      do j=this%JS, this%JE
      do i=this%IS, this%IE
      do k=this%KS-1, this%KE
        flux(k,i,j) = vel(k,i,j) &
          * ( ( F51 * ( q(k+3,i,j)+q(k-2,i,j) ) &
              + F52 * ( q(k+2,i,j)+q(k-1,i,j) ) &
              + F53 * ( q(k+1,i,j)+q(k,i,j) ) ) &
            - ( F51 * ( q(k+3,i,j)-q(k-2,i,j) ) &
              + F54 * ( q(k+2,i,j)-q(k-1,i,j) ) &
              + F55 * ( q(k+1,i,j)-q(k,i,j) ) ) * sign(1.0_RP,vel(k,i,j)) )         
      end do
      end do
      end do
    end select
  end subroutine operator_fvm_C_flux_XYW

  subroutine operator_fvm_C_flux_UYZ( this, flux, vel, q )

    class(operator_fvm), intent(inout) :: this
    real(RP), intent(out) :: flux(this%KA, this%IA, this%JA)
    real(RP), intent(in) :: vel(this%KA, this%IA, this%JA)
    real(RP), intent(in) :: q(this%KA, this%IA, this%JA)

    integer :: i, j, k
    !----------------------------------------

    select case( this%flux_scheme_id )
    case (FLUX_SCHEME_ID_UD1)

    case (FLUX_SCHEME_ID_CD2)
      do j=this%JS, this%JE
      do i=this%IS-1, this%IE
      do k=this%KS, this%KE
        flux(k,i,j) = vel(k,i,j)*( F2*(q(k,i,j) + q(k,i+1,j)) )
      end do
      end do
      end do
    case (FLUX_SCHEME_ID_UD3)
      do j=this%JS, this%JE
      do i=this%IS-1, this%IE
      do k=this%KS, this%KE
        flux(k,i,j) = vel(k,i,j) &
          * ( ( F31 * ( q(k,i+2,j)+q(k,i-1,j) ) + F32 * ( q(k,i+1,j)+q(k,i,j) ) ) &
            - ( F31 * ( q(k,i+2,j)-q(k,i-1,j) ) + F33 * ( q(k,i+1,j)-q(k,i,j) ) ) * sign(1.0_RP,vel(k,i,j)) )
      end do
      end do
      end do 
    case (FLUX_SCHEME_ID_CD4)
      do j=this%JS, this%JE
      do i=this%IS-1, this%IE
      do k=this%KS, this%KE
        flux(k,i,j) = vel(k,i,j)*( F41*(q(k,i,j) + q(k,i+1,j)) + F42*(q(k,i+2,j) + q(k,i-1,j)) )
      end do
      end do
      end do      
    case (FLUX_SCHEME_ID_UD5)
      do j=this%JS, this%JE
      do i=this%IS-1, this%IE
      do k=this%KS, this%KE
        flux(k,i,j) = vel(k,i,j) &
          * ( ( F51 * ( q(k,i+3,j)+q(k,i-2,j) ) &
              + F52 * ( q(k,i+2,j)+q(k,i-1,j) ) &
              + F53 * ( q(k,i+1,j)+q(k,i,j) ) ) &
            - ( F51 * ( q(k,i+3,j)-q(k,i-2,j) ) &
              + F54 * ( q(k,i+2,j)-q(k,i-1,j) ) &
              + F55 * ( q(k,i+1,j)-q(k,i,j) ) ) * sign(1.0_RP,vel(k,i,j)) )         
      end do
      end do
      end do      
    end select
  end subroutine operator_fvm_C_flux_UYZ

  subroutine operator_fvm_C_flux_XVZ( this, flux, vel, q )

    class(operator_fvm), intent(inout) :: this
    real(RP), intent(out) :: flux(this%KA, this%IA, this%JA)
    real(RP), intent(in) :: vel(this%KA, this%IA, this%JA)
    real(RP), intent(in) :: q(this%KA, this%IA, this%JA)

    integer :: i, j, k
    !----------------------------------------

    select case( this%flux_scheme_id )
    case (FLUX_SCHEME_ID_UD1)

    case (FLUX_SCHEME_ID_CD2)
      do j=this%JS-1, this%JE
      do i=this%IS, this%IE
      do k=this%KS, this%KE
        flux(k,i,j) = vel(k,i,j)*( F2*(q(k,i,j) + q(k,i,j+1)) )
      end do
      end do
      end do
    case (FLUX_SCHEME_ID_UD3)

    case (FLUX_SCHEME_ID_CD4)
      do j=this%JS-1, this%JE
      do i=this%IS, this%IE
      do k=this%KS, this%KE
        flux(k,i,j) = vel(k,i,j)*( F41*(q(k,i,j) + q(k,i,j+1)) + F42*(q(k,i,j+2) + q(k,i,j-1)) )
      end do
      end do
      end do      
    case (FLUX_SCHEME_ID_UD5)
    end select
  end subroutine operator_fvm_C_flux_XVZ

end module mod_operator_fvm

