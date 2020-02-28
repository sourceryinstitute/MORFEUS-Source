!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
#ifndef USE_ASSERTIONS
# define USE_ASSERTIONS .false.
#endif
module assertions_interface
  !! summary: Utility for runtime checking of logical assertions.
  !!
  !! Compile with -DNO_ASSERTIONS to turn assertions off
  !!
  !! Use case 1
  !! ----------
  !!    Pass the optional success argument & check for false return value as an indication of assertion failure:
  !!
  !!    use assertions_interface, only : assert,assertions
  !!    if (assertions) call assert( 2 > 1, "always true inequality", success)
  !!    if (error_code/=0) call my_error_handler()
  !!
  !! Use case 2
  !! ----------
  !!    Error-terminate if the assertion fails:
  !!
  !!    use assertions_interface, only : assert,assertions
  !!    if (assertions) call assert( 2 > 1, "always true inequality")
  !!
  implicit none
  private
  public :: assert
  public :: assertions
  public :: max_errmsg_len

! Set the USE_ASSERTIONS constant below using the C preprocessor:
!
!    gfortran -cpp -DUSE_ASSERTIONS=.false. -c assertions_interface.f90
!
! or set the corresponding NO_ASSERTIONS variable defined in the CMakeLists.txt file in this directory:
!
!    FC=caf cmake <path-to-icar-source-dir> -DNO_ASSERTIONS=ON
!
! Conditioning assertion calls on this compile-time constant enables optimizing compilers
! to eliminate assertion calls during a dead-code removal phase of optimization.

  logical, parameter :: assertions=USE_ASSERTIONS
  integer, parameter :: max_errmsg_len = len( &
  "warning (183): FASTMEM allocation is requested but the libmemkind library is not linked in, so using the default allocator.")
  !! longest Intel compiler error messagea (see https://intel.ly/35x84yr).

  interface
#ifndef HAVE_ERROR_STOP_IN_PURE
    impure &
#endif
    elemental module subroutine assert(assertion,description,diagnostic_data,success)
      !! Report on the truth of an assertion or error-terminate on assertion failure
      implicit none
      logical, intent(in) :: assertion
        !! Most assertions will be expressions, e.g., call assert( i>0, "positive i")
      character(len=*), intent(in) :: description
        !! Brief statement of what is being asserted
      class(*), intent(in), optional :: diagnostic_data
        !! Optional assertion result
      logical, intent(out), optional :: success
        !! Optional assertion result
    end subroutine
  end interface
end module
