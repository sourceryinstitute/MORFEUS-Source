!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program test_stopwatch
  !! label: Morfeus-FV
  !!
  !! Test stopwatch class behavior

    use assertions_interface, only : assert
    use class_stopwatch, only : stopwatch
    USE class_psblas

    implicit none

    type(stopwatch) :: watch

    watch = stopwatch_(icontxt_())
    call watch%tic()
    call watch%toc()

    call assert( watch%partial_() == watch%total_(), "test_stopwatch: watch%partial_() == watch%total_()" )

    print *, "Test passed."

end program test_stopwatch
