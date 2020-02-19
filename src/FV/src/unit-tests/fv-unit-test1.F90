program example_fv_unit_test
  ! use: whatever...
  use, intrinsic :: ISO_FORTRAN_ENV, only: stdout => output_unit, stderr => error_unit
  implicit none

  write(stdout,*) "Test passed."
  stop 0
end program
