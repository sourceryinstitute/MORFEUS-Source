!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program test_class_vector
  !! label: Morfeus-FV
  !!
  !! Test vector class behavior

    use assertions_interface, only : assert
    use class_vector, only : vector, operator(.dot.), vector_
    use class_psblas, only : psb_dpk_

    implicit none

    type(vector) :: a(3), b, c
    real(psb_dpk_), parameter :: tolerance=1.0E-06, one=1.0_psb_dpk_, zero=0.0_psb_dpk_

    a = [vector_(zero, zero, zero), vector_(one, zero, zero), vector_(zero, one, zero)]
    b = vector_(zero, zero, one) ! perpindicular to all vectors in a
    c = a .dot. b ! = [a(1) .dot. b,  a(2) .dot. b,  a(3) .dot. b] = [zero, zero, zero]

    test_operators: &
    block
      call assert( c%mag() <= tolerance, "fv-test-class_vector: mag()" )
    end block test_operators

    print *, "Test passed."

end program test_class_vector
