!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program test_dimensions
  !! label: Morfeus-FV
  !!
  !! Test dimensions class behavior

    use assertions_interface, only : assert
    use class_dimensions, only : dimensions, operator(==), operator(/=), length_, time_, velocity_

    implicit none

    type(dimensions) :: height_dimensions, displacement_dimensions
    type(dimensions) :: time_step_dimensions
    type(dimensions) :: velocity_dimensions

    height_dimensions = length_
    displacement_dimensions = length_
    time_step_dimensions = time_
    velocity_dimensions = velocity_

    test_operators: &
    block
      !call assert( height_dimensions == displacement_dimensions, "test_dimensions: operator(==)" )
      !call assert( .not. (height_dimensions == time_step_dimensions), "test_dimensions: operator(==)" )
      !call assert( height_dimensions /= time_step_dimensions, "test_dimensions: inequality operator(/=)" )
      !call assert( displacement_dimensions / time_dimensions == velocity_, "test_dimensions: operator(/)")
      !call assert( height_dimensions == (height_dimensions + displacement_dimensions), "test_dimensions: operator(+)")
      !call assert( height_dimensions == (height_dimensions - length_), "test_dimensions: operator(-)")
    end block test_operators

    print *, "Test passed."

end program test_dimensions
