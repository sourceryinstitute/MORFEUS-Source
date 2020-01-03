!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  use assertions_interface, only : assert
  use problem_discretization_interface, only :  problem_discretization
  use cartesian_grid_interface, only : cartesian_grid
  use inbox_interface, only : inbox, face_index
  implicit none

  type(problem_discretization) block_structured_grid
    !! encapsulate the global grid structure
  type(cartesian_grid) prototype
    !! pass the cartesian_grid type
  integer, parameter :: num_structured_grids(*) = [3,3,3]
    !! number of subdomains in each coordinate direction

  call assert( product(num_structured_grids) >= num_images(), "test-halo-exchange: enough blocks to distribute to images")

  call block_structured_grid%partition( num_structured_grids, prototype )

  print *,"Test passed"
end program
