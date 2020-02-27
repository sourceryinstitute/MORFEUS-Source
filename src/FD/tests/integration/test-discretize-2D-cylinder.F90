!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson and Karla Morris
  !! date: 2/7/2020
  !! Test the initialization of a problem_discretization from cylindrical-geometry json_file input

  use assertions_interface, only : assert
  use cylinder_2D_interface, only : cylinder_2D
  use problem_discretization_interface, only : problem_discretization
  implicit none

  !call create_grid_for_rod(input="2Dcylinder.json", output="2Dcylinder-geometry.vtu")
  call create_grid_for_rod(input="3Dplate-low-resolution-layers.json", output="2Dcylinder-geometry.vtu")

  print *,"Test passed."

contains

  subroutine create_grid_for_rod( input, output)
    implicit none
    character(len=*), intent(in) :: input, output
    type(cylinder_2D) :: rod_geometry
    type(problem_discretization) :: global_grid

    call rod_geometry%build( input ) !! read geometrical information
    call global_grid%initialize_from_geometry( rod_geometry ) !! partition block-structured grid & define grid vertex locations
    call global_grid%write_output (output)

  end subroutine

end program
