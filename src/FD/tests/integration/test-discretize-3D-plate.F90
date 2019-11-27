!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson and Karla Morris
  !! date: 9/9/2019
  !! Test the initialization of a problem_discretization from json_file input

  use assertions_interface, only : assert
  use plate_3D_interface, only : plate_3D
  use problem_discretization_interface, only : problem_discretization
  implicit none

  call create_grid_for_plate(input="3Dplate-low-resolution-layers.json", output="3Dplate-low-resolution-layers")
  call create_grid_for_plate(input="3Dplate-high-resolution-layers.json", output="3Dplate-high-resolution-layers")

  print *,"Test passed."

contains

  subroutine create_grid_for_plate( input, output)
    implicit none
    character(len=*), intent(in) :: input, output
    type(plate_3D) :: plate_geometry
    type(problem_discretization) :: global_grid

    call plate_geometry%build( input ) !! read geometrical information
    call global_grid%initialize_from_geometry( plate_geometry ) !! partition block-structured grid & define grid vertex locations

    call global_grid%write_output (output, 'vtk') !! TODO. Make more sophisticated to allow calling of other output types

  end subroutine

end program
