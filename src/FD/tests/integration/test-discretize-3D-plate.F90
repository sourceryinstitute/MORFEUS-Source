!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson
  !!
  !! Test the initialization of a problem_discretization from json_file input

  use assertions_interface, only : assert, assertions
  use plate_3D_interface, only : plate_3D
  use problem_discretization_interface, only : problem_discretization
  implicit none

  call minimal_resolution

  ! Outer surface BC: 500 K for all external surfaces
  ! Initial condition: 293 K everywhere, 10 kW/m

  print *,"Test passed."

contains

  subroutine minimal_resolution()
    integer, parameter :: z_dir=3, upper_boundary=2, lower_boundary=1
    integer, parameter :: minimum_z_grid_pts = 100
    type(plate_3D) :: plate_geometry
    type(problem_discretization) :: global_grid
    integer file_unit, open_status, io_status
    INTEGER, DIMENSION(0) :: v_list
    character(len=132)  :: io_message
    integer, parameter :: success=0
    character(len=*), parameter :: input_file="problem_description.json", output_file="3Dplate.vtk"

    call plate_geometry%build( input_file ) !! read geometrical information
    call global_grid%initialize_from( plate_geometry ) !! partition the block-structured grid and define the grid vertex locations

    open(newunit=file_unit,file=output_file,iostat=open_status)
    call assert(open_status==success, output_file//" opened succesfully")
#ifdef HAVE_UDDTIO
    write(file_unit,*) global_grid
#else
    call global_grid%write_formatted (file_unit, 'DT', v_list, io_status, io_message)
#endif

  end subroutine minimal_resolution

end program
