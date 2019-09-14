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

  use assertions_interface, only : assert
  use plate_3D_interface, only : plate_3D
  use problem_discretization_interface, only : problem_discretization
  implicit none

  call create_grid_for_plate(input="3Dplate-low-resolution.json", output="3Dplate-low-resolution.vtk")
  call create_grid_for_plate(input="3Dplate-high-resolution.json", output="3Dplate-high-resolution.vtk")

  ! Outer surface BC: 500 K for all external surfaces
  ! Initial condition: 293 K everywhere, 10 kW/m

  print *,"Test passed."

contains

  subroutine output_grid( mesh, file_unit )
    type(problem_discretization) :: mesh
    integer file_unit

#ifdef HAVE_UDDTIO
    write(file_unit,*) mesh
#else
    block
      integer, dimension(0) :: v_list
      character(len=132) io_message
      integer io_status
      call mesh%write_formatted (file_unit, 'DT', v_list, io_status, io_message)
    end block
#endif
  end subroutine

  subroutine create_grid_for_plate( input, output)
    character(len=*), intent(in) :: input, output
    type(plate_3D) :: plate_geometry
    type(problem_discretization) :: global_grid
    integer file_unit, open_status
    integer, parameter :: success=0

    call plate_geometry%build( input ) !! read geometrical information
    call global_grid%initialize_from_geometry( plate_geometry ) !! partition block-structured grid & define grid vertex locations

    open(newunit=file_unit, file=output, iostat=open_status)
    call assert(open_status==success, output//" opened succesfully")

    call output_grid( global_grid, file_unit)

  end subroutine

end program
