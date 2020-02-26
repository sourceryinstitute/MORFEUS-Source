!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson and Karla Morris
  !! date: 2/24/2020
  !!
  !! Test solution of the governing equationsn in a 1D spherical geometry

  use sphere_1D_interface, only : sphere_1D
  use problem_discretization_interface, only : problem_discretization
  use kind_parameters, only : r8k
  implicit none

  integer, parameter :: max_digits=9
  character(len=max_digits) image_number
  character(len=*), parameter:: input_file = "1Dsphere.json"
  character(len=*), parameter:: base_name = "1Dsphere-derivatives"
  character(len=:), allocatable :: output_file
  type(problem_discretization) global_grid
  type(sphere_1D) sphere_geometry

  associate( me => this_image() )
    write(image_number,'(i4)') me
    output_file = base_name //"-image-"// trim(adjustl(image_number)) // ".csv"

    call sphere_geometry%build( input_file ) !! read geometrical information
    call global_grid%initialize_from_geometry( sphere_geometry ) !! partition block-structured grid & define grid vertex locations
    call global_grid%solve_governing_equations( duration = sphere_geometry%get_end_time() )
    call global_grid%write_output (output_file)

    sync all
    if (me==1) print *,"Test passed."
  end associate

end program
