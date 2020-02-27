!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson and Karla Morris
  !! date: 11/15/2019
  !! Test the computation of spatial derivatives on a structured_grid

  use plate_3D_interface, only : plate_3D
  use problem_discretization_interface, only : problem_discretization
  use ellipsoidal_field_interface, only : ellipsoidal_field
  implicit none

  integer, parameter :: max_digits=9
  character(len=max_digits) image_number
  character(len=*), parameter:: input_file = "3Dplate-high-resolution-layers.json"
  character(len=*), parameter:: base_name = "3Dplate-high-resolution-layers-derivatives"
  character(len=:), allocatable :: output_file
  type(problem_discretization) global_grid
  type(plate_3D) plate_geometry
  type(ellipsoidal_field) ellipsoidal_function

  associate( me => this_image() )
    write(image_number,'(i4)') me
    output_file = base_name //"-image-"// trim(adjustl(image_number)) // ".vtu"

    call plate_geometry%build( input_file ) !! read geometrical information
    call global_grid%initialize_from_geometry( plate_geometry ) !! partition block-structured grid & define grid vertex locations
    call global_grid%set_scalars( [ellipsoidal_function] )
    call global_grid%set_scalar_flux_divergence( exact_result=[ellipsoidal_function] )
    call global_grid%write_output (output_file)

    sync all
    if (me==1) print *,"Test passed."
  end associate

end program
