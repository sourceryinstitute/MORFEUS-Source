!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  use assertions_interface, only : assert
  use problem_discretization_interface, only :  problem_discretization
  use surfaces_interface, only : surfaces
  use plate_3D_interface, only : plate_3D
  use ellipsoidal_field_interface, only : ellipsoidal_field
  implicit none

  type(plate_3D) plate_geometry
  type(problem_discretization) global_grid
  integer, parameter :: max_digits=9
  character(len=*), parameter :: input = "3Dplate-low-resolution-halo.json"
  character(len=*), parameter:: base_name = "3Dplate-low-resolution-halo"
  character(len=:), allocatable :: output
  character(len=max_digits) image_number
  type(ellipsoidal_field) ellipsoidal_function

  associate( me => this_image() )
    write(image_number,'(i4)') me
    output = base_name //"-image-"// trim(adjustl(image_number)) // ".vtu"
  end associate

  call plate_geometry%build(input)

  call global_grid%initialize_from_geometry(plate_geometry)
  call global_grid%set_scalars( [ellipsoidal_function] )
  call global_grid%set_scalar_flux_divergence( exact_result=[ellipsoidal_function] )
  call global_grid%write_output(output)

  sync all
  print *,"Test passed"
end program
