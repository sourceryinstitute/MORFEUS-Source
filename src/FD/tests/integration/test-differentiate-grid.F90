!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson and Karla Morris
  !! date: 11/15/2019
  !! Test the computation of spatial derivatives on a structured_grid

  use assertions_interface, only : assert, assertions
  use plate_3D_interface, only : plate_3D
  use problem_discretization_interface, only : problem_discretization, setter
  use ellipsoidal_field_interface, only : ellipsoidal_field
  use kind_parameters, only : i4k, r8k
  implicit none

  integer file_unit, open_status
  integer(i4k), parameter :: success=0, max_digits=9, num_scalars=1
  character(len=max_digits) image_number
  character(len=*), parameter:: input_file = "3Dplate-high-resolution-layers.json"
  character(len=*), parameter:: base_name = "3Dplate-high-resolution-layers-derivatives"
  character(len=:), allocatable :: output_file
  type(problem_discretization) :: global_grid
  type(plate_3D) plate_geometry
  type(setter) scalar_setter
  type(ellipsoidal_field) ellipsoidal_function

  scalar_setter%define_scalar => f

  associate( me => this_image() )
    write(image_number,'(i4)') me
    output_file = base_name //"-image-"// trim(adjustl(image_number)) // ".vtu"

    call plate_geometry%build( input_file ) !! read geometrical information
    call global_grid%initialize_from_geometry( plate_geometry ) !! partition block-structured grid & define grid vertex locations
    call global_grid%set_scalars( [scalar_setter] )
    call global_grid%set_scalar_flux_divergence( exact_result=ellipsoidal_function )
    call global_grid%write_output (output_file)

    sync all
    if (me==1) print *,"Test passed."
  end associate

contains

  pure function f(position_vectors) result(f_value)
    real(r8k), intent(in) :: position_vectors(:,:,:,:)
    real(r8k), allocatable :: f_value(:,:,:)
    real(r8k), parameter ::  x_center = 3*(0.25E-01) + 1.E-01/2., x_max = 2*x_center
    real(r8k), parameter ::  y_center = 0.5E-01 + 2*(0.25E-01) + 3.E-01/2., y_max = 2*y_center
    real(r8k), parameter ::  z_center = 20.E-01/2., z_max=2*z_center
    if (assertions) call assert(lbound(position_vectors,4)==1 .and. ubound(position_vectors,4)==3, "field dimension == 3")
    associate( x=>position_vectors(:,:,:,1), y=>position_vectors(:,:,:,2), z=>position_vectors(:,:,:,3) )
      associate(r_sq=>((x-x_center)/(x_max-x_center))**2 + ((y-y_center)/(y_max-y_center))**2 + ((z-z_center)/(z_max-z_center))**2)
        f_value = 400. - r_sq
      end associate
    end associate
  end function

end program
