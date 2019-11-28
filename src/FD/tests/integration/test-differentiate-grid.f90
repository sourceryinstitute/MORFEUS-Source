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

  scalar_setter%define_scalar => f

  associate( me => this_image() )
    write(image_number,'(i4)') me
    output_file = base_name //"-image-"// trim(adjustl(image_number))! // ".vtk"

    call plate_geometry%build( input_file ) !! read geometrical information
    call global_grid%initialize_from_geometry( plate_geometry ) !! partition block-structured grid & define grid vertex locations
    call global_grid%set_scalars( [scalar_setter] )
    call global_grid%set_div_scalar_flux()

!    open(newunit=file_unit, file=output_file, iostat=open_status)
!    call assert(open_status==success, output_file//" opened succesfully")
    call output_grid(global_grid, output_file)

    sync all
    if (me==1) print *,"Test passed."
  end associate

contains

  subroutine output_grid( mesh, file_name )
    type(problem_discretization), intent(in) :: mesh
    character(len=*), intent(in) :: file_name

    call mesh%write_output (file_name, 'vtk')

  end subroutine

  pure function f(position_vectors) result(f_value)
    real(r8k), intent(in) :: position_vectors(:,:,:,:)
    real(r8k), allocatable :: f_value(:,:,:)
    if (assertions) call assert(lbound(position_vectors,4)==1 .and. ubound(position_vectors,4)>=2, "field dimension >= 2")
    associate( x=>position_vectors(:,:,:,1), y=>position_vectors(:,:,:,2) )
      f_value = x**2 * y
    end associate
  end function

end program
