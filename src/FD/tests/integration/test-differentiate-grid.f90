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

  use assertions_interface, only : assert
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
    output_file = base_name //"-image-"// trim(adjustl(image_number)) // ".vtk"

    call plate_geometry%build( input_file ) !! read geometrical information
    call global_grid%initialize_from_geometry( plate_geometry ) !! partition block-structured grid & define grid vertex locations
    call global_grid%set_scalars( [scalar_setter] )

    open(newunit=file_unit, file=output_file, iostat=open_status)
    call assert(open_status==success, output_file//" opened succesfully")
    call output_grid( global_grid, file_unit)

    sync all
    if (me==1) print *,"Test passed."
  end associate

contains

  subroutine output_grid( mesh, file_unit )
    type(problem_discretization), intent(in) :: mesh
    integer, intent(in) :: file_unit
    integer, dimension(0) :: v_list
    character(len=132) io_message
    integer io_status

    call mesh%write_formatted (file_unit, 'DT', v_list, io_status, io_message)
  end subroutine

  pure function f(x,y,z) result(f_value)
    real(r8k), intent(in), dimension(:,:,:) :: x,y,z
    real(r8k), allocatable :: f_value(:,:,:)
    f_value = x**2 * y
  end function

end program
