!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module problem_discretization_interface
  use object_interface, only : object
  use structured_grid_interface, only : structured_grid
  use geometry_interface, only : geometry
  use kind_parameters, only : r8k, i4k
  use plate_3D_interface, only : plate_3D
  implicit none

  private
  public :: problem_discretization, setter

  integer, parameter :: space_dimension=3

  abstract interface
    pure function field_function(x, y, z) result(f_xyz)
      import r8k
      real(r8k), intent(in), dimension(:,:,:) :: x, y, z
      real(r8k), allocatable :: f_xyz(:,:,:)
    end function
  end interface

  type setter
    procedure(field_function), pointer, nopass :: define_scalar=>null()
  end type

  type, extends(object) :: problem_discretization
    private
    integer global_block_shape_(space_dimension)
      !! global shape of the structured_grid blocks
    type(structured_grid), allocatable :: vertices(:)
      !! grid nodal locations: size(vertices) == number of blocks owned by the executing image
    type(structured_grid), allocatable :: scalar_fields(:,:)
      !! scalar values at the grid nodes: size(scalar_fields,1)==size(vertices), size(scalar_fields,2) == number of scalar fields
    class(geometry), allocatable :: problem_geometry
  contains
    procedure partition
    procedure my_subdomains
    procedure block_indicial_coordinates
    procedure block_identifier
    procedure block_load
    procedure user_defined_vertices
    generic :: set_vertices => user_defined_vertices
    procedure set_analytical_scalars
    generic :: set_scalars => set_analytical_scalars
    procedure initialize_from_plate_3D
    generic :: initialize_from_geometry => initialize_from_plate_3D
    procedure write_formatted
#ifdef HAVE_UDDTIO
    generic :: write(formatted) => write_formatted
#endif
  end type

  interface

    module subroutine write_formatted (this,unit,iotype, v_list, iostat, iomsg)
      !! Formatted user-defined derived type output
      implicit none
      class(problem_discretization), intent(in) ::this
      integer, intent(in) :: unit, v_list(:)
      character (len=*), intent(in) :: iotype
      integer, intent(out) :: iostat
      character(len=*), intent(inout) :: iomsg
    end subroutine

    module subroutine initialize_from_plate_3D(this,plate_3D_geometry)
      !! Define a grid with points only at the corners of each structured-grid block subdomain
      implicit none
      class(problem_discretization), intent(inout) :: this
      type(plate_3D), intent(in) :: plate_3D_geometry
    end subroutine

    module subroutine partition(this,global_block_shape)
      !! Define the distribution of subdomains across images for the given shape of the block-structured partitions
      !! (impure because of image-control statement in emulated co_sum -- may be pure when replaced by intrinsic co_sum)
      implicit none
      class(problem_discretization), intent(inout) :: this
      integer, intent(in) :: global_block_shape(:)
    end subroutine

    module subroutine user_defined_vertices(this,x_nodes,y_nodes,z_nodes,block_identifier)
      !! Define the vertex locations within the specified structured_grid block
      implicit none
      class(problem_discretization), intent(inout) :: this
      real(r8k), intent(in) :: x_nodes(:,:,:), y_nodes(:,:,:), z_nodes(:,:,:)
      integer, intent(in) :: block_identifier
    end subroutine

    module subroutine set_analytical_scalars(this, setters)
      !! Use functions to define scalar values at vertex locations
      implicit none
      class(problem_discretization), intent(inout) :: this
      type(setter), intent(in) :: setters(:)
    end subroutine

    pure module function my_subdomains(this) result(block_identifier_range)
      !! Result contains the first & last identifiers for blocks owned by this image
      implicit none
      class(problem_discretization), intent(in) :: this
      integer, parameter :: num_bounds=2
      integer block_identifier_range(num_bounds)
    end function

    pure module function block_indicial_coordinates(this,n) result(ijk)
      !! Calculate the 3D location of the block that has the provided 1D block identifer
      implicit none
      class(problem_discretization), intent(in) :: this
      integer, intent(in) :: n
      integer :: ijk(space_dimension)
    end function

    pure module function block_identifier(this,ijk) result(n)
      !! Calculate the 1D block identifer associated with the provided 3D block location
      implicit none
      class(problem_discretization), intent(in) :: this
      integer, intent(in) :: ijk(space_dimension)
      integer :: n
    end function

    pure module function block_load(this) result(num_blocks)
      !! the result is a measure of the workload partitioned to the executing image
      implicit none
      class(problem_discretization), intent(in) :: this
      integer num_blocks
    end function

  end interface

end module
