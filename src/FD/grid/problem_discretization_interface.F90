!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module problem_discretization_interface
  !! author: Damian Rouson
  !! date: 9/9/2019
  !!
  !! Encapsulate the problem geometry, grid points, and dependent variables
  use object_interface, only : object
  use structured_grid_interface, only : structured_grid
  use geometry_interface, only : geometry
  use kind_parameters, only : r8k, i4k
  use plate_3D_interface, only : plate_3D
  use differentiable_field_interface, only : differentiable_field
  use surfaces_interface, only : surfaces
  use package_interface, only : package
  implicit none

  private
  public :: problem_discretization

  type, extends(object) :: problem_discretization
    private
    class(structured_grid), allocatable :: block_map
      !! hook for invoking block_indicial_coordinates and block_identifier
    integer, allocatable, dimension(:) :: block_partitions
      !! block ownership specification
    type(surfaces) block_surfaces
      !! global surface information with singleton communication buffers
    class(structured_grid), allocatable :: vertices(:)
      !! grid nodal locations: size(vertices) == number of blocks owned by the executing image
    class(structured_grid), allocatable :: scalar_fields(:,:)
      !! scalar values at the grid nodes: size(scalar_fields,1)==size(vertices), size(scalar_fields,2)==number of scalar fields
    class(structured_grid), allocatable :: scalar_flux_divergence(:,:)
      !! div( D grad(s)): same dimensions as scalar_fields
    class(geometry), allocatable :: problem_geometry
      !! description of problem domain and material identities; child types: plate_3D, cylinder_2D, sphere_1D
  contains
    procedure get_surface_packages
    procedure get_block_surfaces
    procedure get_global_block_shape
    procedure partition
    procedure my_blocks
    procedure block_indicial_coordinates
    procedure block_identifier
    procedure block_load
    procedure user_defined_vertices
    procedure set_analytical_scalars
    procedure num_scalars
    procedure num_scalar_flux_divergences
    procedure initialize_from_plate_3D
    procedure set_scalar_flux_divergence
#ifndef FORD
    generic :: set_vertices => user_defined_vertices
    generic :: set_scalars => set_analytical_scalars
    generic :: initialize_from_geometry => initialize_from_plate_3D
#endif
    procedure :: write_output
  end type

  interface

    module subroutine get_surface_packages(this, this_surface_packages)
      !! get surface data packages of this image
      !! This cannot be a function because the function result would not preserve the desired bounds.
      implicit none
      class(problem_discretization), intent(in) :: this
      type(package), allocatable, dimension(:,:,:), intent(out) :: this_surface_packages
    end subroutine

    module function get_block_surfaces(this) result(this_block_surfaces)
      !! get surfaces owned by this image
      implicit none
      class(problem_discretization), intent(in) :: this
      type(surfaces) this_block_surfaces
    end function

    module function get_global_block_shape(this) result(block_shape)
      !! result is shape of array of structured_grid blocks
      implicit none
      class(problem_discretization), intent(in) :: this
      integer, dimension(:), allocatable :: block_shape
    end function

    module subroutine write_output (this, filename)
      !! Generic write output interface
      implicit none
      class(problem_discretization), intent(in) ::this
      character (len=*), intent(in) :: filename
    end subroutine

    module subroutine initialize_from_plate_3D(this, plate_3D_geometry)
      !! Define a grid with points only at the corners of each structured-grid block block
      implicit none
      class(problem_discretization), intent(inout) :: this
      type(plate_3D), intent(in) :: plate_3D_geometry
    end subroutine

    module subroutine partition(this,global_block_shape,prototype)
      !! Define the distribution of blocks across images for the given shape of the block-structured partitions
      !! (impure because of image-control statement in emulated co_sum -- may be pure when replaced by intrinsic co_sum)
      implicit none
      class(problem_discretization), intent(inout) :: this
      integer, intent(in) :: global_block_shape(:)
      class(structured_grid), intent(in) :: prototype
    end subroutine

    module subroutine user_defined_vertices(this,x_nodes,y_nodes,z_nodes,block_identifier)
      !! Define the vertex locations within the specified structured_grid block
      implicit none
      class(problem_discretization), intent(inout) :: this
      real(r8k), intent(in) :: x_nodes(:,:,:), y_nodes(:,:,:), z_nodes(:,:,:)
      integer, intent(in) :: block_identifier
    end subroutine

    module subroutine set_analytical_scalars(this, scalar_setters )
      !! Use functions to define scalar values at vertex locations
      implicit none
      class(problem_discretization), intent(inout) :: this
      class(differentiable_field), intent(in), dimension(:) :: scalar_setters
    end subroutine

    module subroutine set_scalar_flux_divergence(this, exact_result)
      !! Compute and store div( D grad( S )) for each scalar S and diffusion coefficient D
      implicit none
      class(problem_discretization), intent(inout) :: this
      class(differentiable_field), intent(in), dimension(:), optional :: exact_result
    end subroutine

    pure module function num_scalars(this) result(num_scalar_fields)
      !! Result contains the total number of scalar_field components
      implicit none
      class(problem_discretization), intent(in) :: this
      integer num_scalar_fields
    end function

    pure module function num_scalar_flux_divergences(this) result(num_divergences)
      !! Result contains the total number of scalar_field components
      implicit none
      class(problem_discretization), intent(in) :: this
      integer num_divergences
    end function

   pure module function my_blocks(this) result(block_identifier_range)
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
      integer, dimension(:), allocatable :: ijk(:)
    end function

    pure module function block_identifier(this,ijk) result(n)
      !! Calculate the 1D block identifer associated with the provided 3D block location
      implicit none
      class(problem_discretization), intent(in) :: this
      integer, intent(in) :: ijk(:)
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
