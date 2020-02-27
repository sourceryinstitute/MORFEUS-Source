!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module cartesian_grid_interface
  !! author: Damian Rouson and Karla Morris
  !! date: 9/9/2019
  !! subject: perform coordinate-specific tasks on cartesian structured_grid objects
  use structured_grid_interface, only : structured_grid
  use differentiable_field_interface, only : differentiable_field
  use geometry_interface, only : geometry
  use surfaces_interface, only : surfaces
  implicit none

  private
  public :: cartesian_grid

  type, extends(structured_grid) :: cartesian_grid
  contains
    procedure :: set_up_div_scalar_flux
    procedure :: div_scalar_flux
    procedure :: block_indicial_coordinates
    procedure :: block_identifier
    procedure :: build_surfaces
    procedure :: block_identifier_in_bounds
    procedure :: block_coordinates_in_bounds
  end type

  interface

    module subroutine set_up_div_scalar_flux(this, vertices, block_surfaces, div_flux_internal_points)
      !! define the scalar flux divergence at points internal to grid blocks grid; define block-surface data on halo blocks
      implicit none
      class(cartesian_grid), intent(in) :: this
      class(structured_grid), intent(in) :: vertices
      type(surfaces), intent(inout) :: block_surfaces
      class(structured_grid), intent(inout) :: div_flux_internal_points
    end subroutine

    pure module subroutine div_scalar_flux(this, vertices, block_surfaces, div_flux)
      !! comunicate scalar fluxes between block neighbors in a halo exchange; compute scalar flux divergence at block boundaries
      implicit none
      class(cartesian_grid), intent(in) :: this
      class(structured_grid), intent(in) :: vertices
      type(surfaces), intent(in) :: block_surfaces
      class(structured_grid), intent(inout) :: div_flux
    end subroutine

    module subroutine build_surfaces(this, problem_geometry, vertices, block_faces, block_partitions, num_scalars)
      !! allocate the surfaces array for use in halo exchanges and boundary conditions
      implicit none
      class(cartesian_grid), intent(in) :: this
      class(geometry), intent(in) :: problem_geometry
      class(structured_grid), intent(in), dimension(:), allocatable :: vertices
      type(surfaces), intent(inout) :: block_faces
      integer, intent(in), dimension(:) :: block_partitions
      integer, intent(in) :: num_scalars
    end subroutine

    pure module function block_indicial_coordinates(this,n) result(ijk)
      !! calculate the 3D location of the block that has the provided 1D block identifer
      implicit none
      class(cartesian_grid), intent(in) :: this
      integer, intent(in) :: n
      integer, dimension(:), allocatable :: ijk
    end function

    pure module function block_identifier(this,ijk) result(n)
      !! calculate the 1D block identifer associated with the provided 3D block location
      implicit none
      class(cartesian_grid), intent(in) :: this
      integer, intent(in), dimension(:) :: ijk
      integer :: n
    end function

    pure module function block_identifier_in_bounds(this, id) result(in_bounds)
      implicit none
      class(cartesian_grid), intent(in) :: this
      integer, intent(in) :: id
      logical in_bounds
    end function

    pure module function block_coordinates_in_bounds(this, ijk) result(in_bounds)
      implicit none
      class(cartesian_grid), intent(in) :: this
      integer, intent(in), dimension(:) :: ijk
      logical in_bounds
    end function

  end interface

end module
