!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module package_interface
  !! author: Damian Rouson and Karla Morris
  !! date: 1/2/2019
  !! Encapsulate halo data to be communicated across structured_grid block boundaries
  use kind_parameters, only : r8k
  implicit none

  private
  public :: package, null_neighbor_id

  integer, parameter :: null_neighbor_id=-1

  type flux_planes
    real(r8k), allocatable, dimension(:,:,:) :: fluxes
      !! surface-normal scalar flux components: shape = shape(positions)
  end type

  type package
    !! basic transmission data. extend this type to add coordinate-specific data
    private
    integer :: neighbor_block_id
      !! block id corresponding to the destination block for this outbound package
    integer :: step
      !! time step
    type(flux_planes), allocatable, dimension(:) :: surface_normal_fluxes
      !! size = number of scalars; using a derived type allows for setting the number of scalars without knowing surface
      !! grid resolution and orientation (both of which are determined by shape(positions)).
    real(r8k), allocatable, dimension(:,:,:) :: positions
      !! flux planar locations: shape = [Nx, Ny, Nz], where one of Nx|Ny|Nz is 1, denoting the surface-normal direction
  contains
    procedure set_neighbor_block_id
    procedure set_step
    procedure set_surface_flux_positions
    procedure set_num_scalars
    procedure set_normal_scalar_fluxes
    procedure get_neighbor_block_id
    procedure neighbor_block_id_null
    procedure copy
    generic :: assignment(=) => copy
  end type

  interface

    pure module subroutine set_surface_flux_positions(this, positions)
      !! define grid locations at this structured_grid block's surface
      implicit none
      class(package), intent(inout) :: this
      real(r8k), dimension(:,:,:), intent(in) :: positions
        !! flux planar locations: shape = [Nx, Ny, Nz], where one of Nx|Ny|Nz is 1
    end subroutine

    elemental module subroutine set_num_scalars(this, num_scalars)
      !! establish the number of scalars
      implicit none
      class(package), intent(inout) :: this
      integer, intent(in) :: num_scalars
        !! set size of surface_normal_fluxes array
    end subroutine

    pure module function get_neighbor_block_id(this) result(this_neighbor_block_id)
      !! result is block_id correspdonding to the destination structured_grid block for this package
      implicit none
      class(package), intent(in) :: this
      integer :: this_neighbor_block_id
    end function

    elemental module subroutine set_neighbor_block_id(this, neighbor_block_id)
      !! set the block_id correspdonding to the destination structured_grid block for this package
      implicit none
      class(package), intent(inout) :: this
      integer, intent(in) :: neighbor_block_id
    end subroutine

    elemental module subroutine set_step(this, step)
      !! record the time step in order to verify the correct synchronization across images
      implicit none
      class(package), intent(inout) :: this
      integer, intent(in) :: step
    end subroutine

    pure module subroutine set_surface_positions(this, positions)
      !! define grid point locations for block surfaces
      implicit none
      class(package), intent(inout) :: this
      real(r8k), intent(in), dimension(:,:,:) :: positions
        !! flux locations: shape = [Nx, Ny, Nz], where one of Nx|Ny|Nz is 1
    end subroutine

    pure module subroutine set_normal_scalar_fluxes(this, fluxes, scalar_id)
      !! set datum to be communicated across structured_grid block internal surfaces
      implicit none
      class(package), intent(inout) :: this
      real(r8k), intent(in), dimension(:,:,:) :: fluxes
        !! surface-normal scalar flux components: shape = [Nx, Ny, Nz], where one of Nx|Ny|Nz is 1
      integer, intent(in) :: scalar_id
    end subroutine

    elemental module function neighbor_block_id_null(this) result(is_null)
      !! result is true if for external boundaries (no block sends halo data to a boundary)
      implicit none
      class(package), intent(in) :: this
      logical is_null
    end function

    pure module subroutine copy(this, rhs)
      !! copy rhs package components into this package
      class(package), intent(inout) :: this
      type(package), intent(in) :: rhs
    end subroutine

  end interface

end module package_interface
