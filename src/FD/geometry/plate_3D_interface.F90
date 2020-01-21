!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module plate_3D_interface
  !! author: Damian Rouson and Karla Morris
  !!
  !! Encapsulate a 3D plate geometry and grid-verification method
  use geometry_interface, only : geometry
  use json_module, only : json_file
  use block_metadata_interface, only : block_metadata, subdomain_t, max_name_length, space_dimension, num_end_points
  use kind_parameters, only : r8k

  implicit none

  private
  public :: plate_3D

  type, extends(geometry) :: plate_3D
    !! encapsulate the grid specification for a plate_3D object
    private
    type(json_file) :: grid_specification
    type(block_metadata), dimension(:,:,:), allocatable :: metadata
    character(len=:), allocatable :: units_system
  contains
    procedure set_grid_specification
    procedure set_block_metadata
    procedure get_block_metadata_shape
    procedure get_block_domain
    procedure get_block_metadatum
    procedure get_block_metadata
  end type

  interface

    module subroutine set_grid_specification(this, grid_description_file)
      !! define json_file problem description
      implicit none
      class(plate_3D), intent(out) :: this
      character(len=*), intent(in) :: grid_description_file
    end subroutine

    module subroutine set_block_metadata(this)
      !! read grid metadata from json_file stored in child class
      implicit none
      class(plate_3D), intent(inout) :: this
    end subroutine

    module function get_block_metadata_shape(this) result(shape_)
      !! define the shape of the array of grid blocks
      implicit none
      class(plate_3D), intent(in) :: this
      integer, dimension(space_dimension) :: shape_
    end function

    module function get_block_domain(this, indicial_coordinates) result(this_domain)
      !! result is the spatial subdomain for the grid block with the given indicial_coordinates
      implicit none
      class(plate_3D), intent(in) :: this
      integer, dimension(space_dimension) :: indicial_coordinates
      integer, parameter :: num_end_points=2
      real(r8k), dimension(space_dimension,num_end_points) :: this_domain
    end function

    module function get_block_metadatum(this, indicial_coordinates) result(this_metadata_xyz)
      !! result is the block_metadata component for the block with the given indicial_coordinates
      implicit none
      class(plate_3D), intent(in) :: this
      integer, dimension(space_dimension) :: indicial_coordinates
      type(block_metadata) :: this_metadata_xyz
    end function

    module function get_block_metadata(this) result(this_metadata)
      !! result is the block_metadata component for all blocks
      implicit none
      class(plate_3D), intent(in) :: this
      type(block_metadata), dimension(:,:,:), allocatable :: this_metadata
    end function

  end interface

end module plate_3D_interface
