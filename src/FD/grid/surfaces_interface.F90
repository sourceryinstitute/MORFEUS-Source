!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module surfaces_interface
  !! author: Damian Rouson
  !! date: 12/26/2019
  !!
  !! Encapsulate block boundary data for structured_grid halo exchanges
  use package_interface, only : package
  use kind_parameters, only : r8k
  use iso_c_binding, only : enumeration=>c_int
  implicit none

  private
  public :: surfaces, enumeration, backward, forward, face_name, coordinate_name, x_dir, y_dir, z_dir

  enum, bind(C)
    !! array indices
    enumerator :: backward=1, forward
      !! surface outward-normal direction: 'forward'/'backward' -> direction of increasing/decreasing coordinate values
    enumerator :: x_dir=1, y_dir=2, z_dir=3
  end enum

  character(len=*), parameter, dimension(*) :: face_name = ["backward","forward "]
  character(len=*), parameter, dimension(*) :: coordinate_name = ["x","y","z"]

  type surfaces
    !! hexahedral structured_grid block surface data
    private
    type(package), allocatable, dimension(:,:,:) :: halo_outbox
      !! allocate to dimensions [number of blocks, space_dimension, size([forward, backward]))
  contains
    procedure, nopass :: set_halo_outbox
    procedure, nopass :: set_num_scalars
    procedure, nopass :: set_normal_scalar_fluxes
    procedure, nopass :: get_surface_normal_spacing
    procedure, nopass :: get_halo_outbox
    procedure, nopass :: get_block_image
    procedure, nopass :: get_global_block_partitions
    procedure, nopass :: get_neighbor_block_id
    procedure, nopass :: get_surface_positions
    procedure, nopass :: get_normal_scalar_fluxes
    procedure, nopass :: is_external_boundary
  end type

  interface

    module subroutine set_halo_outbox(my_halo_outbox, block_partitions)
      !! define halo_outbox component array
      implicit none
      type(package), intent(in), dimension(:,:,:) :: my_halo_outbox
      integer, intent(in), dimension(:) :: block_partitions
    end subroutine

    module subroutine set_num_scalars(num_scalars)
      !! allocate surface-normal scalar flux array
      implicit none
      integer, intent(in) :: num_scalars
    end subroutine

    module subroutine set_normal_scalar_fluxes( block_id, coordinate_direction, face, s_flux_normal, scalar_id)
      !! define halo outbox for a specific surface
      implicit none
      integer, intent(in) :: block_id, coordinate_direction, scalar_id
      integer(enumeration), intent(in) :: face
      real(r8k), intent(in), dimension(:,:) :: s_flux_normal
        !! surface-normal scalar flux components: shape = [Ny, Nz] or [Nx, Nz] or [Nx, Ny]
    end subroutine

    pure module function get_surface_normal_spacing(image, block_id, coordinate_direction, face_direction) result(dx_normal)
      !! result is the distance to the nearest plane inside the specified block at the specified boundary
      implicit none
      integer, intent(in) :: image, block_id, coordinate_direction
      integer(enumeration), intent(in) :: face_direction
      real(r8k) dx_normal
    end function

    module function get_halo_outbox() result(singleton_halo_outbox)
      !! output singleton_halo_outbox of shape [number of blocks, space_dimension, size([backward,forward])]
      implicit none
      type(package), dimension(:,:,:), allocatable :: singleton_halo_outbox
    end function

    pure module function get_block_image(block_id) result(image)
      !! result is the image that owns the given structured_grid block
      implicit none
      integer, intent(in) :: block_id
      integer image
    end function

    pure module function get_global_block_partitions() result(block_partitions)
      !! result contains the block identifiers that start the subrange of blocks owned by each image such that
      !! [block_partitions(me), block_partitions(me+1)-1] = spans the block identifiers owned by image me
      implicit none
      integer, allocatable, dimension(:) :: block_partitions
    end function

    pure module function get_neighbor_block_id(my_block_id, coordinate_direction, face_direction) result(neighbor_block_id)
      !! result is the block_id of the neighbor adjacent to the image that owns the given structured_grid block
      implicit none
      integer, intent(in) :: my_block_id, coordinate_direction, face_direction
      integer neighbor_block_id
    end function

    pure module function get_surface_positions(image, block_id, coordinate_direction, face_direction) result(positions)
      !! result contains the vertices inside the designated grid block surface
      implicit none
      integer, intent(in) :: image, block_id, coordinate_direction, face_direction
      real(r8k), allocatable, dimension(:,:,:,:) :: positions
        !! surface vertices: shape=[Nx,Ny,Nz,space_dim] where findloc(shape(positions, value=1)) designates surface-normal direction
    end function

    pure module function get_normal_scalar_fluxes(image, block_id, coordinate_direction, face_direction, scalar_id) &
      result(fluxes)
      !! result contains the the surface-normal component of the designated scalar flus at the designated grid block surface
      implicit none
      integer, intent(in) :: image, block_id, coordinate_direction, face_direction, scalar_id
      real(r8k), allocatable, dimension(:,:) :: fluxes
        !! surface-normal scalar flux components: shape = [Ny, Nz] or [Nx, Nz] or [Nx, Ny] depending on orientation
    end function

    elemental module function is_external_boundary(block_id, coordinate_direction, face) result(is_external)
      !! result is .true. if the identified structured_grid block surface corresponds to a problem domain boundary
      implicit none
      integer, intent(in) :: block_id, coordinate_direction
      integer(enumeration), intent(in) :: face
      logical is_external
    end function

  end interface

end module
