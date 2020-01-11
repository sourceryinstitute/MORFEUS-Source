!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module surfaces_interface
  !! author: Damian Rouson
  !! date: 12/26/2019
  !!
  !! Encapsulate block boundary data for structured_grid halo exchanges
  use package_interface, only : package
  use iso_c_binding, only : enumeration=>c_int
  implicit none

  private
  public :: surfaces, enumeration, face_normal, backward, forward

  enum, bind(C)
    enumerator :: backward=1, forward
  end enum

  integer(enumeration), parameter, dimension(*) :: face_normal = [backward, forward]
    !! surface outward-normal direction: 'forward' for the direction of increasing coordinate; 'backward' for decreasing

  type surfaces
    !! hexahedral structured_grid block surface data: all components will be allocated to have the
    !! dimensions [my_blocks(first):my_blocks(last), space_dimension, size([forward, backward]))
    private
    class(package), allocatable, dimension(:,:,:) :: halo_data
  contains
    procedure, nopass :: is_external_boundary
    procedure, nopass :: set_halo_data
    procedure, nopass :: get_halo_data
  end type

  interface

    elemental module function is_external_boundary(block_id, coordinate_direction, face) result(is_external)
      !! result is .true. if the identified structured_grid block surface corresponds to a problem domain boundary
      implicit none
      integer, intent(in) :: block_id, coordinate_direction
      integer(enumeration), intent(in) :: face
      logical is_external
    end function

    module subroutine set_halo_data(my_halo_data, my_blocks)
      !! define halo_data component array
      implicit none
      class(package), intent(in), dimension(:,:,:) :: my_halo_data
      integer, intent(in), dimension(:) :: my_blocks
    end subroutine

    module subroutine get_halo_data(singleton_halo_data)
      !! output singleton_halo_data with the following bounds:
      !! lbounds=[my_blocks(first), 1, backward], ubounds=[my_blocks(last), space_dimension, forward]
      !! This can't be a function because the function result would not preserve the desired bounds.
      implicit none
      class(package), dimension(:,:,:), allocatable, intent(out) :: singleton_halo_data
    end subroutine

  end interface

end module
