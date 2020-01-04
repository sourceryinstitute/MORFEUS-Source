!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module surfaces_interface
  !! author: Damian Rouson
  !! date: 12/26/2019
  !! Encapsulate information and procedures for structured_grid block halo exchanges
  use package_interface, only : package
  use iso_c_binding, only : enumeration=>c_int
  implicit none

  private
  public :: surfaces, enumeration, face_normal, backward, forward

  enum, bind(C)
    enumerator :: backward=1, forward
  end enum

  integer(enumeration), parameter, dimension(*) :: face_normal = [backward, forward]
    !! surface outward-normal direction for a given block

  type surfaces
    !! hexahedral structured_grid block surface data: all components will be allocated to have the
    !! bounds [my_blocks(first):my_blocks(last), space_dimension, size([forward, backward]))
    private
    class(package), allocatable, dimension(:,:,:) :: halo_data
  contains
    procedure, nopass :: is_external_boundary
    procedure, nopass :: set_halo_data
  end type

  type(surfaces), allocatable :: singleton[:]
    !! Singleton pattern (one instance per image).
    !!
    !! Design: making the object a coarray rather than the components coarrays frees us to use different component array bounds on
    !! different images, which in turn facilitates using the global block ID as the first index.
    !!
    !! Gfortran 8.3 workarounds:
    !! 1. Moving this to the submodule to eliminate the compiler warning about its not being used generates an ICE.
    !! 2. Declaring this as singleton[*] also generates an ICE.

  interface

    module function is_external_boundary(this, block_id, coordinate_direction, face) result(is_external)
      implicit none
      class(surfaces), intent(in) :: this
      integer, intent(in) :: block_id, coordinate_direction, face
      logical is_external
    end function

    module subroutine set_halo_data(my_halo_data)
      !! set halo_data component array
      implicit none
      class(package), intent(in), dimension(:,:,:) :: my_halo_data
    end subroutine

  end interface

end module
