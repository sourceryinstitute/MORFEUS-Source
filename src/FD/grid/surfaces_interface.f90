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
  use iso_c_binding, only : enumeration=>c_int
  implicit none

  private
  public :: surfaces, envelope, enumeration
  public :: face_index, east, west, north, south, up, down

  enum, bind(C)
    enumerator :: east=1, west, north, south, up, down
      !! boundary direction: rename to coordinate-specific labels in (sub)modules that use this module
  end enum

  integer(enumeration), parameter, dimension(*) :: face_index = [east, west, north, south, up, down]

  type envelope
    !! basic transmission data. extend this type to add coordinate-specific data
    integer sender, step
  end type

  type surfaces
    !! hexahedral structured_grid block surface data: all components will be allocated to have the
    !! bounds [my_blocks(first):my_blocks(last), space_dimension, size([east,west]))
    private
    integer(enumeration), allocatable, dimension(:,:,:) :: direction
    logical, allocatable, dimension(:,:,:) :: internal
    type(envelope), allocatable, dimension(:,:,:) :: halo_data
  contains
    procedure, nopass :: set_direction
    procedure, nopass :: set_internal
    procedure, nopass :: set_halo_data
  end type

  interface

    module subroutine set_direction(direction)
      !! set surface direction component array
      implicit none
      integer(enumeration), intent(in), dimension(:,:,:) :: direction
    end subroutine

    module subroutine set_halo_data(halo_data)
      !! set halo_data component array
      implicit none
      class(envelope), intent(in), dimension(:,:,:) :: halo_data
    end subroutine

    module subroutine set_internal(internal)
      !! set logical internal component array
      implicit none
      logical, intent(in), dimension(:,:,:) :: internal
    end subroutine

  end interface

end module
