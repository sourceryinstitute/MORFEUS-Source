!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
include "surfaces_interface.f90"
  !! required to work around a gfortran 8.3 bug 93158 (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93158)

submodule(surfaces_interface) surfaces_implementation
  !! author: Damian Rouson and Karla Morris
  !! date: 12/27/2019
  !! Implement procedures for exchanging information with halo blocks in block-structured grid
  use assertions_interface, only : assert, max_errmsg_len, assertions
  implicit none

  type(surfaces) singleton[*]
    !! Singleton pattern: one instance per image.
    !! Design: using a derived-type coarray instead of component coarrays is both simpler and facilitates setting different
    !! component array bounds on different images, which facilitates using the global block_identifier as the first index.

contains

  module procedure is_external_boundary
    if (assertions) then
      call assert(any(face==[forward,backward]), "surfaces%is_external_boundary: any(face==[forward,backward])")
      call assert(any(coordinate_direction==[1,2,3]), "surfaces%is_external_boundary: any(coordinate_direction==[1,2,3])")
    end if
    is_external = singleton%halo_data(block_id, coordinate_direction, face)%sender_block_id_unset()
  end procedure

  module procedure set_halo_data
    !! With gfortran versions > 8.2, this function can be reduced to a single intrinsic assignment:
    !! singleton%halo_data = my_halo_data
    integer alloc_stat
    integer, parameter :: success=0
    character(len=max_errmsg_len) error_message

    if(allocated(singleton%halo_data)) deallocate(singleton%halo_data)

    associate( lower => lbound(my_halo_data), upper => ubound(my_halo_data) )
      allocate(singleton%halo_data( lower(1):upper(1), lower(2):upper(2), lower(3):upper(3) ), source = my_halo_data, &
        stat = alloc_stat, errmsg = error_message)
      if (assertions) call assert(alloc_stat == success, "surfaces%set_halo_data: allocate(singleton%halo_data)", error_message)
    end associate
  end procedure

  module procedure get_halo_data
    !! With gfortran versions > 8.2, this function can be reduced to a single intrinsic assignment:
    !! singleton_halo_data = singleton%halo_data
    integer alloc_stat
    integer, parameter :: success=0
    character(len=max_errmsg_len) error_message

    if (assertions) call assert(allocated(singleton%halo_data),"surfaces%get_halo_data: allocated(singleton%halo_data)")

    associate( lower => lbound(singleton%halo_data), upper => ubound(singleton%halo_data) )
      allocate(singleton_halo_data( lower(1):upper(1), lower(2):upper(2), lower(3):upper(3) ), source = singleton%halo_data, &
        stat = alloc_stat, errmsg = error_message)
      if (assertions) call assert(alloc_stat == success, "surfaces%get_halo_data: allocate(singleton_halo_data)", error_message)
    end associate
  end procedure

end submodule
