!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
include "surfaces_interface.f90"
  !! required to work around a gfortran 8.3 bug 93158 (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93158)

submodule(surfaces_interface) surfaces_implementation
  !! author: Damian Rouson and Karla Morris
  !! date: 12/27/2019
  !! Implement procedures for exchanging information with halo blocks in block-structured grid
  use assertions_interface, only : assert, max_errmsg_len, assertions
#ifndef HAVE_FINDLOC
  use emulated_intrinsics_interface, only : findloc
#endif
  implicit none

  type(surfaces) singleton[*]
    !! Singleton pattern: one instance per image.
    !! Design: using a derived-type coarray instead of component coarrays is both simpler and facilitates setting different
    !! component array bounds on different images, which facilitates using the global block_identifier as the first index.

  integer, allocatable, dimension(:) :: global_block_partitions
  integer, parameter :: success=0

contains

  module procedure is_external_boundary
    if (assertions) then
      call assert(any(face==[forward,backward]), "surfaces%is_external_boundary: any(face==[forward,backward])")
      call assert(any(coordinate_direction==[1,2,3]), "surfaces%is_external_boundary: any(coordinate_direction==[1,2,3])")
    end if
    is_external = singleton%halo_inbox(block_id, coordinate_direction, face)%sender_block_id_null()
  end procedure

  module procedure set_halo_inbox
    integer alloc_stat
    character(len=max_errmsg_len) error_message

    call assert( size(my_halo_inbox,1)==my_blocks(2)-my_blocks(1)+1, &
       "surfaces%set_halo_inbox: size(my_halo_inbox,1)==my_blocks(2)-my_blocks(1)+1")

    if(allocated(singleton%halo_inbox)) deallocate(singleton%halo_inbox)

    allocate(singleton%halo_inbox, source = my_halo_inbox, stat = alloc_stat, errmsg = error_message)

    if (assertions) then
       call assert(allocated(singleton%halo_inbox), "surfaces% allocated(singleton%halo_inbox)", error_message)
       call assert( &
         all([lbound(singleton%halo_inbox), ubound(singleton%halo_inbox)] == [lbound(my_halo_inbox), ubound(my_halo_inbox)]), &
         "surfaces%is_external_boundary: halo_index bounds copied correctly" )
     end if

    global_block_partitions = block_partitions
  end procedure

  module procedure get_halo_inbox
    integer alloc_stat
    character(len=max_errmsg_len) error_message

    if (assertions) call assert(allocated(singleton%halo_inbox),"surfaces%get_halo_inbox: allocated(singleton%halo_inbox)")

    associate( lower => lbound(singleton%halo_inbox), upper => ubound(singleton%halo_inbox) )
      allocate(singleton_halo_inbox( lower(1):upper(1), lower(2):upper(2), lower(3):upper(3) ), source = singleton%halo_inbox, &
        stat = alloc_stat, errmsg = error_message)
      if (assertions) call assert(alloc_stat == success, "surfaces%get_halo_inbox: allocate(singleton_halo_inbox)", error_message)
    end associate
  end procedure

end submodule
