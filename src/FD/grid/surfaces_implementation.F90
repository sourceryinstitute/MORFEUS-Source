!! category: Morfeus-FD
!!
!! ## Copyright Notice
!!
!!
!!     (c) 2019-2020 Guide Star Engineering, LLC
!!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!!     contract # NRC-HQ-60-17-C-0007
!!

#ifndef FORD
include "surfaces_interface.F90"
  ! required to work around a gfortran 8.3 bug 93158 (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93158)
#endif

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
    is_external = singleton%halo_outbox(block_id, coordinate_direction, face)%neighbor_block_id_null()
  end procedure

  module procedure set_halo_outbox
    !! a shorter implementation of this procedure would simply assign my_halo_outbox to singleton%halo_outbox
    !! GCC 8 compiler bugs necessitate the source allocation and instead
    integer alloc_stat
    character(len=max_errmsg_len) error_message

    if (assertions) then
      call assert(allocated(my_halo_outbox), "surfaces%set_halo_outbox: allocated(my_halo_outbox)")
      associate( me => this_image() )
        associate( my_blocks => [block_partitions(me), block_partitions(me+1)-1] )
          call assert( all( my_blocks == [lbound(my_halo_outbox,1), ubound(my_halo_outbox,1)] ), &
            "my_blocks == all([lbound(my_halo_outbox,1), ubound(my_halo_outbox,1)])")
        end associate
      end associate
    end if

    if(allocated(singleton%halo_outbox)) deallocate(singleton%halo_outbox)

    allocate(singleton%halo_outbox, source = my_halo_outbox, stat = alloc_stat, errmsg = error_message)
    if (assertions) call assert(alloc_stat==success, "surfaces%set_halo_outbox: allocate(singleton%halo_outbox)", error_message)

    global_block_partitions = block_partitions
  end procedure

  module procedure get_halo_outbox
    integer alloc_stat
    character(len=max_errmsg_len) error_message

    if (assertions) call assert(allocated(singleton%halo_outbox),"surfaces%get_halo_outbox: allocated(singleton%halo_outbox)")

    associate( lower => lbound(singleton%halo_outbox), upper => ubound(singleton%halo_outbox) )
      allocate(singleton_halo_outbox( lower(1):upper(1), lower(2):upper(2), lower(3):upper(3) ), source = singleton%halo_outbox, &
        stat = alloc_stat, errmsg = error_message)
      if (assertions) call assert(alloc_stat == success, "surfaces%get_halo_outbox: allocate(singleton_halo_outbox)", error_message)
    end associate
  end procedure

  module procedure set_normal_scalar_fluxes
    if (assertions) &
      call assert(allocated(singleton%halo_outbox), "surfaces%set_normal_scalar_fluxes: allocated(singleton%halo_outbox)")
    call singleton%halo_outbox( block_id, coordinate_direction, face)%set_normal_scalar_fluxes(s_flux_normal, scalar_id)
  end procedure

  module procedure get_block_image
    if (assertions) &
      call assert(allocated(global_block_partitions), "surfaces%set_surface_package: allocated(global_block_partitions)")
    image = findloc( block_id >= global_block_partitions, value=.true., dim=1, back=.true.)
  end procedure

  module procedure set_num_scalars
    call singleton%halo_outbox%set_num_scalars(num_scalars)
  end procedure

end submodule
