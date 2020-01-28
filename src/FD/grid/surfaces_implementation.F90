!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
include "surfaces_interface.F90"
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


  pure function first_index( block_id ) result(outbox_1st_index)
    !! result is the halo_outbox first index corresponding to the block identifier block_id
    integer, intent(in) :: block_id
    integer outbox_1st_index

    if (assertions) call assert(allocated(global_block_partitions), "surfaces(first_index): allocated(global_block_partitions)")

    outbox_1st_index = block_id - global_block_partitions(this_image()) + 1

    if (assertions) then
      call assert(allocated(singleton%halo_outbox), "surfaces(first_index): allocated(singleton%halo_outbox)")
      call assert(lbound(singleton%halo_outbox, 1) <= outbox_1st_index .and. outbox_1st_index <= ubound(singleton%halo_outbox, 1), &
        "surfaces(first_index): outbox_1st_index in bounds")
    end if
  end function

  module procedure get_surface_normal_spacing
    associate( nearest_plane => singleton%halo_outbox(block_id, coordinate_direction, face_direction)%get_positions() )
      if (assertions) then
        call assert(  &
          assertion = count(shape(nearest_plane)==1)==1, &
          description = "surfaces%get_surface_normal_spacing: count(shape(nearest_plane)==1)")
      end if
      associate( normal_direction => findloc(shape(nearest_plane), value=1, dim=1, back=.false.))
         dx_normal = nearest_plane(1,1,1,normal_direction)
      end associate
    end associate
  end procedure

  module procedure set_halo_outbox
    global_block_partitions = block_partitions
    singleton%halo_outbox = ( my_halo_outbox )
     !! parentheses prevent GCC 8.3 internal compiler error for actual arguments for which any lbound /= 1
  end procedure

  module procedure set_num_scalars
    call singleton%halo_outbox%set_num_scalars(num_scalars)
  end procedure

  module procedure set_normal_scalar_fluxes
    if (assertions) then
      call assert(allocated(singleton%halo_outbox), "surfaces%set_normal_scalar_fluxes: allocated(singleton%halo_outbox)")
    end if
    call singleton%halo_outbox(first_index(block_id), coordinate_direction, face)%set_normal_scalar_fluxes(s_flux_normal, scalar_id)
  end procedure

  module procedure get_halo_outbox

    if (assertions) call assert(allocated(singleton%halo_outbox),"surfaces%get_halo_outbox: allocated(singleton%halo_outbox)")

    singleton_halo_outbox = singleton%halo_outbox
  end procedure

  module procedure get_block_image
    if (assertions) then
      call assert(allocated(global_block_partitions), "surfaces%set_surface_package: allocated(global_block_partitions)")
    end if
    image = findloc( block_id >= global_block_partitions, value=.true., dim=1, back=.true.)
  end procedure

  module procedure get_global_block_partitions
    block_partitions = global_block_partitions
  end procedure

  module procedure get_neighbor_block_id
    neighbor_block_id= singleton%halo_outbox(first_index(my_block_id), coordinate_direction, face_direction)%get_neighbor_block_id()
  end procedure

  module procedure get_surface_positions
    type(package) neighbor_package

    neighbor_package = singleton[image]%halo_outbox(first_index(block_id), coordinate_direction, face_direction)
    positions = neighbor_package%get_positions()
  end procedure

  module procedure get_normal_scalar_fluxes
    type(package) neighbor_package

    neighbor_package = singleton[image]%halo_outbox(first_index(block_id), coordinate_direction, face_direction)
    fluxes = neighbor_package%get_fluxes(scalar_id)
  end procedure

  module procedure is_external_boundary
    is_external = singleton%halo_outbox(first_index(block_id), coordinate_direction, face)%neighbor_block_id_null()
  end procedure

end submodule
