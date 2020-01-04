!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(surfaces_interface) surfaces_implementation
  !! author: Damian Rouson and Karla Morris
  !! date: 12/27/2019
  !! Implement procedures for exchanging information with halo blocks in block-structured grid
  use assertions_interface, only : assert, max_errmsg_len, assertions
  implicit none

contains

  module procedure is_external_boundary
    if (assertions) then
      call assert(any(face==[forward,backward]), "surfaces%is_external_boundary: any(face==[forward,backward])")
      call assert(any(coordinate_direction==[1,2,3]), "surfaces%is_external_boundary: any(coordinate_direction==[1,2,3])")
    end if
    is_external = singleton%halo_data(block_id, coordinate_direction, face)%sender_block_id_unset()
  end procedure

  module procedure set_halo_data
    integer alloc_stat
    character(len=max_errmsg_len) error_message
    integer, parameter :: success=0

    if (.not. allocated(singleton)) allocate(singleton[*], stat=alloc_stat, errmsg=error_message)
    call assert( alloc_stat==success, "surfaces%set_halo_data: allocate(singleton[*])", error_message)

    singleton%halo_data = my_halo_data
  end procedure

end submodule
