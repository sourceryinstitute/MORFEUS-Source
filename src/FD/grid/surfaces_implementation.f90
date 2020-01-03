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
  implicit none

  type(surfaces) singleton[*]
    !! single instance per image with codimension for inter-image communication of halo_data

contains

  module procedure set_halo_data
    singleton%halo_data = halo_data
  end procedure

  module procedure set_internal
    singleton%internal = internal
  end procedure

  module procedure set_direction
    singleton%direction = direction
  end procedure

end submodule
    !associate( neighbor_id => block_structured_grid%block_neighbor(block_, face_index(face)) )
    !  associate( neighbor_image => block_structured_grid%neighbor_image(neighbor_id) )
    !     associate( my_blocks_array_b_index=> findloc(block_, my_blocks_array) )
    !        call halo(my_blocks_array_b_index, face_index(face))%exchange() ! atomically decrement counter
    !     end associate
    !   end associate
    ! end associate

   !sync all
