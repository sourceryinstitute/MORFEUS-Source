!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(geometry_interface) geometry_implementation
  !! author: Damian Rouson
  !! date: 8/20/2019
  !! summary: geometry procedures
  implicit none

contains

    module procedure build
      call this%set_grid_specification( grid_description_file )
      call this%set_block_metadata()
    end procedure

end submodule geometry_implementation
