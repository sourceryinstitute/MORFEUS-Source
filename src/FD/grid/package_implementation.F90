!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
submodule(package_interface) package_implementation
  !! author: Damian Rouson and Karla Morris
  !! date: 1/2/2019
  !! Encapsulate information and procedures for structured_grid block halo exchanges
  use assertions_interface, only: assert, assertions
  implicit none

contains

  module procedure get_neighbor_block_id
    this_neighbor_block_id = this%neighbor_block_id
  end procedure

  module procedure set_neighbor_block_id
    this%neighbor_block_id = neighbor_block_id
  end procedure

  module procedure set_step
    this%step = step
  end procedure

  module procedure set_normal_scalar_fluxes
    this%s_flux_normal = s_flux_normal
    this%positions = positions
  end procedure

  module procedure neighbor_block_id_null
    is_null = (this%neighbor_block_id == null_neighbor_id)
  end procedure

  module procedure copy
    if (assertions) then
      call assert( all([allocated(rhs%s_flux_normal), allocated(rhs%positions)]), &
        " all([allocated(rhs%s_flux_normal), allocated(rhs%positions)])" )
    end if
    this%neighbor_block_id = rhs%neighbor_block_id
    this%step = rhs%step
    this%s_flux_normal = rhs%s_flux_normal
    this%positions = rhs%positions
  end procedure

end submodule package_implementation
