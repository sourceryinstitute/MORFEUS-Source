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
  implicit none

contains

  module procedure get_sender_block_id
    this_sender_block_id = this%sender_block_id
  end procedure

  module procedure set_sender_block_id
    this%sender_block_id = sender_block_id
  end procedure

  module procedure set_step
    this%step = step
  end procedure

  module procedure set_message
    this%x_f = this_x_f
    this%x_b = this_x_b
    this%s_flux_f = this_s_flux_f
    this%s_flux_b = this_s_flux_b
  end procedure

  module procedure sender_block_id_null
    is_null = (this%sender_block_id == null_sender_id)
  end procedure

end submodule package_implementation
