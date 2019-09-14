!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(block_metadata_interface) block_metadata_implementation
  implicit none
  !! author: Damian Rouson
  !! date: August 8, 2019
  !! summary: procedure implementations for encapsulatingd metadata for blocks in block-structured grids

contains

  module procedure set_tag
    this%tag_ = tag
  end procedure

  module procedure set_label
    this%label_ = label
  end procedure

  module procedure set_subdomain
  !call assert( shape(subdomain)==[space_dimension, num_end_points], "shape(subdomain)==[space_dimension, num_end_points]" )
    this%subdomain_ = subdomain
  end procedure

  module procedure set_max_spacing
    this%max_spacing_ = max_spacing
  end procedure

  module procedure get_tag
    this_tag = this%tag_
  end procedure

  module procedure get_label
    this_label = this%label_
  end procedure

  module procedure get_subdomain
    this_subdomain = this%subdomain_
  end procedure

  module procedure get_max_spacing
    this_max_spacing = this%max_spacing_
  end procedure

end submodule block_metadata_implementation
