!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(string_functions_interface) string_functions_implementation
  implicit none

contains

  module procedure file_extension
    character(len=:), allocatable :: name_

    name_ = trim(file_name)
    associate( dot_location => index(name_, '.', back=.true.) )
      extension = merge( name_(dot_location+1:), "", dot_location<len(name_))
    end associate
  end procedure

  module procedure base_name
    character(len=:), allocatable :: name_

    name_ = trim(file_name)
    associate( dot_location => index(name_, '.', back=.true.) )
      base = merge( name_(1:dot_location-1), "", dot_location<len(name_))
    end associate
  end procedure

end submodule
