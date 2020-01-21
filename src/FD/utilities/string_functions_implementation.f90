!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
submodule(string_functions_interface) string_functions_implementation
  implicit none

contains

  module procedure file_extension
    character(len=:), allocatable :: name_

    name_ = trim(file_name)
    associate( dot_location => index(name_, '.', back=.true.) )
      if (dot_location < len(name_)) then
        extension = name_(dot_location+1:)
      else
        extension = ""
      end if
    end associate
  end procedure

  module procedure base_name
    character(len=:), allocatable :: name_

    name_ = trim(file_name)
    associate( dot_location => index(name_, '.', back=.true.) )
      if (dot_location < len(name_)) then
        base = name_(1:dot_location-1)
      else
        base = ""
      end if
    end associate
  end procedure

end submodule
