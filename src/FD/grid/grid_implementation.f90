!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(grid_interface) grid_implementation
    use assertions_interface, only : assertions,assert
    implicit none

    contains

        module procedure set_units
            implicit none
            !! define grid units
            this%units_ = units_obj
            call this%mark_as_defined
        end procedure

        module procedure get_units
            implicit none
            if (assertions) call assert(this%user_defined(),"grid%get_units: operand defined")
            !! Require user-defined passed-object dummy argument
            this_units = this%units_
        end procedure

end submodule grid_implementation
