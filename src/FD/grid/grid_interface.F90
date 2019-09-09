!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module grid_interface
    use object_interface, only : object
    use units_interface, only : units
    implicit none

    private
    public :: grid

    type, extends(object) :: grid
        !! Morfeus universal base type for all grids
        private
        type(units) :: units_
    contains
        procedure set_units
        procedure get_units
    end type

    interface

#ifndef HAVE_ERROR_STOP_IN_PURE
        impure &
#endif
        elemental module subroutine set_units(this,units_obj)
            implicit none
            class(grid), intent(inout) :: this
            type(units), intent(in) :: units_obj
        end subroutine

#ifndef HAVE_ERROR_STOP_IN_PURE
        impure &
#endif
        elemental module function get_units(this) result(this_units)
            implicit none
            class(grid), intent(in) :: this
            type(units) :: this_units
        end function

    end interface

end module
