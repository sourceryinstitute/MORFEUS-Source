!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module differentiable_field_interface
  use grid_interface, only : grid
  implicit none

  private
  public :: differentiable_field

  type, abstract :: differentiable_field
  contains
    procedure(field_interface), deferred :: evaluate
    procedure(field_interface), deferred :: laplacian
  end type

  abstract interface
    function field_interface(this, grid_points) result(f)
      use grid_interface, only : grid
      import differentiable_field
      implicit none
      class(differentiable_field), intent(in)  :: this
      class(grid), intent(in)  :: grid_points
      class(grid), allocatable :: f
    end function
  end interface

end module differentiable_field_interface
