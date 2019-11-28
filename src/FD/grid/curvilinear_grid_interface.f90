!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module curvilinear_grid_interface
  use structured_grid_interface, only : structured_grid
  implicit none

  private
  public :: curvilinear_grid

  type, extends(structured_grid) :: curvilinear_grid
    class(structured_grid), allocatable :: metric_tensor(:,:)
  end type

end module
