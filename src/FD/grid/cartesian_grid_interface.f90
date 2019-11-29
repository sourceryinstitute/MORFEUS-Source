!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module cartesian_grid_interface
  !! author: Damian Rouson and Karla Morris
  use structured_grid_interface, only : structured_grid
  implicit none

  private
  public :: cartesian_grid

  type, extends(structured_grid) :: cartesian_grid
  contains
    procedure div_scalar_flux
  end type

  interface
    module function div_scalar_flux( this, vertices, diffusion_coefficient) result(div_flux)
     class(cartesian_grid), intent(in) :: this
     class(structured_grid), intent(in) :: vertices, diffusion_coefficient
     class(structured_grid), allocatable :: div_flux
    end function
  end interface

end module
