!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module ellipsoidal_field_interface
  !! author: Damian Rouson
  !! date: 12/19/2019
  !!
  !! Define a 3D scalar field with ellipsoidal isosurfaces and provide differential operators
  use differentiable_field_interface, only : differentiable_field
  use grid_interface, only : grid
  implicit none

  private
  public :: ellipsoidal_field

  type, extends(differentiable_field) :: ellipsoidal_field
  contains
    procedure evaluate
    procedure laplacian
  end type

  interface

    module function evaluate(this, grid_points) result(f)
      !! Evaluate the function at the provided grid points
      implicit none
      class(ellipsoidal_field), intent(in)  :: this
      class(grid), intent(in)  :: grid_points
      class(grid), allocatable :: f
    end function

    module function laplacian(this, grid_points) result(laplacian_f)
      !! Compute the Laplacian of the ellipsoidal function employed in "evaluate" above
      implicit none
      class(ellipsoidal_field), intent(in)  :: this
      class(grid), intent(in)  :: grid_points
      class(grid), allocatable :: laplacian_f
    end function

  end interface

end module ellipsoidal_field_interface
