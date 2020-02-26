!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module spherical_1D_solver_interface
  !! author: Xiaofeng Xu and Damian Rouson
  !! date: 2/24/2020
  !!
  !! Solve the 1D heat equation in spherically symmetric radial coordinates
  use kind_parameters, only : r8k, i4k
  implicit none

  private
  public :: grid_block

  type grid_block
    !! encapsulate all grid data
    private
    real(r8k), allocatable :: v(:,:)          !! v(:,1) = r, v(:,2) = T, shape = [nr,2]
    real(r8k), allocatable :: rho(:), cp(:)   !! density and specific heat (size = nr)
    real(r8k), allocatable :: T_analytical(:) !! expected solution (size = nr)
  contains
    procedure :: set_v
    procedure :: set_material_properties_size
    procedure :: set_expected_solution_size
    procedure :: set_rho
    procedure :: set_cp
    procedure :: time_advance_heat_equation
  end type grid_block

  interface

    module subroutine set_v( this, nr, constants )
      implicit none
      class(grid_block), intent(inout) :: this
      integer, intent(in) :: nr
      real(r8k), intent(in) :: constants(:)
    end subroutine

    module subroutine set_material_properties_size(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine set_expected_solution_size(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine set_rho(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine set_cp(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine time_advance_heat_equation(this, duration)
      implicit none
      class(grid_block), intent(inout) :: this
      real(r8k), intent(in) :: duration
    end subroutine

  end interface

end module spherical_1D_solver_interface
