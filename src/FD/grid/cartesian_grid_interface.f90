!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module cartesian_grid_interface
  !! author: Damian Rouson and Karla Morris
  use structured_grid_interface, only : structured_grid
  use differentiable_field_interface, only : differentiable_field
  implicit none

  private
  public :: cartesian_grid

  type, extends(structured_grid) :: cartesian_grid
  contains
    procedure div_scalar_flux
    procedure assign_structured_grid
    procedure block_indicial_coordinates
    procedure block_identifier
  end type

  interface

    module function div_scalar_flux(this, vertices, exact_result) result(div_flux)
     implicit none
     class(cartesian_grid), intent(in) :: this
     class(structured_grid), intent(in) :: vertices
     class(differentiable_field), intent(in), optional :: exact_result
     class(structured_grid), allocatable :: div_flux
    end function

    module subroutine assign_structured_grid(this, rhs)
     implicit none
     class(cartesian_grid), intent(inout) :: this
     class(structured_grid), intent(in) :: rhs
    end subroutine

    pure module function block_indicial_coordinates(this,n) result(ijk)
      !! calculate the 3D location of the block that has the provided 1D block identifer
      implicit none
      class(cartesian_grid), intent(in) :: this
      integer, intent(in) :: n
      integer, dimension(:), allocatable :: ijk(:)
    end function

    pure module function block_identifier(this,ijk) result(n)
      !! calculate the 1D block identifer associated with the provided 3D block location
      implicit none
      class(cartesian_grid), intent(in) :: this
      integer, intent(in), dimension(:) :: ijk
      integer :: n
    end function

  end interface

end module
