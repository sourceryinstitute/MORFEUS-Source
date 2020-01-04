!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module package_interface
  !! author: Damian Rouson and Karla Morris
  !! date: 1/2/2019
  !! Encapsulate halo data to be communicated across structured_grid block boundaries
  use kind_parameters, only : r8k
  implicit none

  private
  public :: package

  integer, parameter :: unset=0

  type package
    !! basic transmission data. extend this type to add coordinate-specific data
    private
    integer :: sender_block_id = unset
    integer :: step = unset
    real(r8k) datum
  contains
    procedure set_sender_block_id
    procedure set_step
    procedure set_datum
    procedure sender_block_id_unset
  end type

  interface

    elemental module subroutine set_sender_block_id(this, sender_block_id)
      !! set sender_bock_id in order to verify sender identity
      implicit none
      class(package), intent(inout) :: this
      integer, intent(in) :: sender_block_id
    end subroutine

    elemental module subroutine set_step(this, step)
      !! set time step in order to verify time step on package data
      implicit none
      class(package), intent(inout) :: this
      integer, intent(in) :: step
    end subroutine

    elemental module subroutine set_datum(this, datum)
      !! set datum to be communicated across structured_grid block internal surfaces
      implicit none
      class(package), intent(inout) :: this
      real(r8k), intent(in) :: datum
    end subroutine

    elemental module function sender_block_id_unset(this) result(is_unset)
      implicit none
      class(package), intent(in) :: this
      logical is_unset
    end function

  end interface

end module package_interface
