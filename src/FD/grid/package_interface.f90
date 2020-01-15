!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module package_interface
  !! author: Damian Rouson and Karla Morris
  !! date: 1/2/2019
  !! Encapsulate halo data to be communicated across structured_grid block boundaries
  use kind_parameters, only : r8k
  implicit none

  private
  public :: package, null_sender_id

  integer, parameter :: null_sender_id=-1

  type package
    !! basic transmission data. extend this type to add coordinate-specific data
    private
    integer :: sender_block_id
    integer :: step
    real(r8k) datum
  contains
    procedure assign
    generic :: assignment(=)=>assign
    procedure get_sender_block_id
    procedure set_sender_block_id
    procedure set_step
    procedure set_datum
    procedure sender_block_id_null
  end type

  interface

    module subroutine assign(this, rhs)
      !! copy rhs package components into this package
      implicit none
      class(package), intent(inout) :: this
      type(package), intent(in) :: rhs
    end subroutine

    module function get_sender_block_id(this) result(this_sender_block_id)
      !! result is sender_block_id for this package
      implicit none
      class(package), intent(in) :: this
      integer :: this_sender_block_id
    end function

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

    elemental module function sender_block_id_null(this) result(is_null)
      implicit none
      class(package), intent(in) :: this
      logical is_null
    end function

  end interface

end module package_interface
