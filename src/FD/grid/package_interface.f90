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
  implicit none

  private
  public :: package

  type package
    !! basic transmission data. extend this type to add coordinate-specific data
    private
    integer sender_block_id, step
  contains
    procedure set_sender_block_id
    procedure set_step
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

  end interface

end module package_interface
