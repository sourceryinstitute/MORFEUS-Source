!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module inbox_interface
  !! author: Damian Rouson
  !! date:
  !! Encapsulate information and procedures for structured_grid block halo exchanges
  use iso_fortran_env, only : event_type
  implicit none

  private
  public :: inbox
  public :: face_index, x_forward, x_backward, y_forward, y_backward, z_forward, z_backward

  enum, bind(C)
    enumerator :: x_forward=1, x_backward, y_forward, y_backward, z_forward, z_backward
  end enum

  integer, parameter, dimension(*) :: face_index = [x_forward, x_backward, y_forward, y_backward, z_forward, z_backward]

  type inbox
    private
    type(event_type) message_ready, message_read
    integer sender_block_id, step
  contains
    procedure  set_sender_block_id
    procedure  set_step
  end type

  interface

    pure module subroutine set_sender_block_id(this, id)
      !! set sender_block_id component
      implicit none
      class(inbox), intent(inout) :: this
      integer, intent(in) :: id
    end subroutine

    pure module subroutine set_step(this, step)
      !! set step component
      implicit none
      class(inbox), intent(inout) :: this
      integer, intent(in) :: step
    end subroutine

  end interface

end module
