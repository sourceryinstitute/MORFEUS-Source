!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module block_metadata_interface
  !! author: Damian Rouson
  !! date: August 2, 2019
  !! summary: encapsulate metadata describing structured-grid blocks
  use iso_c_binding, only : c_int
  use Kinds, only : r8k
  implicit none

  private
  public :: block_metadata, tag_kind, untagged, lower, upper

  integer, parameter :: space_dimension=3, num_end_points=2, tag_kind=c_int, lower=1, upper=2

  enum, bind(C)
    enumerator :: untagged = -huge(1_tag_kind)
  end enum

  type block_metadata
    !! structured-grid block descriptor
    private
    real(r8k) subdomain_(space_dimension,num_end_points)
    real(r8k) max_spacing_
    integer(tag_kind) :: tag_ = untagged
    character(len=len('unlabeled')) :: label_='unlabeled'
  contains
    procedure set_tag
    procedure set_label
    procedure set_subdomain
    procedure set_max_spacing
    procedure get_tag
    procedure get_label
    procedure get_subdomain
    procedure get_max_spacing
  end type

  interface

    elemental module subroutine set_tag( this, tag )
      !! Define the identification tag of a this block_metadata object
      implicit none
      class(block_metadata), intent(inout) :: this
      integer, intent(in) :: tag
    end subroutine

    elemental module subroutine set_label( this, label )
      !! Define the label of a this block_metadata object
      implicit none
      class(block_metadata), intent(inout) :: this
      character(len=*), intent(in) :: label
    end subroutine

#ifdef HAVE_ERROR_STOP_IN_PURE
    pure &
#endif
    module subroutine set_subdomain(this, subdomain)
      !! Define the end point of a block-structured grid coordinate direction
      implicit none
      class(block_metadata), intent(inout) :: this
      real(r8k), dimension(:,:), intent(in) :: subdomain
    end subroutine

#ifdef HAVE_ERROR_STOP_IN_PURE
    pure &
#endif
    module subroutine set_max_spacing(this, max_spacing)
      !! Define the maximum allowable grid spacing
      implicit none
      class(block_metadata), intent(inout) :: this
      real(r8k), intent(in) :: max_spacing
    end subroutine

    elemental module function get_tag( this ) result(this_tag)
      !! Result is the identification tag of a this block_metadata object
      implicit none
      class(block_metadata), intent(in) :: this
      integer(tag_kind) this_tag
    end function

    pure module function get_label( this ) result(this_label)
      !! Result is the label of a this block_metadata object
      implicit none
      class(block_metadata), intent(in) :: this
      character(len=:), allocatable :: this_label
    end function

    pure module function get_subdomain( this ) result(this_subdomain)
      !! Result contains the coordinate intervals delimiting this block_metadata object
      implicit none
      class(block_metadata), intent(in) :: this
      real(r8k), dimension(space_dimension, num_end_points) :: this_subdomain
    end function

    pure module function get_max_spacing( this ) result(this_max_spacing)
      !! Result contains the maximum allowable grid spacing for this block_metadata object
      implicit none
      class(block_metadata), intent(in) :: this
      real this_max_spacing
    end function

  end interface

end module block_metadata_interface
