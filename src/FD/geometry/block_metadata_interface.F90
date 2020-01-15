!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module block_metadata_interface
  !! author: Damian Rouson
  !! date: 8/2/2019
  !! summary: encapsulate metadata describing structured-grid blocks
  use iso_c_binding, only : c_int
  use kind_parameters, only : r8k
  implicit none

  private
  public :: block_metadata, tag_kind, untagged, lower, upper, subdomain_t, max_name_length, space_dimension, num_end_points

  integer, parameter :: space_dimension=3, num_end_points=2, tag_kind=c_int, lower=1, upper=2, max_name_length=32

  enum, bind(C)
    enumerator :: untagged = -huge(1_tag_kind)
  end enum

  type subdomain_t !! scalar argument for elemental set_subdomain procedure
    real(r8k), dimension(space_dimension,num_end_points) :: edges
  end type

  type block_metadata
    !! structured-grid block descriptor
    private
    type(subdomain_t) subdomain
    real(r8k) max_spacing_
    integer(tag_kind) :: tag_ = untagged
    character(len=max_name_length) :: label_='unlabeled'
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

    pure module subroutine set_label( this, label )
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
      type(subdomain_t), intent(in) :: subdomain
    end subroutine

#ifndef HAVE_ERROR_STOP_IN_PURE
    !impure &
#endif
    elemental module subroutine set_max_spacing(this, max_spacing)
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

    pure module function get_subdomain( this ) result(edges)
      !! Result contains the coordinate intervals delimiting this block_metadata object
      implicit none
      class(block_metadata), intent(in) :: this
      real(r8k), dimension(space_dimension,num_end_points) :: edges
    end function

    pure module function get_max_spacing( this ) result(this_max_spacing)
      !! Result contains the maximum allowable grid spacing for this block_metadata object
      implicit none
      class(block_metadata), intent(in) :: this
      real this_max_spacing
    end function

  end interface

end module block_metadata_interface
