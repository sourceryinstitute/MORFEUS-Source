!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module co_object_interface
  implicit none

  private
  public :: co_object

  ! Define an abstract parent type to ensure basic functionality expected to be provided by all non-abstract types.
  ! Each non-abstract type provides the functionality by extending this type and implementing its deferred binding(s).  This
  ! type resembles java's Object class in the sense that it is intended to be the ultimate ancester of every other type.
  type, abstract :: co_object
    private
    logical :: defined=.false.
      !! Default initialization indicates not yet user-defined
    logical, allocatable :: facilitate_type_extension[:]
  contains
    procedure :: mark_as_defined
    procedure :: user_defined
  end type

  interface

    pure module subroutine mark_as_defined(this)
      !! Mark the co_object as user-defined
      implicit none
      class(co_object), intent(inout) :: this
    end subroutine

    pure module function user_defined(this) result(is_defined)
      !! Return a boolean result indicating whether this co_object has been initialized since its declaration
      implicit none
      class(co_object), intent(in) :: this
      logical :: is_defined
    end function

  end interface

end module co_object_interface
