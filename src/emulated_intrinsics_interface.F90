!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
module emulated_intrinsics_interface
  !! author: Damian Rouson
  !!
  !! Fortran 2008 coarray emulations of Fortran 2018 intrinsic collective subroutines
  !! and Fortran 2003 emulation of Fortran 2008 intrinsic procedures (e.g, findloc)
  implicit none

#ifndef HAVE_FINDLOC
  interface findloc
    !! result is the last occurrence of a value in an array or zero if not found
    module procedure findloc_integer_dim1, findloc_logical_dim1, findloc_character_dim1
  end interface
#endif

#ifndef HAVE_COLLECTIVE_SUBROUTINES
  interface co_sum
    !! parallel computation of the sum of the first argument
    module procedure co_sum_integer
  end interface

  interface co_broadcast
    !! parallel one-to-all communication of the value of first argument
    module procedure co_broadcast_integer
  end interface
#endif

  interface

#ifndef HAVE_COLLECTIVE_SUBROUTINES
    module subroutine co_sum_integer(a,result_image,stat,errmsg)
      implicit none
      integer, intent(inout) :: a
      integer, intent(in), optional :: result_image
      integer, intent(out), optional ::  stat
      character(len=*), intent(inout), optional :: errmsg
    end subroutine

    module subroutine co_broadcast_integer(a,source_image,stat,errmsg)
      implicit none
      integer, intent(inout) :: a
      integer, intent(in) :: source_image
      integer, intent(out), optional ::  stat
      character(len=*), intent(inout), optional :: errmsg
    end subroutine
#endif

#ifndef HAVE_FINDLOC
    pure module function findloc_integer_dim1(array, value, dim, back) result(location)
      implicit none
      integer, intent(in) :: array(:), value, dim
      logical, intent(in), optional :: back
      integer location
    end function

    pure module function findloc_logical_dim1(array, value, dim, back) result(location)
      implicit none
      logical, intent(in) :: array(:), value, back
      integer, intent(in) :: dim
      integer location
    end function

    pure module function findloc_character_dim1(array, value, dim, back) result(location)
      implicit none
      character(len=*), intent(in) :: array(:), value
      integer, intent(in) :: dim
      logical, intent(in) :: back
      integer location
    end function
#endif

  end interface

end module emulated_intrinsics_interface
