!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
MODULE kind_parameters
    USE iso_fortran_env, ONLY : i4k => int32, i8k => int64, r4k => real32, r8k =>real64
    IMPLICIT NONE
    !! author: Damian Rouson
    !! date: 09/27/2019
    !!
    !! This module contains the kinds used for specifying the precision of variables
    !!
    PRIVATE
    PUBLIC :: i4k, i8k, r4k, r8k

    ! r4k - Single precision for reals
    ! r8k - Double precision for reals
    ! i4k - Single precision for integers
    ! i8k - Double precision for integers

END MODULE kind_parameters
