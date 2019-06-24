!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
!
!    NEMO - Numerical Engine (for) Multiphysics Operators
! Copyright (c) 2007, Stefano Toninel
!                     Gian Marco Bianchi  University of Bologna
!              David P. Schmidt    University of Massachusetts - Amherst
!              Salvatore Filippone University of Rome Tor Vergata
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!
!     1. Redistributions of source code must retain the above copyright notice,
!        this list of conditions and the following disclaimer.
!     2. Redistributions in binary form must reproduce the above copyright notice,
!        this list of conditions and the following disclaimer in the documentation
!        and/or other materials provided with the distribution.
!     3. Neither the name of the NEMO project nor the names of its contributors
!        may be used to endorse or promote products derived from this software
!        without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!---------------------------------------------------------------------------------
!
! $Id: renum.F90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Renumbering
!
MODULE renum

    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: start_renum, stop_renum, apply_renum, print_renum
    PUBLIC :: to_a_, to_b_, to_a_and_b_


    ! Direct and Inverse permutation array (public)
    INTEGER, ALLOCATABLE, SAVE :: perm(:), pinv(:)

    ! ----- Named Constants -----
    INTEGER, PARAMETER :: to_a_       = 1
    INTEGER, PARAMETER :: to_b_       = 2
    INTEGER, PARAMETER :: to_a_and_b_ = 3

 INTERFACE

    ! ----- Constructor -----

    MODULE SUBROUTINE start_renum(irenum,c2c)
        USE class_connectivity
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: irenum
        TYPE(connectivity), INTENT(IN) :: c2c
    END SUBROUTINE start_renum


    ! ----- Destuctor -----

    MODULE SUBROUTINE stop_renum
        IMPLICIT NONE
    END SUBROUTINE stop_renum


    ! ----- Broadcast -----

    MODULE SUBROUTINE print_renum(iout)
        IMPLICIT NONE
        INTEGER :: iout
    END SUBROUTINE print_renum
    ! ----- Broadcast -----

    MODULE SUBROUTINE build_pinv
        IMPLICIT NONE
    END SUBROUTINE build_pinv


    ! ----- Interfaces To Renumbering Methods -----

    ! Gibbs-Poole-Stockmeyer renumbering algorithm
    MODULE SUBROUTINE cmp_gps(c2c)
        USE class_connectivity, ONLY : connectivity
        IMPLICIT NONE
        TYPE(connectivity) :: c2c
    END SUBROUTINE cmp_gps

  END INTERFACE

  ! ----- Routines For The Application Of The Renumbering -----

  ! ----- Generic Interfaces -----

  INTERFACE apply_renum

    ! Integer array
    MODULE SUBROUTINE apply_renum_array(a)
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: a(:)
    END SUBROUTINE apply_renum_array

    ! CELL objects array
    MODULE SUBROUTINE apply_renum_cell(c)
        USE class_cell, ONLY : cell
        IMPLICIT NONE
        TYPE(cell), INTENT(INOUT) :: c(:)
    END SUBROUTINE apply_renum_cell

    ! FACE objects array
    MODULE SUBROUTINE apply_renum_face(f)
        USE class_face, ONLY : face
        IMPLICIT NONE
        TYPE(face), INTENT(INOUT) :: f(:)
    END SUBROUTINE apply_renum_face

    ! CONNECTIVITY object
    MODULE SUBROUTINE apply_renum_conn(a2b,apply)
        USE class_connectivity, ONLY : connectivity
        IMPLICIT NONE
        TYPE(connectivity), INTENT(INOUT) :: a2b
        INTEGER, INTENT(IN) :: apply
    END SUBROUTINE apply_renum_conn

  END INTERFACE apply_renum

END MODULE renum
