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
! $Id: class_iterating.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_iterating
    USE class_psblas, ONLY : psb_dpk_
    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: iterating                           ! Class
    PUBLIC :: create_iterating                    ! Constructor
    PUBLIC :: nmax_, delta_, tol_, current_iteration, & ! Getters
        &    next_iteration, previous_iteration  !    "
    PUBLIC :: increment, reset                    ! Setters
    PUBLIC :: stop_iterating                      ! Utilities
  
    TYPE iterating
        PRIVATE
        INTEGER :: nmax
        INTEGER :: icurrent
        REAL(psb_dpk_) :: delta
        REAL(psb_dpk_) :: tol
    CONTAINS
        PROCEDURE, PRIVATE :: nemo_iterating_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_iterating_sizeof
    END TYPE iterating


    ! ----- Generic Interfaces -----

    INTERFACE

      ! ----- Constructor -----

      MODULE SUBROUTINE create_iterating(iter,input_file,sec,itype)
          TYPE(iterating),  INTENT(OUT) :: iter
          CHARACTER(len=*), INTENT(IN) :: input_file
          CHARACTER(len=*), INTENT(IN) :: sec
          INTEGER         , INTENT(IN) :: itype
      END SUBROUTINE create_iterating


    ! ----- Getters -----

      MODULE FUNCTION nmax_(iter)
        IMPLICIT NONE
        INTEGER :: nmax_
        TYPE(iterating), INTENT(IN) :: iter
      END FUNCTION nmax_


      MODULE FUNCTION delta_(iter)
        IMPLICIT NONE
        REAL(psb_dpk_) :: delta_
        TYPE(iterating), INTENT(IN) :: iter
      END FUNCTION delta_

      MODULE FUNCTION tol_(iter)
        IMPLICIT NONE
        REAL(psb_dpk_) :: tol_
        TYPE(iterating), INTENT(IN) :: iter
      END FUNCTION tol_

      MODULE FUNCTION current_iteration(iter)
        IMPLICIT NONE
        INTEGER :: current_iteration
        TYPE(iterating), INTENT(IN) :: iter
      END FUNCTION current_iteration


      MODULE FUNCTION next_iteration(iter)
        IMPLICIT NONE
        INTEGER :: next_iteration
        TYPE(iterating) :: iter
      END FUNCTION next_iteration

      MODULE FUNCTION previous_iteration(iter)
        IMPLICIT NONE
        INTEGER :: previous_iteration
        TYPE(iterating) :: iter
      END FUNCTION previous_iteration


      ! ----- Utilities -----

      MODULE FUNCTION stop_iterating(iter,eps)
        IMPLICIT NONE
        LOGICAL :: stop_iterating
        TYPE(iterating),  INTENT(IN) :: iter
        REAL(psb_dpk_), INTENT(IN), OPTIONAL :: eps
      END FUNCTION stop_iterating

    END INTERFACE


    ! ----- Setters -----

    INTERFACE increment
      MODULE SUBROUTINE increment_iterating(iter)
        IMPLICIT NONE
        TYPE(iterating), INTENT(INOUT) :: iter
      END SUBROUTINE increment_iterating
    END INTERFACE increment


    INTERFACE reset
      MODULE SUBROUTINE reset_iterating(iter)
        IMPLICIT NONE
        TYPE(iterating), INTENT(INOUT) :: iter
      END SUBROUTINE reset_iterating
    END INTERFACE reset


    INTERFACE
      MODULE FUNCTION nemo_iterating_sizeof(iter)
        USE class_psblas, ONLY : nemo_int_long_
        IMPLICIT NONE
        INTEGER(kind=nemo_int_long_)   :: nemo_iterating_sizeof
        CLASS(iterating), INTENT(IN) :: iter
      END FUNCTION nemo_iterating_sizeof
    END INTERFACE


END MODULE class_iterating
