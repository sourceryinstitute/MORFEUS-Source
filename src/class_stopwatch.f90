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
! $Id: class_stopwatch.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    Defines a parallel stopwatch for timing (``tic-toc'') a code section.
!
MODULE class_stopwatch

    USE psb_base_mod
    USE tools_psblas, ONLY : psb_dpk_
    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: stopwatch                ! Class
    PUBLIC :: stopwatch_               ! Constructor

    TYPE stopwatch
        PRIVATE
        INTEGER :: icontxt
        REAL(psb_dpk_) :: total
        REAL(psb_dpk_) :: split
        REAL(psb_dpk_) :: partial
    CONTAINS
        PROCEDURE :: partial_, total_   ! Getters
        PROCEDURE, PRIVATE :: synchro_stopwatch, tic_stopwatch, toc_stopwatch  !Setters
        GENERIC, PUBLIC :: synchro => synchro_stopwatch
        GENERIC, PUBLIC :: tic => tic_stopwatch
        GENERIC, PUBLIC :: toc => toc_stopwatch
        PROCEDURE, PRIVATE :: reset_stopwatch
        GENERIC, PUBLIC:: reset => reset_stopwatch
    END TYPE stopwatch


    ! ----- Generic Interfaces -----

    INTERFACE
        MODULE SUBROUTINE reset_stopwatch(sw)
            IMPLICIT NONE
            CLASS(stopwatch), INTENT(INOUT) :: sw
        END SUBROUTINE reset_stopwatch

        MODULE SUBROUTINE tic_stopwatch(sw)
            IMPLICIT NONE
            CLASS(stopwatch), INTENT(INOUT) :: sw
        END SUBROUTINE tic_stopwatch

        MODULE SUBROUTINE toc_stopwatch(sw)
            IMPLICIT NONE
            CLASS(stopwatch), INTENT(INOUT) :: sw
        END SUBROUTINE toc_stopwatch

        MODULE SUBROUTINE synchro_stopwatch(sw)
            IMPLICIT NONE
            CLASS(stopwatch), INTENT(INOUT) :: sw
        END SUBROUTINE synchro_stopwatch

    ! ----- Constructor -----

        MODULE FUNCTION stopwatch_(icontxt)
            IMPLICIT NONE
            TYPE(stopwatch) :: stopwatch_
            INTEGER, INTENT(IN) :: icontxt
        END FUNCTION stopwatch_

    ! ----- Getters -----
        MODULE FUNCTION partial_(sw)
            IMPLICIT NONE
            REAL(psb_dpk_) :: partial_
            CLASS(stopwatch), INTENT(IN) :: sw
        END FUNCTION partial_

        MODULE FUNCTION total_(sw)
            IMPLICIT NONE
            REAL(psb_dpk_) :: total_
            CLASS(stopwatch), INTENT(IN) :: sw
        END FUNCTION total_
    END INTERFACE
END MODULE class_stopwatch
