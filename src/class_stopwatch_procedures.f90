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
SUBMODULE(class_stopwatch) class_stopwatch_procedures

    USE psb_base_mod
    USE tools_psblas, ONLY : psb_dpk_
    IMPLICIT NONE

CONTAINS

    ! ----- Constructor -----

    MODULE PROCEDURE stopwatch_

        stopwatch_%icontxt = icontxt
        stopwatch_%partial = 0.d0
        stopwatch_%split   = 0.d0
        stopwatch_%total   = 0.d0

    END PROCEDURE stopwatch_


    ! ----- Getters -----

    MODULE PROCEDURE partial_
        partial_ = sw%partial
    END PROCEDURE partial_


    MODULE PROCEDURE total_
        total_ = sw%total
    END PROCEDURE total_


    ! ----- Setters -----

    MODULE PROCEDURE reset_stopwatch

        sw%split = 0.d0
        sw%partial = 0.d0

    END PROCEDURE reset_stopwatch


    MODULE PROCEDURE tic_stopwatch

        sw%split   = psb_wtime()
        sw%partial = 0.d0

    END PROCEDURE tic_stopwatch


    MODULE PROCEDURE toc_stopwatch

        sw%partial = psb_wtime() - sw%split
        sw%total = sw%total + sw%partial

    END PROCEDURE toc_stopwatch


    MODULE PROCEDURE synchro_stopwatch

        CALL psb_amx(sw%icontxt,sw%partial)
        CALL psb_amx(sw%icontxt,sw%total)

    END PROCEDURE synchro_stopwatch


END SUBMODULE class_stopwatch_procedures
