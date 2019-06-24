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
SUBMODULE(class_iterating) class_iterating_procedures

    IMPLICIT NONE

CONTAINS

    ! ----- Constructor -----

    MODULE PROCEDURE create_iterating
        USE class_psblas, ONLY : abort_psblas
        USE tools_input, ONLY : mandatory_d_, mandatory_i_, read_par
        USE tools_math, ONLY : it_counter_, it_convergence_, it_time_

        ! Sets common members
        iter%icurrent = 0
        iter%nmax     = read_par(input_file,sec,'nmax',mandatory_i_)

        ! Sets specific members
        SELECT CASE(itype)
        CASE(it_time_)
            iter%delta = read_par(input_file,sec,'delta',mandatory_d_)
            iter%tol   = 0.d0
        CASE(it_convergence_)
            iter%delta = 0.d0
            iter%tol   = read_par(input_file,sec,'tol',mandatory_d_)
        CASE(it_counter_)
            iter%delta = 0.d0
            iter%tol   = 0.d0
        CASE default
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

100     FORMAT(' ERROR! Illegal ITYPE value in CREATE_ITERATING')

    END PROCEDURE create_iterating


    ! ----- Getters -----

    MODULE PROCEDURE nmax_

        nmax_ = iter%nmax

    END PROCEDURE nmax_


    MODULE PROCEDURE delta_

        delta_ = iter%delta

    END PROCEDURE delta_

    MODULE PROCEDURE tol_

        tol_ = iter%tol

    END PROCEDURE tol_

    MODULE PROCEDURE current_iteration

        current_iteration = iter%icurrent

    END PROCEDURE current_iteration


    MODULE PROCEDURE next_iteration

        next_iteration = iter%icurrent + 1

    END PROCEDURE next_iteration


    MODULE PROCEDURE previous_iteration

        previous_iteration = iter%icurrent - 1

    END PROCEDURE previous_iteration


    ! ----- Setters -----

    MODULE PROCEDURE increment_iterating

        iter%icurrent = iter%icurrent + 1

    END PROCEDURE increment_iterating


    MODULE PROCEDURE reset_iterating

        iter%icurrent = 0

    END PROCEDURE reset_iterating


    ! ----- Utilities -----

    MODULE PROCEDURE stop_iterating

        ! Initialization
        stop_iterating = .FALSE.

        IF(iter%icurrent >= iter%nmax) THEN
            stop_iterating = .TRUE.
            RETURN
        END IF

        IF(PRESENT(eps) .AND. eps <= iter%tol) THEN
            stop_iterating = .TRUE.
        END IF

    END PROCEDURE stop_iterating

    MODULE PROCEDURE nemo_iterating_sizeof
        USE class_psblas, ONLY : nemo_int_long_, nemo_sizeof_dp, nemo_sizeof_int

        nemo_iterating_sizeof =  2 * nemo_sizeof_int + 2* nemo_sizeof_dp

    END PROCEDURE nemo_iterating_sizeof

END SUBMODULE class_iterating_procedures
