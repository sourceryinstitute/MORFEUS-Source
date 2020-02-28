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
! $Id: get_par_l.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE (tools_input) get_par_l_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE get_par_l
            USE class_psblas
            IMPLICIT NONE
            !
            LOGICAL, PARAMETER :: debug = .FALSE.
            !
            LOGICAL :: found
            CHARACTER(len=15) :: str
            INTEGER :: i, k

            ! File pointer is supposed to at the beginning of SEC section.

            found = .FALSE.

            k = 0
            reading: DO
                READ(inp,'(a)') str
                k = k + 1
                IF(str == par) THEN
                    BACKSPACE(inp)
                    READ(inp,*) str, get_par_l
                    IF(debug) WRITE(*,100) str, get_par_l
                    found = .TRUE.
                    EXIT reading
                ELSEIF(str == 'END OF SECTION') THEN
                    EXIT reading
                END IF
            END DO reading

            ! Rewinds the section
            DO i = 1, k
                BACKSPACE(inp)
            END DO

            IF(found) RETURN

            ! Parameter not found in input file.
            WRITE(*,200) TRIM(par), TRIM(sec), default
            get_par_l = default

            ! REMARK. A logical variable can assume only two values. No PAR_L_MANDATORY_
            ! named constant has been set. Using such a parameter could conflict with
            ! the user supplied DEFAULT value.

100         FORMAT(1x,a15,1x,l5)
200         FORMAT(' WARNING! Parameter "',a,'" in section ',a,'.',&
                &    ' Set to default = ',l5)

        END PROCEDURE get_par_l

END SUBMODULE get_par_l_implementation
