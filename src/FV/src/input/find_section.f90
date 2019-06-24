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
! $Id: find_section.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
MODULE PROCEDURE find_section
    USE class_psblas

    IMPLICIT NONE
    !
    CHARACTER(len=80) :: input_file
    CHARACTER(len=80) :: str
    INTEGER :: i, ir, n


    INQUIRE(unit=inp,name=input_file)

    REWIND(inp)

    i = 0
    n = 0
    findSec: DO
        READ(inp,'(a)',END=999) str
        i = i + 1
        IF(str == sec) THEN
            n = n + 1
            IF(n == 1) ir = i
        END IF
    END DO findSec

999 CONTINUE

    IF(n == 0) THEN
        WRITE(*,100) TRIM(sec), TRIM(input_file)
        CALL abort_psblas
    ELSEIF(n > 1) THEN
        WRITE(*,200) TRIM(sec), TRIM(input_file)
        CALL abort_psblas
    END IF

    REWIND(inp)

    ! Repositions file pointer at section beginning
    DO i = 1, ir
        READ(inp,'()')
    END DO

100 FORMAT(' ERROR! "',a,'" section not found in "',a,'" file.')
200 FORMAT(' ERROR! Multiple ',a,' sections found in ',a,' file.')

END PROCEDURE find_section
