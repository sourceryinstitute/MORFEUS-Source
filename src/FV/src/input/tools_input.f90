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
! $Id: tools_input.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    To be added...
!
MODULE tools_input

    USE class_vector, ONLY : vector
    USE class_psblas, ONLY : psb_dpk_
    IMPLICIT NONE

    PUBLIC
    PRIVATE :: vector, psb_dpk_

    INTERFACE get_par
        MODULE FUNCTION get_par_l(inp,sec,par,default)
            IMPLICIT NONE
            LOGICAL :: get_par_l
            INTEGER, INTENT(IN) :: inp
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            LOGICAL, INTENT(IN) :: default
        END FUNCTION get_par_l

        MODULE FUNCTION get_par_i(inp,sec,par,default)
            IMPLICIT NONE
            INTEGER :: get_par_i
            INTEGER, INTENT(IN) :: inp
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            INTEGER, INTENT(IN) :: default
        END FUNCTION get_par_i

        MODULE FUNCTION get_par_d(inp,sec,par,default)
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            REAL(psb_dpk_) :: get_par_d
            INTEGER, INTENT(IN) :: inp
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            REAL(psb_dpk_), INTENT(IN) :: default
        END FUNCTION get_par_d

        MODULE FUNCTION get_par_h(inp,sec,par,default)
            IMPLICIT NONE
            CHARACTER(len=32) :: get_par_h
            INTEGER, INTENT(IN) :: inp
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            CHARACTER(len=*), INTENT(IN) :: default
        END FUNCTION get_par_h

        MODULE FUNCTION get_par_v(inp,sec,par,default)
            USE class_vector
            IMPLICIT NONE
            TYPE(vector) :: get_par_v
            INTEGER,          INTENT(IN) :: inp
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            TYPE(vector),     INTENT(IN) :: default
        END FUNCTION get_par_v
    END INTERFACE get_par


    INTERFACE read_par
        MODULE FUNCTION read_par_l(input_file,sec,par,default)
            IMPLICIT NONE
            LOGICAL :: read_par_l
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            LOGICAL, INTENT(IN) :: default
        END FUNCTION read_par_l

        MODULE FUNCTION read_par_i(input_file,sec,par,default)
            IMPLICIT NONE
            INTEGER :: read_par_i
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            INTEGER, INTENT(IN) :: default
        END FUNCTION read_par_i
     
        MODULE FUNCTION read_par_d(input_file,sec,par,default)
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            REAL(psb_dpk_) :: read_par_d
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            REAL(psb_dpk_), INTENT(IN) :: default
        END FUNCTION read_par_d

        MODULE FUNCTION read_par_h(input_file,sec,par,default)
            IMPLICIT NONE
            CHARACTER(len=32) :: read_par_h
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            CHARACTER(len=*), INTENT(IN) :: default
        END FUNCTION read_par_h

        MODULE FUNCTION read_par_v(input_file,sec,par,default)
            USE class_vector
            IMPLICIT NONE
            TYPE(vector) :: read_par_v
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(IN) :: par
            TYPE(vector),     INTENT(IN) :: default
        END FUNCTION read_par_v
    END INTERFACE read_par


    INTERFACE
        MODULE SUBROUTINE find_section(sec,inp)
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN) :: sec
            INTEGER, INTENT(IN) :: inp
        END SUBROUTINE find_section
    END INTERFACE


    INTERFACE
        MODULE SUBROUTINE open_file(input_file,inp)
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN) :: input_file
            INTEGER, INTENT(OUT) :: inp
        END SUBROUTINE open_file
    END INTERFACE


    ! ----- Named Constants -----

    LOGICAL,           PARAMETER :: mandatory_l_ = .FALSE.
    INTEGER,           PARAMETER :: mandatory_i_ = -999
    REAL(psb_dpk_),    PARAMETER :: mandatory_d_ = -9999.d0
    CHARACTER(len=10), PARAMETER :: mandatory_h_ = 'mandatory'

    ! MANDATORY_V_ is defined inside the VECTOR class

END MODULE tools_input
