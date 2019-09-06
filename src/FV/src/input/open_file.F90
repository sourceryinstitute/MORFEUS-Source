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
! $Id: open_file.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(tools_input) open_file_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE open_file
        USE class_psblas, ONLY : abort_psblas
        USE json_module, ONLY : json_file
        IMPLICIT NONE
        !! Opens a file

        LOGICAL :: ex
        INTEGER :: path_length, var_status, inp
        CHARACTER(LEN=:), ALLOCATABLE :: dirname
#ifdef WIN32
        CHARACTER(LEN=*), PARAMETER :: cwd_var = 'CD'
#else
        CHARACTER(LEN=*), PARAMETER :: cwd_var = 'PWD'
#endif

        INQUIRE(file=input_file,exist=ex,number=inp)
        IF(.NOT. ex) THEN
            CALL get_environment_variable (cwd_var, length=path_length, status=var_status)
            SELECT CASE (var_status)
            CASE (1)
                !! Undefined
                ERROR STOP 'Undefined cwd_var in open_file_implementation, open_file.F90'
            CASE (0)
                !! Good
                ALLOCATE(CHARACTER(LEN=path_length)::dirname)
                CALL get_environment_variable (cwd_var, value=dirname)
            CASE DEFAULT
                !! Something went wrong
                ERROR STOP 'Environment variable does not exist or processor does not support environment variables, open_file.F90'
            END SELECT

            WRITE(*,100) TRIM(input_file), dirname

            CALL abort_psblas
        END IF

        IF(inp == -1) THEN
            CALL nemo_json%initialize()
            CALL nemo_json%load_file(filename=input_file)
        END IF

100     FORMAT(' ERROR! Input file "',a,'" doesn''t exist',/,' in directory ',a)

        END PROCEDURE open_file

END SUBMODULE open_file_implementation
