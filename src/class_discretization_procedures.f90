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
! $Id: $
!
! Description:
!    To be added...
!
SUBMODULE(class_discretization) class_discretization_procedures

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_discretization_sizeof
        USE psb_base_mod
        USE class_psblas

        nemo_discretization_sizeof = nemo_sizeof_int + nemo_sizeof_dp

    END PROCEDURE nemo_discretization_sizeof

    ! ----- Constructors -----

    MODULE PROCEDURE read_par_discretization
        USE class_psblas
        USE json_module, ONLY : json_file
        USE tools_input

        LOGICAL :: found
        TYPE(json_file) :: nemo_json
        INTEGER :: icontxt, mypnum
        INTEGER :: inp
        CHARACTER(len=15) :: par_

        icontxt = icontxt_()
        mypnum  = mypnum_()

        found = .FALSE.

        IF(mypnum == 0) THEN

            CALL open_file(input_file,nemo_json)
            CALL find_section(sec,nemo_json)

            !       read(inp,'()')

            seek_par: DO
                READ(inp,'(a)') par_

                !          backspace(inp)

                IF(TRIM(par_) == TRIM(par)) THEN
                    BACKSPACE(inp)
                    found = .TRUE.

                    READ(inp,*) par_, r%id
                    SELECT CASE(r%id)
                    CASE(id_cd_, id_up_)
                        !                read(inp,'()')
                        r%blend = 1.d0
!!$             case(id_tvd_) ! Any scheme requiring a blending factor too
!!$                read(inp,*) r%blend
                    CASE DEFAULT
                        WRITE(*,200)
                        CALL abort_psblas
                    END SELECT

                ELSEIF(TRIM(par_) == 'END OF SECTION') THEN
                    EXIT seek_par
                ELSE
                    READ(inp,'()')
                END IF

            END DO seek_par

            CLOSE(inp)
        END IF

        CALL psb_bcast(icontxt,found)

        IF(found) THEN
            CALL psb_bcast(icontxt,r%id)
            CALL psb_bcast(icontxt,r%blend)
        ELSE
            r = default
        END IF

100     FORMAT(a,i2)
200     FORMAT(' ERROR! Unsupported differencing scheme in',&
            & ' READ_PAR_DISCRETIZATION')

    END PROCEDURE read_par_discretization


    ! ----- Getters -----

    MODULE PROCEDURE id_

        id_ = ds%id

    END PROCEDURE id_


END SUBMODULE class_discretization_procedures
