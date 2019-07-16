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
! $Id$
!
! Description:
!    To be added...
!
SUBMODULE(class_scalar_source) class_scalar_source_procedures

    USE class_psblas
    USE class_dimensions

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_scalar_source_sizeof
        USE psb_base_mod
        USE class_psblas

        nemo_scalar_source_sizeof = 2 * nemo_sizeof_dp + src%dim%nemo_sizeof()

    END PROCEDURE nemo_scalar_source_sizeof

    ! ----- Constructor -----

    MODULE PROCEDURE create_scalar_source
        USE tools_input
        USE json_module, ONLY : json_file

        TYPE(json_file) :: nemo_json
        CHARACTER(LEN=80) :: src_sec
        LOGICAL :: found

        src%dim = dim

        src_sec = 'MORFEUS_FV.Source-terms.'//TRIM(sec)

        CALL open_file(input_file,nemo_json)
        IF(mypnum_() == 0) WRITE(*,*) 'Reading ', TRIM(src_sec), &
            &             ' section from ', TRIM(input_file)

        CALL nemo_json%get(TRIM(src_sec)//'.sc.value', src%sc, found)
        IF (.NOT.found) THEN
            src%sc = 0.0d0
        END IF
        CALL nemo_json%get(TRIM(src_sec)//'.sp.value', src%sc, found)
        IF (.NOT.found) THEN
            src%sp = 0.0d0
        END IF

        IF(mypnum_() == 0) WRITE(*,*)

    END PROCEDURE create_scalar_source


    ! ----- Getter -----

    MODULE PROCEDURE get_scalar_source_dim
        get_scalar_source_dim = src%dim
    END PROCEDURE get_scalar_source_dim


    MODULE PROCEDURE get_scalar_source_sc
        get_scalar_source_sc = src%sc
    END PROCEDURE get_scalar_source_sc


    MODULE PROCEDURE get_scalar_source_sp
        get_scalar_source_sp = src%sp
    END PROCEDURE get_scalar_source_sp

END SUBMODULE class_scalar_source_procedures
