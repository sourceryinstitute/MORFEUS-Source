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
! $Id: class_output.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_output) class_output_procedures
    use class_iterating
    IMPLICIT NONE

CONTAINS


    MODULE PROCEDURE nemo_output_sizeof
        USE class_psblas, ONLY : nemo_sizeof_int

        nemo_output_sizeof = nemo_sizeof_int + LEN(obj%basepath)  + LEN(obj%path)

    END PROCEDURE nemo_output_sizeof

    ! ----- Constructor -----

    MODULE PROCEDURE create_output
        USE tools_input
        USE tools_output_basics

        CHARACTER(len=32) :: basepath_

        ! Gets format
        out%fmt  = read_par(input_file,sec,'format',mandatory_i_)

        ! Set default basepath
        SELECT CASE(out%fmt)
        CASE(vtk_)
            basepath_ = 'out_nemo'
        END SELECT

        ! Gets path basename
        out%basepath = TRIM(read_par(input_file,sec,'basepath',basepath_))

        ! Sets initial path
        out%path = out%basepath

    END PROCEDURE create_output


    ! ----- Getters -----

    MODULE PROCEDURE fmt_

        fmt_ = out%fmt

    END PROCEDURE fmt_


    MODULE PROCEDURE path_
        USE tools_output_basics
        !
        INTEGER :: l
        CHARACTER(len=1) :: path_end

        SELECT CASE(out%fmt)
        CASE(vtk_)
            path_ = TRIM(out%path) // '.vtk'

        END SELECT

    END PROCEDURE path_


    ! ----- Setters -----

    MODULE PROCEDURE set_output_path_h

        out%path = TRIM(path)

    END PROCEDURE set_output_path_h


    MODULE PROCEDURE set_output_path_iter
        USE tools_output_basics

        !
        INTEGER :: it, ndigits

        ndigits = INT(LOG10(REAL(iter%nmax_()))) + 1
        it = iter%current_iteration()

        out%path = TRIM(out%basepath)//itoh(it,ndigits)

    END PROCEDURE set_output_path_iter

END SUBMODULE class_output_procedures
