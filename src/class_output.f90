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

MODULE class_output
    USE class_psblas, ONLY : psb_dpk_, nemo_int_long_
    USE class_mesh,   ONLY : mesh
    !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_scalar_field, ONLY : scalar_field
    !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_vector_field, ONLY : vector_field
    !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: output, create_output

    TYPE output
        PRIVATE
        INTEGER :: fmt
        CHARACTER(len=:), ALLOCATABLE :: basepath
        CHARACTER(len=:), ALLOCATABLE :: path
    CONTAINS
        PROCEDURE, NOPASS :: create_output  ! Constructor
        PROCEDURE :: fmt_, path_    ! Getters
        PROCEDURE, PRIVATE :: set_output_path_h, set_output_path_iter
        GENERIC, PUBLIC :: set_output_path => set_output_path_h, set_output_path_iter ! Setters
        PROCEDURE, PRIVATE :: nemo_output_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_output_sizeof
        PROCEDURE, PRIVATE :: write_output
        GENERIC, PUBLIC :: write => write_output
        PROCEDURE, NOPASS :: get_scalar_field
        PROCEDURE, NOPASS :: get_vector_field
    END TYPE output

    ! ----- Generic Interface -----

    INTERFACE

        MODULE FUNCTION nemo_output_sizeof(obj)
            IMPLICIT NONE
            CLASS(output), INTENT(IN) :: obj
            INTEGER(kind=nemo_int_long_) :: nemo_output_sizeof
        END FUNCTION nemo_output_sizeof

        ! ----- Setters -----

        MODULE SUBROUTINE set_output_path_h(out,path)
            IMPLICIT NONE
            CLASS(output),    INTENT(INOUT) :: out
            CHARACTER(LEN=*), INTENT(IN)    :: path
        END SUBROUTINE set_output_path_h

        MODULE SUBROUTINE set_output_path_iter(out,iter)
            USE class_iterating, ONLY : iterating
            IMPLICIT NONE
            CLASS(output),   INTENT(INOUT) :: out
            TYPE(iterating), INTENT(IN)    :: iter
        END SUBROUTINE set_output_path_iter

        ! ----- Constructor -----

        MODULE FUNCTION create_output(input_file,sec) RESULT(out)
            IMPLICIT NONE
            CLASS(output),   ALLOCATABLE :: out
            CHARACTER(LEN=*), INTENT(IN) :: input_file
            CHARACTER(LEN=*), INTENT(IN) :: sec
        END FUNCTION create_output

        ! ----- Getters -----

        MODULE FUNCTION fmt_(out)
            IMPLICIT NONE
            INTEGER :: fmt_
            CLASS(output), INTENT(IN) :: out
        END FUNCTION fmt_

        MODULE FUNCTION path_(out) RESULT(path)
            IMPLICIT NONE
            CLASS(output),    INTENT(IN)  :: out
            CHARACTER(len=:), ALLOCATABLE :: path
        END FUNCTION path_

        MODULE SUBROUTINE write_output(out, msh, sfield, vfield, iter)
            USE class_iterating, ONLY : iterating
            IMPLICIT NONE
            !! author: Ian Porter, GSE
            !! date: 01/04/2020
            !!
            !! This subroutine is a generic writer
            !!
            CLASS(output),                    INTENT(INOUT)        :: out          !! DT of output file info
            TYPE(mesh),                       INTENT(IN)           :: msh          !! DT of mesh info
            TYPE(scalar_field), DIMENSION(:), INTENT(IN), OPTIONAL :: sfield       !! DT of scalar info
            TYPE(vector_field), DIMENSION(:), INTENT(IN), OPTIONAL :: vfield       !! DT of vector info
            TYPE(iterating),                  INTENT(IN), OPTIONAL :: iter         !! DT of iteration info
        END SUBROUTINE write_output

        MODULE FUNCTION get_scalar_field(fld) RESULT (x_glob)
            IMPLICIT NONE
            !! Returns a scalar field
            TYPE(scalar_field), INTENT(IN) :: fld
            REAL(psb_dpk_), ALLOCATABLE :: x_glob(:)

        END FUNCTION get_scalar_field

        MODULE FUNCTION get_vector_field(fld) RESULT (x_glob)
            IMPLICIT NONE
            !! Returns a vector field
            TYPE(vector_field), INTENT(IN) :: fld
            REAL(psb_dpk_), ALLOCATABLE :: x_glob(:,:)

        END FUNCTION get_vector_field

    END INTERFACE

END MODULE class_output
