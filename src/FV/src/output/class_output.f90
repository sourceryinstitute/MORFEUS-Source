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
! $Id: class_output.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_output

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: output

    TYPE output
        PRIVATE
        INTEGER :: fmt
        CHARACTER(len=32) :: basepath
        CHARACTER(len=32) :: path
    CONTAINS
        PROCEDURE :: create_output  ! Constructor
        PROCEDURE :: fmt_, path_    ! Getters
        PROCEDURE, PRIVATE :: set_output_path_h, set_output_path_iter
        GENERIC, PUBLIC :: set_output_path => set_output_path_h, set_output_path_iter ! Setters
        PROCEDURE, PRIVATE :: nemo_output_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_output_sizeof
    END TYPE output


    ! ----- Generic Interface -----


  INTERFACE

    MODULE FUNCTION nemo_output_sizeof(obj)
        USE psb_base_mod
        USE class_psblas
        CLASS(output), INTENT(IN) :: obj
        INTEGER(kind=nemo_int_long_)   :: nemo_output_sizeof
    END FUNCTION nemo_output_sizeof

  END INTERFACE

    ! ----- Setters -----

  INTERFACE

    MODULE SUBROUTINE set_output_path_h(out,path)
        CLASS(output),      INTENT(INOUT) :: out
        CHARACTER(len=32), INTENT(IN) :: path
    END SUBROUTINE set_output_path_h


    MODULE SUBROUTINE set_output_path_iter(out,iter)
        USE class_iterating
        USE tools_output_basics
        CLASS(output),    INTENT(INOUT) :: out
        TYPE(iterating), INTENT(IN) :: iter
    END SUBROUTINE set_output_path_iter

    ! ----- Constructor -----

    MODULE SUBROUTINE create_output(out,input_file,sec)
        USE tools_input
        USE tools_output_basics
        CLASS(output),      INTENT(INOUT) :: out
        CHARACTER(len=*), INTENT(IN) :: input_file
        CHARACTER(len=*), INTENT(IN) :: sec
    END SUBROUTINE create_output


    ! ----- Getters -----

    MODULE FUNCTION fmt_(out)
        INTEGER :: fmt_
        CLASS(output), INTENT(IN) :: out
    END FUNCTION fmt_


    MODULE FUNCTION path_(out)
        USE tools_output_basics
        CHARACTER(len=32) :: path_
        CLASS(output), INTENT(IN) :: out
        !
        INTEGER :: l
        CHARACTER(len=1) :: path_end
    END FUNCTION path_

  END INTERFACE

END MODULE class_output
