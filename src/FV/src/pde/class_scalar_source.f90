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
! $Id$
!
! Description:
!    To be added...
!
MODULE class_scalar_source

    USE class_psblas
    USE class_dimensions

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: scalar_source         ! Class

    TYPE scalar_source
        PRIVATE
        TYPE(dimensions) :: dim
        REAL(psb_dpk_) :: sc
        REAL(psb_dpk_) :: sp
    CONTAINS
        PROCEDURE, PRIVATE :: create_scalar_source
        GENERIC, PUBLIC :: create_source => create_scalar_source
        PROCEDURE, PRIVATE :: get_scalar_source_sc, get_scalar_source_sp  ! Getters
        GENERIC, PUBLIC :: sc_ => get_scalar_source_sc
        GENERIC, PUBLIC :: sp_ => get_scalar_source_sp
        PROCEDURE, PRIVATE :: get_scalar_source_dim
        GENERIC, PUBLIC :: dim_ => get_scalar_source_dim   ! Getter
        PROCEDURE, PRIVATE :: nemo_scalar_source_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_scalar_source_sizeof
    END TYPE scalar_source

    ! ----- Generic Interfaces -----

    INTERFACE

        MODULE FUNCTION nemo_scalar_source_sizeof(src)
            USE psb_base_mod
            USE class_psblas
            IMPLICIT NONE
            CLASS(scalar_source), INTENT(IN) :: src
            INTEGER(kind=nemo_int_long_)   :: nemo_scalar_source_sizeof
        END FUNCTION nemo_scalar_source_sizeof


        ! Constructor
        MODULE SUBROUTINE create_scalar_source(src,input_file,sec,dim)
            USE tools_input
            IMPLICIT NONE
            CLASS(scalar_source), INTENT(INOUT) :: src
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            TYPE(dimensions), INTENT(IN) :: dim
        END SUBROUTINE create_scalar_source


        ! Getters

        MODULE FUNCTION get_scalar_source_dim(src)
            IMPLICIT NONE
            TYPE(dimensions) :: get_scalar_source_dim
            CLASS(scalar_source), INTENT(IN) :: src
        END FUNCTION get_scalar_source_dim


        MODULE FUNCTION get_scalar_source_sc(src)
            IMPLICIT NONE
            REAL(psb_dpk_) :: get_scalar_source_sc
            CLASS(scalar_source), INTENT(IN) :: src
        END FUNCTION get_scalar_source_sc

        MODULE FUNCTION get_scalar_source_sp(src)
            REAL(psb_dpk_) :: get_scalar_source_sp
            CLASS(scalar_source), INTENT(IN) :: src
        END FUNCTION get_scalar_source_sp

    END INTERFACE

END MODULE class_scalar_source
