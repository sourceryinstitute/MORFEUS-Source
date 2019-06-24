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
MODULE class_discretization

    USE class_psblas
    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: discretization        ! Class
    PUBLIC :: read_par              ! Other constructor
    PUBLIC :: cd_, up_              ! Named constants

    TYPE discretization
        PRIVATE
        INTEGER :: id
        REAL(psb_dpk_) :: blend
    CONTAINS
        PROCEDURE :: id_            ! Getters
        PROCEDURE, PRIVATE :: nemo_discretization_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_discretization_sizeof
    END TYPE discretization


  ! ----- Generic Interfaces -----

  INTERFACE
    MODULE FUNCTION nemo_discretization_sizeof(dsc)
        USE psb_base_mod
        USE class_psblas
        IMPLICIT NONE
        CLASS(discretization), INTENT(IN) :: dsc
        INTEGER(kind=nemo_int_long_)   :: nemo_discretization_sizeof
    END FUNCTION nemo_discretization_sizeof
  END INTERFACE

  ! ----- Constructors -----

  INTERFACE read_par
    MODULE FUNCTION read_par_discretization(input_file,sec,par,default)RESULT(r)
        USE class_psblas
        USE tools_input
        IMPLICIT NONE
        TYPE(discretization) :: r
        CHARACTER(len=*),     INTENT(IN) :: input_file
        CHARACTER(len=*),     INTENT(IN) :: sec
        CHARACTER(len=*),     INTENT(IN) :: par
        TYPE(discretization), INTENT(IN) :: default
    END FUNCTION read_par_discretization
  END INTERFACE read_par


    ! ----- Named Constants -----

    ! Private
    INTEGER, PARAMETER :: id_cd_ = 1
    INTEGER, PARAMETER :: id_up_ = 2

    ! Public
    TYPE(discretization), PARAMETER :: cd_ = discretization(id_cd_,1.d0)
    TYPE(discretization), PARAMETER :: up_ = discretization(id_up_,1.d0)


    ! ----- Getters -----
  INTERFACE

    MODULE FUNCTION id_(ds)
        INTEGER :: id_
        CLASS(discretization), INTENT(IN) :: ds
    END FUNCTION id_

  END INTERFACE

END MODULE class_discretization
