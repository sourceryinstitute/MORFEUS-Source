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
! $Id: class_least_squares.f90 3175 2008-06-13 12:59:07Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_least_squares
    USE class_psblas, ONLY : nemo_int_long_, psb_dpk_
    USE class_connectivity, ONLY : connectivity
    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: least_squares           !! Class
    PUBLIC :: free_least_squares      !! Procedure
    PUBLIC :: set_least_squares       !! Procedure

    TYPE least_squares
        REAL(psb_dpk_), ALLOCATABLE :: A(:)
    CONTAINS
        PROCEDURE, PRIVATE :: nemo_least_squares_sizeof
        GENERIC,   PUBLIC  :: nemo_sizeof => nemo_least_squares_sizeof
        PROCEDURE, PUBLIC  :: solve_least_squares         !! Math Operations
    END TYPE least_squares

  INTERFACE

    ELEMENTAL MODULE FUNCTION nemo_least_squares_sizeof(lsq)
        IMPLICIT NONE
        CLASS(least_squares), INTENT(IN) :: lsq
        INTEGER(kind=nemo_int_long_)     :: nemo_least_squares_sizeof
    END FUNCTION nemo_least_squares_sizeof

    ! ----- Constructor -----

    MODULE SUBROUTINE alloc_least_squares(lsr,n,ncd)
        IMPLICIT NONE
        TYPE(least_squares), ALLOCATABLE, INTENT(OUT) :: lsr(:)
        INTEGER, INTENT(IN) :: n, ncd
    END SUBROUTINE alloc_least_squares

    ! ----- Destructor -----

    MODULE SUBROUTINE free_least_squares(lsr)
        IMPLICIT NONE
        TYPE(least_squares), ALLOCATABLE  :: lsr(:)
    END SUBROUTINE free_least_squares

    ! ----- Setter -----

    MODULE SUBROUTINE set_least_squares(lsr,ncd,desc,c2c,f2b,faces,cell_cntr,face_cntr)
        USE class_face, ONLY : face
        USE class_vector, ONLY : vector
        USE tools_math
        USE psb_base_mod
        IMPLICIT NONE
        TYPE(least_squares), ALLOCATABLE :: lsr(:)
        INTEGER, INTENT(IN) :: ncd
        TYPE(psb_desc_type), INTENT(IN) :: desc
        TYPE(connectivity), INTENT(IN) :: c2c, f2b
        TYPE(face), INTENT(IN) :: faces(:)
        TYPE(vector), INTENT(IN) :: cell_cntr(:), face_cntr(:)
    END SUBROUTINE set_least_squares

    ! ----- Computational Routines -----

    MODULE SUBROUTINE solve_least_squares(lsr,rhs)
        IMPLICIT NONE
        CLASS(least_squares), INTENT(IN) :: lsr
        REAL(psb_dpk_), INTENT(INOUT) :: rhs(:)
    END SUBROUTINE solve_least_squares

  END INTERFACE

END MODULE class_least_squares
