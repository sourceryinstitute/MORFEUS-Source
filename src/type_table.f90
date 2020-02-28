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
!    It provides the derived data type TABLE, used in READ_*_MESH, supplying
!    the same functionalities of the CONNECTIVITY class, but with attribute
!    PUBLIC access and CSR-like implementation.
!
MODULE type_table

    USE class_psblas, ONLY : nemo_int_long_

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: table

    TYPE table
        INTEGER, ALLOCATABLE :: lookup(:)
        INTEGER, ALLOCATABLE :: tab(:)
    CONTAINS
        PROCEDURE:: alloc_table, free_table, get_dual_table
        PROCEDURE, PRIVATE :: nemo_table_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_table_sizeof
    END TYPE table

    INTERFACE
        ELEMENTAL MODULE FUNCTION nemo_table_sizeof(tab)
            USE psb_base_mod
            IMPLICIT NONE
            CLASS(table), INTENT(IN) :: tab
            INTEGER(kind=nemo_int_long_)   :: nemo_table_sizeof
        END FUNCTION nemo_table_sizeof

    ! ----- Constructor -----

        MODULE SUBROUTINE alloc_table(a2b, nel, ntab)
            IMPLICIT NONE
            CLASS(table), INTENT(INOUT) :: a2b
            INTEGER, INTENT(IN), OPTIONAL :: nel, ntab
        END SUBROUTINE alloc_table

    ! ----- Destructor -----

        MODULE SUBROUTINE free_table(a2b)
            IMPLICIT NONE
            CLASS(table), INTENT(INOUT) :: a2b
        END SUBROUTINE free_table

    ! ----- Other -----
        MODULE SUBROUTINE get_dual_table(a2b,b2a)
        ! The B2A dual TABLE can be thought as a set of stacks which are
        ! eventually sticked together sequentially.
        ! There is one stack for every element of B which appears at least
        ! once in A2B%tab. The I-th stack is built with the indices of
        ! the elements of A (a2b%lookup superscripts) which are connected to
        ! the I-th element of B.
        ! A practical application is extracting the face2vert connectivity,
        ! starting from the vert2face one.
        ! WARNING: we are not referring to a STACK according to the common
        ! vocabulary of the programming science, i.e. as a special container
        ! implemented with pointers etc.
            IMPLICIT NONE
            CLASS(table), INTENT(IN) :: a2b
            TYPE(table), INTENT(OUT) :: b2a
        END SUBROUTINE get_dual_table
    END INTERFACE

END MODULE type_table
