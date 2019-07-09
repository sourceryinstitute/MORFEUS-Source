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
SUBMODULE(type_table) type_table_procedures

    USE class_psblas

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_table_sizeof
        USE psb_base_mod

        INTEGER(kind=nemo_int_long_)   :: val

        val = 0
        IF (ALLOCATED(tab%lookup)) &
            & val = val + nemo_sizeof_int * SIZE(tab%lookup)
        IF (ALLOCATED(tab%tab)) &
            &  val = val + nemo_sizeof_int * SIZE(tab%tab)
        nemo_table_sizeof = val

    END PROCEDURE nemo_table_sizeof


    ! ----- Constuctor -----

    MODULE PROCEDURE alloc_table
        !
        INTEGER :: info


        IF ( COUNT((/PRESENT(nel),PRESENT(ntab)/)) == 0 ) THEN
            WRITE(*,400)
            CALL abort_psblas
        END IF

        IF(PRESENT(nel)) THEN
            IF (ALLOCATED(a2b%lookup)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
            ALLOCATE(a2b%lookup(nel+1),stat=info)
            IF(info/=0) THEN
                WRITE(*,300)
                CALL abort_psblas
            END IF

            ! Initialization
            a2b%lookup(1) = 1
        END IF

        IF(PRESENT(ntab)) THEN
            IF (ALLOCATED(a2b%tab)) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF
            ALLOCATE(a2b%tab(ntab),stat=info)
            IF(info /= 0) THEN
                WRITE(*,300)
                CALL abort_psblas
            END IF
        END IF

100     FORMAT(' ERROR! A2B%LOOKUP pointer is already allocated')
200     FORMAT(' ERROR! A2B%TAB pointer is already allocated')
300     FORMAT(' ERROR! Memory allocation failure in ALLOC_TABLE')
400     FORMAT(' ERROR! Unsufficient actual arguments in ALLOC_TABLE')

    END PROCEDURE alloc_table


    ! ----- Destructor -----

    MODULE PROCEDURE free_table
        !
        INTEGER :: info

        DEALLOCATE(a2b%lookup,a2b%tab,stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory deallocation failure in FREE_TABLE')

    END PROCEDURE free_table


    ! ----- Other -----

    MODULE PROCEDURE get_dual_table
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

        !
        INTEGER :: i, i1, j, info, ir, is
        INTEGER :: nel, nrd, nrt, ntot
        INTEGER, ALLOCATABLE :: level(:)

        ! Number of rows in the original A2B TABLE
        nrt = SIZE(a2b%lookup) - 1

        ! Number of rows of the B2A TABLE, i.e. number of stacks
        nrd = MAXVAL(a2b%tab)

        ! The total amount of elements in the stacks is equal to:
        ! - the number of elements in A2B%tab
        ! - the number of elements in B2A%tab
        ! I.e. the number of links expressed by a TABLE is equal
        ! in the DUAL one.
        ntot = SIZE(a2b%tab)

        ! LEVEL(I) = level of I-th stack container in B2A dual TABLE.
        ALLOCATE(level(nrd),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! 1st swap of A2B%tab -> Evaluates the level of each stack
        level(:) = 0
        DO ir = 1, nrt
            ! Offset for the elements contained in the IR-th row of A2B
            i1 = a2b%lookup(ir) - 1
            !
            ! Number of elements contained in the I-th row of A2B
            nel = a2b%lookup(ir+1) - a2b%lookup(ir)
            !
            DO i = 1, nel
                ! Index of the stack
                is = a2b%tab(i1+i)

                ! Increase the filling level of the IS-th stack
                level(is) = level(is) + 1
            END DO
        END DO

        CALL b2a%alloc_table(nel=nrd,ntab=ntot)

        ! The following operation is equivalent to size the empty stacks
        b2a%lookup(1) = 1
        DO is = 1, nrd
            b2a%lookup(is+1) = b2a%lookup(is) + level(is)
        END DO

        ! 2nd swap of A2B%tab --> Building of the stacks set B2A%tab
        !
        ! Resets to zero the level counter of the stacks set
        level(:) = 0
        DO ir = 1, nrt
            i1 = a2b%lookup(ir) -1
            nel = a2b%lookup(ir+1) - a2b%lookup(ir)
            DO i = 1, nel
                is = a2b%tab(i1+i)
                level(is) = level(is) + 1
                !
                ! Sets the position on the B2A TABLE of the new element found
                ! to belong to the IS-th stack
                j = b2a%lookup(is) -1 + level(is)
                !
                ! Stores index IR in the IS-th stack
                b2a%tab(j) = ir
            END DO
        END DO

        DEALLOCATE(level)

100     FORMAT(' ERROR! Memory allocation failure in GET_DUAL_TABLE')

    END PROCEDURE get_dual_table


END SUBMODULE type_table_procedures
