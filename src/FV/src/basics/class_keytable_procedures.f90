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
! $Id: class_keytable.f90 3206 2008-07-02 12:39:02Z sfilippo $
!
! Description:
!    Provides storage of an irregular 2d integer information, where rows would all
!    have different length.  Similar to the connectivity class, except that
!    keytables cannot currently be inverted.  However, keytables provide convenient
!    storage for non-contiguous data and provide bounds checking on calls
!
! Provides:
!    KEYTABLE             class.
!    NEW_KEYTABLE         constructor for KEYTABLE class.
!    FREE_KEYTABLE        destructor for KEYTABLE class.
!    GET_KT_ROW           returns ith row of the keytable
!    SET_KT_ROW           sets, allocating if necessary, the key table row
!    GET_ROW_UB           returns the upper bound of rows in a keytable
!    GET_ROW_LB           returns the lower bound of rows in a keytable
!
SUBMODULE(class_keytable) class_keytable_procedures

    ! use class_psblas ! Move it to functions

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_a_row_sizeof
        USE class_psblas, ONLY : nemo_int_long_, nemo_sizeof_int ! Move it to functions

        IF (ALLOCATED(row%entries)) THEN
            nemo_a_row_sizeof = nemo_sizeof_int * SIZE(row%entries)
        ELSE
            nemo_a_row_sizeof = 0
        END IF

    END PROCEDURE nemo_a_row_sizeof

    MODULE PROCEDURE nemo_keytable_sizeof
        USE class_psblas, ONLY : nemo_int_long_ ! Move it to functions

        INTEGER(kind=nemo_int_long_)   :: val
        INTEGER :: i

        val = 0
        IF (ALLOCATED(table%row)) &
            & val = val + SUM(nemo_a_row_sizeof(table%row))
        nemo_keytable_sizeof = val

    END PROCEDURE nemo_keytable_sizeof

    ! ----- Constructor -----

    MODULE PROCEDURE alloc_keytable
        USE class_psblas, ONLY : abort_psblas ! Move it to functions
        !
        INTEGER :: info

        IF (ALLOCATED(table%row)) THEN
            WRITE(6,100)
            CALL abort_psblas
        ENDIF

        ALLOCATE(table%row(lb:ub), stat = info)

        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Keytable already exists...cannot allocate memory.')
200     FORMAT(' ERROR! Memory allocation failure in ALLOC_KEYTABLE')

    END PROCEDURE alloc_keytable

    ! copy constructor or "reconstructor"
    MODULE PROCEDURE copy_keytable
        USE class_psblas, ONLY : abort_psblas ! Move it to functions
        INTEGER          :: lb,ub
        INTEGER, POINTER :: row(:)

        INTEGER :: info,i

        ! This statement should be redundant.....
        IF ( ALLOCATED(table2%row) ) CALL free_keytable(table2)

        lb  = get_row_lb(table1)
        ub  = get_row_ub(table1)

        ALLOCATE(table2%row(lb:ub),stat=info)

        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        DO i = lb, ub

            CALL get_kt_row(table1,i,row)
            CALL set_kt_row(table2,i,row)
        ENDDO

100     FORMAT(' ERROR! Memory allocation failure in COPY_KEYTABLE')

    END PROCEDURE copy_keytable


    ! ----- Destructor -----

    MODULE PROCEDURE free_keytable
        USE class_psblas, ONLY : abort_psblas, mypnum_ ! Move it to functions
        !
        INTEGER :: i

        IF (.NOT.ALLOCATED(table%row)) THEN
            IF(mypnum_() == 0) WRITE(*,100)
            RETURN
        ENDIF

        DO i = get_row_lb(table), get_row_ub(table)
            IF ( ALLOCATED(table%row(i)%entries) ) DEALLOCATE(table%row(i)%entries)
        ENDDO

        DEALLOCATE(table%row)

100     FORMAT(' WARNING! Attempt to free unallocated keytable in FREE_KEYTABLE')

    END PROCEDURE free_keytable


    ! ----- Getters -----

    MODULE PROCEDURE get_kt_row
        USE class_psblas, ONLY : abort_psblas ! Move it to functions

        IF ( i < get_row_lb(table) ) THEN
            WRITE(6,100)
            CALL abort_psblas
        ENDIF

        IF ( i > get_row_ub(table) ) THEN
            WRITE(6,200)
            CALL abort_psblas
        ENDIF

        irow => table%row(i)%entries

100     FORMAT(' ERROR! Attempt to get below lower bound of KEYTABLE')
200     FORMAT(' ERROR! Attempt to get beyond upper bound of KEYTABLE')

    END PROCEDURE get_kt_row

    ! indicates if the keytable memory has been allocated
    MODULE PROCEDURE exists
        !use class_psblas ! Move it to functions

        IF (ALLOCATED(table%row)) THEN
            exists = .TRUE.
        ELSE
            exists = .FALSE.
        ENDIF

    END PROCEDURE exists

    ! ----- Get info about table  -----

    MODULE PROCEDURE get_row_lb
        !use class_psblas ! Move it to functions

        get_row_lb = LBOUND(table%row,1)

    END PROCEDURE get_row_lb


    MODULE PROCEDURE get_row_ub

        get_row_ub = UBOUND(table%row,1)

    END PROCEDURE get_row_ub

    MODULE PROCEDURE get_rows

        get_rows = UBOUND(table%row,1) - LBOUND(table%row,1) + 1

    END PROCEDURE get_rows


    MODULE PROCEDURE get_row_size
        USE class_psblas, ONLY : abort_psblas ! Move it to functions

        IF ( i < get_row_lb(table) ) THEN
            WRITE(6,100)
            CALL abort_psblas
        ENDIF

        IF ( i > get_row_ub(table) ) THEN
            WRITE(6,200)
            CALL abort_psblas
        ENDIF

        get_row_size = SIZE(table%row(i)%entries)

100     FORMAT("Error:  Attempt to size below lower bound in GET_ROW_SIZE")
200     FORMAT("Error:  Attempt to size beyond upper bound of GET_ROW_SIZE")

    END PROCEDURE get_row_size


    ! ----- Setters -----

    MODULE PROCEDURE set_kt_row
        USE class_psblas, ONLY : abort_psblas ! Move it to functions

        !
        INTEGER :: info

        IF ( i < get_row_lb(table) ) THEN
            WRITE(6,200)
            CALL abort_psblas
        ENDIF

        IF ( i > get_row_ub(table) ) THEN
            WRITE(6,300)
            CALL abort_psblas
        ENDIF

        ! wipe out existing entries

        IF  ( ALLOCATED(table%row(i)%entries) ) THEN
            DEALLOCATE(table%row(i)%entries)
        ENDIF

        ALLOCATE(table%row(i)%entries(SIZE(irow)),stat = info)

        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        table%row(i)%entries = irow

        !else leave row nullified

100     FORMAT(' ERROR! Memory allocation failure in SET_KT_ROW')
200     FORMAT("Error:  Attempt to set below lower bound of KEYTABLE")
300     FORMAT("Error:  Attempt to set beyond upper bound of KEYTABLE")
    END PROCEDURE set_kt_row

END SUBMODULE class_keytable_procedures
