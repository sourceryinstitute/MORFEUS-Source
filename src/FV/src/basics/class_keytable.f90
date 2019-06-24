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
MODULE class_keytable

    ! use class_psblas ! Move it to functions

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: keytable ! Class
    PUBLIC :: alloc_keytable, free_keytable                 ! Constructor/Destructor
    PUBLIC :: get_kt_row, get_rows, get_row_ub, get_row_lb  ! Getters
    PUBLIC :: get_row_size                                  ! Getters, cont.
    PUBLIC :: set_kt_row, ASSIGNMENT(=)                     ! Setters

    !  `public :: assignment(=)'  ==  `public :: copy_keytable'

    ! Abstract data type for handling sets of information where
    ! A two dimensional array would be wasteful, since the length
    ! Of each row in the table varies tremendously.

    TYPE a_row
        INTEGER, ALLOCATABLE :: entries(:)
    END TYPE a_row

    TYPE keytable
        PRIVATE
        TYPE(a_row), ALLOCATABLE :: row(:)
    CONTAINS
        PROCEDURE, PRIVATE :: nemo_keytable_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_keytable_sizeof
        PROCEDURE, PUBLIC :: exists
    END TYPE keytable

    INTERFACE
        ELEMENTAL MODULE FUNCTION nemo_keytable_sizeof(table)
            USE class_psblas, ONLY : nemo_int_long_ ! Move it to functions
            IMPLICIT NONE
            CLASS(keytable), INTENT(IN) :: table
            INTEGER(kind=nemo_int_long_)   :: nemo_keytable_sizeof
        END FUNCTION nemo_keytable_sizeof
    END INTERFACE

    ! Description of class members.
    ! entry(:)  is a_row with varying numbers of elements
    ! columns(:)    is an array of rows


    ! ----- Generic Interfaces -----

    ! Setter
    INTERFACE ASSIGNMENT(=)
        MODULE SUBROUTINE copy_keytable(table2,table1)
            USE class_psblas, ONLY : abort_psblas ! Move it to functions
            IMPLICIT NONE
            TYPE(keytable), INTENT(IN)  :: table1
            TYPE(keytable), INTENT(OUT) :: table2
        END SUBROUTINE copy_keytable
    END INTERFACE ASSIGNMENT(=)

    INTERFACE
        ELEMENTAL MODULE FUNCTION nemo_a_row_sizeof(row)
        USE class_psblas, ONLY : nemo_int_long_, nemo_sizeof_int ! Move it to functions
            IMPLICIT NONE
            TYPE(a_row), INTENT(IN) :: row
            INTEGER(kind=nemo_int_long_)   :: nemo_a_row_sizeof
        END FUNCTION nemo_a_row_sizeof

    ! ----- Constructor -----

        MODULE SUBROUTINE alloc_keytable(table,lb,ub)
            USE class_psblas, ONLY : abort_psblas ! Move it to functions
            IMPLICIT NONE
            TYPE(keytable), INTENT(INOUT) :: table
            INTEGER, INTENT(IN) :: lb, ub
        END SUBROUTINE alloc_keytable

    ! ----- Destructor -----

        MODULE SUBROUTINE free_keytable(table)
            USE class_psblas, ONLY : abort_psblas, mypnum_ ! Move it to functions
            IMPLICIT NONE
            TYPE(keytable), INTENT(INOUT) :: table
        END SUBROUTINE free_keytable

    ! ----- Getters -----

        MODULE SUBROUTINE get_kt_row(table,i,irow)
            USE class_psblas, ONLY : abort_psblas ! Move it to functions
            IMPLICIT NONE
            TYPE(keytable), INTENT(IN), TARGET :: table
            INTEGER, INTENT(IN)        :: i
            INTEGER, POINTER           :: irow(:)
        END SUBROUTINE get_kt_row

    ! indicates if the keytable memory has been allocated
        MODULE FUNCTION exists(table)
            !use class_psblas ! Move it to functions
            IMPLICIT NONE
            CLASS(keytable), INTENT(IN) :: table
            LOGICAL                    :: exists
        END FUNCTION exists

    ! ----- Get info about table  -----

        MODULE FUNCTION get_row_lb (table)
            !use class_psblas ! Move it to functions
            IMPLICIT NONE
            INTEGER :: get_row_lb
            TYPE (keytable), INTENT(IN) :: table
        END FUNCTION get_row_lb


        MODULE FUNCTION get_row_ub(table)
            IMPLICIT NONE
            INTEGER :: get_row_ub
            TYPE (keytable), INTENT(IN) :: table
        END FUNCTION get_row_ub

        MODULE FUNCTION get_rows(table)
            IMPLICIT NONE
            INTEGER :: get_rows
            TYPE (keytable), INTENT(IN) :: table
        END FUNCTION get_rows

        MODULE  FUNCTION get_row_size(table,i)
            USE class_psblas, ONLY : abort_psblas ! Move it to functions
            IMPLICIT NONE
            INTEGER :: get_row_size
            TYPE (keytable), INTENT(IN) :: table
            INTEGER,INTENT(IN) :: i
        END FUNCTION get_row_size

    ! ----- Setters -----

        MODULE SUBROUTINE set_kt_row(table,i,irow)
            USE class_psblas, ONLY : abort_psblas ! Move it to functions
            IMPLICIT NONE
            TYPE(keytable), INTENT(INOUT) :: table
            INTEGER, INTENT(IN)        :: i
            INTEGER, INTENT(IN)        :: irow(:)
        END SUBROUTINE set_kt_row

    END INTERFACE

END MODULE class_keytable
