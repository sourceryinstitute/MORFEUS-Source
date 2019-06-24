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
! $Id: renum.F90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Renumbering
!
SUBMODULE(renum) renum_procedures

    USE class_psblas

    IMPLICIT NONE

CONTAINS

    ! ----- Constructor -----

    MODULE PROCEDURE start_renum
        USE class_connectivity

        INTEGER, ALLOCATABLE  :: temp(:)
        INTEGER :: info, mypnum, ncells

        mypnum = mypnum_()
        ncells = nel_(c2c)

        ! Builds permutation array PERM on P0, according to IRENUM

        SELECT CASE(irenum)
        CASE(1)
            ! Gibbs-Poole-Stockmeyer
            WRITE(*,*) 'Matrix renumbering: Gibbs-Poole-Stockmeyer'
            CALL cmp_gps(c2c)
        END SELECT
        ! OTHER METHODS: TO BE IMPLEMENTED ...

        ! Extends permutation array backward to PERM(0) = 0.
        ! 0 = index of GHOST cells at domain boundaries.
        ALLOCATE(temp(0:ncells),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        temp(0) = 0
        temp(1:ncells) = perm(1:ncells)

#if defined(HAVE_MOVE_ALLOC)
        CALL move_ALLOC(temp,perm)
#else
        DEALLOCATE(perm)
        ALLOCATE(perm(LBOUND(temp,1):UBOUND(temp,1)))
        perm = temp
        DEALLOCATE(temp)
#endif


        ! Builds inverse permutation array PINV
        CALL build_pinv

100     FORMAT(' ERROR! Memory allocation failure in START_RENUM')

    END PROCEDURE start_renum


    ! ----- Destuctor -----

    MODULE PROCEDURE stop_renum
        DEALLOCATE(perm,pinv)
    END PROCEDURE stop_renum


    ! ----- Broadcast -----

    MODULE PROCEDURE print_renum
        INTEGER :: i,n
        CHARACTER(len=*),PARAMETER :: fmt='(3(i8,1x))'
        n=SIZE(perm) - 1
        WRITE(iout,*) 'Renumbering permutation and inverse: ',n
        DO i=1,n
            WRITE(iout,fmt) i,perm(i), pinv(i)
        END DO

    END PROCEDURE print_renum
    ! ----- Broadcast -----

    MODULE PROCEDURE build_pinv

        INTEGER :: icontxt, mypnum
        INTEGER :: i, info, j, n

        icontxt = icontxt_()
        mypnum  = mypnum_()

        ! PERM(:) is supposed to be associated only in P0
        ! => Check and allocation on P > 0
        IF(mypnum == 0) THEN
            IF(.NOT.ALLOCATED(perm)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            n=SIZE(perm)


        END IF

        IF(mypnum == 0) THEN

            ! Build inverse permutation array PINV(:)
            !  * new = perm(old)
            !  * old = pinv(new)
            n = SIZE(perm)-1 ! = ncells
            ALLOCATE(pinv(0:n),stat=info)
            IF(info /= 0) THEN
                WRITE(*,300)
                CALL abort_psblas
            END IF

            pinv(0) = 0
            DO i = 1, n
                j = perm(i)
                pinv(j) = i
            END DO

        ENDIF

100     FORMAT(' ERROR! Illegal Send: PERM array not associated')
200     FORMAT(' ERROR! Illegal Recv: PERM array already associated')
300     FORMAT(' ERROR! Memory allocation failure in BCAST_RENUM')

    END PROCEDURE build_pinv


    ! ----- Interfaces To Renumbering Methods -----

    ! Gibbs-Poole-Stockmeyer renumbering algorithm
    MODULE PROCEDURE cmp_gps
        USE class_connectivity
        USE psb_gps_mod, ONLY : psb_gps_reduce

        INTEGER :: i, info, j
        !
        ! GPS variables
        INTEGER :: ideg, nodes
        INTEGER :: ibw,  idpth, ipf
        INTEGER, ALLOCATABLE :: ndstk(:,:), iold(:), ndeg(:)
        INTEGER, POINTER :: iconn(:) => NULL()
        INTEGER :: nconn

        !------------------------------------------------------------------

        ! Order of matrix, number of cells
        nodes = nel_(c2c)

        ! Evaluation of graph maximum degree --> maximum connectivity
        ideg = max_conn(c2c)

        ! Array allocation
        ALLOCATE(ndstk(nodes,ideg), iold(nodes), perm(nodes+1), &
            & ndeg(nodes),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! table of connectivities --> ndtsk
        ndstk = 0
        DO i = 1, nodes
            CALL get_ith_conn(iconn,c2c,i)
            nconn = SIZE(iconn)
            DO j = 1, nconn
                ndstk(i,j) = iconn(j)
            END DO
        END DO
        ! REMARK: bad do loop! The inversion of the stride would require
        ! a nested if check, and repeated calls to get_ith_conn.

        ! Old node-numbering
        DO i = 1, nodes
            iold(i) = i
        END DO
        perm = 0

        ! Call routine for Gibbs-Poole-Stockmeyer algorithm (ACM-TOMS no. 508)
        CALL psb_gps_reduce(ndstk,nodes,ideg,iold,perm,ndeg,ibw,ipf,idpth)

        DEALLOCATE(ndstk,iold,ndeg)
        NULLIFY(iconn)

100     FORMAT(' ERROR! Memory allocation failure in CMP_GPS')

    END PROCEDURE cmp_gps


    ! ----- Routines For The Application Of The Renumbering -----

    ! Integer array
    MODULE PROCEDURE apply_renum_array
        !
        INTEGER :: i, j, n

        n = SIZE(a)

        DO i = 1, n
            j = a(i)
            a(i) = perm(j)
        END DO

    END PROCEDURE apply_renum_array


    ! CELL objects array
    MODULE PROCEDURE apply_renum_cell
        USE class_cell

        INTEGER :: i, j, n
        TYPE(cell), ALLOCATABLE  :: work(:)

        n = SIZE(c)

        CALL alloc_cell(work,n)

        work(:) = c(:)

        DO i = 1, n
            j = pinv(i)
            c(i) = work(j)
        END DO

        CALL free_cell(work)

    END PROCEDURE apply_renum_cell


    ! FACE objects array
    MODULE PROCEDURE apply_renum_face
        USE class_face

        INTEGER :: i, j, n

        n = SIZE(f)

        DO i = 1, n
            j = master_(f(i))
            CALL set_face(f(i),master=perm(j))
            j = slave_(f(i))
            CALL set_face(f(i),slave=perm(j))
        END DO

    END PROCEDURE apply_renum_face


    ! CONNECTIVITY object
    MODULE PROCEDURE apply_renum_conn
        USE class_connectivity

        INTEGER :: i, iold, imax, info, j, nb, nconn
        INTEGER, POINTER :: iconn_old(:) => NULL()
        INTEGER, POINTER :: iconn_new(:) => NULL()
        TYPE(connectivity) :: work

        IF(.NOT.(apply == to_a_ .OR. &
            &   apply == to_b_ .OR. &
            &   apply == to_a_and_b_)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        work = a2b     ! That's a copy.
        nb = nel_(a2b) ! Number of B elements

        ! Applies permutation to A elements of A2B connectivity
        IF(apply == to_a_ .OR. apply == to_a_and_b_) THEN

            imax = max_conn(a2b) ! Maximum connectivity degree

            ALLOCATE(iconn_new(imax),stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            DO i = 1, nb
                CALL get_ith_conn(iconn_old,a2b,i) ! Pointing
                nconn = SIZE(iconn_old)

                DO j = 1, nconn
                    iold = iconn_old(j)
                    iconn_new(j) = perm(iold)
                END DO

                CALL set_ith_conn(a2b,i,iconn_new(1:nconn)) ! Copy
            END DO

            DEALLOCATE(iconn_new)
            CALL free_conn(work)
        END IF

        ! The result of the permutation of A elements is stored in A2B.
        IF(apply == to_a_and_b_) work = a2b

        ! Applies permutation to B elements of A2B connectivity
        IF(apply == to_b_ .OR. apply == to_a_and_b_) THEN

            DO i = 1, nb
                j = pinv(i)

                ! That's a pointing.
                CALL get_ith_conn(iconn_old,work,j)

                ! That's a copy of the section of work pointed by ICONN_OLD.
                CALL set_ith_conn(a2b,i,iconn_old)
            END DO

            CALL free_conn(work)
        END IF

        NULLIFY(iconn_old)


100     FORMAT(' ERROR! Illegal value of APPLY argument in APPLY_RENUM_CONN')
200     FORMAT(' ERROR! Memory allocation failure in APPLY_RENUM_CONN')

    END PROCEDURE apply_renum_conn


END SUBMODULE renum_procedures
