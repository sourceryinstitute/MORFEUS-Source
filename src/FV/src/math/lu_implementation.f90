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
SUBMODULE (tools_math) lu_implementation
    IMPLICIT NONE

    CONTAINS

MODULE PROCEDURE lu_solve
    USE class_psblas, ONLY : psb_dpk_
    IMPLICIT NONE
    !! $Id: lu_solve.f90 2469 2007-10-08 10:34:43Z sfilippo $
    !!
    !! Description:
    !!    Solve a system Ax=b where A has already been factored.
    !!    Adapted from LAPACK, this is supposed to be called with
    !!    small sizes (4 or thereabout)
    !!
    INTEGER :: i, j
    REAL(psb_dpk_) :: tmp

    info = 0

    DO i=1, n
        j = ipiv(i)
        IF (j /= i) THEN
            tmp = b(i)
            b(i) = b(j)
            b(j) = tmp
        END IF
    END DO

    DO i=1,n
        b(i+1:n) = b(i+1:n) - b(i) * a(i+1:n,i)
    ENDDO
    ! Upper triangular solve
    DO i=n,1,-1
        b(i) = b(i) - dot_PRODUCT(a(i,i+1:n),b(i+1:n))
        b(i) = b(i)/a(i,i)
    END DO

END PROCEDURE lu_solve


MODULE PROCEDURE lu_fact
    USE class_psblas, ONLY : psb_dpk_
    IMPLICIT NONE
    !!
    !! $Id: lu_fact.f90 2469 2007-10-08 10:34:43Z sfilippo $
    !!
    !! Description:
    !!    Adapted from LAPACK and unrolled for small sizes.
    !!
    REAL(psb_dpk_) :: tmp
    INTEGER :: j, k, jp
    REAL(psb_dpk_), PARAMETER :: tol = 1.d-20

    info = 0


    SELECT CASE(n)

    CASE(4)
        CALL lu_fact_4(a,ipiv,info)

    CASE(3)
        CALL lu_fact_3(a,ipiv,info)

    CASE(2)
        CALL lu_fact_2(a,ipiv,info)

    CASE default

        info = 0
        DO j = 1, n
            jp = j - 1 + MAXLOC(ABS(a(j:,j)),dim=1)
            ipiv(j) = jp
            IF (jp /= j) THEN
                DO k=1,n
                    tmp     = a(j,k)
                    a(j,k)  = a(jp,k)
                    a(jp,k) = tmp
                END DO
            END IF

            IF(ABS(a(j,j)) < tol) THEN
                info = j
            ELSE
                a(j+1:,j) = a(j+1:,j) * (1.d0/a(j,j))
            END IF


            DO k = j + 1, n
                a(j+1:,k) = a(j+1:,k) - a(j,k) * a(j+1:,j)
            END DO

        END DO
    END SELECT

    RETURN

CONTAINS

    SUBROUTINE  lu_fact_4(a,ipiv,info)
        USE class_psblas, ONLY : psb_dpk_

        REAL(psb_dpk_), INTENT(INOUT) :: A(4,4)
        INTEGER, INTENT(OUT) :: ipiv(4)
        INTEGER, INTENT(OUT) :: info
        !
        REAL(psb_dpk_) :: tmp, swp(4)
        INTEGER :: j, k, jp
        REAL(psb_dpk_), PARAMETER :: tol = 1.d-20

        info = 0

        IF (.FALSE.) THEN
            DO j = 1, 4
                jp = j - 1 + MAXLOC(ABS(a(j:,j)),dim=1)
                ipiv(j) = jp
                IF (jp /= j) THEN
                    DO k=1,4
                        tmp     = a(j,k)
                        a(j,k)  = a(jp,k)
                        a(jp,k) = tmp
                    END DO
                END IF

                IF(ABS(a(j,j)) < tol) THEN
                    info = j
                ELSE
                    a(j+1:,j) = a(j+1:,j) * (1.d0/a(j,j))
                END IF

                DO k = j + 1, 4
                    a(j+1:,k) = a(j+1:,k) - a(j,k) * a(j+1:,j)
                END DO

            END DO
        ELSE

            jp = 1 - 1 + MAXLOC(ABS(a(1:,1)),dim=1)
            ipiv(1) = jp
            IF (jp /= 1) THEN
                swp(1:4) = a(1,1:4)
                a(1,1:4) = a(jp,1:4)
                a(jp,1:4) = swp(1:4)
            END IF

            IF(ABS(a(1,1)) < tol) THEN
                info = 1
            ELSE
                a(1+1:,1) = a(1+1:,1) * (1.d0/a(1,1))
            END IF

            DO k = 1 + 1, 4
                a(1+1:,k) = a(1+1:,k) - a(1,k) * a(1+1:,1)
            END DO

            jp = 2 - 1 + MAXLOC(ABS(a(2:,2)),dim=1)
            ipiv(2) = jp
            IF (jp /= 2) THEN
                swp(1:4) = a(2,1:4)
                a(2,1:4) = a(jp,1:4)
                a(jp,1:4) = swp(1:4)
            END IF

            IF(ABS(a(2,2)) < tol) THEN
                info = 2
            ELSE
                a(2+1:,2) = a(2+1:,2) * (1.d0/a(2,2))
            END IF

            DO k = 2 + 1, 4
                a(2+1:,k) = a(2+1:,k) - a(2,k) * a(2+1:,2)
            END DO

            jp = 3 - 1 + MAXLOC(ABS(a(3:,3)),dim=1)
            ipiv(3) = jp
            IF (jp /= 3) THEN
                swp(1:4) = a(3,1:4)
                a(3,1:4) = a(jp,1:4)
                a(jp,1:4) = swp(1:4)
            END IF

            IF(ABS(a(3,3)) < tol) THEN
                info = 3
            ELSE
                a(3+1:,3) = a(3+1:,3) * (1.d0/a(3,3))
            END IF

            DO k = 3 + 1, 4
                a(3+1:,k) = a(3+1:,k) - a(3,k) * a(3+1:,3)
            END DO


            ipiv(4) = 4

            IF(ABS(a(4,4)) < tol) THEN
                info = 4
            END IF

        ENDIF
    END SUBROUTINE lu_fact_4

    SUBROUTINE  lu_fact_3(a,ipiv,info)
        USE class_psblas, ONLY : psb_dpk_
        REAL(psb_dpk_), INTENT(INOUT) :: A(3,3)
        INTEGER, INTENT(OUT) :: ipiv(3)
        INTEGER, INTENT(OUT) :: info
        !
        REAL(psb_dpk_) ::tmp, swp(3)
        INTEGER :: j, k, jp
        REAL(psb_dpk_), PARAMETER :: tol = 1.d-20
        info = 0
        IF (.FALSE.) THEN
            DO j = 1, 3
                jp = j - 1 + MAXLOC(ABS(a(j:,j)),dim=1)
                ipiv(j) = jp
                IF (jp /= j) THEN
                    DO k=1,3
                        tmp     = a(j,k)
                        a(j,k)  = a(jp,k)
                        a(jp,k) = tmp
                    END DO
                END IF

                IF(ABS(a(j,j)) < tol) THEN
                    info = j
                ELSE
                    a(j+1:,j) = a(j+1:,j) * (1.d0/a(j,j))
                END IF


                DO k = j + 1, 3
                    a(j+1:,k) = a(j+1:,k) - a(j,k) * a(j+1:,j)
                END DO

            END DO
        ELSE

            jp = 1 - 1 + MAXLOC(ABS(a(1:,1)),dim=1)
            ipiv(1) = jp
            IF (jp /= 1) THEN
                swp(1:3) = a(1,1:3)
                a(1,1:3) = a(jp,1:3)
                a(jp,1:3) = swp(1:3)
            END IF

            IF(ABS(a(1,1)) < tol) THEN
                info = 1
            ELSE
                a(1+1:,1) = a(1+1:,1) * (1.d0/a(1,1))
            END IF

            DO k = 1 + 1, 3
                a(1+1:,k) = a(1+1:,k) - a(1,k) * a(1+1:,1)
            END DO

            jp = 2 - 1 + MAXLOC(ABS(a(2:,2)),dim=1)
            ipiv(2) = jp
            IF (jp /= 2) THEN
                swp(1:3) = a(2,1:3)
                a(2,1:3) = a(jp,1:3)
                a(jp,1:3) = swp(1:3)
            END IF

            IF(ABS(a(2,2)) < tol) THEN
                info = 2
            ELSE
                a(2+1:,2) = a(2+1:,2) * (1.d0/a(2,2))
            END IF

            DO k = 2 + 1, 3
                a(2+1:,k) = a(2+1:,k) - a(2,k) * a(2+1:,2)
            END DO

            ipiv(3) = 3

            IF(ABS(a(3,3)) < tol) THEN
                info = 3
            END IF

        END IF
    END SUBROUTINE lu_fact_3


    SUBROUTINE  lu_fact_2(a,ipiv,info)
        USE class_psblas, ONLY : psb_dpk_
        REAL(psb_dpk_), INTENT(INOUT) :: A(2,2)
        INTEGER, INTENT(OUT) :: ipiv(2)
        INTEGER, INTENT(OUT) :: info
        !
        REAL(psb_dpk_) ::tmp, swp(2)
        INTEGER :: j, k, jp
        REAL(psb_dpk_), PARAMETER :: tol = 1.d-20
        info = 0

        IF (.FALSE.) THEN
            DO j = 1, 2
                jp = j - 1 + MAXLOC(ABS(a(j:,j)),dim=1)
                ipiv(j) = jp
                IF (jp /= j) THEN
                    DO k=1,2
                        tmp     = a(j,k)
                        a(j,k)  = a(jp,k)
                        a(jp,k) = tmp
                    END DO
                END IF

                IF(ABS(a(j,j)) < tol) THEN
                    info = j
                ELSE
                    a(j+1:,j) = a(j+1:,j) * (1.d0/a(j,j))
                END IF


                DO k = j + 1, 2
                    a(j+1:,k) = a(j+1:,k) - a(j,k) * a(j+1:,j)
                END DO

            END DO

        ELSE

            jp = 1 - 1 + MAXLOC(ABS(a(1:,1)),dim=1)
            ipiv(1) = jp
            IF (jp /= 1) THEN
                swp(1:2) = a(1,1:2)
                a(1,1:2) = a(jp,1:2)
                a(jp,1:2) = swp(1:2)
            END IF

            IF(ABS(a(1,1)) < tol) THEN
                info = 1
            ELSE
                a(1+1:,1) = a(1+1:,1) * (1.d0/a(1,1))
            END IF

            DO k = 1 + 1, 2
                a(1+1:,k) = a(1+1:,k) - a(1,k) * a(1+1:,1)
            END DO

            ipiv(2) = 2

            IF(ABS(a(2,2)) < tol) THEN
                info = 2
            END IF

        ENDIF
    END SUBROUTINE lu_fact_2

END PROCEDURE lu_fact

END SUBMODULE lu_implementation