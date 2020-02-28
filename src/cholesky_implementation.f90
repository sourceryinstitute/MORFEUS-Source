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
SUBMODULE (tools_math) cholesky_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE cholesky_solve_v
        USE class_psblas, ONLY : psb_dpk_
        IMPLICIT NONE
        !! $Id$
        !!
        !! Description:
        !!    Solves a linear system Ax=b, where A has been previously decomposed
        !!    according to the Cholesky factorization.
        !!    A: type REAL, rank 1.
        !!
        INTEGER :: i, is, j, ii, ij, ji, n, n2, nn
        REAL(psb_dpk_) :: sum, y(SIZE(b))

        is = SIZE(A)
        n  = (-1 + SQRT(REAL(1 + is * 8))) * 0.5
        n2 = 2 * n

        y(1) = b(1) / A(1)
        DO i = 2, n
            sum = b(i)
            DO j = 1, i - 1
                ij = i + (n2 - j) * (j - 1) * 0.5
                sum = sum - A(ij) * y(j)
            END DO

            ii = i + (n2 - i) * (i - 1) * 0.5
            y(i) = sum / A(ii)
        END DO

        nn = SIZE(A)
        b(n) = y(n) / A(nn)
        DO i = n - 1, 1, -1
            sum = y(i)
            DO j = i + 1, n
                ji = j + (n2 - i) * (i - 1) * 0.5
                sum = sum - A(ji) * b(j)
            END DO

            ii = i + (n2 - i) * (i - 1) * 0.5
            b(i) = sum / A(ii)
        END DO


        END PROCEDURE cholesky_solve_v


        MODULE PROCEDURE cholesky_solve_m
        USE class_psblas, ONLY : psb_dpk_
        IMPLICIT NONE
        !!
        !! Description:
        !!    Solves a linear system Ax=b, where A has been previously decomposed
        !!    according to the Cholesky factorization.
        !!    A: type REAL, rank 2.
        !!
        INTEGER :: i, j, n
        REAL(psb_dpk_) :: y(SIZE(b)), sum

        n = SIZE(b)

        y(1) = b(1) / A(1,1)
        DO i = 2, n
            sum = b(i)
            DO j = 1, i - 1
                sum = sum - A(i,j) * y(j)
            END DO

            y(i) = sum / A(i,i)
        END DO

        b(n) = y(n) / A(n,n)
        DO i = n - 1, 1, -1
            sum = y(i)
            DO j = i + 1, n
                sum = sum - A(j,i) * b(j)
            END DO

            b(i) = sum / A(i,i)
        END DO

        END PROCEDURE cholesky_solve_m

        MODULE PROCEDURE choloesky_fact_v
        USE class_psblas, ONLY : psb_dpk_
        IMPLICIT NONE
        !! Description:
        !!    Computes the Cholesky factorization of a symmetric positive definite
        !!    matrix:  A = LL'.
        !!    A: type REAL, rank 1
        !!
        INTEGER :: i, j, k
        INTEGER :: ii, ij, ik, jj, jk, is, n, n2
        REAL(psb_dpk_) :: sum

        is = SIZE(A)
        n  = (-1 + SQRT(REAL(1 + is * 8))) * 0.5
        n2 = 2 * n

        A(1) = SQRT(A(1))
        DO i = 2, n
            DO j =1, i - 1

                ij = i + (n2 - j) * (j - 1) * 0.5
                sum = A(ij)

                DO k = 1, j - 1
                    ik = i + (n2 - k) * (k - 1) * 0.5
                    jk = j + (n2 - k) * (k - 1) * 0.5
                    sum = sum - A(ik) * A(jk)
                END DO

                jj = j + (n2 - j) * (j - 1) * 0.5
                A(ij) = sum / A(jj)
            END DO


            ii = i + (n2 - i) * (i - 1) * 0.5
            sum = A(ii)

            DO k = 1, i - 1
                ik = i + (n2 - k) * (k - 1) * 0.5
                sum = sum - A(ik) ** 2
            END DO

            A(ii) = SQRT(sum)
        END DO

        END PROCEDURE choloesky_fact_v


        MODULE PROCEDURE  cholesky_fact_m
        USE class_psblas, ONLY : psb_dpk_
        IMPLICIT NONE
        !! Description:
        !!    Computes the Cholesky factorization of a symmetric positive definite
        !!    matrix A: = LL'.
        !!    A: type REAL, rank 2
        !!
        INTEGER :: i, j, k, n
        REAL(psb_dpk_) :: sum

        n = SIZE(A,1)

        A(1,1) = SQRT(A(1,1))
        DO i = 2, n
            DO j =1, i - 1
                sum = A(i,j)
                DO k = 1, j - 1
                    sum = sum - A(i,k) * A(j,k)
                END DO
                A(i,j) = sum / A(j,j)
            END DO

            sum = A(i,i)
            DO k = 1, i - 1
                sum = sum - A(i,k) ** 2
            END DO
            A(i,i) = SQRT(sum)
        END DO

        END PROCEDURE cholesky_fact_m

END SUBMODULE cholesky_implementation
