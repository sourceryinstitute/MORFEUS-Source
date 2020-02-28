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
! $Id: tools_math.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    To be added...
!
MODULE tools_math
    USE class_psblas, ONLY : psb_dpk_
    USE class_vector, ONLY : vector
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: cart_to_polar
    PUBLIC :: factorize
    PUBLIC :: solve_sys
    PUBLIC :: build_prec
    PUBLIC :: lin_interp
    PUBLIC :: pwl_nearest
    PUBLIC :: pwl_deriv
    PUBLIC :: pwl_interp
    PUBLIC :: sort
    PUBLIC :: pi, it_time_, it_convergence_, it_counter_

    INTERFACE factorize
        MODULE PROCEDURE :: lu_fact
        MODULE PROCEDURE :: cholesky_fact_m
        MODULE PROCEDURE :: choloesky_fact_v
    END INTERFACE factorize

    INTERFACE solve_sys
        MODULE PROCEDURE :: lu_solve
        MODULE PROCEDURE :: cholesky_solve_m
        MODULE PROCEDURE :: cholesky_solve_v
        MODULE PROCEDURE :: psb_solve_sys
    END INTERFACE solve_sys

    INTERFACE build_prec
        PROCEDURE :: psb_build_prec
    END INTERFACE build_prec

    INTERFACE lin_interp
        MODULE PROCEDURE :: lin_interp_s
        MODULE PROCEDURE :: lin_interp_v
    END INTERFACE lin_interp

    INTERFACE pwl_deriv
        MODULE PROCEDURE :: pwl_deriv_x_s
        MODULE PROCEDURE :: pwl_deriv_x_v
        MODULE PROCEDURE :: pwl_deriv_x_vec
    END INTERFACE pwl_deriv

    INTERFACE pwl_interp
        MODULE PROCEDURE :: pwl_interp_dx_s
        MODULE PROCEDURE :: pwl_interp_dx_v
        MODULE PROCEDURE :: pwl_interp_x_s
        MODULE PROCEDURE :: pwl_interp_x_v
        MODULE PROCEDURE :: pwl_interp_x_vec
    END INTERFACE pwl_interp

    INTERFACE sort
        MODULE PROCEDURE :: isort
    END INTERFACE sort

    INTERFACE

        MODULE SUBROUTINE cart_to_polar(x,y,rho,theta,square)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(IN) :: x, y
            REAL(psb_dpk_), INTENT(OUT) :: rho, theta
            LOGICAL, INTENT(IN), OPTIONAL :: square
        END SUBROUTINE cart_to_polar

        MODULE SUBROUTINE lu_fact(n,A,ipiv,info)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n
            REAL(psb_dpk_), INTENT(INOUT) :: A(n,n)
            INTEGER, INTENT(OUT) :: ipiv(n)
            INTEGER, INTENT(OUT) :: info
        END SUBROUTINE lu_fact

        MODULE SUBROUTINE cholesky_fact_m(A)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(INOUT) :: A(:,:)
        END SUBROUTINE cholesky_fact_m

        MODULE SUBROUTINE choloesky_fact_v(A)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(INOUT) :: A(:)
        END SUBROUTINE choloesky_fact_v

        MODULE SUBROUTINE lu_solve(n,A,ipiv,b,info)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n
            REAL(psb_dpk_), INTENT(IN) :: A(n,n)
            INTEGER, INTENT(IN) :: ipiv(n)
            REAL(psb_dpk_), INTENT(INOUT) :: b(n)
            INTEGER, INTENT(OUT) :: info
        END SUBROUTINE lu_solve

        MODULE SUBROUTINE cholesky_solve_m(A,b)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(IN) :: A(:,:)
            REAL(psb_dpk_), INTENT(INOUT) :: b(:)
        END SUBROUTINE cholesky_solve_m

        MODULE SUBROUTINE cholesky_solve_v(A,b)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(IN) :: A(:)
            REAL(psb_dpk_), INTENT(INOUT) :: b(:)
        END SUBROUTINE cholesky_solve_v

        MODULE SUBROUTINE psb_solve_sys(A,prec,b,x,desc,cmethod,eps,itmax,iter,err)
            USE psb_base_mod, ONLY : psb_dspmat_type, psb_desc_type
            USE psb_prec_mod, ONLY : psb_dprec_type
            IMPLICIT NONE
            TYPE(psb_dspmat_type), INTENT(IN)    :: A
            TYPE(psb_dprec_type),  INTENT(INOUT) :: prec
            REAL(psb_dpk_),      INTENT(IN)    :: b(:)
            REAL(psb_dpk_),      INTENT(INOUT) :: x(:)
            TYPE(psb_desc_type),   INTENT(IN)    :: desc
            CHARACTER(len=*),      INTENT(IN)    :: cmethod
            REAL(psb_dpk_),      INTENT(IN)    :: eps
            INTEGER,               INTENT(IN)    :: itmax
            INTEGER,               INTENT(OUT)   :: iter
            REAL(psb_dpk_),      INTENT(OUT)   :: err
        END SUBROUTINE psb_solve_sys

        SUBROUTINE psb_build_prec(cprec,nlev,cmethod,A,desc_a,prec)
        !! Note: This can't be a MODULE SUBROUTINE b/c of a compiler bug with
        !!       gfortran-8.3. IP 6/6/2019
            USE psb_base_mod, ONLY : psb_dspmat_type, psb_desc_type
            USE psb_prec_mod, ONLY : psb_dprec_type
            IMPLICIT NONE
            CHARACTER(len=*),      INTENT(IN) :: cprec
            INTEGER,               INTENT(IN) :: nlev
            CHARACTER(len=*),      INTENT(IN) :: cmethod
            TYPE(psb_dspmat_type), INTENT(IN) :: A
            TYPE(psb_desc_type),   INTENT(INOUT) :: desc_a
            TYPE(psb_dprec_type),  INTENT(INOUT) :: prec
        END SUBROUTINE psb_build_prec

        MODULE FUNCTION lin_interp_s(f1,f2,fac)
            IMPLICIT NONE
            REAL(psb_dpk_) :: lin_interp_s
            REAL(psb_dpk_), INTENT(IN) :: f1, f2
            REAL(psb_dpk_), INTENT(IN) :: fac
        END FUNCTION lin_interp_s

        MODULE FUNCTION lin_interp_v(f1,f2,fac)
            IMPLICIT NONE
            TYPE(vector) :: lin_interp_v
            TYPE(vector),   INTENT(IN) :: f1, f2
            REAL(psb_dpk_), INTENT(IN) :: fac
        END FUNCTION lin_interp_v

        MODULE SUBROUTINE pwl_nearest(x,x_data,i1,i2)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(IN) :: x
            REAL(psb_dpk_), INTENT(IN) :: x_data(:)
            INTEGER, INTENT(OUT) :: i1, i2
        END SUBROUTINE pwl_nearest

        MODULE SUBROUTINE pwl_deriv_x_s(dydx,x,y_data,x_data)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(OUT) :: dydx
            REAL(psb_dpk_), INTENT(IN) :: x
            REAL(psb_dpk_), INTENT(IN) :: y_data(:)
            REAL(psb_dpk_), INTENT(IN) :: x_data(:)
        END SUBROUTINE pwl_deriv_x_s

        MODULE SUBROUTINE pwl_deriv_x_v(dydx,x,y_data,x_data)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(OUT) :: dydx(:)
            REAL(psb_dpk_), INTENT(IN) :: x
            REAL(psb_dpk_), INTENT(IN) :: y_data(:,:)
            REAL(psb_dpk_), INTENT(IN) :: x_data(:)
        END SUBROUTINE pwl_deriv_x_v

        MODULE SUBROUTINE pwl_deriv_x_vec(dydx,x,y_data,x_data)
            IMPLICIT NONE
            TYPE(vector), INTENT(OUT) :: dydx
            REAL(psb_dpk_), INTENT(IN) :: x
            TYPE(vector), INTENT(IN) :: y_data(:)
            REAL(psb_dpk_), INTENT(IN) :: x_data(:)
        END SUBROUTINE pwl_deriv_x_vec

        MODULE SUBROUTINE pwl_interp_dx_s(y,x,y_data,xmin,dx)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(OUT) :: y
            REAL(psb_dpk_), INTENT(IN) :: x
            REAL(psb_dpk_), INTENT(IN) :: y_data(:)
            REAL(psb_dpk_), INTENT(IN) :: xmin, dx
        END SUBROUTINE pwl_interp_dx_s

        MODULE SUBROUTINE pwl_interp_dx_v(y,x,y_data,xmin,dx)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(OUT) :: y(:)
            REAL(psb_dpk_), INTENT(IN) :: x(:)
            REAL(psb_dpk_), INTENT(IN) :: y_data(:)
            REAL(psb_dpk_), INTENT(IN) :: xmin, dx
        END SUBROUTINE pwl_interp_dx_v

        MODULE SUBROUTINE pwl_interp_x_s(y,x,y_data,x_data)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(OUT) :: y
            REAL(psb_dpk_), INTENT(IN) :: x
            REAL(psb_dpk_), INTENT(IN) :: y_data(:)
            REAL(psb_dpk_), INTENT(IN) :: x_data(:)
        END SUBROUTINE pwl_interp_x_s

        MODULE SUBROUTINE pwl_interp_x_v(y,x,y_data,x_data)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(OUT) :: y(:)
            REAL(psb_dpk_), INTENT(IN) :: x
            REAL(psb_dpk_), INTENT(IN) :: y_data(:,:)
            REAL(psb_dpk_), INTENT(IN) :: x_data(:)
        END SUBROUTINE pwl_interp_x_v

        MODULE SUBROUTINE pwl_interp_x_vec(y,x,y_data,x_data)
            IMPLICIT NONE
            TYPE(vector), INTENT(OUT) :: y
            REAL(psb_dpk_), INTENT(IN) :: x
            TYPE(vector), INTENT(IN) :: y_data(:)
            REAL(psb_dpk_), INTENT(IN) :: x_data(:)
        END SUBROUTINE pwl_interp_x_vec

    ! WARNING! One cannot use an ELEMENTAL function instead of a
    ! subroutine, because the DATA argument is not scalar.
    ! If rank(X) = 1 -> a subroutine is mandatory
    ! If rank(X) = 0 -> a function could be used
    ! We use in both cases a subroutine + generic interface

        MODULE SUBROUTINE isort(ivet)
            IMPLICIT NONE
            INTEGER, INTENT(INOUT) :: ivet(:)
        END SUBROUTINE isort

    END INTERFACE

    ! ----- SLATEC routines -----

    ! INTERFACE
    !     MODULE SUBROUTINE dnls1e(fcn,iopt,m,n,x,fvec,tol,nprint, &
    !         &            info,iw,wa,lwa)
    !         USE class_psblas, ONLY : psb_dpk_
    !         INTEGER :: iopt, m, n, nprint, info, lwa, iw(n)
    !         REAL(psb_dpk_) :: tol, x(n), fvec(m), wa(lwa)
    !         EXTERNAL fcn
    !     END SUBROUTINE dnls1e

    !     MODULE SUBROUTINE xsetf(kontrol)
    !         INTEGER :: kontrol
    !     END SUBROUTINE xsetf
    ! END INTERFACE


    ! ----- Named Constants -----

    ! Greek PI
    REAL(psb_dpk_), PARAMETER :: pi = 3.14159265358979323846264338327950288d0

    ! ITERATING class
    INTEGER, PARAMETER :: it_time_        = 1
    INTEGER, PARAMETER :: it_convergence_ = 2
    INTEGER, PARAMETER :: it_counter_     = 3

END MODULE tools_math
