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

!

SUBMODULE (tools_math) pwl_implementation
    IMPLICIT NONE
    !! author: sfilippo
    !! date: 2008-04-22 14:51:09

    CONTAINS

        MODULE PROCEDURE pwl_interp_x_s
        USE class_psblas, ONLY : psb_dpk_, abort_psblas
        USE tools_math, ONLY : pwl_nearest
        IMPLICIT NONE
        !! Description:
        !!    Piecewise linear interpolation for a x-y data series.
        !!    Type REAL, rank = 0
        !!
        INTEGER :: i1, i2
        REAL(psb_dpk_) :: dx, dy, m

        IF(SIZE(x_data) /= SIZE(y_data)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pwl_nearest(x,x_data,i1,i2)

        IF(x == x_data(i1)) THEN
            y = y_data(i1)
        ELSE
            dx = x_data(i2) - x_data(i1)
            dy = y_data(i2) - y_data(i1)
            m = dy / dx

            y = y_data(i1) + m * (x - x_data(i1))
        END IF

100     FORMAT(' ERROR! Array size mismatch in PWL_INTERP_X_S')

        END PROCEDURE pwl_interp_x_s


        MODULE PROCEDURE pwl_interp_x_v
        USE class_psblas, ONLY : psb_dpk_, abort_psblas
        USE tools_math, ONLY : pwl_nearest
        IMPLICIT NONE
        !! Description:
        !!    Piecewise linear interpolation for a x-y data series.
        !!    Type REAL, rank = 1
        !!
        INTEGER :: i1, i2
        REAL(psb_dpk_) :: dx
        REAL(psb_dpk_), DIMENSION(SIZE(y)) :: dy, m

        IF(    SIZE(x_data) /= SIZE(y_data,2) .OR. &
            & SIZE(y)      /= SIZE(y_data,1)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pwl_nearest(x,x_data,i1,i2)

        IF(x == x_data(i1)) THEN
            y(:) = y_data(:,i1)
        ELSE
            dx = x_data(i2) - x_data(i1)
            dy(:) = y_data(:,i2) - y_data(:,i1)

            m = dy / dx

            y = y_data(:,i1) + m * (x - x_data(i1))
        END IF

100     FORMAT(' ERROR! Array size mismatch in PWL_INTERP_X_V')

        END PROCEDURE pwl_interp_x_v


        MODULE PROCEDURE pwl_interp_x_vec
        USE class_psblas, ONLY : psb_dpk_, abort_psblas
        USE class_vector, ONLY : vector, OPERATOR(*), OPERATOR(+), OPERATOR(-)
        USE tools_math, ONLY : pwl_nearest
        IMPLICIT NONE
        !! Description:
        !!    Piecewise linear interpolation for a x-y data series.
        !!    Type VECTOR, rank = 1
        !!
        INTEGER :: i1, i2
        REAL(psb_dpk_) :: dx
        TYPE(vector) :: dy, m

        IF(SIZE(x_data) /= SIZE(y_data)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pwl_nearest(x,x_data,i1,i2)

        IF(x_data(i1) == x) THEN
            y = y_data(i1)
        ELSE
            dx = x_data(i2) - x_data(i1)
            dy = y_data(i2) - y_data(i1)

            m = (1.d0 / dx) * dy

            y = y_data(i1) + (x - x_data(i1)) * m
        END IF

100     FORMAT(' ERROR! Array size mismatch in PWL_INTERP_X_VEC')

        END PROCEDURE pwl_interp_x_vec


        MODULE PROCEDURE pwl_nearest
        USE class_psblas, ONLY : abort_psblas
        IMPLICIT NONE
        !! $Id: pwl_nearest.f90 3093 2008-04-22 14:51:09Z sfilippo $
        !!
        !! Description:
        !!    Given an array X_DATA and a scalar X returns the subscripts of the two
        !!    nearest elements in the data series.
        !!
        INTEGER :: i, n

        n = SIZE(x_data)

        IF(x < x_data(1)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        i = 1
        DO
            IF(x < x_data(i)) EXIT
            i = i + 1
            IF(i > n) EXIT
        END DO

        IF(i <= n) THEN
            i1 = i - 1
            i2 = i
        ELSE
            i1 = n - 1
            i2 = n
        END IF

100     FORMAT(' ERROR! Illegal entry in PWL_NEAREST')

        END PROCEDURE pwl_nearest

        MODULE PROCEDURE pwl_interp_dx_s
        USE class_psblas, ONLY : psb_dpk_
        IMPLICIT NONE
        !! Description:
        !!    Piecewise linear interpolation for a set of equally spaced data.
        !!    Type REAL, rank = 0
        !!
        INTEGER :: i1, i2
        REAL(psb_dpk_) :: dx_inv, slope, x1, xrel

        dx_inv = 1.d0 / dx

        xrel = (x - xmin) * dx_inv

        i1 = FLOOR(xrel) + 1
        i2 = CEILING(xrel) + 1

        slope = (y_data(i2) - y_data(i1)) * dx_inv

        x1 = xmin + dx * (i1 - 1)

        y = y_data(i1) + slope * (x - x1)

        END PROCEDURE pwl_interp_dx_s


        MODULE PROCEDURE pwl_interp_dx_v
        USE class_psblas, ONLY : psb_dpk_, abort_psblas
        IMPLICIT NONE
        !! Description:
        !!    Piecewise linear interpolation for a set of equally spaced data.
        !!    Type REAL, rank = 1
        !!
        INTEGER :: i, i1, i2, n
        REAL(psb_dpk_) :: dx_inv, slope, x1, xrel

        IF(SIZE(y) /= SIZE(x)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        dx_inv = 1.d0 / dx

        n = SIZE(x)

        DO i = 1, n
            xrel = (x(i) - xmin) * dx_inv

            i1 = FLOOR(xrel) + 1
            i2 = CEILING(xrel) + 1

            slope = (y_data(i2) - y_data(i1)) * dx_inv

            x1 = xmin + dx * (i1 - 1)

            y(i) = y_data(i1) + slope * (x(i) - x1)
        END DO

100     FORMAT(' ERROR! Size mismatch in PWL_INTERP_DX_V')

        END PROCEDURE pwl_interp_dx_v

        MODULE PROCEDURE pwl_deriv_x_s
        USE class_psblas, ONLY : abort_psblas
        USE tools_math, ONLY : pwl_nearest

        IMPLICIT NONE
        !! $Id: pwl_deriv_x.f90 3093 2008-04-22 14:51:09Z sfilippo $
        !!
        !! Description:
        !!    Compute the first derivative of a piecewise linear x-y data series.
        !!    Type REAL, rank = 0
        !!
        INTEGER :: i1, i2

        IF(SIZE(x_data) /= SIZE(y_data)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pwl_nearest(x,x_data,i1,i2)

        dydx = (y_data(i2) - y_data(i1)) / (x_data(i2) - x_data(i1))

100     FORMAT(' ERROR! Array size mismatch in PWL_DERIV_X_S')

        END PROCEDURE pwl_deriv_x_s


        MODULE PROCEDURE pwl_deriv_x_v
        USE class_psblas, ONLY : abort_psblas
        USE tools_math, ONLY : pwl_nearest
        IMPLICIT NONE
        !! Description:
        !!    Compute the first derivative of a piecewise linear x-y data series.
        !!    Type REAL, rank = 1
        !!
        INTEGER :: i1, i2

        IF(    SIZE(x_data) /= SIZE(y_data,2) .OR. &
            & SIZE(dydx)   /= SIZE(y_data,1)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pwl_nearest(x,x_data,i1,i2)

        dydx = (y_data(:,i2) - y_data(:,i1)) / (x_data(i2) - x_data(i1))

100     FORMAT(' ERROR! Array size mismatch in PWL_DERIV_X_V')

        END PROCEDURE pwl_deriv_x_v


        MODULE PROCEDURE pwl_deriv_x_vec
        USE class_psblas, ONLY : psb_dpk_, abort_psblas
        USE tools_math, ONLY : pwl_nearest
        USE class_vector, ONLY : OPERATOR(*), OPERATOR(-)
        IMPLICIT NONE
        !! Description:
        !!    Compute the first derivative of a piecewise linear x-y data series.
        !!    Type VECTOR, rank = 0
        !!
        INTEGER :: i1, i2
        REAL(psb_dpk_) :: dx

        IF(SIZE(x_data) /= SIZE(y_data)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pwl_nearest(x,x_data,i1,i2)

        dx = x_data(i2) - x_data(i1)
        dydx = (1.d0 / dx) * (y_data(i2) - y_data(i1))

100     FORMAT(' ERROR! Array size mismatch in PWL_DERIV_X_VEC')
        ! REMARK! It returns the right derivative!
        END PROCEDURE pwl_deriv_x_vec

END SUBMODULE pwl_implementation
