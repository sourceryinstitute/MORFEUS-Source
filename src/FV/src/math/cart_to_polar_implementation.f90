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
! SOFTWARE, EVEN if ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!---------------------------------------------------------------------------------
!
SUBMODULE (tools_math) cart_to_polar_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE cart_to_polar
        USE tools_math, ONLY : pi
        USE class_psblas, ONLY : psb_dpk_
        IMPLICIT NONE
        !! $Id: cart_to_polar.f90 2469 2007-10-08 10:34:43Z sfilippo $
        !!
        !! Description:
        !!    Converts cartesian into polar coordinates.
        !!             (X,Y) -> (RHO,THETA)
        !!    THETA is in degrees.
        !!
        LOGICAL :: square_
        REAL(psb_dpk_) :: r2

        IF(PRESENT(square)) THEN
            square_ = square
        ELSE
            square_ = .FALSE. ! Default
        END IF

        ! Computes RHO
        r2= x ** 2 + y ** 2
        IF(square_) THEN
            rho = r2
        ELSE
            rho = SQRT(r2)
        END IF

        ! Computes THETA
        theta = ATAN(y/x)
        IF(x < 0) THEN
            theta = theta + pi
        ELSEIF(y < 0) THEN
            theta = theta + 2 * pi
        END IF

        ! Conversion radians -> degrees
        theta = 360.d0 * theta / (2 * pi)

        END PROCEDURE cart_to_polar

END SUBMODULE cart_to_polar_implementation