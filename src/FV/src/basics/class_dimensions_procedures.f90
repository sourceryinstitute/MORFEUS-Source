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
! $Id: class_dimensions.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_dimensions) class_dimensions_procedures

    USE class_psblas

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_dimensions_sizeof
        USE psb_base_mod

        nemo_dimensions_sizeof = 4 * nemo_sizeof_int

    END PROCEDURE nemo_dimensions_sizeof


    ! ----- Broadcast -----

    MODULE PROCEDURE bcast_dim
        !
        LOGICAL :: psblas_on
        INTEGER :: icontxt, mypnum
        INTEGER :: intbuf(4)

        psblas_on = psblas_is_on()
        IF(.NOT. psblas_on) CALL start_psblas

        icontxt = icontxt_()
        mypnum  = mypnum_()

        IF(mypnum == 0) THEN
            intbuf(1) = dim%l
            intbuf(2) = dim%m
            intbuf(3) = dim%t
            intbuf(4) = dim%theta
            CALL psb_bcast(icontxt,intbuf)
        ELSE
            CALL psb_bcast(icontxt,intbuf)
            dim%l = intbuf(1)
            dim%m = intbuf(2)
            dim%t = intbuf(3)
            dim%theta = intbuf(4)
        END IF

    END PROCEDURE bcast_dim


    ! ----- Boolean Operations With DIMENSIONS Objects -----

    MODULE PROCEDURE dim_eq

        dim_eq = (dim1%l == dim2%l .AND.&
            & dim1%m == dim2%m .AND.&
            & dim1%t == dim2%t .AND.&
            & dim1%theta == dim2%theta)

    END PROCEDURE dim_eq


    MODULE PROCEDURE dim_ne
        !

        dim_ne = .NOT. dim_eq(dim1,dim2)

    END PROCEDURE dim_ne


    ! ----- Algebra Operations With DIMENSIONS Objects -----

    MODULE PROCEDURE dim_sum

        IF (dim1 /= dim2)  error stop ' ERROR! Illegal algebraic sum between physical dimensions'
            
        dim_sum = dim1

    END PROCEDURE dim_sum


    MODULE PROCEDURE dim_diff

        dim_diff = dim_sum(dim1,dim2)

    END PROCEDURE dim_diff


    MODULE PROCEDURE dim_mul

        dim_mul%l = dim1%l + dim2%l
        dim_mul%m = dim1%m + dim2%m
        dim_mul%t = dim1%t + dim2%t
        dim_mul%theta = dim1%theta + dim2%theta

    END PROCEDURE dim_mul


    MODULE PROCEDURE dim_div

        dim_div%l = dim1%l - dim2%l
        dim_div%m = dim1%m - dim2%m
        dim_div%t = dim1%t - dim2%t
        dim_div%theta = dim1%theta - dim2%theta

    END PROCEDURE dim_div


    MODULE PROCEDURE dim_pow

        dim_pow%l = dim%l * n
        dim_pow%m = dim%m * n
        dim_pow%t = dim%t * n
        dim_pow%theta = dim%theta * n

    END PROCEDURE dim_pow


    MODULE PROCEDURE dim_sqrt

        IF(MOD(dim%l,2) /= 0 .OR. &
            & MOD(dim%l,2) /= 0 .OR. &
            & MOD(dim%l,2) /= 0 .OR. &
            & MOD(dim%l,2) /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        dim_sqrt%l = dim%l / 2
        dim_sqrt%m = dim%m / 2
        dim_sqrt%t = dim%t / 2
        dim_sqrt%theta = dim%theta / 2

100     FORMAT(' ERROR! Non-integer dimension resulting from DIM_SQRT')

    END PROCEDURE dim_sqrt


    ! ----- Debug & Log -----

    MODULE PROCEDURE debug_dim

        WRITE(*,100) dim%l, dim%m, dim%t, dim%theta

100     FORMAT(' L=',i2,' M=',i2,' T=',i2,' Theta=',i2)

    END PROCEDURE debug_dim


    MODULE PROCEDURE quantity

        IF     (dim == null_dim_) THEN
            quantity = 'pure number'
        ELSE IF(dim == length_) THEN
            quantity = 'length'
        ELSE IF(dim == mass_) THEN
            quantity = 'mass'
        ELSE IF(dim == time_) THEN
            quantity = 'time'
        ELSE IF(dim == temperature_) THEN
            quantity = 'temperature'
        ELSE IF(dim == surface_) THEN
            quantity = 'surface'
        ELSE IF(dim == volume_) THEN
            quantity = 'volume'
        ELSE IF(dim == velocity_) THEN
            quantity = 'velocity'
        ELSE IF(dim == acceleration_) THEN
            quantity = 'acceleration'
        ELSE IF(dim == force_) THEN
            quantity = 'force'
        ELSE IF(dim == pressure_) THEN
            quantity = 'pressure'
        ELSE IF(dim == energy_) THEN
            quantity = 'energy'
        ELSE IF(dim == power_) THEN
            quantity = 'power'
        ELSE IF(dim == density_) THEN
            quantity = 'density'
        ELSE IF(dim == viscosity_) THEN
            quantity = 'viscosity'
        ELSE IF(dim == conductivity_) THEN
            quantity = 'conductivity'
        ELSE IF(dim == specific_heat_) THEN
            quantity = 'specific heat'
        ELSE
            quantity = '"Unsupported quantity"'
        END IF

    END PROCEDURE quantity


END SUBMODULE class_dimensions_procedures
