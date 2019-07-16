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
! $Id: class_dimensions.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_dimensions

    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: dimensions
    PUBLIC :: OPERATOR(==), OPERATOR(/=)

    ! Constants
    PUBLIC :: null_dim_, quantity
    PUBLIC :: length_, mass_, time_, temperature_
    PUBLIC :: surface_, volume_
    PUBLIC :: velocity_, acceleration_
    PUBLIC :: force_, pressure_
    PUBLIC :: energy_, power_
    PUBLIC :: density_, viscosity_, conductivity_, specific_heat_, youngs_modulus_, therm_exp_coeff_

    TYPE dimensions
        PRIVATE
        INTEGER :: l
        INTEGER :: m
        INTEGER :: t
        INTEGER :: theta
    CONTAINS
        PROCEDURE :: bcast_dim, quantity
        PROCEDURE, PRIVATE :: dim_sqrt
        GENERIC, PUBLIC :: sqrt => dim_sqrt
        PROCEDURE :: debug_dim
        PROCEDURE :: dim_sum, dim_diff, dim_mul, dim_div, dim_pow
        GENERIC :: OPERATOR(+)  => dim_sum
        GENERIC :: OPERATOR(-)  => dim_diff
        GENERIC :: OPERATOR(*)  => dim_mul
        GENERIC :: OPERATOR(/)  => dim_div
        GENERIC :: OPERATOR(**) => dim_pow
        PROCEDURE, PRIVATE :: nemo_dimensions_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_dimensions_sizeof
    END TYPE dimensions


    ! ----- Named Constants -----

    TYPE(dimensions), PARAMETER :: null_dim_     = dimensions(0,0,0,0)

    ! Fundamental Quantities
    TYPE(dimensions), PARAMETER :: length_        = dimensions(1,0,0,0)
    TYPE(dimensions), PARAMETER :: mass_          = dimensions(0,1,0,0)
    TYPE(dimensions), PARAMETER :: time_          = dimensions(0,0,1,0)
    TYPE(dimensions), PARAMETER :: temperature_   = dimensions(0,0,0,1)
    !
    ! Spatial Quantities
    TYPE(dimensions), PARAMETER :: surface_       = dimensions(2,0,0,0)
    TYPE(dimensions), PARAMETER :: volume_        = dimensions(3,0,0,0)
    !
    ! Kinematics
    TYPE(dimensions), PARAMETER :: velocity_      = dimensions(1,0,-1,0)
    TYPE(dimensions), PARAMETER :: acceleration_  = dimensions(1,0,-2,0)
    !
    ! Dynamics
    TYPE(dimensions), PARAMETER :: force_         = dimensions(1,1,-2,0)
    TYPE(dimensions), PARAMETER :: pressure_      = dimensions(-1,1,-2,0)
    !
    ! Energy and Power
    TYPE(dimensions), PARAMETER :: energy_        = dimensions(2,1,-2,0)
    TYPE(dimensions), PARAMETER :: power_         = dimensions(2,1,-3,0)
    !
    ! Physical Properties
    TYPE(dimensions), PARAMETER :: density_       = dimensions(-3,1,0,0)
    TYPE(dimensions), PARAMETER :: viscosity_     = dimensions(-1,1,-1,0)
    TYPE(dimensions), PARAMETER :: conductivity_  = dimensions(1,1,-3,-1)
    TYPE(dimensions), PARAMETER :: specific_heat_ = dimensions(2,0,-2,-1)
    TYPE(dimensions), PARAMETER :: youngs_modulus_ = dimensions(-1,1,-2,0)
    TYPE(dimensions), PARAMETER :: therm_exp_coeff_ = dimensions(0,0,0,-1)


    ! ----- Generic Interfaces -----

    ! Computation
    INTERFACE
        MODULE FUNCTION dim_sqrt(dim)
            IMPLICIT NONE
            TYPE(dimensions) :: dim_sqrt
            CLASS(dimensions), INTENT(IN) :: dim
        END FUNCTION dim_sqrt
    END INTERFACE

    ! Debug
    INTERFACE
        MODULE SUBROUTINE debug_dim(dim)
            IMPLICIT NONE
            CLASS(dimensions), INTENT(IN) :: dim
        END SUBROUTINE debug_dim
    END INTERFACE


    ! ----- Operator Overloading -----

    INTERFACE OPERATOR(==)
        PURE MODULE FUNCTION dim_eq(dim1,dim2)
            IMPLICIT NONE
            LOGICAL :: dim_eq
            TYPE(dimensions), INTENT(IN) :: dim1
            TYPE(dimensions), INTENT(IN)  :: dim2
        END FUNCTION dim_eq
    END INTERFACE OPERATOR(==)

    INTERFACE OPERATOR(/=)
        PURE MODULE FUNCTION dim_ne(dim1,dim2)
            IMPLICIT NONE
            LOGICAL :: dim_ne
            TYPE(dimensions), INTENT(IN) :: dim1
            TYPE(dimensions), INTENT(IN)  :: dim2
        END FUNCTION dim_ne
    END INTERFACE OPERATOR(/=)

    INTERFACE
        PURE MODULE FUNCTION dim_sum(dim1,dim2)
            IMPLICIT NONE
            TYPE(dimensions) :: dim_sum
            CLASS(dimensions), INTENT(IN) :: dim1
            TYPE(dimensions), INTENT(IN)  :: dim2
        END FUNCTION dim_sum

        PURE MODULE FUNCTION dim_diff(dim1,dim2)
            IMPLICIT NONE
            TYPE(dimensions) :: dim_diff
            CLASS(dimensions), INTENT(IN) :: dim1
            TYPE(dimensions), INTENT(IN)  :: dim2
        END FUNCTION dim_diff

        PURE MODULE FUNCTION dim_mul(dim1,dim2)
            IMPLICIT NONE
            TYPE(dimensions) :: dim_mul
            CLASS(dimensions), INTENT(IN) :: dim1
            TYPE(dimensions), INTENT(IN)  :: dim2
        END FUNCTION dim_mul

        PURE MODULE FUNCTION dim_div(dim1,dim2)
            IMPLICIT NONE
            TYPE(dimensions) :: dim_div
            CLASS(dimensions), INTENT(IN) :: dim1
            TYPE(dimensions), INTENT(IN) :: dim2
        END FUNCTION dim_div

        PURE MODULE FUNCTION dim_pow(dim,n)
            IMPLICIT NONE
            TYPE(dimensions) :: dim_pow
            CLASS(dimensions), INTENT(IN) :: dim
            INTEGER, INTENT(IN) :: n
        END FUNCTION dim_pow

        ELEMENTAL MODULE FUNCTION nemo_dimensions_sizeof(dims)
            USE psb_base_mod
            IMPLICIT NONE
            CLASS(dimensions), INTENT(IN) :: dims
            INTEGER(kind=nemo_int_long_)   :: nemo_dimensions_sizeof
        END FUNCTION nemo_dimensions_sizeof

    ! ----- Broadcast -----

        MODULE SUBROUTINE bcast_dim(dim)
            IMPLICIT NONE
            CLASS(dimensions), INTENT(INOUT) :: dim
        END SUBROUTINE bcast_dim
    END INTERFACE

    ! ----- Debug & Log -----
    INTERFACE
        MODULE FUNCTION quantity(dim)
            IMPLICIT NONE
            CHARACTER(len=32) :: quantity
            CLASS(Dimensions), INTENT(IN) :: dim
        END FUNCTION quantity
    END INTERFACE

END MODULE class_dimensions
