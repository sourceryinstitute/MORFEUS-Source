!! category: Morfeus-FV
!! summary: Procedures for matrix LU factorization and solution
!!
!! Copyright Notice
!! ----------------
!!
!!     (c) 2019 Guide Star Engineering, LLC
!!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!!     under contract
!!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!!     contract # NRC-HQ-60-17-C-0007

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

MODULE op_field
    USE class_scalar_field, ONLY : scalar_field
    USE class_vector_field, ONLY : vector_field
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: OPERATOR(*), fld_flux, rhie_chow

    INTERFACE fld_flux
        PROCEDURE :: vector_field_flux
    END INTERFACE fld_flux

    INTERFACE OPERATOR(*)
        PROCEDURE :: scalar_vector_fld_mul
        PROCEDURE :: scalar_fld_vector_mul
    END INTERFACE

    INTERFACE

        MODULE FUNCTION vector_field_flux(fld) RESULT(res)
            IMPLICIT NONE
            TYPE(scalar_field) :: res
            TYPE(vector_field), INTENT(IN) :: fld
        END FUNCTION vector_field_flux

        IMPURE MODULE FUNCTION scalar_vector_fld_mul(fld_s,fld_v) RESULT(res)
            IMPLICIT NONE
            TYPE(vector_field) :: res
            TYPE(scalar_field), INTENT(IN) :: fld_s
            TYPE(vector_field), INTENT(IN) :: fld_v
        END FUNCTION scalar_vector_fld_mul

        IMPURE MODULE FUNCTION scalar_fld_vector_mul(fld_s,v) RESULT(res)
            USE class_vector, ONLY : vector
            IMPLICIT NONE
            TYPE(vector_field) :: res
            TYPE(scalar_field), INTENT(IN) :: fld_s
            TYPE(vector), INTENT(IN) :: v
        END FUNCTION scalar_fld_vector_mul

        MODULE FUNCTION rhie_chow(phi,fld) RESULT(res)
            IMPLICIT NONE
            TYPE(scalar_field) :: res
            TYPE(scalar_field), INTENT(IN) :: phi,fld
        END FUNCTION rhie_chow

    END INTERFACE

END MODULE op_field
