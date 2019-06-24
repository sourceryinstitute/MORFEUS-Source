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
! $Id: $
!
! Description:
!    To be added...
!
MODULE op_div
     USE class_scalar_pde, ONLY : scalar_pde
       !! A gfortran 8.3.0 bug prevents moving this into the scalar_pde_div or flux_pde_div interface bodies
     USE class_vector_pde, ONLY : vector_pde
       !! A gfortran 8.3.0 bug prevents moving this into the vector_pde_div interface body
    USE class_vector_field, ONLY : vector_field
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_scalar_field, ONLY : scalar_field
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_discretization, ONLY : discretization
    USE class_psblas, ONLY : psb_dpk_
    USE class_mesh, ONLY : mesh, check_mesh_consistency
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: pde_div

    INTERFACE pde_div
        PROCEDURE :: scalar_pde_div
        PROCEDURE :: vector_pde_div
        PROCEDURE :: flux_pde_div
    END INTERFACE pde_div

    INTERFACE

        MODULE SUBROUTINE scalar_pde_div(sign,pde,phi,flux,side,ds)
            IMPLICIT NONE
            CHARACTER(LEN=1),     INTENT(IN) :: sign
            TYPE(scalar_pde),     INTENT(INOUT) :: pde
            TYPE(scalar_field),   INTENT(IN)    :: phi
            TYPE(scalar_field),   INTENT(INOUT)    :: flux
            REAL(psb_dpk_),     INTENT(IN), OPTIONAL :: side
            TYPE(discretization), INTENT(IN), OPTIONAL :: ds
        END SUBROUTINE scalar_pde_div

        MODULE SUBROUTINE vector_pde_div(sign,pde,phi,flux,side,ds)
            IMPLICIT NONE
            CHARACTER(LEN=1),     INTENT(IN) :: sign
            TYPE(vector_pde),     INTENT(INOUT) :: pde
            TYPE(vector_field),   INTENT(IN) :: phi
            TYPE(scalar_field),   INTENT(INOUT) :: flux
            REAL(psb_dpk_),     INTENT(IN), OPTIONAL :: side
            TYPE(discretization), INTENT(IN), OPTIONAL :: ds
        END SUBROUTINE vector_pde_div

        MODULE SUBROUTINE flux_pde_div(sign,pde,flux,side,ds)
            IMPLICIT NONE
            CHARACTER(LEN=1),     INTENT(IN) :: sign
            TYPE(scalar_pde),     INTENT(INOUT) :: pde
            TYPE(scalar_field),   INTENT(IN) :: flux
            REAL(psb_dpk_),       INTENT(IN), OPTIONAL :: side
            TYPE(discretization), INTENT(IN), OPTIONAL :: ds
        END SUBROUTINE flux_pde_div

    END INTERFACE

END MODULE op_div
