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
! $Id: op_laplacian.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    To be added...
!
MODULE op_laplacian
    USE class_scalar_pde, ONLY : scalar_pde
    !! A gfortran 8.3.0 bug precludes putting this in the scalar_pde_laplacian_phi and
    !! scalar_pde_laplacian_gamma_phi interface bodies
    USE class_vector_pde, ONLY : vector_pde
    !! A gfortran 8.3.0 bug precludes putting this in the vector_pde_laplacian_gamma_phi and
    !! vector_pde_laplacian_gamma_phi interface bodies
    USE class_vector_field, ONLY : vector_field
    !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_scalar_field, ONLY : scalar_field
    !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_psblas, ONLY : psb_dpk_
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: pde_laplacian

    INTERFACE pde_laplacian
        MODULE PROCEDURE :: scalar_pde_laplacian_phi
        MODULE PROCEDURE :: scalar_pde_laplacian_gamma_phi
        MODULE PROCEDURE :: vector_pde_laplacian_phi
        MODULE PROCEDURE :: vector_pde_laplacian_gamma_phi
    END INTERFACE pde_laplacian
    ! ----- Public Generic Interfaces -----

    INTERFACE

        ! ----- Laplacian Operator Wrappers -----

        MODULE SUBROUTINE scalar_pde_laplacian_phi(sign,pde,phi,side)
            IMPLICIT NONE
            CHARACTER(len=1),   INTENT(IN) :: sign
            TYPE(scalar_pde),   INTENT(INOUT) :: pde
            TYPE(scalar_field), INTENT(IN) :: phi
            REAL(psb_dpk_),   INTENT(IN), OPTIONAL :: side
        END SUBROUTINE scalar_pde_laplacian_phi

        MODULE SUBROUTINE scalar_pde_laplacian_gamma_phi(sign,pde,gamma,phi,side)
            IMPLICIT NONE
            CHARACTER(len=1),   INTENT(IN) :: sign
            TYPE(scalar_pde),   INTENT(INOUT) :: pde
            TYPE(scalar_field), INTENT(IN) :: gamma
            TYPE(scalar_field), INTENT(IN) :: phi
            REAL(psb_dpk_),   INTENT(IN), OPTIONAL :: side
        END SUBROUTINE scalar_pde_laplacian_gamma_phi

        MODULE SUBROUTINE vector_pde_laplacian_phi(sign,pde,phi,side)
            IMPLICIT NONE
            CHARACTER(len=1),   INTENT(IN) :: sign
            TYPE(vector_pde),   INTENT(INOUT) :: pde
            TYPE(vector_field), INTENT(IN) :: phi
            REAL(psb_dpk_),   INTENT(IN), OPTIONAL :: side
        END SUBROUTINE vector_pde_laplacian_phi

        MODULE SUBROUTINE vector_pde_laplacian_gamma_phi(sign,pde,gamma,phi,side)
            IMPLICIT NONE
            CHARACTER(len=1),   INTENT(IN) :: sign
            TYPE(vector_pde),   INTENT(INOUT) :: pde
            TYPE(scalar_field), INTENT(IN) :: gamma
            TYPE(vector_field), INTENT(IN) :: phi
            REAL(psb_dpk_),   INTENT(IN), OPTIONAL :: side
        END SUBROUTINE vector_pde_laplacian_gamma_phi

        ! ----- Private Interfaces -----

        MODULE SUBROUTINE scalar_pde_laplacian(sign,pde,phi,gamma,side)
            IMPLICIT NONE
            CHARACTER(len=1),   INTENT(IN) :: sign
            TYPE(scalar_pde),   INTENT(INOUT) :: pde
            TYPE(scalar_field), INTENT(IN) :: phi
            TYPE(scalar_field), INTENT(IN), OPTIONAL :: gamma
            REAL(psb_dpk_),   INTENT(IN), OPTIONAL :: side
        END SUBROUTINE scalar_pde_laplacian

        MODULE SUBROUTINE vector_pde_laplacian(sign,pde,phi,gamma,side)
            IMPLICIT NONE
            CHARACTER(len=1),   INTENT(IN) :: sign
            TYPE(vector_pde),   INTENT(INOUT) :: pde
            TYPE(vector_field), INTENT(IN) :: phi
            TYPE(scalar_field), INTENT(IN), OPTIONAL :: gamma
            REAL(psb_dpk_),   INTENT(IN), OPTIONAL :: side
        END SUBROUTINE vector_pde_laplacian

    END INTERFACE

END MODULE op_laplacian
