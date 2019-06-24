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
! $Id: op_grad.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    To be added...
!
MODULE op_grad
    USE class_vector_pde, ONLY : vector_pde
      !! Needed only by vector_pde_grad; a gfortran bug requires locating it here
    USE class_vector, ONLY : vector
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_vector_field, ONLY : vector_field
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_scalar_field, ONLY : scalar_field
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_psblas, ONLY : psb_dpk_
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: fld_grad, pde_grad

    INTERFACE fld_grad
        PROCEDURE :: scalar_field_grad
        PROCEDURE :: vector_field_grad
    END INTERFACE fld_grad

    INTERFACE pde_grad
        PROCEDURE :: vector_pde_grad
    END INTERFACE pde_grad

    INTERFACE

        MODULE FUNCTION scalar_field_grad(phi)RESULT(grad)
            IMPLICIT NONE
            TYPE(vector), ALLOCATABLE :: grad(:)
            TYPE(scalar_field), INTENT(IN) :: phi
        END FUNCTION scalar_field_grad

        MODULE FUNCTION vector_field_grad(phi)RESULT(grad)
            IMPLICIT NONE
            TYPE(vector), ALLOCATABLE :: grad(:,:)
            TYPE(vector_field), INTENT(IN) :: phi
        END FUNCTION vector_field_grad

        MODULE SUBROUTINE vector_pde_grad(sign,pde,phi,side)
            IMPLICIT NONE
            CHARACTER(len=1),   INTENT(IN)           :: sign
            CLASS(vector_pde),   INTENT(INOUT)        :: pde
            TYPE(scalar_field), INTENT(IN)           :: phi
            REAL(psb_dpk_),   INTENT(IN), OPTIONAL :: side
        END SUBROUTINE vector_pde_grad

    END INTERFACE

END MODULE op_grad
