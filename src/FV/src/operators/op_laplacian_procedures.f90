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
SUBMODULE(op_laplacian) op_laplacian_procedures
    IMPLICIT NONE

CONTAINS

    ! ----- Laplacian Operator Wrappers -----

    MODULE PROCEDURE scalar_pde_laplacian_phi

        CALL scalar_pde_laplacian(sign,pde,phi=phi,side=side)

    END PROCEDURE scalar_pde_laplacian_phi


    MODULE PROCEDURE scalar_pde_laplacian_gamma_phi

        CALL scalar_pde_laplacian(sign,pde,phi=phi,gamma=gamma,side=side)

    END PROCEDURE scalar_pde_laplacian_gamma_phi


    MODULE PROCEDURE vector_pde_laplacian_phi

        CALL vector_pde_laplacian(sign,pde,phi=phi,side=side)

    END PROCEDURE vector_pde_laplacian_phi


    MODULE PROCEDURE vector_pde_laplacian_gamma_phi

        CALL vector_pde_laplacian(sign,pde,phi=phi,gamma=gamma,side=side)

    END PROCEDURE vector_pde_laplacian_gamma_phi

END SUBMODULE op_laplacian_procedures
