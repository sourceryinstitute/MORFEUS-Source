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
! $Id: op_d2dt2.f90 2469 2018-07-18 00:00:00Z hariradh $
!
! Description:
!    To be added...
!
SUBMODULE(op_d2dt2) op_d2dt2_procedures

    IMPLICIT NONE

CONTAINS

    ! ----- Time Derivative Operator Wrappers -----

    MODULE PROCEDURE scalar_pde_d2dt2_phi
        IMPLICIT NONE

        CALL scalar_pde_d2dt2(sign,pde,phi,dt,side=side)

    END PROCEDURE scalar_pde_d2dt2_phi


    MODULE PROCEDURE scalar_pde_d2dt2_fld_phi
        IMPLICIT NONE

        CALL scalar_pde_d2dt2(sign,pde,phi,dt,fld=fld,side=side)

    END PROCEDURE scalar_pde_d2dt2_fld_phi


    MODULE PROCEDURE vector_pde_d2dt2_phi
        IMPLICIT NONE

        CALL vector_pde_d2dt2(sign,pde,phi,dt,side=side)

    END PROCEDURE vector_pde_d2dt2_phi


    MODULE PROCEDURE vector_pde_d2dt2_fld_phi
        IMPLICIT NONE

        CALL vector_pde_d2dt2(sign,pde,phi,dt,fld=fld,side=side)

    END PROCEDURE vector_pde_d2dt2_fld_phi

END SUBMODULE op_d2dt2_procedures
