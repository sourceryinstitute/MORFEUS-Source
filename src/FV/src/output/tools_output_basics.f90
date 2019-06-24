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
! $Id: tools_output_basics.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    To be added...
!
MODULE tools_output_basics

    IMPLICIT NONE

    ! ----- Matrix Market Format -----

    INTERFACE
        SUBROUTINE wr_mtx_matrix(A,desc,name)
            USE psb_base_mod
            TYPE(psb_dspmat_type), INTENT(IN) :: A
            TYPE(psb_desc_type),   INTENT(IN) :: desc
            CHARACTER(len=*),      INTENT(IN) :: name
        END SUBROUTINE wr_mtx_matrix
    END INTERFACE

    INTERFACE
        MODULE SUBROUTINE wr_mtx_pattern(c2c,name)
            USE class_connectivity
            IMPLICIT NONE
            TYPE(connectivity), INTENT(IN) :: c2c
            CHARACTER(len=*),   INTENT(IN) :: name
        END SUBROUTINE wr_mtx_pattern
    END INTERFACE

    INTERFACE
        SUBROUTINE wr_mtx_vector(loc_vect,desc,name)
            USE psb_base_mod
            REAL(psb_dpk_)                :: loc_vect(:)
            TYPE(psb_desc_type), INTENT(IN) :: desc
            CHARACTER(len=*),    INTENT(IN) :: name
        END SUBROUTINE wr_mtx_vector
    END INTERFACE

    ! ----- String Manipulation ----

    INTERFACE
        MODULE FUNCTION itoh(i,max)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i
            INTEGER, INTENT(IN) :: max
            CHARACTER(len=max)  :: itoh
        END FUNCTION itoh
    END INTERFACE


    INTERFACE
        MODULE FUNCTION htoi(h)
            IMPLICIT NONE
            INTEGER :: htoi
            CHARACTER(len=*), INTENT(IN) :: h
        END FUNCTION htoi
    END INTERFACE


    ! ----- Named Constants -----

    ! OUTPUT formats
    INTEGER, PARAMETER :: vtk_  = 2
    INTEGER, PARAMETER :: cgns_ = 3

END MODULE tools_output_basics
