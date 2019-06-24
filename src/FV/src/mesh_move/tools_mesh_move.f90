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
! $Id: tools_mesh_move.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    To be added...
!
MODULE tools_mesh_move

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: rd_inp_motion_law
    PUBLIC :: stationary_, moving_
    PUBLIC :: sticky_, sliding_
    PUBLIC :: ml_position_, ml_velocity_

    ! ----- Generic Interfaces -----

    INTERFACE

        MODULE SUBROUTINE rd_inp_motion_law(ml_file,iml,law_x,law_y)
            USE class_vector, ONLY : vector
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN) :: ml_file
            INTEGER, INTENT(OUT) :: iml
            REAL(psb_dpk_), ALLOCATABLE  :: law_x(:)
            TYPE(vector), ALLOCATABLE    :: law_y(:)
        END SUBROUTINE rd_inp_motion_law

    END INTERFACE

    ! ----- Named Constants -----

    ! Flags for moving surfaces
    INTEGER, PARAMETER :: stationary_ = 1
    INTEGER, PARAMETER :: moving_     = 2

    ! Flags for moving vertices
    ! Vertices may move even on stationary surfaces, so this flag describes
    ! the relative motion of nodes
    INTEGER, PARAMETER :: sticky_  = 1
    INTEGER, PARAMETER :: sliding_ = 2

    ! Flags for motion law
    INTEGER, PARAMETER :: ml_position_ = 1
    INTEGER, PARAMETER :: ml_velocity_ = 2

END MODULE tools_mesh_move
