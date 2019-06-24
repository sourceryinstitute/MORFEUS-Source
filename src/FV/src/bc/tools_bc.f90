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
! $Id: tools_bc.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
MODULE tools_bc

    IMPLICIT NONE

    INTERFACE
        MODULE SUBROUTINE rd_inp_bc(input_file,sec,nbc_msh,id,mot)
            USE class_motion
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            INTEGER, INTENT(IN) :: nbc_msh
            INTEGER, INTENT(INOUT) :: id(nbc_msh)
            TYPE(motion), INTENT(INOUT) :: mot(nbc_msh)
        END SUBROUTINE rd_inp_bc
    END INTERFACE


    INTERFACE
        MODULE SUBROUTINE rd_inp_bc_math(input_file,sec,nbf,id,a,b,c)
            USE class_psblas, ONLY : psb_dpk_
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            INTEGER, INTENT(IN) :: nbf
            INTEGER, INTENT(OUT) :: id
            REAL(psb_dpk_), ALLOCATABLE, INTENT(OUT)  :: a(:), b(:), c(:)
        END SUBROUTINE rd_inp_bc_math
    END INTERFACE


    ! ----- Named Constants -----

    ! BC IDs
    INTEGER, PARAMETER :: bc_math_ = 1
    INTEGER, PARAMETER :: bc_wall_ = 2

    ! Math BC IDs (also used in input files)
    INTEGER, PARAMETER :: bc_dirichlet_        = 1
    INTEGER, PARAMETER :: bc_neumann_          = 2
    INTEGER, PARAMETER :: bc_robin_            = 3
    INTEGER, PARAMETER :: bc_dirichlet_map_    = 4
    INTEGER, PARAMETER :: bc_neumann_map_      = 5
    INTEGER, PARAMETER :: bc_robin_map_        = 6
    INTEGER, PARAMETER :: bc_neumann_flux_     = 7
    INTEGER, PARAMETER :: bc_robin_convection_ = 8

    ! Temperature BC IDs (also used in input files)
    INTEGER, PARAMETER :: bc_temp_ = 1
    INTEGER, PARAMETER :: bc_temp_fixed_           = 1
    INTEGER, PARAMETER :: bc_temp_adiabatic_       = 2
    INTEGER, PARAMETER :: bc_temp_flux_            = 3
    INTEGER, PARAMETER :: bc_temp_convection_      = 4
    INTEGER, PARAMETER :: bc_temp_convection_map_  = 5

    ! Concentration BC IDs (also used in input files)
    INTEGER, PARAMETER :: bc_conc_ = 2
    INTEGER, PARAMETER :: bc_conc_fixed_           = 1
    INTEGER, PARAMETER :: bc_conc_adiabatic_       = 2

    ! Velocity BC IDs (also used in input files)
    INTEGER, PARAMETER :: bc_vel_ = 3
    INTEGER, PARAMETER :: bc_vel_no_slip_        = 1
    INTEGER, PARAMETER :: bc_vel_free_slip_      = 2
    INTEGER, PARAMETER :: bc_vel_sliding_        = 3
    INTEGER, PARAMETER :: bc_vel_moving_         = 4
    INTEGER, PARAMETER :: bc_vel_free_sliding_   = 5

    ! Position BC IDs (also used in input files)
    INTEGER, PARAMETER :: bc_pos_ = 4
    INTEGER, PARAMETER :: bc_pos_stationary_     = 1
    INTEGER, PARAMETER :: bc_pos_moving_         = 2

    ! Stress BC IDs (also used in input files)
    INTEGER, PARAMETER :: bc_stress_ = 5
    INTEGER, PARAMETER :: bc_stress_free_           = 1
    INTEGER, PARAMETER :: bc_stress_prescribed_     = 2

    ! Pressure BC IDs (also used in input files)
    INTEGER, PARAMETER :: bc_pressure_ = 6
    INTEGER, PARAMETER :: bc_pressure_free_           = 1
    INTEGER, PARAMETER :: bc_pressure_prescribed_     = 2

END MODULE tools_bc
