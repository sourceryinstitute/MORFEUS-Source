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
! $Id: tools_material.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    To be added...
!
MODULE tools_material

    USE class_psblas, ONLY : psb_dpk_
    IMPLICIT NONE

    ! Physical properties IDs
    INTEGER, PARAMETER :: irho    = 1 ! Density
    INTEGER, PARAMETER :: imu     = 2 ! Viscosity
    INTEGER, PARAMETER :: ilambda = 3 ! Thermal Conductivity
    INTEGER, PARAMETER :: ish     = 4 ! Specific heat

    ! Standard ambient temperature and pressure
    REAL(psb_dpk_), PARAMETER :: std_temp_ = 298.15d0    ! Temperature [K]
    REAL(psb_dpk_), PARAMETER :: std_pres_ = 1.01325d+05 ! Pressure [Pa]


    ! ----- Interfaces -----

    INTERFACE
        MODULE SUBROUTINE rd_inp_material(input_file,sec,name,ilaw,TYPE,id)
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=*), INTENT(OUT) :: name,TYPE
            INTEGER, INTENT(OUT) :: ilaw(:),id
        END SUBROUTINE rd_inp_material
    END INTERFACE


    INTERFACE
        MODULE SUBROUTINE load_material(name,state,dtemp,tmin,tmax,rho,mu,lambda,sh)
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN)  :: name
            CHARACTER(len=1), INTENT(OUT) :: state
            REAL(psb_dpk_), INTENT(OUT) :: dtemp, tmin, tmax
            REAL(psb_dpk_), ALLOCATABLE, INTENT(OUT) :: rho(:), mu(:), lambda(:), sh(:)
        END SUBROUTINE load_material
    END INTERFACE


END MODULE tools_material
