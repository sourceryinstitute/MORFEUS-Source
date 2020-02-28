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
! $Id: class_psblas.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Class with static members (i.e. global variables) shared by all other
!    units: timing, etc.
!
MODULE class_psblas

    USE class_stopwatch, ONLY: stopwatch, stopwatch_
    USE tools_psblas

    IMPLICIT NONE

    INTEGER, SAVE       :: nemo_sizeof_dp, nemo_sizeof_int, nemo_sizeof_long_int

    ! The default access has been left PUBLIC so that the usage of this
    ! module also implies public access for the PSB_SPARSE_MOD module
    ! loaded by TOOLS_PSBLAS

    LOGICAL, PARAMETER :: debug_mat_bld=.FALSE.

    ! variables kept PRIVATE in order to avoid accidental access
    LOGICAL, PRIVATE, SAVE :: psblas_on = .FALSE.
    INTEGER, PRIVATE, SAVE :: icontxt, mypnum, nprocs
    INTEGER, PRIVATE, SAVE :: myprow, mypcol, nprows, npcols

    ! stopwatches
    TYPE(stopwatch), SAVE :: sw_msh, sw_ord, sw_par, sw_dsc, sw_g2l, &
        &                   sw_geo, sw_lsr, sw_ins, sw_asb, sw_pre, &
        &                   sw_sol, sw_pde, sw_out, sw_tot

    ! ----- Constructor -----

    INTERFACE
        MODULE SUBROUTINE start_psblas
            IMPLICIT NONE
        END SUBROUTINE start_psblas
    END INTERFACE


    ! ----- Destructors -----

    INTERFACE

    ! Normal termination
        MODULE SUBROUTINE stop_psblas(time)
            IMPLICIT NONE
            LOGICAL, OPTIONAL :: time
        END SUBROUTINE stop_psblas

    ! Immediate abortion
        MODULE SUBROUTINE abort_psblas
            IMPLICIT NONE
        END SUBROUTINE abort_psblas
    END INTERFACE

    ! ----- Getters -----

    INTERFACE
        MODULE FUNCTION psblas_is_on()
            IMPLICIT NONE
            LOGICAL :: psblas_is_on
        END FUNCTION psblas_is_on


        MODULE FUNCTION icontxt_()
            IMPLICIT NONE
            INTEGER :: icontxt_
        END FUNCTION icontxt_

        MODULE FUNCTION mypnum_()
            IMPLICIT NONE
            INTEGER :: mypnum_
        END FUNCTION mypnum_

        MODULE FUNCTION nprocs_()
            IMPLICIT NONE
            INTEGER :: nprocs_
        END FUNCTION nprocs_
    END INTERFACE

    ! ----- Utilities -----

    INTERFACE
        MODULE SUBROUTINE stop_timing
            IMPLICIT NONE
        END SUBROUTINE stop_timing
    END INTERFACE

END MODULE class_psblas
