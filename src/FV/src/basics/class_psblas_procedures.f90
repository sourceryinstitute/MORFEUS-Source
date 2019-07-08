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
SUBMODULE(class_psblas) class_psblas_procedures

    USE class_stopwatch
    USE tools_psblas

    IMPLICIT NONE

CONTAINS

    ! ----- Constructor -----

    MODULE PROCEDURE start_psblas
        INTEGER :: ip, ipnum

        IF (psblas_on) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_set_errverbosity(2)
        CALL psb_set_erraction(1)

        ! Initializes
        CALL psb_init(icontxt)
        CALL psb_info(icontxt,mypnum,nprocs)

        ! Sets psblas_on flag on TRUE
        psblas_on = .TRUE.
        nemo_sizeof_dp       = psb_sizeof_dp
        nemo_sizeof_int      = psb_sizeof_int
        ! WARNING: to be reviewed
        nemo_sizeof_long_int = 2*psb_sizeof_int
        ! Initializes all stopwatches
        sw_msh = stopwatch_(icontxt)
        sw_ord = stopwatch_(icontxt)
        sw_par = stopwatch_(icontxt)
        sw_dsc = stopwatch_(icontxt)
        sw_g2l = stopwatch_(icontxt)
        sw_geo = stopwatch_(icontxt)
        sw_ins = stopwatch_(icontxt)
        sw_asb = stopwatch_(icontxt)
        sw_pre = stopwatch_(icontxt)
        sw_sol = stopwatch_(icontxt)
        sw_lsr = stopwatch_(icontxt)
        sw_pde = stopwatch_(icontxt)
        sw_out = stopwatch_(icontxt)
        sw_tot = stopwatch_(icontxt)

        CALL psb_barrier(icontxt)

        ! Starts timing of parallel job
        CALL sw_tot%tic()

        ! Echoes status of processes
        IF(mypnum == 0) THEN
            WRITE(*,'()')
            WRITE(*,200) mypnum, nprocs
        END IF
        DO ip = 1, nprocs - 1
            IF(mypnum == ip) THEN
                CALL psb_snd(icontxt,mypnum,0)
            ELSEIF(mypnum == 0) THEN
                CALL psb_rcv(icontxt,ipnum,ip)
                WRITE(*,200) ipnum, nprocs
            END IF
        END DO

        IF(mypnum == 0) WRITE(*,'()')

100     FORMAT(' ERROR! The PSBLAS job has already been started')
200     FORMAT(' Process',i3,' on',i3,' active')

    END PROCEDURE start_psblas


    ! ----- Destructors -----

    ! Normal termination
    MODULE PROCEDURE stop_psblas
        !
        INTEGER :: ip, ipnum
        LOGICAL :: time_

        IF(PRESENT(time)) THEN
            time_ = time
        ELSE
            time_ = .TRUE.
        END IF

        IF(.NOT. psblas_on) THEN
            WRITE(*,100)
            STOP
        END IF

        IF(time_) CALL stop_timing

        ! Echoes Normal Termination
        IF(mypnum == 0) WRITE(*,200) mypnum, nprocs
        DO ip = 1, nprocs - 1
            IF(mypnum == ip) THEN
                CALL psb_snd(icontxt,mypnum,0)
            ELSEIF(mypnum == 0) THEN
                CALL psb_rcv(icontxt,ipnum,ip)
                WRITE(*,200) ipnum, nprocs
            END IF
        END DO
        IF(mypnum == 0) WRITE(*,'()')

        ! Exits from PSBLAS environment
        CALL psb_exit(icontxt,CLOSE=.TRUE.)

100     FORMAT(' ERROR! The PSBLAS job is not running: it can''t be stopped')
200     FORMAT(' Process',i3,' of',i3,': Normal Termination')

    END PROCEDURE stop_psblas


    ! Immediate abortion
    MODULE PROCEDURE abort_psblas
        CALL psb_abort(icontxt)
    END PROCEDURE abort_psblas


    ! ----- Getters -----

    MODULE PROCEDURE psblas_is_on
        psblas_is_on = psblas_on
    END PROCEDURE psblas_is_on


    MODULE PROCEDURE icontxt_
        icontxt_ = icontxt
    END PROCEDURE icontxt_


    MODULE PROCEDURE mypnum_
        mypnum_ = mypnum
    END PROCEDURE mypnum_


    MODULE PROCEDURE nprocs_
        nprocs_ = nprocs
    END PROCEDURE nprocs_


    ! ----- Utilities -----

    MODULE PROCEDURE stop_timing
        ! Stops timing of parallel job
        CALL sw_tot%toc()

        ! Synchronizes stopwatches
        CALL sw_msh%synchro()
        CALL sw_ord%synchro()
        CALL sw_par%synchro()
        CALL sw_dsc%synchro()
        CALL sw_g2l%synchro()
        CALL sw_geo%synchro()
        CALL sw_lsr%synchro()

        CALL sw_ins%synchro()
        CALL sw_asb%synchro()
        CALL sw_pre%synchro()
        CALL sw_sol%synchro()

        CALL sw_pde%synchro()
        CALL sw_out%synchro()
        CALL sw_tot%synchro()

        ! Time consumption log message
        IF(mypnum == 0) THEN
            WRITE(*,*)   '---------- PHASES TIMINGS ----------'
            WRITE(*,*)
            WRITE(*,*)   '* mesh:'
            WRITE(*,100) '  Importing:        ', sw_msh%total_(), ' s'
            WRITE(*,100) '  Reordering:       ', sw_ord%total_(), ' s'
            WRITE(*,100) '  Partitioning:     ', sw_par%total_(), ' s'
            WRITE(*,100) '  Descriptor bld.:  ', sw_dsc%total_(), ' s'
            WRITE(*,100) '  G2L reallocation: ', sw_g2l%total_(), ' s'
            WRITE(*,100) '  Geometry comp.:   ', sw_geo%total_(), ' s'
            WRITE(*,*)
            WRITE(*,*)   '* PSBLAS:'
            WRITE(*,100) '  Inserting:        ', sw_ins%total_(), ' s'
            WRITE(*,100) '  Assembling:       ', sw_asb%total_(), ' s'
            WRITE(*,100) '  Preconditioning:  ', sw_pre%total_(), ' s'
            WRITE(*,100) '  System solving:   ', sw_sol%total_(), ' s'
            WRITE(*,*)
            WRITE(*,*)   '* PDE operators:'
            WRITE(*,100) '  w   asb/ins:      ', sw_pde%total_(), ' s'
            WRITE(*,100) '  w/o asb/ins:      ', sw_pde%total_() - &
                &                               sw_ins%total_() - &
                &                               sw_asb%total_(), ' s'
            WRITE(*,100) '  Least squares:    ', sw_lsr%total_(), ' s'
            WRITE(*,*)
            WRITE(*,100) '* Results writing:  ', sw_out%total_(), ' s'
            WRITE(*,*)
            WRITE(*,100) '* Total running:    ', sw_tot%total_(), ' s'
            WRITE(*,*)
            WRITE(*,*)   '------------------------------------'
            WRITE(*,*)
        END IF

100     FORMAT(1x,a,es13.6,a)

    END PROCEDURE stop_timing


END SUBMODULE class_psblas_procedures
