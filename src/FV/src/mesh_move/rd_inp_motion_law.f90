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
! $Id: rd_inp_motion_law.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(tools_mesh_move) rd_inp_motion_law_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE rd_inp_motion_law
            USE class_psblas, ONLY : psb_dpk_, abort_psblas, mypnum_, icontxt_, psb_bcast
            USE class_vector
            USE tools_input
            IMPLICIT NONE
            !
            CHARACTER(len=32), PARAMETER :: sec = 'MOTION LAW'
            INTEGER :: i, info, inp, icontxt, mypnum, nsteps
            REAL(psb_dpk_) :: law_y1, law_y2, law_y3

            mypnum  = mypnum_()
            icontxt = icontxt_()

            IF(mypnum == 0) THEN
                CALL open_file(ml_file,inp)

                iml    = get_par(inp,sec=TRIM(ml_file),par='iml',default=ml_position_)
                nsteps = get_par(inp,sec=TRIM(ml_file),par='nsteps',default=mandatory_i_)

                IF(.NOT.(iml == ml_position_ .OR. &
                    &   iml == ml_velocity_)) THEN
                    WRITE(*,100)
                    CALL abort_psblas
                END IF

            END IF

            CALL psb_bcast(icontxt,iml)
            CALL psb_bcast(icontxt,nsteps)

            CALL alloc_vector(law_y,nsteps)
            ALLOCATE(law_x(nsteps),stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            IF(mypnum == 0) THEN
                CALL find_section(TRIM(sec),inp)

                ! Skip header
                READ(inp,'()')

                ! Reads motion law table
                DO i = 1, nsteps
                    READ(inp,*) law_x(i), law_y1, law_y2, law_y3
                    law_y(i) = vector_(law_y1,law_y2,law_y3)
                END DO

                CLOSE(inp)
            END IF

            CALL psb_bcast(icontxt,law_x)
            CALL bcast_vector(law_y)

100         FORMAT(' ERROR! Unsupported type of motion law in RD_INP_MOTION_LAW')
200         FORMAT(' ERROR! Memory allocation failure in RD_INP_MOTION_LAW')

        END PROCEDURE rd_inp_motion_law

END SUBMODULE rd_inp_motion_law_implementation
