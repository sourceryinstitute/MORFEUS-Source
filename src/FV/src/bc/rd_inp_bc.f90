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
! $Id: rd_inp_bc.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(tools_bc) rd_inp_bc_implementation
    USE class_motion, only : motion
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE rd_inp_bc
        USE class_psblas, ONLY : psb_bcast, icontxt_, mypnum_, abort_psblas
        USE json_module
        USE tools_bc, ONLY: bc_math_, bc_wall_, bc_vel_moving_
        USE tools_input, ONLY : get_par, find_section, open_file
        USE tools_mesh_move, ONLY : moving_, stationary_, sticky_

        IMPLICIT NONE
        !
        TYPE(json_file) :: nemo_json
        TYPE(json_core) :: core
        CHARACTER(KIND=json_CK, LEN=:), DIMENSION(:), ALLOCATABLE :: bcs
        CHARACTER(KIND=json_CK,LEN=:),ALLOCATABLE :: cval, name
        CHARACTER(len=15) :: par
        CHARACTER(LEN=2) :: ib_string
        CHARACTER(len=32) :: bc_sec, bctype
        CHARACTER(len=13)  :: aformat
        LOGICAL :: found
        INTEGER(json_IK) :: vartype, nr
        INTEGER :: icontxt, mypnum
        INTEGER :: ib, inp, sec_len, nbc_inp
        INTEGER :: id_vel
        INTEGER :: digits
        INTEGER, DIMENSION(nbc_msh) :: surface_motion, vertex_motion
        CHARACTER(len=80), DIMENSION(nbc_msh) :: motion_law_file


        icontxt = icontxt_()
        mypnum  = mypnum_()

        ! Preliminary check
        IF(mypnum == 0) THEN

            CALL open_file(input_file,nemo_json)

            sec_len = LEN(TRIM(sec))

            ! Counts BC sections in external unit INP


            ! Checks NBC_MSH vs. NBC_INP
            ! IF(nbc_msh /= nbc_inp) THEN
            !     WRITE(*,100)
            !     CALL abort_psblas
            ! END IF
        END IF

        ! First reads IDs and MOTION stuff on P0...
        IF(mypnum == 0) THEN
            bc_loop: DO ib = 1, nbc_msh

                WRITE(ib_string, '(i0)') ib
                bc_sec = "MORFEUS_FV.BCS("//trim(ib_string)//")"

                !CALL find_section(bc_sec,nemo_json)
                CALL nemo_json%get(trim(bc_sec)//'.type', cval, found)
                IF (.NOT.found) THEN
                    bctype = 'wall'
                ELSE
                    bctype = TRIM(cval)
                END IF

                ! ----- POLYMORPHISM -----
                SELECT CASE(bctype)
                CASE('math')
                    id(ib) = bc_math_
                CASE('wall')
                    id(ib) = bc_wall_
                CASE default
                    WRITE(*,200)
                    CALL abort_psblas
                END SELECT
                ! ------------------------


                ! Gets motion condition for boundary vertices
                !vertex_motion(ib) = get_par(inp,bc_sec,par='vertex_motion',&
                !    default=sticky_)

                ! ----- Gets motion condition for boundary surface -----

                ! Defaults
                surface_motion(ib) = stationary_
                motion_law_file(ib) = ''
                !
                ! Wall BC
                IF(id(ib) == bc_wall_) THEN

                    !scan_for_vel: DO
                        ! READ(inp,'(a15)') par
                        ! par = TRIM(par)

                        ! BACKSPACE(inp)
                        ! IF(par == 'velocity') THEN
                        !     READ(inp,'(a,i1)',advance='no') par, id_vel

                        !     IF(id_vel == bc_vel_moving_) THEN
                        !         surface_motion(ib) = moving_
                        !         READ(inp,*) motion_law_file(ib)
                        !         EXIT scan_for_vel
                        !     END IF
                        ! ELSE IF(par == 'END OF SECTION') THEN
                        !     EXIT scan_for_vel
                        ! ELSE
                        !     READ(inp,'()')
                        ! END IF

                    !END DO scan_for_vel
                END IF

            END DO bc_loop

            CLOSE(inp)
        END IF


        ! ... then bcast_ it to all other processes.
        CALL psb_bcast(icontxt,id)
        CALL psb_bcast(icontxt,surface_motion)
        CALL psb_bcast(icontxt,vertex_motion)
        CALL psb_bcast(icontxt,motion_law_file)


        ! Eventually create the motion object according to file inputs
        DO ib = 1, nbc_msh
            CALL mot(ib)%create_motion(surface_motion(ib),vertex_motion(ib),&
                & TRIM(motion_law_file(ib)))
        END DO


100     FORMAT(' ERROR! NBC_MSH doesn''t match NBC_INP in RD_INP_BC')
200     FORMAT(' ERROR! Unsupported BC in RD_INP_BC')

        END PROCEDURE rd_inp_bc

END SUBMODULE rd_inp_bc_implementation
