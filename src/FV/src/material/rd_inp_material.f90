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
! $Id: rd_inp_material.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(tools_material) rd_inp_material_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE rd_inp_material
            USE class_psblas
            USE json_module
            USE tools_input
            USE tools_material, ONLY: irho, imu, ilambda, ish

            IMPLICIT NONE
            !
            LOGICAL, PARAMETER :: debug = .FALSE.
            !
            TYPE(json_file) :: nemo_json
            CHARACTER(LEN=100) :: str
            LOGICAL :: found
            INTEGER :: mypnum, icontxt
            INTEGER :: inp


            icontxt = icontxt_()
            mypnum  = mypnum_()


            ! Reads parameters on P0 and performs consistency checks
            IF(mypnum == 0) THEN

                WRITE(*,*) 'Reading MATERIAL section from ',TRIM(input_file)," ", sec

                CALL open_file(input_file,nemo_json)
                str = TRIM('MORFEUS_FV.MATERIALS.'//TRIM(sec))
                CALL nemo_json%get(str, id, found)
                ! Gets material name from input file
                name = 'copper'
                type = 'default'

                ! material%ilaw(irho) -> density
                ilaw(irho) = 1

                ! material%ilaw(imu) -> viscosity
                ilaw(imu) = 1

                ! material%ilaw(ilambda) -> thermal conductivity
                ilaw(ilambda) = 1

                ! material%ilaw(ish) -> specific heat
                ilaw(ish) = 1

                WRITE(*,*)
            END IF


            ! Broadcasting
            IF(mypnum == 0) THEN
                ! Sends CHARACTER type members
                CALL psb_bcast(icontxt,name)
                CALL psb_bcast(icontxt,type)

                ! Sends INTEGER type members
                CALL psb_bcast(icontxt,ilaw)
                CALL psb_bcast(icontxt,id)
            ELSE
                ! Receives CHARACTER type members
                CALL psb_bcast(icontxt,name)
                CALL psb_bcast(icontxt,type)

                ! Receives INTEGER type members
                CALL psb_bcast(icontxt,ilaw)
                CALL psb_bcast(icontxt,id)
            END IF


            IF(debug) THEN
                WRITE(*,*)
                WRITE(*,800) mypnum_()
                WRITE(*,810) 'material%name  = ', name
                WRITE(*,810) 'material%mat_type  = ', type
                WRITE(*,820) 'material%ilaw  = ', ilaw(:)
                WRITE(*,830) 'material%mat_id  = ', id
                WRITE(*,*)
            END IF

        100 FORMAT(' ERROR! Illegal value of MATERIAL%ILAW(IRHO)')
        200 FORMAT(' ERROR! Illegal value of MATERIAL%ILAW(IMU)')
        300 FORMAT(' ERROR! Illegal value of MATERIAL%ILAW(ILAMBDA)')
        400 FORMAT(' ERROR! Illegal value of MATERIAL%ILAW(ISH)')

        800 FORMAT(' ----- Process ID = ',i2,' -----')
        810 FORMAT(1x,a,a)
        820 FORMAT(1x,a,4(1x,i1))
        830 FORMAT(1x,a,i5)

        END PROCEDURE rd_inp_material

END SUBMODULE rd_inp_material_implementation
