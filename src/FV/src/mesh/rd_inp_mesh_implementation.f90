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
SUBMODULE (tools_mesh) rd_inp_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE rd_inp_mesh
            USE class_psblas
            USE json_module
            USE tools_input
            IMPLICIT NONE
            !! $Id: rd_inp_mesh.f90 3093 2008-04-22 14:51:09Z sfilippo $
            !!
            !! Description:
            !!    To be added...
            !!
            TYPE(json_file) :: nemo_json
            INTEGER, PARAMETER :: nlen = 80
            !
            LOGICAL, PARAMETER :: debug = .FALSE.
            !
            CHARACTER(len=nlen) :: mesh_file_1, mesh_file_2
            CHARACTER(KIND=json_CK,LEN=:),ALLOCATABLE :: cval
            LOGICAL :: found
            INTEGER :: icontxt, mypnum
            INTEGER :: inp, intbuf(3)

            icontxt = icontxt_()
            mypnum = mypnum_()

            ! First reads parameters on P0...
            IF(mypnum == 0) THEN

                WRITE(*,*) 'Reading MESH section from ', TRIM(input_file)

                CALL open_file(input_file,nemo_json)
                CALL nemo_json%get('MORFEUS_FV.MESH.mesh-dir',cval,found)
                mesh_file_1 = trim(cval)
                CALL nemo_json%get('MORFEUS_FV.MESH.mesh-file',cval,found)
                mesh_file_2 = trim(cval)
                mesh_file=TRIM(mesh_file_1)//TRIM(mesh_file_2)

                CALL nemo_json%get('MORFEUS_FV.MESH.scale',scale,found)
                CALL nemo_json%get('MORFEUS_FV.MESH.irenum',irenum,found)
                CALL nemo_json%get('MORFEUS_FV.MESH.ipart',ipart,found)
                CALL nemo_json%get('MORFEUS_FV.MESH.nswpref',nswpref,found)
                CALL nemo_json%get('MORFEUS_FV.MESH.mtx_pat',mtx_pat,found)

                WRITE(*,*)
            END IF


            ! ... then bcast_ them to all other processes.
            IF(mypnum == 0) THEN
                ! Sends CHARACTER kind variables
                CALL psb_bcast(icontxt,mesh_file)

                ! Sends INTEGER type variables
                intbuf(1) = irenum
                intbuf(2) = ipart
                intbuf(3) = nswpref
                CALL psb_bcast(icontxt,intbuf)

                ! Sends LOGICAL type variables
                CALL psb_bcast(icontxt,mtx_pat)

                ! Sends REAL(8) type variables
                CALL psb_bcast(icontxt,scale)

            ELSE

                ! Receives CHARACTER type variables
                CALL psb_bcast(icontxt,mesh_file)

                ! Receives INTEGER type variables
                CALL psb_bcast(icontxt,intbuf)
                irenum     = intbuf(1)
                ipart      = intbuf(2)
                nswpref    = intbuf(3)

                ! Receives LOGICAL type variables
                CALL psb_bcast(icontxt,mtx_pat)

                ! Receives REAL(8) type variables
                CALL psb_bcast(icontxt,scale)
            END IF


            ! ----- Debug -----
            IF(debug) THEN
                WRITE(*,*)
                WRITE(*,100) mypnum_()
                WRITE(*,200) 'mesh%mesh_file   = ', mesh_file
                WRITE(*,300) 'mesh%scale      = ', scale
                WRITE(*,400) 'mesh%irenum     = ', irenum
                WRITE(*,400) 'mesh%ipart      = ', ipart
                WRITE(*,400) 'mesh%nswpref    = ', nswpref
                WRITE(*,500) 'mesh%mtx_pat     = ', mtx_pat
                WRITE(*,*)
            END IF
            ! -----------------


100         FORMAT(' ----- Process ID = ',i2,' -----')
200         FORMAT(1x,a,a)
300         FORMAT(1x,a,es10.3)
400         FORMAT(1x,a,i5)
500         FORMAT(1x,a,l5)

        END PROCEDURE rd_inp_mesh

END SUBMODULE rd_inp_implementation
