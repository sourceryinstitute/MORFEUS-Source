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
! $Id: tools_part.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!
MODULE tools_part

    USE class_connectivity
    USE part_block
    USE part_graph
    USE part_random

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: part_cells, c2v, c2f
    PUBLIC :: part_verts, part_faces

    INTEGER, ALLOCATABLE, SAVE :: part_cells(:)
    TYPE(connectivity), SAVE :: c2v, c2f

    ! PART_CELLS is created in CMP_MESH_PART
    ! C2V and C2F are created in CMP_MESH_DESC

    ! PART_CELLS, C2V and C2F are destroyed in CMP_MESH_DESC

    INTERFACE

        MODULE SUBROUTINE part_verts(iv,nverts,nprocs,pv,nv)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: iv, nverts, nprocs
            INTEGER, INTENT(OUT) :: pv(*)
            INTEGER, INTENT(OUT) :: nv
        END SUBROUTINE part_verts


        MODULE SUBROUTINE part_faces(IF,nfaces,nprocs,pv,nv)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: IF, nfaces, nprocs
            INTEGER, INTENT(OUT) :: pv(*)
            INTEGER, INTENT(OUT) :: nv
        END SUBROUTINE part_faces

        ! REMARK-1: the interface of PART_VERTS and PART_FACES, needed for
        ! allocating the PSBLAS descriptors, is fixed by PARTS.FH, in the
        ! PSBLAS library. Hence, the only way for carrying the C2V and C2F
        ! variables inside the two procedures is to handle them as public
        ! shared variables.

        ! REMARK-2: PV could be initialized equal to any negative value.
        ! The choice of NVERTS and NFACES is only for avoiding the detection
        ! of an unused variable.

    END INTERFACE

END MODULE tools_part
