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
SUBMODULE(tools_part) tools_part_procedures


    IMPLICIT NONE

    !INTEGER, ALLOCATABLE, SAVE :: part_cells(:)
    !TYPE(connectivity), SAVE :: c2v, c2f

    ! PART_CELLS is created in CMP_MESH_PART
    ! C2V and C2F are created in CMP_MESH_DESC

    ! PART_CELLS, C2V and C2F are destroyed in CMP_MESH_DESC

CONTAINS

    MODULE PROCEDURE part_verts
        USE tools_psblas
        !
        INTEGER :: i, ic, n, pid
        INTEGER :: ppv(nprocs)
        INTEGER, POINTER :: ic2v(:) => NULL()

        ppv(:) = -nverts

        CALL c2v%get_ith_conn(ic2v,iv)
        n = SIZE(ic2v)

        DO i = 1, n
            ic = ic2v(i)
            pid = part_cells(ic)
            ppv(pid+1) = pid
        END DO

        nv       = COUNT(ppv >= 0)
        ppv(1:nv) = PACK(ppv,ppv >= 0)
        IF (nv <= 0) THEN
            WRITE(0,*) 'part_verts: orphan vertex',iv,nverts,n
        END IF
        NULLIFY(ic2v)

        pv(1:nv) = ppv(1:nv)

    END PROCEDURE part_verts


    MODULE PROCEDURE part_faces
        !
        INTEGER :: i, ic, n, pid
        INTEGER, POINTER :: ic2f(:) => NULL()
        INTEGER :: ppv(nprocs)

        ppv(:) = -nfaces

        CALL c2f%get_ith_conn(ic2f,IF)
        n = SIZE(ic2f)

        DO i = 1, n
            ic = ic2f(i)
            pid = part_cells(ic)
            ppv(pid+1) = pid
        END DO

        nv        = COUNT(ppv >= 0)
        ppv(1:nv) = PACK(ppv,ppv >= 0)

        NULLIFY(ic2f)
        pv(1:nv) = ppv(1:nv)

    END PROCEDURE part_faces


    ! REMARK-1: the interface of PART_VERTS and PART_FACES, needed for
    ! allocating the PSBLAS descriptors, is fixed by PARTS.FH, in the
    ! PSBLAS library. Hence, the only way for carrying the C2V and C2F
    ! variables inside the two procedures is to handle them as public
    ! shared variables.

    ! REMARK-2: PV could be initialized equal to any negative value.
    ! The choice of NVERTS and NFACES is only for avoiding the detection
    ! of an unused variable.

END SUBMODULE tools_part_procedures
