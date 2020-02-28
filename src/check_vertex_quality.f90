!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
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
!  Description:
!     Provides the quality of a single vertex.
!     Looks at the cells that use a vertex and returns the quality of the worst cell.
!     This routine is designed for speed and does not use allocatable variables
!
!
SUBMODULE (tools_mesh_check) tools_mesh_check_vertex
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE check_vertex_quality
        USE class_psblas
        USE class_connectivity
        USE class_mesh
        USE tools_mesh_basics, ONLY : geom_tet_quality
        IMPLICIT NONE
        !
        INTEGER, PARAMETER :: cmax = 50 ! the max number of cells connected to a vertex
        INTEGER, POINTER :: iv2c(:) => NULL()
        INTEGER :: ic, iv1, iv2, iv3, iv4, neighbor
        INTEGER :: ncells               ! number of cells connected to this pointer
        REAL(psb_dpk_) :: qlist(cmax) ! list of qualities

        ncells = SIZE(ic2v)

        IF (ncells>cmax) THEN
            WRITE(*,100)
            CALL abort_psblas
        ENDIF


        DO neighbor = 1, ncells
            ic = ic2v(neighbor)
            CALL msh%v2c%get_ith_conn(iv2c,ic)
            iv1 = iv2c(1)
            iv2 = iv2c(2)
            iv3 = iv2c(3)
            iv4 = iv2c(4)
            qlist(neighbor) = geom_tet_quality(&
                & msh%verts(iv1),msh%verts(iv2),msh%verts(iv3),msh%verts(iv4), &
                & msh%vol(ic))
        ENDDO

        quality = MINVAL(qlist(1:ncells))

        NULLIFY(iv2c)

100     FORMAT(' ERROR! Too many cells connected to a vertex in CHECK_VERTEX_QUALITY.',&
            &    ' Check value of CMAX in CHECK_VERTEX_QUALITY.' )
        END PROCEDURE check_vertex_quality

END SUBMODULE tools_mesh_check_vertex
