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
! $Id$
!
! Description:
!    Checks that a tet is valid using a scheme that is consistent with how OptMS does it.  Can only be used for strictly local tetrahedral cells
!
FUNCTION tet_valid(msh,ic)

    USE class_psblas
    USE class_connectivity
    USE class_mesh
    USE class_vector
    USE class_vertex
    USE tools_mesh_optimize, ONLY: right_handed


    IMPLICIT NONE
    !
    LOGICAL                :: tet_valid
    TYPE(mesh), INTENT(IN) :: msh            ! the mesh structure
    INTEGER, INTENT(IN)    :: ic             ! the cell ID number
    !
    REAL(psb_dpk_)      :: vtx1(3),vtx2(3),vtx3(3),vtx4(3)
    INTEGER ::  iv1, iv2, iv3, iv4
    INTEGER, POINTER :: iv2c(:) => NULL()
    INTEGER :: valid_flag,ierr

    ! Get vertex indices
    CALL get_ith_conn(iv2c,msh%v2c,ic)
    iv1 = iv2c(1)
    iv2 = iv2c(2)
    iv3 = iv2c(3)
    iv4 = iv2c(4)

    vtx1 = position_( msh%verts(iv1) )
    vtx2 = position_( msh%verts(iv2) )
    vtx3 = position_( msh%verts(iv3) )
    vtx4 = position_( msh%verts(iv4) )

    valid_flag = right_handed(vtx1, vtx2, vtx3, vtx4)

    IF ( valid_flag == 1) THEN
        tet_valid = .TRUE.
    ELSE
        tet_valid = .FALSE.
    ENDIF

END FUNCTION tet_valid
