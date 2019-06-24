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
!  $Id: mobile_verts.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
!  Description: Constructs a list of
!     (1) which vertices are free to move
!     (2) which vertices can move, but are constrained (sliding on a 2D surface)
!
SUBMODULE (tools_mesh_optimize) mobile_verts_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE mobile_verts
        USE class_psblas
        USE class_bc
        USE class_connectivity
        USE class_cell
        USE class_face
        USE class_keytable
        USE class_mesh
        USE class_vertex
        USE tools_mesh_basics
        USE tools_mesh_move
        IMPLICIT NONE

        ! we potentially waste some memory by allocating unconstrained, constrained,
        ! and shared to the size of nverts, but these are integer and logical
        ! arrays, so it is worth the savings in CPU time and coding

        ! ----- Local Variables -----

        TYPE(connectivity) :: f2v,b2v
        ! given vertices, find faces, and given vertices, find cells
        INTEGER,POINTER :: ic2v(:)   ! list of cells connected to a vertex
        INTEGER,POINTER :: iv2c(:)   ! list of cells connected to a vertex
        INTEGER,POINTER :: ib2v(:)   ! list of bndrys to which a vertex is connected
        INTEGER :: nverts            ! total number of vertices
        INTEGER :: iv,ic     ! looping indices for vertices, faces, or cells
        INTEGER :: i         ! loop index
        LOGICAL :: is_slip
        INTEGER :: bc_id        ! the id number of a bc used by this vertex
        INTEGER :: glob_iv      ! global index number
        INTEGER :: nfaces
        INTEGER :: size_first   ! size of first cell used by vertex

        nfaces=SIZE(msh%faces(:))

        ! 0. set up data structures
        n_unconstrained = 0
        n_constrained = 0

        ! we only consider local vertices for movement (incl. overlap)
        nverts = psb_cd_get_local_rows(msh%desc_v)

        CALL get_dual_conn(msh%v2f,f2v)
        CALL get_dual_conn(msh%v2b,b2v)

        ! 1.  Loop over all vertices and check their status.
        DO iv=1,nverts

            ! are strictly local vertices connected only to tets? or more than one kind of cell?

            IF (.NOT. shared_flag(iv)) THEN
                CALL get_ith_conn(ic2v,c2v,iv) ! list of cells connected to this vertex

                all_tets(iv) = .TRUE.
                mixed(iv)    = .FALSE.

                DO i = 1,SIZE(ic2v)
                    ic = ic2v(i)
                    CALL get_ith_conn(iv2c,msh%v2c,ic)
                    IF (SIZE(iv2c) /= 4) all_tets(iv)=.FALSE.

                    IF ( i == 1) size_first = SIZE(iv2c)
                    IF ( SIZE(iv2c) /= size_first ) mixed(iv) = .TRUE.
                ENDDO
            ELSE
                ! are shared vertices connected only to tets?

                ! the overlap vertex supplemental information is based on global numbering
                glob_iv = loc_to_glob_(msh%desc_v, iv)

                ! get cells attached to the vertex
                CALL get_kt_row(msh%c2ov_sup, glob_iv, ic2v)

                all_tets(iv)=.TRUE.
                mixed(iv)    = .FALSE.


                DO i=1,SIZE(ic2v)  ! if the cell has 4 vertices, we presume it is a tet
                    ic=ic2v(i)
                    CALL get_kt_row(msh%ov2c_sup, ic, iv2c)
                    IF (SIZE(iv2c) /= 4) all_tets(iv) = .FALSE.

                    IF ( i == 1) size_first = SIZE(iv2c)
                    IF ( SIZE(iv2c) /= size_first) mixed(iv) = .TRUE.

                ENDDO

            ENDIF

            !Check to see if the vertex is on a boundary

            IF ( on_boundary_(msh%verts(iv)) ) THEN

                ! we are on a physical boundary, so then check if it is stick or slip.
                ! if it is a slip boundary, then the vertex is mobile, but constrained.
                ! if we are used by more than one boundary, then we treat this as a stick
                ! boundary, since motion with two or more constraints would be unworkable.

                CALL get_ith_conn(ib2v,b2v,iv)

                is_slip=.FALSE. ! the default for vertices on multiple boundaries

                ! for vertex on only 1 boundary get the actual value from the boundary

                IF (SIZE(ib2v) == 1) THEN ! this vertex is only used by one bc
                    bc_id=ib2v(1)

                    !DEBUG
                    IF ((bc_id<0) .OR. (bc_id>msh%nbc)) THEN
                        WRITE(6,*)"Error:  Boundary condition index out of bounds."
                        WRITE(6,*)"Value for vertex ",iv," is equal to ",bc_id
                        WRITE(6,*)"but the mesh only has ",msh%nbc," boundaries."
                        WRITE(6,*)"This vertex is located at ",x_(msh%verts(iv)),", ", &
                            & y_(msh%verts(iv)),", ",z_(msh%verts(iv))
                    ENDIF
                    !end DEBUG

                    IF (vertex_motion_(bc(bc_id)) == sliding_) is_slip=.TRUE.

                ENDIF ! end of if-check for a vertex shared by multiple boundaries

                IF (is_slip) THEN ! this vtx, connected to tets, can slide on the boundary

                    n_constrained=n_constrained+1
                    constrained(n_constrained)=iv

                ENDIF

            ELSE ! this is an interior vertex that can move
                IF (.NOT. (mixed(iv)) ) THEN ! unless shared by mixed cells
                    n_unconstrained=n_unconstrained+1  ! increment count
                    unconstrained(n_unconstrained)=iv  !add to the unconstrained list
                ENDIF
            ENDIF

        ENDDO ! end of loop over all vertices

        ! clear memory
        NULLIFY(ic2v); NULLIFY(iv2c) ; NULLIFY(ib2v)
        CALL free_conn(f2v) ; CALL free_conn(b2v)

        END PROCEDURE mobile_verts

END SUBMODULE mobile_verts_implementation
