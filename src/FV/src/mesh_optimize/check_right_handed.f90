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
! $Id: check_right_handed.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description: sweeps through the mesh to see if there are tangled cells
!
SUBMODULE (tools_mesh_optimize) check_right_handed_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE check_right_handed
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_keytable
        USE class_mesh
        USE class_vector
        USE class_vertex
        USE tools_mesh_basics
        USE tools_mesh_check
        USE tools_mesh_optimize, ONLY: initoptms,freeoptms,right_handed, &
            &                          smooth_interior_vtx
        IMPLICIT NONE

        ! Local variables

        INTEGER :: iv,ic,i,j          ! vertex number, cell id, and loop index
        INTEGER :: nverts           !number of local vertices in mesh
        INTEGER :: info
        REAL(psb_dpk_) :: pos(3,4)  !x,y,z position for a tet
        INTEGER,PARAMETER :: max_incident_vtx=50   ! sanity check on connectivity
        ! for interior tets, this should hold
        INTEGER,PARAMETER :: max_incident_tet=2*max_incident_vtx-4

        ! parameters for calling Opt-MS library
        INTEGER,POINTER :: iv2c(:) ! given a cell, finds connected vertices
        INTEGER,POINTER :: ic2v(:) ! given a vert, finds connected cells
        INTEGER :: index_copy(4)

        !check if tet cells are tangled (assumes connectivity is right-handed)
        tangled=0

        nverts = psb_cd_get_local_rows(msh%desc_v) !local and shared vertices

        DO i=1,msh%v2c%nel_()  ! check local cells

            CALL msh%v2c%get_ith_conn(iv2c,i)  !list of verts connected to ith cell

            IF ( msh%cells(i)%geo_() == 'tet' ) THEN ! is a tet

                ! get positions of the cell's vertices
                DO iv=1,SIZE(iv2c)
                    pos(1,iv) = msh%verts(iv2c(iv))%x_()
                    pos(2,iv) = msh%verts(iv2c(iv))%y_()
                    pos(3,iv) = msh%verts(iv2c(iv))%z_()
                ENDDO

                !if info is 0 then the cell is degenerate, -1 is left-handed
                info = right_handed(pos(1,1),pos(1,2),pos(1,3),pos(1,4))

                IF (info <= 0) tangled = tangled + 1 ! there is at least 1 inverted cell
            ENDIF
        ENDDO

        IF (ALLOCATED(shared)) THEN
            ! check if cells attached to overlap vertices are tangled
            DO i  = 1, nverts
                IF ( shared_flag(i) ) THEN

                    iv = loc_to_glob_(msh%desc_v,i)

                    ! get list of global id of cells attached to a vertex
                    CALL msh%c2ov_sup%get_kt_row(iv, ic2v)

                    ! loop over cells connected to vertex..even if we've checked
                    !  them previously

                    DO j = 1,SIZE(ic2v)

                        !loop over cells
                        ic=ic2v(j)

                        !get vertices of the cell (with global numbering)
                        CALL msh%ov2c_sup%get_kt_row(ic, iv2c)

                        IF ( all_tets(i) ) THEN ! then we can use OptMS code for tet

                            index_copy = iv2c

                            ! get local vertex numbers
                            CALL psb_glob_to_loc(index_copy,msh%desc_v,info)

                            IF(info /= 0) THEN
                                WRITE(*,200)
                                CALL abort_psblas
                            END IF

                            DO iv=1,SIZE(iv2c) ! get position of vertices
                                pos(1,iv) = msh%verts(index_copy(iv))%x_()
                                pos(2,iv) = msh%verts(index_copy(iv))%y_()
                                pos(3,iv) = msh%verts(index_copy(iv))%z_()
                            ENDDO

                            !if info is 0 then the cell is degenerate, -1 is left-handed
                            info = right_handed(pos(1,1),pos(1,2),pos(1,3),pos(1,4))

                            IF (info < 1) tangled = tangled + 1 ! there is at least 1 inverted cell
                        ENDIF ! end of "4 vertices" check

                    ENDDO ! end of loop over cells connected to shared vertex

                ENDIF ! end of shared check

            ENDDO ! end of loop over all vertices
        ENDIF

200     FORMAT(' ERROR! Failed to get overlap points in CHECK_RIGHT_HANDED')

        END PROCEDURE check_right_handed

END SUBMODULE check_right_handed_implementation
