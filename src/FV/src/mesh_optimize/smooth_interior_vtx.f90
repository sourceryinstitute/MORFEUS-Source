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
! $Id: smooth_interior_vtx.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description: moves an interior point in the mesh to optimize adjoining cell quality
!
SUBMODULE (tools_mesh_optimize) smooth_interior_vtx_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE smooth_interior_vtx
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_keytable
        USE class_mesh
        USE class_vector
        USE class_vertex
        USE tools_mesh_basics
        USE tools_mesh_check
        USE tools_mesh_optimize, ONLY: call_smooth,optimize_vertex_rand,right_handed
        IMPLICIT NONE

        ! Local Variables
        INTEGER,POINTER :: iv2c(:) ! given a cell, finds connected vertices
        INTEGER,POINTER :: iv2v(:) ! given a vert, finds connected vertices
        INTEGER,POINTER :: ic2v(:) ! given a vert, finds connected cells
        REAL(psb_dpk_) :: free_pos(3)   ! position of free vertex
        INTEGER :: num_incident_vtx,num_incident_tet
        INTEGER,PARAMETER :: max_incident_vtx=50   ! sanity check on connectivity
        INTEGER,PARAMETER :: max_incident_tet=2*max_incident_vtx-4 !for interior tets, this should hold
        INTEGER :: ic,j,jv,k,itet             ! vertex number, cell id, and loop index
        INTEGER :: itmp,ifree_vtx ! used for creating a right-handed list of cell vertices
        INTEGER :: info
        REAL(psb_dpk_),ALLOCATABLE :: tet_pos(:,:)
        INTEGER,ALLOCATABLE :: tet_verts(:,:)
        INTEGER :: index_copy(4)
        INTEGER :: iv1, iv2, iv3, iv4  !numbering of four tet cells
        INTEGER,ALLOCATABLE :: relative_numbering(:)
        INTEGER,ALLOCATABLE :: idloc(:),idglob(:)
        TYPE(vector) :: new_pos  ! the new vertex position, after smoothing
        INTEGER :: local_tangled, valid_flag
        REAL(psb_dpk_) :: vtx1(3),vtx2(3),vtx3(3),vtx4(3)

        IF ( .NOT. (all_tets(iv) ) ) THEN

            ! simply average with neighbor vertex positions

            CALL get_ith_conn(iv2v, msh%v2v,iv) ! get neighbors

            new_pos = vector_(0.0d0,0.0d0,0.0d0)

            DO j = 1,SIZE(iv2v)
                jv = iv2v(j)
                new_pos = new_pos + position_(msh%verts(jv))
            ENDDO

            msh%verts(iv) = (1.0d0/SIZE(iv2v)) * new_pos

            RETURN

        ENDIF

        !else do optimization

        ! allocate storage for point-wise connectivity (should probably use regular array here)
        ALLOCATE(relative_numbering(SIZE(msh%verts)),stat=info)

        ALLOCATE(tet_pos(3,max_incident_vtx),tet_verts(3,max_incident_tet),stat=info)
        IF(info /= 0) THEN
            WRITE(*,600)
            CALL abort_psblas
        END IF

        ! prepare info for calling OptMS
        CALL get_ith_conn(iv2v,msh%v2v,iv)  !list of verts connected to ith vert

        num_incident_vtx = SIZE(iv2v)

        IF (num_incident_vtx > max_incident_vtx) THEN
            WRITE(*,400) ! too many tets connected to this vertex--something must be wrong
            CALL abort_psblas
        ENDIF

        IF ( shared_flag(iv) ) THEN ! we have to use supplemental connectivity

            jv = loc_to_glob_(msh%desc_v,iv) ! j is now global index number

            CALL get_kt_row(msh%c2ov_sup, jv, ic2v) ! get global cell numbers

        ELSE

            CALL get_ith_conn(ic2v, c2v, iv)  !list of tets connected to ith vert

        ENDIF

        num_incident_tet = SIZE(ic2v)

        IF (num_incident_tet > max_incident_tet) THEN
            WRITE(*,300)  ! too many tets use this vertex--something must be wrong
            CALL abort_psblas
        ENDIF

        ! note location of the free vertex, whoose position will be optimized
        free_pos(1) = x_(msh%verts(iv))
        free_pos(2) = y_(msh%verts(iv))
        free_pos(3) = z_(msh%verts(iv))

        ! load up an array with the positions of the other vertices
        DO j=1,num_incident_vtx

            jv=iv2v(j)  ! look up vertex id number from ith conn
            relative_numbering(jv)=j ! store array position so that we can go backward

            tet_pos(1,j) = x_(msh%verts(jv))
            tet_pos(2,j) = y_(msh%verts(jv))
            tet_pos(3,j) = z_(msh%verts(jv))
        ENDDO

        ! set up connectivity array and check if any cells are tangled
        local_tangled = 0

        DO itet = 1,num_incident_tet

            IF ( shared_flag(iv) ) THEN ! use global id's set above
                ic = ic2v(itet)

                !lookup vertex numbers, global cell and vertex id's
                CALL get_kt_row (msh%ov2c_sup, ic, iv2c)

                ! make a copy to dereference iv2c
                index_copy(1:4) = iv2c(1:4)

                ! convert to local vertex id's
                CALL psb_glob_to_loc(index_copy,msh%desc_v,info)

            ELSE
                ic = ic2v(itet)  ! look up cell id number from ith conn

                !get the vertices for this cell
                CALL get_ith_conn(iv2c,msh%v2c,ic)

                ! make a copy to dereference iv2c
                index_copy(1:4) = iv2c(1:4)
            ENDIF

            ! check to see if the cell is inverted
            iv1 = index_copy(1)
            iv2 = index_copy(2)
            iv3 = index_copy(3)
            iv4 = index_copy(4)

            vtx1 = position_( msh%verts(iv1) )
            vtx2 = position_( msh%verts(iv2) )
            vtx3 = position_( msh%verts(iv3) )
            vtx4 = position_( msh%verts(iv4) )

            ! if any cells are invalid, the local mesh is tangled
            valid_flag =  right_handed(vtx1, vtx2, vtx3, vtx4)
            IF ( valid_flag /= 1 ) local_tangled = 1

            ! begin process of shifting vertices so that the free vertex is first in the array
            ! this is required by OptMS
            DO k=1,4  !find where the free vertex is located...we want it listed first
                IF (index_copy(k) == iv) THEN

                    ifree_vtx = k
                ENDIF
            ENDDO

            IF (ifree_vtx /= 1 ) THEN ! if it is 1, we are fine
                itmp = index_copy(ifree_vtx)  ! else, swap with first of the 4 vertices
                index_copy(ifree_vtx) = index_copy(1)
                index_copy(1) = itmp

                ! the listing of vertices is now left-handed.  Switch elements 2 & 3
                ! to make this a right-handed listing
                itmp = index_copy(3)
                index_copy(3) = index_copy(2)
                index_copy(2) = itmp
            ENDIF

            !store the connected vertices for each tet (except the vertex of interest)
            !in the 2d array

            DO k = 2,4 ! loop over verts in this tet skipping the free vertex

                !get relative numbering (in the tet_pos array above)
                j = relative_numbering(index_copy(k))  ! index copy is the vertex number

                tet_verts(k-1,itet) = j-1 !shift numbering by 1 to start with 0
            ENDDO ! end of k loop

        ENDDO ! end of loop over connected tets

        !         --------- call smoother ------
        IF (num_incident_tet == num_incident_vtx*2-4) THEN

            info = call_smooth (num_incident_vtx,num_incident_tet, free_pos,  &
                & tet_pos, tet_verts, local_tangled)
            !returns the optimal vertex location, free_pos

            new_pos = vector_(free_pos(1),free_pos(2),free_pos(3))
            msh%verts(iv) = new_pos

        ELSE ! something is wrong, but not fatal, so produce helpful information

            WRITE(6,'(a,i2,a,i2,a,i2)') &
                & "Abnormal vertex with: ",num_incident_vtx," incident vertices and ", &
                & num_incident_tet," incident tets on processor ",mypnum_()
            ALLOCATE(idloc(1),idglob(1))
            idloc(1) = iv
            CALL psb_loc_to_glob(idloc,idglob,msh%desc_v,info)
            WRITE(6,'(a,i5,a,i5)')"Local vertex number: ",idloc(1), &
                &                " and global number ",idglob(1)

            DEALLOCATE(idloc,idglob)

            WRITE(6,*)
            WRITE(6,*)" Index   Vertex"
            DO k=1,num_incident_vtx
                WRITE(6,'(i5,3x,i5)')k,iv2v(k)
            ENDDO

            WRITE(6,*)
            WRITE(6,*)" Index    Tet Id   "
            DO k=1,num_incident_tet
                WRITE(6,'(i5,3x,i5,3x,i5)')k,ic2v(k)
            ENDDO

        ENDIF

        IF(info /= 0) THEN  !OptMS failed
            WRITE(*,500)
            WRITE(6,*)
            WRITE(6,*)"The problem vertex is:"
            WRITE(6,*)"ID #",iv
            WRITE(6,*)"Located at: (",x_(msh%verts(iv)),",",y_(msh%verts(iv)),",", &
                & z_(msh%verts(iv)),")"
            WRITE(6,*)"With ",num_incident_vtx," incident vertices and ",num_incident_tet, &
                &" incident tets."
            WRITE(6,*)
            CALL abort_psblas
        END IF


!!$  j = 0
!!$  random_iteration: do
!!$     j = j + 1
!!$
!!$     ! set up connectivity array and check if any cells are tangled
!!$     local_tangled = 0
!!$
!!$     update_volume: do itet = 1,num_incident_tet
!!$
!!$        if ( shared_flag(iv) ) then ! use global id's set above
!!$           ic = ic2v(itet)
!!$
!!$           !lookup vertex numbers, global cell and vertex id's
!!$           call get_kt_row (msh%ov2c_sup, ic, iv2c)
!!$
!!$           ! make a copy to dereference iv2c
!!$           index_copy(1:4) = iv2c(1:4)
!!$
    ! convert to local vertex id's
!!$           call psb_glob_to_loc(index_copy,msh%desc_v,info)
!!$
!!$        else
!!$           ic = ic2v(itet)  ! look up cell id number from ith conn
!!$
!!$           !get the vertices for this cell
!!$           call get_ith_conn(iv2c,msh%v2c,ic)
!!$
!!$           ! make a copy to dereference iv2c
!!$           index_copy(1:4) = iv2c(1:4)
!!$        endif
!!$
!!$        ! check to see if the cell is inverted
!!$        iv1 = index_copy(1)
!!$        iv2 = index_copy(2)
!!$        iv3 = index_copy(3)
!!$        iv4 = index_copy(4)
!!$
!!$        valid_flag =  right_handed(vtx1, vtx2, vtx3, vtx4)
!!$
!!$        if ( valid_flag /= 1 ) local_tangled = 1
!!$     end do update_volume
!!$
!!$     if ( ( local_tangled == 0 ) .or. (j >=50 ) ) exit
!!$
!!$     if ( local_tangled == 1 )  call optimize_vertex_rand(msh,c2v,iv)
!!$

!!$  end do random_iteration
!!$
!!$  deallocate(relative_numbering)
!!$
!!$  deallocate( tet_pos,tet_verts)
!!$
!!$  if ( allocated(idloc) )  deallocate(idloc)
!!$
!!$  if ( allocated(idglob) ) deallocate(idglob)
!!$

300     FORMAT(' ERROR! Too many tets connected to a single vertex in SMOOTH_INTERIOR_VTX')
400     FORMAT(' ERROR! Too many verts connected to a single vertex in SMOOTH_INTERIOR_VTX')
500     FORMAT(' ERROR! Call to smoother from SMOOTH_INTERIOR_VTX returned error.')
600     FORMAT(' ERROR! Memory allocation failure in SMOOTH_INTERIOR_VTX')

        END PROCEDURE smooth_interior_vtx

END SUBMODULE smooth_interior_vtx_implementation