
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
! $Id: smooth_surf_vtx.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description: moves an surf point in the mesh to optimize adjoining cell quality
!
SUBMODULE (tools_mesh_optimize) smooth_surf_vtx_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE smooth_surf_vtx
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_face
        USE class_keytable
        USE class_mesh
        USE class_surface
        USE class_vector
        USE class_vertex
        USE tools_mesh_basics
        USE tools_mesh_check
        USE tools_mesh_optimize, ONLY: call_smooth2d,right_handed2d
        USE iso_c_binding, ONLY : c_loc
        IMPLICIT NONE

        ! Local Variables
        INTEGER,POINTER :: iv2f(:)=>null() ! given a cell, finds connected vertices
        INTEGER,POINTER :: if2v(:)=>null() ! given a vert, finds connected cells

        TYPE(vector)    ::   free_pos      ! position of free vertex
        REAL(psb_dpk_) :: free_p2d(2)   ! position of free vertex in 2d local coord. sys
        INTEGER, TARGET :: num_incident_vtx, num_incident_tri
        INTEGER,PARAMETER :: max_incident_vtx=150   ! sanity check on connectivity
        INTEGER,PARAMETER :: max_incident_tri=max_incident_vtx !for tris, this should hold
        INTEGER :: iface,j,jv,k,itri,num_f       ! vertex number, cell id, and loop index
        INTEGER :: ifree_vtx       ! used for creating a right-handed list of cell vertices
        INTEGER :: info
        TYPE(vector)     ,ALLOCATABLE :: vtx_pos(:)   ! 3d position
        REAL(psb_dpk_),ALLOCATABLE :: vtx_p2d(:,:) ! 2d position in local coord sys

        INTEGER,ALLOCATABLE :: tri_verts(:,:)
        INTEGER,ALLOCATABLE, TARGET      :: index_copy(:)
        INTEGER,ALLOCATABLE :: relative_numbering(:)
        INTEGER,ALLOCATABLE :: idloc(:),idglob(:)
        TYPE(vector) :: new_pos  ! the new vertex position, after smoothing
        TYPE(vector) :: ahat,bhat! unit vectors in the local coordinate system
        TYPE(vector) :: nhat     ! normal unit vector to the surface
        TYPE(vector) :: prel     ! relative position between neighbor point and free point
        REAL(psb_dpk_) :: xlocal,ylocal  ! posititions of points in the local coord. system
        INTEGER :: ibv2v(max_incident_vtx)  ! connected vertices that lie on boundaries
        INTEGER :: ibf2v(max_incident_tri)  ! connected triangles that lie on boundaries
        REAL:: cfac                         ! correction factor for non-planar surface dist.
        INTEGER, TARGET :: local_tangled ! for checking if the mesh is locally twisted
        INTEGER :: valid_flag  ! for checking if the mesh is locally twisted
        REAL(psb_dpk_) :: vtx1(2),vtx2(2),vtx3(2)
        TYPE(vector) :: vtx_pos_temp

        ! allocate storage for point-wise connectivity
        ALLOCATE(relative_numbering(SIZE(msh%verts)),stat=info)

        ALLOCATE(vtx_pos(  max_incident_vtx)  , tri_verts(2,max_incident_tri), &
            &   vtx_p2d(2,max_incident_vtx)  , stat=info)

        IF(info /= 0) THEN
            WRITE(*,600)
            CALL abort_psblas
        END IF

        ! prepare info for calling OptMS

        ! first, construct list of boundary faces that are attached to the vertex
        IF ( shared_flag(iv) ) THEN ! we have to use supplemental connectivity

            jv = loc_to_glob_(msh%desc_v,iv) ! jv is now global vertex index number

            CALL msh%f2ov_sup%get_kt_row(jv, if2v) ! get global face numbers

            ! assume that we have culled faces that do not lie on boundaries
            ! when we created the supplemental info
            num_incident_tri = SIZE(if2v)

            ibf2v(1:num_incident_tri) = if2v

        ELSE

            CALL f2v%get_ith_conn(if2v, iv)  !list of tris connected to ith vert

            ! cull faces that do not lie on boundaries

            num_incident_tri = 0
            DO j = 1, SIZE(if2v)
                iface = if2v(j)

                IF ( msh%faces(iface)%flag_() > 0 ) THEN
                    num_incident_tri = num_incident_tri + 1
                    ibf2v(num_incident_tri) = iface
                ENDIF

            ENDDO

        ENDIF

        IF (num_incident_tri > max_incident_tri) THEN
            WRITE(*,300)  ! too many tris use this vertex--something must be wrong
            CALL abort_psblas
        ENDIF

        ! construct vertex list from triangle's vertices, avoiding repeats
        num_incident_vtx = 0  ! this will eventually tally to equal num tri
        ibv2v(:) = 0          ! this will contain the vertices that belong

        DO k=1,num_incident_tri
            itri=ibf2v(k)      ! itri is the triangle ID extracted from tri list

            IF ( shared_flag(iv) ) THEN ! we have global face ID's

                !lookup vertex numbers using global face and vertex id's
                CALL msh%ov2f_sup%get_kt_row (itri, iv2f)

                ! make a copy to dereference the iv2f pointer
                num_f = SIZE(iv2f)

                ALLOCATE(index_copy(num_f))

                index_copy(1:num_f) = iv2f(1:num_f)

                ! convert global to local vertex id's
                CALL psb_glob_to_loc(index_copy,msh%desc_v,info)

                ! save local vert ids in the same storage as we would use in the else-branch
                iv2f => index_copy(1:num_f)

            ELSE
                ! get list of vertices that belong to this triangle
                CALL msh%v2f%get_ith_conn(iv2f,itri)
            ENDIF

            DO j = 1, SIZE(iv2f)
                jv = iv2f(j)  ! one of the vertices used by the tri
                ! we don't want to count our free vertex or already-counted vertiex
                IF ( ( jv /= iv ) .AND. (.NOT. ( ANY( ibv2v == jv) ) ) ) THEN
                    num_incident_vtx = num_incident_vtx + 1
                    ibv2v( num_incident_vtx ) = jv
                ENDIF
            ENDDO

            NULLIFY(iv2f)  ! end the association with index_copy
            IF ( ALLOCATED(index_copy) ) DEALLOCATE(index_copy)

        ENDDO

        IF (num_incident_vtx > max_incident_vtx) THEN
            WRITE(*,400) ! too many tris connected to this vertex--something must be wrong
            CALL abort_psblas
        ENDIF

        ! the position to be optimized
        free_pos = msh%verts(iv)%position_()

        ! load up an array with the positions of the non-free surface vertices
        DO j=1,num_incident_vtx

            jv=ibv2v(j)  ! look up vertex id number from ith conn
            relative_numbering(jv)=j ! store array position so that we can go backward

            vtx_pos(j) = msh%verts(jv)%position_()
        ENDDO

        ! establish 2d basis vectors so that we can operate in a local x,y coord sys

        free_p2d(:) = 0.0d0 ! the free vertex is at the origin

        ! normal of conceptual surface
        nhat = msh%surf(ib)%normal_(free_pos)

        ! the conceptual surface normal should be constructed so that the normal points outside
        !  of the domain.  The local processor may not have ownership of all the faces, if this
        ! is an overlap vertex, so find the first local face.  Surface smoothing requires
        ! at least one local face
        IF ( shared_flag(iv) ) THEN
            k = 0
            DO
                k = k + 1
                itri=ibf2v(k)      ! itri is the triangle ID extracted from tri list
                iface = glob_to_loc_(msh%desc_f,itri) ! get the local ID

                IF ( ( iface >= 0 ) .OR. ( k == num_incident_tri) ) EXIT
            ENDDO

            IF ( iface <= 0 ) THEN  ! none of these faces are local, so no smoothing possible
                ! exit and let other processors negotiate the vertex location

                DEALLOCATE(relative_numbering, vtx_pos, vtx_p2d, tri_verts)

                IF ( ALLOCATED(idloc) )  DEALLOCATE(idloc)

                IF ( ALLOCATED(idglob) ) DEALLOCATE(idglob)

                RETURN

            ELSE
                ! found local face ID
                IF ( ( nhat .dot. msh%af( iface ) ) < 0.0d0 ) nhat = (-1.0d0) * nhat
            ENDIF
        ELSE
            ! for local vertices, any face normal will do
            IF ( ( nhat .dot. msh%af( ibf2v(1) ) ) < 0.0d0) nhat = (-1.0d0) * nhat
        ENDIF

        ! the x_local unit vector is arbitrarily constructed from the first vtx position
        jv = ibv2v(1)
        vtx_pos_temp = vtx_pos(1) - free_pos
        !ahat = unit( vtx_pos(1) - free_pos )
        ahat = vtx_pos_temp %unit()
        !because the points do not generally lie in a plane, ahat is not yet orthogonal to nhat
        !so subtract off any component of nhat

        ahat = ahat - (nhat .dot. ahat) * ahat
        ahat = ahat%unit()

        ! bhat is constructed so that ahat cross bhat gives nhat
        bhat = ( -1.0d0 ) * ( ahat .cross. nhat )

        ! calculate xlocal and ylocal for each point and store in 2d coord sys
        DO j=1,num_incident_vtx

            jv=ibv2v(j)  ! look up vertex id number from ith conn

            prel = msh%verts(jv)-msh%verts(iv) !position relative to free vtx

            !calculate projection into ahat/bhat plane
            xlocal = ahat .dot. prel
            ylocal = bhat .dot. prel

            !correct projected lengths to restore original prel
            cfac = prel%mag()/SQRT(xlocal**2 + ylocal **2)

            vtx_p2d(1,j) = xlocal*cfac
            vtx_p2d(2,j) = ylocal*cfac

        ENDDO

        IF ( all_tets(iv) ) THEN ! use optimization on triangular faces

            ALLOCATE(index_copy(3))
            local_tangled = 0

            ! set up connectivity array
            DO itri = 1,num_incident_tri

                IF ( shared_flag(iv) ) THEN ! use global id's set above

                    iface = ibf2v(itri)

                    !lookup vertex numbers using global face and vertex id's
                    CALL msh%ov2f_sup%get_kt_row (iface, iv2f)

                    ! make a copy to dereference iv2f
                    index_copy(1:3) = iv2f(1:3)

                    ! convert to local vertex id's
                    CALL psb_glob_to_loc(index_copy,msh%desc_v,info)

                    !DEBUG--to be deleted once mesh adaptation coding is complete
                    IF (ANY(index_copy < 1)) THEN
                        WRITE(6,*) "Processor ",mypnum_(),"....vertices out of range on face",iface
                        DO k =1,3
                            WRITE(6,'(a,i3,a,i1,a,i5,a,i6)') &
                                & "Processor ",mypnum_()," index_copy(",k,") =",index_copy(k), &
                                & "   global vertex ID =",iv2f(k)
                        ENDDO

                        CALL abort_psblas
                    ENDIF

                    ! end DEBUG

                ELSE
                    iface = ibf2v(itri)  ! look up cell id number from ith conn

                    !get the vertices for this face and shift the list until
                    !the free vertex is first.  Since this is a bndry face, all
                    ! vertices should be boundary vertices
                    CALL msh%v2f%get_ith_conn(iv2f,iface)

                    ! make a copy to dereference iv2f
                    index_copy(1:3) = iv2f(1:3)
                ENDIF

                DO k=1,3  !find where the free vertex is located...we want it listed first
                    IF (index_copy(k) == iv) THEN

                        ifree_vtx = k
                    ENDIF
                ENDDO

                ! circular shift to make the free vertex first yet keep right-handed
                ! orientation
                IF (ifree_vtx == 2 ) index_copy = CSHIFT(index_copy, 1)
                IF (ifree_vtx == 3)  index_copy = CSHIFT(index_copy,-1)

                !store the connected vertices for each tri (except the vertex of interest)
                !in the 2d array

                DO k = 2,3 ! loop over verts in this tri skipping the free vertex

                    !get relative numbering (in the vtx_pos array above)
                    j = relative_numbering(index_copy(k))  ! index copy is the vertex number

                    tri_verts(k-1,itri) = j-1 !shift numbering by 1 to start with 0
                ENDDO ! end of k loop

                ! check that the triangle is not twisted or invalid
                vtx1 = (/ 0.0d0, 0.0d0 /)
                j= relative_numbering(index_copy(2))
                vtx2 = vtx_p2d(1:2,j)
                j= relative_numbering(index_copy(3))
                vtx3 = vtx_p2d(1:2,j)

                ! we use 3 element vectors and the 3rd element is ignored in right_handed2d
                valid_flag =  right_handed2d (vtx1, vtx2, vtx3)
                IF ( valid_flag == -1 ) local_tangled = 1

            ENDDO ! end of loop over connected tris

            !        ------------ call smoother -----------------
            IF (num_incident_tri == num_incident_vtx) THEN

                info = call_smooth2d (c_loc(num_incident_vtx),c_loc(num_incident_tri), free_p2d,  &
                    & vtx_p2d, tri_verts, c_loc(local_tangled))

                !returns the optimal vertex location, free_p2d, in the local coord. system

                IF(info /= 0) THEN
                    WRITE(*,500)
                    CALL abort_psblas
                END IF

            ELSE ! something is wrong, but not fatal, so produce helpful information

                WRITE(6,*)
                WRITE(6,*)
                WRITE(6,'(a,i2,a,i2,a,i2)') &
                    & "Abnormal vertex with: ",num_incident_vtx," incident vertices and ", &
                    & num_incident_tri," incident tris on processor ",mypnum_()

                IF ( msh%verts(iv)%on_boundary_() ) THEN
                    WRITE(6,'(a,i4,a,3(f8.3,1x))')"Vertex on boundary",ib," at position:", &
                        & free_pos%x_(),free_pos%y_(),free_pos%z_()
                ELSE
                    WRITE(6,'(a)')"Vertex not on boundary"
                ENDIF

                ALLOCATE(idloc(1),idglob(1))
                idloc(1) = iv
                CALL psb_loc_to_glob(idloc,idglob,msh%desc_v,info)
                WRITE(6,'(a,i5,a,i5)')"Local vertex number: ",idloc(1), &
                    &                " and global number ",idglob(1)

                DEALLOCATE(idloc,idglob)

                WRITE(6,*)
                WRITE(6,*)" Index   Vertex   On Boundary"
                DO k=1,num_incident_vtx
                    jv = ibv2v(k)
                    WRITE(6,'(i5,3x,i5,9x,l1)')k,jv,msh%verts(jv)%on_boundary_()
                ENDDO

                WRITE(6,*)
                WRITE(6,*)" Index   Tri ID  "

                DO k=1,num_incident_tri
                    itri=ibf2v(k)
                    WRITE(6,'(2(i5,5x))') k,itri
                ENDDO

            ENDIF

            DEALLOCATE(index_copy)

        ELSE  ! if not all tets, find free_pos by simple averaging

            free_p2d(1) = SUM(vtx_p2d(1,1:num_incident_vtx))/num_incident_vtx
            free_p2d(2) = SUM(vtx_p2d(2,1:num_incident_vtx))/num_incident_vtx

        ENDIF ! end of non-tet branch

        ! convert to 3d coord system.  Note that this conversion is not perfect
        ! since we are working in a projected 2d coord. system, but as the movement
        ! of free_pos becomes small, it becomes a very good approximation (2nd order)
        new_pos = free_pos + free_p2d(1) * ahat + free_p2d(2) * bhat

        ! for curved surfaces, smoothing introduces small errors that require the
        ! point to be forced back onto the precise surface
        new_pos = msh%surf(ib)%get_closest_point(new_pos)

        msh%verts(iv) = new_pos

        IF(info /= 0) THEN  !OptMS failed
            WRITE(*,500)
            WRITE(6,*)
            WRITE(6,*)"The problem vertex is:"
            WRITE(6,*)"ID #",iv
            WRITE(6,*)"Located at: (",msh%verts(iv)%x_(),",",msh%verts(iv)%y_(),",", &
                & msh%verts(iv)%z_(),")"
            WRITE(6,*)"With ",num_incident_vtx," incident vertices and ",num_incident_tri, &
                &" incident tris."
            WRITE(6,*)
            CALL abort_psblas
        END IF

        DEALLOCATE(relative_numbering)

        DEALLOCATE( vtx_pos,vtx_p2d,tri_verts)

        IF ( ALLOCATED(idloc) )  DEALLOCATE(idloc)

        IF ( ALLOCATED(idglob) ) DEALLOCATE(idglob)

300     FORMAT(' ERROR! Too many triangles connected to a single vertex in SMOOTH_SURF_VTX')
400     FORMAT(' ERROR! Too many verts connected to a single vertex in SMOOTH_SURF_VTX')
500     FORMAT(' ERROR! Call to smoother from SMOOTH_SURF_VTX returned error.')
600     FORMAT(' ERROR! Memory allocation failure in SMOOTH_SURF_VTX')

        END PROCEDURE smooth_surf_vtx

END SUBMODULE smooth_surf_vtx_implementation
