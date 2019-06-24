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
SUBMODULE (tools_mesh) cmp_mesh_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE cmp_moving_surf
        USE class_psblas!, ONLY : mypnum_
        USE class_connectivity
        USE class_surface!, ONLY : alloc_surface
        USE class_vertex
        IMPLICIT NONE
        !! $Id: cmp_moving_surf.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description: Allocates and calculates shapes of moving surfaces
        !!
        INTEGER :: ib, info

        IF(mypnum_() == 0) THEN
            WRITE(*,*) 'Autodetecting virtual surfaces:'

            ! Not necessary since INTENT(OUT)
    !!$  if (associated (surf) ) then
    !!$     call abort_psblas
    !!$  endif

            ! Allocates pointer of SURFACE objects
            ALLOCATE(surf(nbc),stat=info)
            IF (info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            ENDIF

            ! Identifies surfaces
            DO ib = 1, nbc
                CALL alloc_surface(v2b,ib,verts,surf(ib))
            END DO

            WRITE(*,*)

        ENDIF

100     FORMAT(' ERROR! Failure to allocate memory in CMP_MOVING_SURF')

        END PROCEDURE cmp_moving_surf


        MODULE PROCEDURE cmp_mesh_c2c
        USE class_psblas
        USE class_connectivity
        USE class_face
        IMPLICIT NONE
        !! $Id: cmp_mesh_c2c.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:
        !!
        INTEGER :: ic, IF, im, is, info
        INTEGER :: n, ncells, nfl_faces
        INTEGER, ALLOCATABLE :: itab(:,:), kconn(:)


        IF(mypnum_() == 0) THEN
            WRITE(*,*) 'Computing cells adjacency graph C2C'

            n         = max_conn(f2c)            ! Maximum connectivity degree
            ncells    = nel_(f2c)                ! Number of cells
            nfl_faces = COUNT(flag_(faces) == 0) ! Number of fluid faces

            ALLOCATE(itab(n,ncells),kconn(ncells),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            itab = 0
            kconn = 0

            ! WARNING! CMP_MESH_C2C is called before the global-to-local reallocation.
            ! Thus faces are still numbered in such a way that fluid ones come first

            DO IF = 1, nfl_faces
                im = master_(faces(IF))
                is = slave_(faces(IF))
                kconn(im) = kconn(im) + 1
                kconn(is) = kconn(is) + 1
                itab(kconn(im),im) = is
                itab(kconn(is),is) = im
            END DO

            ! Counts total amount of connectivities
            n = 0
            DO ic = 1, ncells
                n = n + kconn(ic)
            END DO

            ! Creates C2C connectivity object
            CALL alloc_conn(c2c,nel=ncells,nconn=n)

            ! Sets C2C connectivity
            DO ic = 1, ncells
                n = kconn(ic)
                CALL set_ith_conn(c2c,ic,itab(1:n,ic))
            END DO

            DEALLOCATE(itab,kconn,stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

        ENDIF

100     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_C2C')
200     FORMAT(' ERROR! Memory deallocation failure in CMP_MESH_C2C')

        END PROCEDURE cmp_mesh_c2c


        MODULE PROCEDURE cmp_mesh_f2b
        USE class_psblas
        USE class_connectivity
        USE class_face
        IMPLICIT NONE
        !! $Id: cmp_mesh_f2b.f90 3175 2008-06-13 12:59:07Z sfilippo $
        !!
        !! Description:
        !!
        INTEGER :: ib, IF, info, maxconn, n, nfaces,j
        INTEGER :: flag(SIZE(faces)), iface(SIZE(faces)), kf(-1:nbc)
        INTEGER, ALLOCATABLE :: iconn(:)


        IF(mypnum_() == 0) &
            & WRITE(*,*) 'Computing face to boundary connectivity F2B'

        nfaces = SIZE(faces)
        iface = (/(IF, IF = 1, nfaces)/)

        flag = flag_(faces)

        ! Computes maximum connectivity degree
        IF (.FALSE.) THEN
            maxconn = 0
            DO ib = -1, nbc
                kf(ib) = COUNT(flag == ib)
                maxconn = MAX(maxconn,kf(ib))
            END DO
        ELSE
            CALL psb_msort(flag(1:nfaces),ix=iface(1:nfaces),flag=psb_sort_keep_idx_)
            IF (nfaces > 0) THEN
                IF ((flag(1) < -1).OR.(flag(nfaces)>nbc)) THEN
                    WRITE(0,*) 'Error in cmp_mesh_f2b: flags out of bounds ',flag(1),flag(nfaces)
                    CALL abort_psblas
                END IF
            END IF
            maxconn = 0
            kf(:)   = 0
            DO IF = 1, nfaces
                ib  = flag(IF)
                kf(ib) = kf(ib) + 1
                maxconn = MAX(maxconn,kf(ib))
            END DO
        END IF

        CALL alloc_conn(f2b,nel=nbc+2,nconn=nfaces,lb=-1)

        ALLOCATE(iconn(maxconn),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        j = 0
        DO ib = -1, nbc
            n = kf(ib)
            iconn(1:n) = iface(j+1:j+n)
    !!$     iconn(1:n) = pack(iface,flag == ib)
            CALL set_ith_conn(f2b,ib,iconn(1:n))
            j = j + n
        END DO

        DEALLOCATE(iconn)

100     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_F2B')

        END PROCEDURE cmp_mesh_f2b


        MODULE PROCEDURE cmp_mesh_f2f
        USE class_psblas
        USE class_connectivity
        IMPLICIT NONE
        !! $Id: cmp_mesh_f2f.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:
        !!    Provides the Face To Face adjacency graph.
        !!    The criterium for associating the graph nodes is arbitrarly defined by the
        !!    user. In this case, using a co-located cell-centered strategy for variable
        !!    displacement, the F2F connectivity is obtained by linking each face to all
        !!    the ones belonging to its master and slave cells.
        !!
        INTEGER :: i, ic, IF, if1, if2, info, j, n, ncells, ncf
        INTEGER, POINTER :: if2c(:) => NULL()
        INTEGER, ALLOCATABLE :: kconn(:), itab(:,:)


        IF(mypnum_() == 0) THEN
            WRITE(*,*) 'Computing faces adjacency graph F2F'

            n = (max_conn(f2c) - 1) * 2 ! F2F maximum connectivity degree
            ncells = nel_(f2c)          ! Number of cells

            ALLOCATE(itab(n,nfaces),kconn(nfaces),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            kconn = 0
            DO ic = 1, ncells
                CALL get_ith_conn(if2c,f2c,ic)
                ncf = SIZE(if2c)
                DO i = 1, ncf - 1
                    if1 = if2c(i)
                    DO j = i + 1, ncf
                        if2 = if2c(j)

                        ! Add direct arc
                        kconn(if1) = kconn(if1) + 1
                        itab(kconn(if1),if1) = if2

                        ! Add indirect arc
                        kconn(if2) = kconn(if2) + 1
                        itab(kconn(if2),if2) = if1
                    END DO
                END DO
            END DO


            ! Computes total number of connectivities
            n = 0
            DO IF = 1, nfaces
                n = n + kconn(IF)
            END DO
            CALL alloc_conn(f2f,nel=nfaces,nconn=n)


            ! Sets F2F connectivity
            DO IF = 1, nfaces
                n = kconn(IF)
                CALL set_ith_conn(f2f,IF,itab(1:n,IF))
            END DO


            ! Frees memory
            NULLIFY(if2c)
            DEALLOCATE(itab,kconn)

        ENDIF

100     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_F2F')

        END PROCEDURE cmp_mesh_f2f


        MODULE PROCEDURE cmp_mesh_part
        USE class_psblas
        USE class_connectivity
        USE tools_part
        USE part_graph, ONLY : bld_part_graph, get_part_graph
        USE part_block, ONLY : bld_part_block
        USE part_random, ONLY : bld_part_random
        IMPLICIT NONE
        !! $Id: cmp_mesh_part.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:
        !!
        INTEGER :: info, nprocs, ncells
        INTEGER, ALLOCATABLE :: xadj_glob(:)
        INTEGER, ALLOCATABLE :: adjncy_glob(:)

        IF (mypnum_() == 0) THEN
            CALL tic(sw_par)

            nprocs = nprocs_()

            IF(nprocs > 1 ) WRITE(*,*)

            ncells = nel_(c2c)

            ! Checks status of PART_CELLS vector
            IF(ALLOCATED(part_cells)) DEALLOCATE(part_cells)
            ALLOCATE(part_cells(ncells),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            ! Retrieves connectivity data in CSR-like format
            IF(ipart > 0) THEN
                ! ParMetis requires global connectivity in CSR format, therefore
                ! the call to GET_CONN_CSR is done here, outside BLD_PART_GRAPH.
                CALL get_conn_csr(c2c,xadj_glob,adjncy_glob)
            END IF

            SELECT CASE(ipart)
            CASE(-1)
                ! Random partitioning
                CALL bld_part_random(ncells,nprocs,part_cells)
            CASE(0)
                ! HPF Block partitioning
                CALL bld_part_block(ncells,nprocs,part_cells)
            CASE(1)
                ! METIS Partitioning
                CALL bld_part_graph(xadj_glob,adjncy_glob,ipart=1)
                CALL get_part_graph(part_cells)
            CASE default
                WRITE(*,200)
                CALL abort_psblas
            END SELECT
            WRITE(*,*)


            IF(ipart > 0) THEN
                DEALLOCATE(xadj_glob,adjncy_glob,stat=info)
                IF(info /= 0) THEN
                    WRITE(*,300)
                    CALL abort_psblas
                END IF
            END IF

            CALL toc(sw_par)
        ENDIF

100     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_PART')
200     FORMAT(' ERROR! Unknown partitioning method in CMP_MESH_PART')
300     FORMAT(' ERROR! Memory deallocation failure in CMP_MESH_PART')

        END PROCEDURE cmp_mesh_part


        MODULE PROCEDURE cmp_mesh_renum
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_face
        USE renum
        IMPLICIT NONE
        !! $Id: cmp_mesh_renum.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:
        !!
        INTEGER :: icontxt, mypnum


        icontxt = icontxt_()
        mypnum = mypnum_()

        CALL tic(sw_ord)

        SELECT CASE(irenum)
        CASE(0)
            IF(mypnum == 0) THEN
                WRITE(*,*) 'Matrix renumbering: none'
                WRITE(*,*)
            END IF
        CASE(1)
            ! Builds permutation array on P0, according to IRENUM value.
            IF(mypnum == 0) THEN
                CALL start_renum(irenum,c2c)
                !!     call print_renum(0)
                ! Applies permutation to cell-related structures

                CALL apply_renum(c2c,to_a_and_b_)
                CALL apply_renum(v2c,to_b_)
                CALL apply_renum(f2c,to_b_)
                CALL apply_renum(c2g,to_a_)
                CALL apply_renum(cells)
                CALL apply_renum(faces)

                ! Deallocates private storage of renum module
                CALL stop_renum
            ENDIF
        CASE default
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

        CALL toc(sw_ord)

100     FORMAT(' ERROR! Unknown renumbering method in CMP_MESH_RENUM')

        END PROCEDURE cmp_mesh_renum


        MODULE PROCEDURE cmp_mesh_v2b
        USE class_psblas
        USE class_connectivity
        USE class_face
        IMPLICIT NONE
        !! $Id: cmp_mesh_v2b.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:  Sets up vertex to boundary connectivity, so that given a boundary, the
        !! GET_ITH_CONN method can be used to obtain the list of vertices that belong to the boundary.
        !! Boundary ID's correspond to those used to flag faces:
        !!
        !!---------------------------|
        !! Vertex Flag | Vertex Type |
        !!---------------------------|
        !!      0      | interior    |
        !!    < 0      | proc. bndry |
        !!    > 0      | phys. bndry |
        !!---------------------------|
        !!
        LOGICAL :: inspected(0:nbc)
        INTEGER :: flag, IF, info, iv, k, n, nverts
        INTEGER :: kconn(0:nbc)
        INTEGER, ALLOCATABLE :: itab(:,:)
        INTEGER, POINTER :: if2v(:) => NULL()
        TYPE(connectivity) :: f2v

        IF(mypnum_() == 0) THEN
            WRITE(*,*) 'Computing vertex to boundary connectivity V2B'


            CALL get_dual_conn(v2f,f2v)

            nverts = nel_(f2v)

            ALLOCATE(itab(nverts,0:nbc),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            itab = 0
            kconn = 0

            DO iv = 1, nverts
                inspected(:) = .FALSE.

                CALL get_ith_conn(if2v,f2v,iv)
                n = SIZE(if2v)

                DO k = 1, n
                    IF = if2v(k)
                    flag = flag_(faces(IF))
                    IF(flag > 0 .AND. .NOT.inspected(flag)) THEN
                        inspected(flag) = .TRUE.
                        kconn(flag) = kconn(flag) + 1
                        itab(kconn(flag),flag) = iv
                    END IF
                END DO

                IF(ALL(inspected .EQV. .FALSE.)) THEN
                    kconn(0) = kconn(0) + 1
                    itab(kconn(0),0) = iv
                END IF
            END DO

            ! Computes number of V2B connectivity
            n = 0
            DO k = 0, nbc
                n = n + kconn(k)
            END DO

            CALL alloc_conn(v2b,nel=nbc+1,nconn=n,lb=0)

            DO k = 0, nbc
                n = kconn(k)
                CALL set_ith_conn(v2b,k,itab(1:n,k))
            END DO

            DEALLOCATE(itab)
            NULLIFY(if2v)
            CALL free_conn(f2v)

        ENDIF

100     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_V2B')

        END PROCEDURE cmp_mesh_v2b


        MODULE PROCEDURE  cmp_mesh_v2e
        USE class_psblas
        USE class_connectivity
        IMPLICIT NONE
        !! $Id: cmp_mesh_v2e.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:
        !!    Provides the Vertex To Edge connectivity.
        !!
        LOGICAL, ALLOCATABLE :: inspected(:), twice(:)
        INTEGER :: ie, ie1, ie2
        INTEGER :: iv, iv1, iv2, iv3, iv4
        INTEGER :: i, IF, info, j, n, nedges, nfaces, nverts
        INTEGER, POINTER :: iv2f(:) => NULL()
        INTEGER, POINTER :: iv2e(:) => NULL(), ie2v(:) => NULL()
        TYPE(connectivity) :: work, e2v
        TYPE(stopwatch) :: sw_v2e

        INTEGER :: d1, d2, s1, s2

        ! V2E makes sense only for 3D meshes!
        IF(ncd == 2) RETURN

        IF(mypnum_() == 0) &
            WRITE(*,*) 'Computing vertex to edge connectivity V2E'

        ! Constructs stopwatch
        sw_v2e = stopwatch_(icontxt_())

        ! Start timing
        CALL tic(sw_v2e)

        IF(mypnum_() == 0) THEN

            nfaces = nel_(v2f)

            nedges = 0
            DO IF = 1, nfaces
                CALL get_ith_conn(iv2f,v2f,IF)
                nedges = nedges + SIZE(iv2f)
            END DO


            CALL alloc_conn(work,nel=nedges,nconn=2 * nedges)

            ie = 1
            DO IF = 1, nfaces

                CALL get_ith_conn(iv2f, v2f,IF)

                n = SIZE(iv2f)
                ! # of face vertices is equal to # of face edges N

                DO i = 1, n - 1
                    CALL set_ith_conn(work,ie,(/ iv2f(i), iv2f(i+1) /))
                    ie = ie + 1
                END DO
                CALL set_ith_conn(work,ie,(/ iv2f(n), iv2f(1) /))
                ie = ie + 1
            END DO

            CALL get_dual_conn(work,e2v)

            ALLOCATE(inspected(nedges),twice(nedges),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            nverts = nel_(e2v)

            inspected = .FALSE.
            twice     = .FALSE.
            DO iv = 1,  nverts
                CALL get_ith_conn(ie2v,e2v,iv)
                n = SIZE(ie2v)

                edge1: DO i = 1, n - 1
                    ie1 = ie2v(i)
                    IF(twice(ie1) .OR. inspected(ie1)) CYCLE edge1
                    inspected(ie1) = .TRUE.

                    CALL get_ith_conn(iv2e,work,ie1)
                    iv1 = iv2e(1)
                    iv2 = iv2e(2)
                    s1 = iv1 + iv2
                    d1 = ABS(iv1 - iv2)

                    edge2: DO j = i + 1, n
                        ie2 = ie2v(j)
                        IF(twice(ie2) .OR. inspected(ie2)) CYCLE edge2

                        CALL get_ith_conn(iv2e,work,ie2)
                        iv3 = iv2e(1)
                        iv4 = iv2e(2)
                        s2 = iv3 + iv4
                        d2 = ABS(iv3 - iv4)

                        IF(s1 == s2 .AND. d1 == d2) twice(ie2) = .TRUE.
    !!$           if((iv4 == iv1 .and. iv3 == iv2) .or.&
    !!$               iv4 == iv2 .and. iv3 == iv1 ) then
    !!$              twice(ie2) = .true.
    !!$           end if
                    END DO edge2

                END DO edge1
            END DO

            nedges = COUNT(.NOT.twice)
            CALL alloc_conn(v2e,nel=nedges,nconn=nedges*2)

            ie = 0
            n = nel_(work)
            DO i = 1, n
                IF(twice(i)) CYCLE
                ie = ie + 1
                CALL get_ith_conn(iv2e,work,i)
                CALL set_ith_conn(v2e,ie,iv2e)
            END DO


            ! Frees storage
            NULLIFY(iv2e,ie2v,iv2f)
            DEALLOCATE(inspected,twice)
            CALL free_conn(e2v)
            CALL free_conn(work)

        ENDIF

        ! Stop timing
        CALL toc(sw_v2e)

        ! Synchronization
        CALL synchro(sw_v2e)

        IF(mypnum_() == 0) WRITE(*,'(a,es13.6)') &
            &'  - Elapsed time for V2E building:', partial_(sw_v2e)

100     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_V2E')

        END PROCEDURE cmp_mesh_v2e


        MODULE PROCEDURE cmp_mesh_v2v
        USE class_psblas
        USE class_connectivity
        IMPLICIT NONE
        !! $Id: cmp_mesh_v2v.f90 3093 2008-04-22 14:51:09Z sfilippo $
        !!
        !! Description:
        !!    Provides the Vertex To Vertex adjacency graph.
        !!    The criterium for associating the graph nodes is arbitrarly defined by the
        !!    user. In this case, using a co-located cell-centered strategy for variable
        !!    displacement, the V2V connectivity is obtained by linking each vertex to all
        !!    the ones belonging to the cells which share it.
        !!
        INTEGER, PARAMETER :: ncv_max = 50 ! Max # of cells sharing a vertex
        INTEGER :: i, ic, inb, info, iv, iv1, iv2, j, n, ncv, ncells
        INTEGER, POINTER :: iv2c(:) => NULL()
        INTEGER, ALLOCATABLE :: kconn(:), itab(:,:)
        LOGICAL :: inspected(nverts)
        TYPE(stopwatch) :: sw_v2v


        IF(mypnum_() == 0) WRITE(*,*) 'Computing vertices adjacency graph V2V'

        ! Constructs stopwatch
        sw_v2v = stopwatch_(icontxt_())

        ! Start timing
        CALL tic(sw_v2v)

        n = (max_conn(v2c) - 1) * ncv_max ! Maximum conn. degree (overestimated)
        ncells = nel_(v2c)

        ALLOCATE(itab(n,nverts),kconn(nverts),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        kconn = 0
        DO ic = 1, ncells
            CALL get_ith_conn(iv2c,v2c,ic)
            ncv = SIZE(iv2c)
            DO i = 1, ncv - 1
                iv1 = iv2c(i)
                DO j = i + 1, ncv
                    iv2 = iv2c(j)

                    ! Add direct arc
                    kconn(iv1) = kconn(iv1) + 1
                    itab(kconn(iv1),iv1) = iv2

                    ! Add indirect arc
                    kconn(iv2) = kconn(iv2) + 1
                    itab(kconn(iv2),iv2) = iv1
                END DO
            END DO
        END DO

        ! Eliminates duplicates
        inspected = .FALSE.
        DO iv = 1, nverts
            j = 0

            n = kconn(iv)
            DO i = 1, n
                inb = itab(i,iv)
                IF(.NOT.inspected(inb)) THEN
                    j = j + 1
                    itab(j,iv) = inb
                    inspected(inb) = .TRUE.
                END IF
            END DO
            kconn(iv) = j

            ! Reset inspected
            DO i = 1, j
                inb = itab(i,iv)
                inspected(inb) = .FALSE.
            END DO
        END DO

        ! Computes total number of connectivities
        n = 0
        DO iv = 1, nverts
            n = n + kconn(iv)
        END DO
        CALL alloc_conn(v2v,nel=nverts,nconn=n)


        ! Sets V2V connectivity
        DO iv = 1, nverts
            n = kconn(iv)
            CALL set_ith_conn(v2v,iv,itab(1:n,iv))
        END DO


        ! Frees memory
        NULLIFY(iv2c)
        DEALLOCATE(kconn,itab)


        ! Stop timing
        CALL toc(sw_v2v)

        ! Synchronization
        CALL synchro(sw_v2v)

        IF(mypnum_() == 0) WRITE(*,'(a,es13.6)') &
            & '  - Elapsed time for V2V building:', partial_(sw_v2v)

100     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_V2V')

        END PROCEDURE cmp_mesh_v2v


        MODULE PROCEDURE cmp_mesh_v2ve
        USE class_psblas
        USE class_connectivity
        USE tools_mesh, ONLY : cmp_mesh_v2e
        IMPLICIT NONE
        !! $Id: cmp_mesh_v2ve.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:
        !!    Builds a vertex to vertex connectivity based on edge connections.
        !!
        INTEGER :: i, ie, info, iv, iv_nb, nconn, nmax, nverts
        INTEGER, POINTER :: ie2v(:) => NULL(), iv2e(:) => NULL()
        INTEGER, ALLOCATABLE :: iv2v(:,:), kconn(:)
        TYPE(connectivity) :: v2e, e2v
        TYPE(stopwatch) :: sw_v2v


        IF(mypnum_() == 0) WRITE(*,*) 'Computing vertices adjacency graph V2VE'


        ! Constructs stopwatch
        sw_v2v = stopwatch_(icontxt_())

        ! Start timing
        CALL tic(sw_v2v)

        SELECT CASE(ncd)
        CASE(2)
            IF(mypnum_() == 0)    CALL get_dual_conn(v2f,e2v)

        CASE(3)
            CALL cmp_mesh_v2e(ncd,v2f,v2e)
            IF(mypnum_() == 0)    CALL get_dual_conn(v2e,e2v)
        END SELECT

        IF(mypnum_() == 0) THEN

            nverts = nel_(e2v)
            nmax   = max_conn(e2v)

            ALLOCATE(iv2v(nverts,nmax),kconn(nverts),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF


            SELECT CASE(ncd)
            CASE(2) ! 2D

                DO iv = 1, nverts

                    CALL get_ith_conn(ie2v,e2v,iv)
                    nconn = SIZE(ie2v)

                    DO i = 1, nconn
                        ie = ie2v(i)

                        ! In 2D edges = faces => use V2F
                        CALL get_ith_conn(iv2e,v2f,ie)
                        iv_nb = iv2e(1) + iv2e(2) - iv

                        iv2v(iv,i) = iv_nb
                    END DO

                    kconn(iv) = nconn
                END DO

            CASE(3) ! 3D

                DO iv = 1, nverts

                    CALL get_ith_conn(ie2v,e2v,iv)
                    nconn = SIZE(ie2v)

                    DO i = 1, nconn
                        ie = ie2v(i)

                        ! In 3D edges /= faces => use V2E
                        CALL get_ith_conn(iv2e,v2e,ie)
                        iv_nb = iv2e(1) + iv2e(2) - iv

                        iv2v(iv,i) = iv_nb
                    END DO

                    kconn(iv) = nconn
                END DO

            END SELECT

            ! Computes total number of connectivities
            nconn = 0
            DO iv = 1, nverts
                nconn = nconn + kconn(iv)
            END DO

            ! Allocate and sets V2V
            CALL alloc_conn(v2v,nel=nverts,nconn=nconn)
            DO iv = 1, nverts
                nconn = kconn(iv)
                CALL set_ith_conn(v2v,iv,iv2v(iv,1:nconn))
            END DO


            ! Free memory storage
            NULLIFY(iv2e,ie2v)
            DEALLOCATE(iv2v,kconn)
            CALL free_conn(e2v)
            IF(ncd == 3) CALL free_conn(v2e)

        ENDIF

        ! Stop timing
        CALL toc(sw_v2v)

        ! Synchronization
        CALL synchro(sw_v2v)

        IF(mypnum_() == 0) WRITE(*,'(a,es13.6)') &
            &'  - Elapsed time for V2VE building:', partial_(sw_v2v)

100     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_V2VE')

        END PROCEDURE cmp_mesh_v2ve

END SUBMODULE cmp_mesh_implementation
