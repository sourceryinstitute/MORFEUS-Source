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
! $Id: rd_gambit_mesh.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE (tools_mesh) rd_gmsh_implementation
    USE type_table, ONLY : table
    USE class_connectivity
    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE rd_gmsh_mesh
        USE class_psblas
        USE class_cell
        !USE class_connectivity
        USE class_face
        USE class_vertex
        IMPLICIT NONE

        ! --- Local variables ---
        !
        ! Parameters
        LOGICAL, PARAMETER :: debug = .TRUE.
        INTEGER, PARAMETER :: mesh = 10
        !
        ! Local copies of V2F_, V2C_, F2C_ objects (connectivity class)
        TYPE(table) :: v2f_, v2c_, f2c_
        !
        ! Vertex-related variables
        INTEGER :: nverts
        LOGICAL, ALLOCATABLE :: on_boundary(:)
        REAL(psb_dpk_), ALLOCATABLE :: xv(:), yv(:), zv(:)
        !
        ! face-related variables
        INTEGER :: nfaces
        INTEGER, ALLOCATABLE :: fnv(:), imaster(:), islave(:), iflag(:)
        INTEGER, ALLOCATABLE :: facenv(:), facemaster(:), faceslave(:), faceflag(:)
        !
        ! cell-related variables
        INTEGER :: ncells
        INTEGER, ALLOCATABLE :: cellnv(:), cellnf(:), cellgroup(:)
        CHARACTER(len=3), ALLOCATABLE :: cellgeo(:)
        !
        ! Group-related variables
        INTEGER :: ngroups
        TYPE(table) :: c2g_
        !
        ! BC-related variables
        INTEGER :: nbcfaces
        TYPE(table) :: f2b
        !
        ! Work variables
        LOGICAL :: found
        INTEGER :: i, i1, i2, info, j, j1, j2, k, l, m, n
        INTEGER :: ib, ic, ig
        INTEGER :: IF, if1, if2
        INTEGER :: iv, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8, cellid
        INTEGER, ALLOCATABLE :: perm(:), pinv(:), nvs(:), nbcf(:), faceid(:), ftype(:), facev(:,:), cellv(:,:), gcells(:)
        INTEGER, ALLOCATABLE :: aux(:,:), buf(:), work(:), pftags(:,:), pctags(:,:), facebc(:)
        TYPE(table) :: f2v, dmy
        INTEGER :: status


        CHARACTER(LEN=32) :: pname
        CHARACTER(LEN=80) :: adummy
        REAL(psb_dpk_) :: vers
        INTEGER :: ftype, dsize, nnames, pdim, pid, itemp, i1, i2, i3, i4, ngroups, &
                      dim_entt, tag_entt, nv_in_entt, edim, etag, etype, ibc, icell, ielem, &
                      ientt, iface, igroup, inv, ncmax, ncmin, nelems, nentts, npt, nv, nvmax, &
                      nvmin, tag, ibc_temp, vi(4), vj(4), jf
        REAL(psb_dpk_) :: xn, yn, zn, xx, yx, zx
        ! ------------------------------------------------------------------

        WRITE(*,*) 'Importing mesh in .MSH format on P0'
        WRITE(6,*) "Reading mesh from ", TRIM(ADJUSTL(mesh_file))

        ! Opens file *.msh
        OPEN(unit=mesh,file=ADJUSTL(mesh_file),status='old',iostat=status)
        IF(status /=0 ) THEN
            WRITE(*,050) ADJUSTL(mesh_file)
            CALL abort_psblas
        ENDIF

        READ (mesh,'()') ! $MeshFormat statement

        READ (mesh, *) vers, ftype, dsize

        ! Abort of mesh file version is less than 4.1 or if filetype is not ASCII
        IF (vers < 4.0) THEN
          WRITE(*, 100)
          CALL abort_psblas
        END IF

        DO i = 1,2
          READ (mesh, '()') ! Read section footer and header
        END DO

        ! Read all the named physical entities (materials [3D] and BCs [2D])
        ! and count the number of BCs and volumes
        nbc = 0
        ngroups = 0
        ncd = 0
        READ (mesh, '(I2)') nnames

        DO i = 1, nnames
          READ (mesh, *) pdim, pid, pname
          IF (pdim == 2) THEN
            ncd =2
            nbc = nbc + 1
          END IF
          IF (pdim == 3) THEN
            ncd = 3
            ngroups = ngroups + 1
          END IF
        END DO

        IF (ncd < 2) THEN
          WRITE(6,*) "Mesh is less than 2 dimensions"
          CALL abort_psblas
        END IF

        DO i = 1,2
          READ (mesh, '()') ! Read section footer and header
        END DO

        ! Read the entities section. This matches the numerical entities to the physical entities.
        ! The numerical entities correspond to the vertices, lines, faces, volumes that were
        ! constructed to create the geometry.
        READ (mesh, *) i1, i2, i3, i4
        DO i = 1, (i1+i2)
          READ (mesh, '()')
          ! We ignore the nodes and the lines.
        END DO
        ALLOCATE(pftags(i3, 2), pctags(i4, 2), stat=info)
        IF(info /= 0) THEN
          WRITE (*,100)
          CALL abort_psblas
        END IF
        DO i = 1, i3
          READ(mesh, *) ientt, xn, yn, zn, xx, yx, zx, npt, tag, adummy
          IF (npt == 1) THEN
            ! key-value store (ientt, tag) for the geometric and physical entity tags
            pftags(i, 1) = ientt
            pftags(i, 2) = tag
          END IF
        END DO
        DO i = 1, i4
          READ(mesh, *) ientt, xn, yn, zn, xx, yx, zx, npt, tag, adummy
          IF (npt == 1) THEN
            ! key-value store (ientt, tag) for the geometric and physical entity tags
            pctags(i, 1) = ientt
            pctags(i, 2) = tag - nbc
          END IF
        END DO

        DO i = 1,2
          READ (mesh, '()') ! Read section footer and header
        END DO

        ! Read the nodes information and their coordinates
        READ (mesh, *) nentts, nverts, nvmin, nvmax
        IF (nverts /= nvmax) THEN
          WRITE (*,100)
          write (*, *) "Node numbering is not continuous"
          CALL abort_psblas
        END IF
        ALLOCATE (xv(nverts), yv(nverts), zv(nverts),stat=info)
        IF(info /= 0) THEN
            WRITE (*,100)
            CALL abort_psblas
        END IF
        DO ientt = 1, nentts
          READ (mesh, *) dim_entt, tag_entt, itemp, nv_in_entt
          ALLOCATE (nvs(nv_in_entt))
          DO inv = 1, nv_in_entt
            READ (mesh, *) nvs(inv)
          END DO
          DO inv = 1, nv_in_entt
            nv = nvs(inv)
            READ (mesh, *) xv(nv), yv(nv), zv(nv)
          END DO
          DEALLOCATE (nvs)
        END DO

        DO i = 1,2
          READ (mesh, '()') ! Read section footer and header
        END DO

        ! Read the cells information
        READ (mesh, *) nentts, ncells, ncmin, ncmax
        IF (ncells /= ncmax) THEN
          WRITE (*,100)
          write (*, *) "Cell and face numbering is not continuous"
          CALL abort_psblas
        END IF
        ALLOCATE (nbcf(nbc), faceid(ncells), facev(ncells,4), facebc(ncells), ftype(ncells), cellgroup(ncells), &
                      cellnv(ncells), cellnf(ncells), cellgeo(ncells), cellv(ncells, 8), gcells(ngroups), stat=info)
        IF(info /= 0) THEN
            WRITE (*,100)
            CALL abort_psblas
        END IF
        icell = 0
        iface = 0
        nbcf = 0
        gcells = 0
        DO ientt = 1, nentts
          READ(mesh, *) edim, etag, etype, nelems
          IF (edim == 2) THEN
            ibc_temp = MINLOC(pftags(:,1), DIM=1, MASK=(pftags(:,1) == etag))
            ! Above call should be replaced by FINDLOC when gfortran 9.0 bugs are resolved
            ibc = pftags(ibc_temp, 2)
            nbcf(ibc) = nbcf(ibc) + nelems
            DO ielem = 1, nelems
              iface = iface + 1
              READ (mesh, *) faceid(iface), facev(iface, :)
              facebc(iface) = ibc
              ftype(iface) = etype
            END DO
          ELSE
            igroup = MINLOC(pctags(:,1), DIM=1, MASK=(pctags(:,1) == etag))
            ! Above call should be replaced by FINDLOC when gfortran 9.0 bugs are resolved
            ig = pctags(igroup, 2)
            gcells(ig) = gcells(ig) + nelems
            DO ielem = 1, nelems
              icell = icell + 1
              READ (mesh, *) cellid, cellv(icell,:)
              cellgroup(icell) = ig
              SELECT CASE(etype)
              CASE (5)
                cellnv(icell) = 8
                cellnf(icell) = 6
                cellgeo(icell) = 'hex'
              CASE default
                WRITE(*,*) 'WARNING! Unsupported celltype: ', etype, ' IGNORING'
                CONTINUE
              END SELECT
            END DO
          END IF
        END DO

        CLOSE(mesh)

        nfaces = iface
        ncells = icell
        ALLOCATE(work(ncells*8))
        CALL v2c_%alloc_table(nel=ncells)

        i1 = 1; i2 = 0; v2c_%lookup(1) = 1
        DO ic = 1, ncells
            i2 = i1 + cellnv(ic) - 1
            work(i1:i2) = cellv(ic,1:cellnv(ic))
            i1 = i2 + 1
            v2c_%lookup(ic+1) = i1
        END DO
        CALL v2c_%alloc_table(ntab=i2)

        v2c_%tab = work(1:i2)
        DEALLOCATE(work)

        ! Reads groups definition.
        CALL c2g_%alloc_table(nel=ngroups,ntab=ncells)

        i1 = 1
        DO ig = 1, ngroups
          c2g_%lookup(ig) = i1
          i2 = i1 + gcells(ig) - 1
          DO ic = i1, i2
            c2g_%tab(ic) = ic
          END DO
          i1 = i2 + 1
        END DO

        ! Creates face-cell connectivity
        !
        ! Estimation of elements number for f2c_%tab and v2f_%tab
        nfaces = 0; n = 0
        DO ic = 1, ncells
            j = cellnf(ic)
            nfaces = nfaces + j

            SELECT CASE(cellgeo(ic))
            CASE('qua')
                k = 8
            CASE('tri')
                k = 6
            CASE('hex')
                k = 24
            CASE('pri')
                k = 18
            CASE('tet')
                k = 12
            CASE('pyr')
                k = 16
            END SELECT
            n = n + k
        END DO

        ! Allocation of face-related connectivity tables
        CALL f2c_%alloc_table(nel=ncells,ntab=nfaces)
        CALL v2f_%alloc_table(nel=nfaces,ntab=n)

        ! Allocation of temporary face-related arrays
        ALLOCATE(iflag(nfaces), fnv(nfaces),  &
            &   imaster(nfaces), islave(nfaces), stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        i1 = 1; i2 = 0; f2c_%lookup(1) = 1
        j1 = 1; j2 = 0; v2f_%lookup(1) = 1
        DO ic = 1, ncells
            i2 = i1 + cellnf(ic) - 1

            iv  = v2c_%lookup(ic) - 1
            iv1 = v2c_%tab(iv+1)
            iv2 = v2c_%tab(iv+2)
            iv3 = v2c_%tab(iv+3)
            IF(cellnv(ic) > 3) iv4 = v2c_%tab(iv+4) ! quad, tet
            IF(cellnv(ic) > 4) iv5 = v2c_%tab(iv+5) ! pyr
            IF(cellnv(ic) > 5) iv6 = v2c_%tab(iv+6) ! pri
            IF(cellnv(ic) > 7) THEN                ! hex
                iv7 = v2c_%tab(iv+7)
                iv8 = v2c_%tab(iv+8)
            END IF

            DO IF = i1, i2
                f2c_%tab(IF) = IF
                iflag(IF) =0
                imaster(IF) = ic
                islave(IF) = 0 ! Further defined
            END DO

            SELECT CASE(cellgeo(ic))
            CASE('qua')
                j2 = j1 + 8 - 1
                v2f_%tab(j1:j2) = (/ iv1,iv2, iv2,iv3, iv3,iv4, iv4,iv1 /)
                v2f_%lookup(i1:i2) = (/ j1, (j1+2), (j1+4), (j1+6) /)
                fnv(i1:i2) = (/ 2, 2, 2, 2 /)
            CASE('tri')
                j2 = j1 + 6 - 1
                v2f_%tab(j1:j2) = (/ iv1,iv2, iv2,iv3, iv3,iv1 /)
                v2f_%lookup(i1:i2) = (/ j1, (j1+2), (j1+4) /)
                fnv(i1:i2) = (/ 2, 2, 2 /)
            CASE('hex')
                j2 = j1 + 24 - 1
                v2f_%tab(j1:j2) =  (/ iv1,iv2,iv6,iv5, iv2,iv3,iv7,iv6, iv3,iv4,iv8,iv7, &
                    &             iv4,iv1,iv5,iv8, iv2,iv1,iv4,iv3, iv5,iv6,iv7,iv8 /)
                v2f_%lookup(i1:i2) = (/ j1, (j1+4), (j1+8), (j1+12), (j1+16), (j1+20) /)
                fnv(i1:i2) = (/ 4, 4, 4, 4, 4, 4 /)
            CASE('pri')
                j2 = j1 + 18 - 1
                v2f_%tab(j1:j2) = (/ iv1,iv2,iv5,iv4, iv2,iv3,iv6,iv5, iv3,iv1,iv4,iv6, &
                    &             iv1,iv3,iv2, iv4,iv5,iv6 /)
                v2f_%lookup(i1:i2) = (/ j1, (j1+4), (j1+8), (j1+12), (j1+15) /)
                fnv(i1:i2) = (/ 4, 4, 4, 3, 3  /)
            CASE('tet')
                j2 = j1 + 12 - 1
                v2f_%tab(j1:j2) = (/ iv2,iv1,iv3, iv1,iv2,iv4, iv2,iv3,iv4, iv3,iv1,iv4 /)
                v2f_%lookup(i1:i2) = (/ j1, (j1+3), (j1+6), (j1+9) /)
                fnv(i1:i2) = (/ 3, 3, 3, 3  /)
            CASE('pyr')
                j2 = j1 + 16 - 1
                v2f_%tab(j1:j2) = (/ iv1,iv4,iv3,iv2, iv1,iv2,iv5, iv2,iv3,iv5, &
                    &            iv3,iv4,iv5, iv4,iv1,iv5 /)
                v2f_%lookup(i1:i2) = (/ j1, (j1+4), (j1+7), (j1+10), (j1+13) /)
                fnv(i1:i2) = (/ 4, 3, 3, 3, 3  /)
            END SELECT

            i1 = i2 + 1
            j1 = j2 + 1
            f2c_%lookup(ic+1) = i1
            v2f_%lookup(i2+1) = j1
        END DO

        ! #############################################################################
        IF(debug) THEN
            DO ic = 1, ncells
                i1 = f2c_%lookup(ic)
                i2 = f2c_%lookup(ic+1) - 1
                j1 = f2c_%tab(i1)
                j2 = f2c_%tab(i2)
                DO IF = j1, j2
                    iv1 = v2f_%lookup(IF)
                    iv2 = v2f_%lookup(IF+1) - 1
                    WRITE(*,300) ic, cellgeo(ic), IF, v2f_%tab(iv1:iv2);
                END DO
            END DO
            WRITE(*,*)
    300     FORMAT(' cell: ',i7,' type: ',a3, ' face: ',i6,' vertices:', 8(1x,i7))
        END IF
        ! ###########################################################################

        ! Builds connection table f2v, dual of v2f
        CALL v2f_%get_dual_table(f2v)

        ! Initializes face permutation array
        ALLOCATE(perm(nfaces), stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        DO IF = 1, nfaces
            perm(IF) = IF
        END DO

        ! Locates repeated faces and sets their flags = -1
        IF(ncd == 2) THEN
            pile_2D:  DO iv = 1, nverts
                i = f2v%lookup(iv) - 1! Offset
                n = f2v%lookup(iv+1) - f2v%lookup(iv) ! # of elements
                if1_2D: DO j = 1, n
                    if1 = f2v%tab(i+j)
                    IF(iflag(if1) == -1 .OR. islave(if1) /= 0) CYCLE
                    i1 = v2f_%lookup(if1) - 1
                    iv1 = v2f_%tab(i1+1) ! 1st iv of if1
                    iv2 = v2f_%tab(i1+2) ! 2nd iv of if1
                    if2_2D: DO k = j+1, n
                        if2 = f2v%tab(i+k)
                        i2 = v2f_%lookup(if2)-1
                        iv3 = v2f_%tab(i2+1) ! 1st iv of if2
                        iv4 = v2f_%tab(i2+2) ! 2nd iv of if2
                        !---
                        IF(iv4 == iv1 .AND. iv3 == iv2) THEN
                            iflag(if2) = -1
                            islave(if1) = imaster(if2)
                            imaster(if2) = 0
                            perm(if2) = if1
                            EXIT if2_2D
                        END IF
                    END DO if2_2D
                END DO if1_2D
            END DO pile_2D
        ELSE IF(ncd == 3) THEN
            ALLOCATE(buf(6),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
            ! ---
            pile_3D:  DO iv = 1, nverts
                i = f2v%lookup(iv) - 1! Offset
                n = f2v%lookup(iv+1) - f2v%lookup(iv) ! # of elements
                if1_3D: DO j = 1, n-1
                    if1=f2v%tab(i+j)
                    IF(iflag(if1) == -1 .OR. islave(if1) /= 0) CYCLE
                    ! iflag(if1) = -1  --> repeated face
                    ! islave(if1) /= 0 --> already examinated face
                    i1 = v2f_%lookup(if1) - 1
                    iv1 = v2f_%tab(i1+1)
                    iv2 = v2f_%tab(i1+2)
                    iv3 = v2f_%tab(i1+3)
                    found = .FALSE.
                    if2_3D: DO k = j+1, n
                        if2 = f2v%tab(i+k)
                        i2 = v2f_%lookup(if2)
                        j2 = v2f_%lookup(if2+1) - 1
                        l = fnv(if2)
                        ! Copies if2 iv-sequence in buf
                        buf(1) = v2f_%tab(j2)
                        buf(2:(l+1)) = v2f_%tab(i2:j2)
                        buf(l+2) = v2f_%tab(i2)
                        IF(l == 3) buf(6) = 0
                        !---
                        ! Searches for iv2 in buf
                        match: DO m = 2, (l+1)
                            IF(buf(m) /= iv2) CYCLE
                            IF(buf(m+1) == iv1 .AND. buf(m-1) == iv3) THEN
                                iflag(if2) = -1
                                islave(if1) = imaster(if2)
                                imaster(if2) = 0
                                perm(if2) = if1
                                found = .TRUE.
                            END IF
                        END DO match
                        IF(found) EXIT if2_3D
                    END DO if2_3D
                END DO if1_3D
            END DO pile_3D
            DEALLOCATE(buf)
        END IF
        CALL f2v%free_table()

        ! ###########################################################################
        IF (debug) THEN
            DO IF = 1, nfaces
                WRITE(*,310) IF, iflag(IF), imaster(IF), islave(IF)
            END DO
            WRITE(*,*)
    310     FORMAT(' face: ',i4,' flag: ',i2,' master: ',i4,' slave: ',i4)
        END IF
        ! ###########################################################################

        ! Completes the definition of permutation array
        ! new=perm(old)
        ! old=pinv(new)
        !
        ! By first faces with iflag(if)=0 ...
        j = 0
        DO IF = 1, nfaces
            IF(iflag(IF) == 0) THEN
                j = j + 1
                perm(IF) = j
            END IF
        END DO
        ! ... then faces with iflag(if)=-1
        DO IF = 1, nfaces
            IF(iflag(IF) == -1) THEN
                j = perm(IF)
                perm(IF) = perm(j)
            END IF
        END DO

        ! ###########################################################################
        IF(debug) THEN
            n = COUNT(iflag == 0)
            DO IF = 1, nfaces
                i1 = v2f_%lookup(IF)
                k = fnv(IF)
                WRITE(*,320) IF, iflag(IF), perm(IF), (v2f_%tab(j), j = i1, i1+k-1)
            END DO
            WRITE(*,*)
            DO IF = 1, nfaces
                i1 = v2f_%lookup(IF)
                k = fnv(IF)
            END DO
    320     FORMAT(' face: ',i7,' flag: ',i2,' perm: ',i7,' vertices',4i7)
        END IF
        ! ###########################################################################

        ! Reordering of face-based arrays

        ! Updates NFACES and stores OLD NFACES in N
        n = nfaces
        nfaces = COUNT(iflag == 0)

        ! Allocation of PINV and definitive face-related arrays
        ALLOCATE(pinv(nfaces), facenv(nfaces), faceflag(nfaces),&
            & facemaster(nfaces), faceslave(nfaces), stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Builds inverse permutation array for the NEW first NFACES faces
        DO IF = 1, n
            IF (iflag(IF) == 0) THEN
                j = perm(IF)
                pinv(j) = IF
            END IF
        END DO

        ! Reordering and assignement of definitive face-related arrays
        DO IF = 1, nfaces
            j = pinv(IF)

            facenv(IF)     = fnv(j)
            faceflag(IF)   = 0
            facemaster(IF) = imaster(j)
            faceslave(IF)  = islave(j)
        END DO
        DEALLOCATE(iflag,imaster,islave,fnv)

        ! Reordering and reallocation of v2f_
        k = 0
        DO IF = 1, nfaces
            j = k
            k = j + facenv(IF)
        END DO
        !
        CALL dmy%alloc_table(nel=nfaces,ntab=k)
        i1 = 1; i2 = 0; dmy%lookup(1) = 1
        DO IF = 1, nfaces
            j = pinv(IF)
            i2 = i1 + facenv(IF) - 1
            j1 = v2f_%lookup(j)
            j2 = j1 + facenv(IF) - 1
            dmy%tab(i1:i2) = v2f_%tab(j1:j2)
            i1 = i2 + 1
            dmy%lookup(IF+1) = i1
        END DO
        CALL v2f_%free_table()

        CALL v2f_%alloc_table(nel=nfaces,ntab=k)
        v2f_%lookup = dmy%lookup(1:nfaces+1)
        v2f_%tab = dmy%tab(1:k)

        ! Reordering of f2c_%tab
        DO ic = 1, ncells
            i1 = f2c_%lookup(ic)
            i2 = f2c_%lookup(ic+1) - 1
            DO i = i1, i2
                j = f2c_%tab(i)
                f2c_%tab(i) = perm(j)
            END DO
        END DO

        CALL dmy%free_table()
        DEALLOCATE(perm,pinv)

        ! ###########################################################################
        IF(debug) THEN
            DO IF = 1, nfaces
                i1 = v2f_%lookup(IF)
                i2 = v2f_%lookup(IF+1) - 1
                WRITE(*,330) IF, facemaster(IF), faceslave(IF), v2f_%tab(i1:i2)
            END DO
            WRITE(*,*)
            DO ic = 1, ncells
                i1 = f2c_%lookup(ic)
                i2 = f2c_%lookup(ic+1) - 1
                WRITE(*,340) ic, f2c_%tab(i1:i2)
            END DO
            WRITE(*,*)
    330     FORMAT(' face: ',i7, ' master: ',i7,' slave: ',i7,' vertices: ',4i7)
    340     FORMAT(' cell: ',i7, ' faces:',6(1x,i7))
        END IF
        ! ###########################################################################

        ! - Reads and assigns boundary conditions. Boundary flags are integer
        ! in the range [1-nbc]. If NBC = 0 from import, automatically NBC = 1 and where
        ! faceslave = 0 a boundary face is asssumed.
        ! - Builds list of faces, grouping them by their flag value.
        !
        ! faceflag(if)=0 --> fluid face
        ! faceflag(if)=i  --> i-th boundary face
        !

        ! Counts boundary faces (faceslave = 0)
        k = COUNT(faceslave == 0)

        IF (nbc > 0) THEN
          ALLOCATE(f2b%lookup(nbc), f2b%tab(k))
          i1 = 1
          DO ib = 1, nbc
            f2b%lookup(ib) = i1
            i2 = i1 + nbcf(ib) - 1
            DO j = 1, nbcf(ib)
              if = i1 + j - 1
              vi = facev(if, :)
              jf_loop: DO jf = 1, nfaces
                IF (faceslave(jf) /= 0) CYCLE
                j1 = v2f_%lookup(jf)
                j2 = v2f_%lookup(jf+1) - 1
                vj = v2f_%tab(j1:j2)
                IF (ANY(vi == vj(1)) .AND. ANY(vi == vj(2)) .AND. ANY(vi == vj(3))) THEN
                  f2b%tab(if) = jf
                  faceflag(jf) = ib
                  EXIT jf_loop
                END IF
              END DO jf_loop
            END DO
            i1 = i2 + 1
          END DO

        ELSE
          WRITE( *, *) "ERROR! No boundary conditions in GMSH file."
          CALL abort_psblas
        END IF

        ! Reordering of if-indexing following the flag sequence: 0,1,2,...,nbc
        ! f2b%tab is used for building the inverse permutation array pinv
        ! old = pinv(new)
        ! new = perm(old)

        k=SIZE(v2f_%tab)
        ALLOCATE(perm(nfaces),pinv(nfaces),aux(4,nfaces),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        CALL dmy%alloc_table(nel=nfaces,ntab=k)

        ! Builds pinv
        k = 0
        ! First fluid faces ...
        DO IF = 1, nfaces
            IF(faceflag(IF) == 0) THEN
                k = k + 1
                pinv(k) = IF
            END IF
        END DO
        ! ... then bc faces.

        IF ( nfaces-k /= SIZE(f2b%tab) ) THEN
            WRITE(6,*) nfaces, k, nfaces-k, size(f2b%tab)
            WRITE(6,400)
            WRITE(6,*)"Check that the CAD model exported all faces."
            CALL abort_psblas
        ENDIF

        pinv(k+1:nfaces) = f2b%tab

        DO IF = 1, nfaces
            aux(1,IF) = facenv(IF)
            aux(2,IF) = faceflag(IF)
            aux(3,IF) = facemaster(IF)
            aux(4,IF) = faceslave(IF)
        END DO

        DO IF = 1, nfaces
            j = pinv(IF)
            perm(j) = IF ! Permutation array
            facenv(IF) = aux(1,j)
            faceflag(IF) = aux(2,j)
            facemaster(IF) = aux(3,j)
            faceslave(IF) = aux(4,j)
        END DO

        ! Reordering of v2f_%tab and v2f_%lookup
        k  = SIZE(v2f_%tab)
        dmy%lookup  = v2f_%lookup(1:nfaces+1)
        dmy%tab = v2f_%tab(1:k)
        i1 = 1; i2 = 0; v2f_%tab(i1) = 1
        DO IF = 1, nfaces
            i2 = i1 + facenv(IF) - 1

            j = pinv(IF)
            j1 = dmy%lookup(j)
            j2 = j1 + facenv(IF) - 1

            v2f_%tab(i1:i2) = dmy%tab(j1:j2)

            i1 = i2 + 1
            v2f_%lookup(IF+1) = i1
        END DO

        ! Reordering of f2c_%tab
        k = SIZE(f2c_%tab)
        DO i = 1, k
            j = f2c_%tab(i)
            f2c_%tab(i) = perm(j)
        END DO

        ! Counts total number of bc faces
        nbcfaces = COUNT(faceslave == 0)

        ! Redifinition of f2b%tab
        k = COUNT(faceflag == 0) ! Fluid faces
        DO IF = 1, nbcfaces
            f2b%tab(IF) = IF + k
        END DO

        ! ###########################################################################
        IF (debug) THEN
            WRITE(*,350) (f2b%tab(i), i = 1, f2b%lookup(1) -1)
            WRITE(*,*)

            DO ib = 1, nbc
                i1 = f2b%lookup(ib)
                i2 = i1 + nbcf(ib) - 1
                WRITE(*,360) ib, f2b%tab(i1:i2)
            END DO
            WRITE(*,*)

            DO IF = 1, nfaces
                WRITE(*,370) IF, faceflag(IF), facemaster(IF), faceslave(IF)
            END DO
            WRITE(*,*)

            DO i = 0, nbc
                WRITE(*,380) i, COUNT(faceflag == i)
            END DO

            WRITE(*,*)
    350     FORMAT(' Fluid faces: ',5i8:/(14x,5i8))
    360     FORMAT(' Bcset: ',i2,' List of faces: ', 5i8:/(40x,5i8:))
    370     FORMAT(' faces: ',i7,' Flag: ',i2,' Master: ',i7,' Slave: ',i7)
    380     FORMAT(' Total amount of faces with flag',i3,':',i8)
    400     FORMAT(' Fatal Error in RD_GAMBIT_MESH:  inconsistent face, boundary lists.')
        END IF
        ! ###########################################################################

        ! f2b is no longer useful. It can be deallocated
        CALL f2b%free_table()
        CALL dmy%free_table()
        DEALLOCATE(perm,pinv,aux,nbcf)

        ! Defines ON_BOUNDARY array
        ALLOCATE(on_boundary(nverts),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        on_boundary = .FALSE.

        k = COUNT(faceflag == 0) + 1 ! Offset

        ! Loop on boundary faces
        DO IF = k, nfaces
            i1 = v2f_%lookup(IF)
            i2 = v2f_%lookup(IF+1) - 1
            DO i = i1, i2
                iv = v2f_%tab(i)
                on_boundary(iv) = .TRUE.
            END DO
        END DO

        ! Copies vertex-related data in to the object VERTS of VERTEX class
        CALL alloc_vertex(verts,nverts)
        verts = vertex_(xv,yv,zv,on_boundary)
        DEALLOCATE(xv,yv,zv,on_boundary)

        ! Copies face-related data in to the object FACES of FACE class
        CALL alloc_face(faces,nfaces)
        faces = face_(facenv,facemaster,faceslave,faceflag)
        DEALLOCATE(facenv,facemaster,faceslave,faceflag)

        ! Copies cell-related data in to the object CELLS of CELL class
        CALL alloc_cell(cells,ncells)
        cells = cell_(cellnv,cellnf,cellgroup,cellgeo)
        DEALLOCATE(cellnv,cellnf,cellgroup,cellgeo)

        ! Copies connectivities  V2F, V2C, F2C to dummy arguments V2F_, V2C_, F2C_
        ! (objects of CONNECTIVITY class)
        CALL alloc_conn(v2f,nel=nfaces,nconn=SIZE(v2f_%tab))
        CALL alloc_conn(v2c,nel=ncells,nconn=SIZE(v2c_%tab))
        CALL alloc_conn(f2c,nel=ncells,nconn=SIZE(f2c_%tab))
        CALL alloc_conn(c2g,nel=ngroups,nconn=SIZE(c2g_%tab))

        DO IF = 1, nfaces
            i1 = v2f_%lookup(IF)
            i2 = v2f_%lookup(IF+1) - 1
            CALL v2f%set_ith_conn(IF,v2f_%tab(i1:i2))
        END DO

        DO ic = 1, ncells
            i1 = v2c_%lookup(ic)
            i2 = v2c_%lookup(ic+1) - 1
            CALL v2c%set_ith_conn(ic,v2c_%tab(i1:i2))
        END DO

        DO ic = 1, ncells
            i1 = f2c_%lookup(ic)
            i2 = f2c_%lookup(ic+1) - 1
            CALL f2c%set_ith_conn(ic,f2c_%tab(i1:i2))
        END DO

        DO ig = 1, ngroups
            i1 = c2g_%lookup(ig)
            i2 = c2g_%lookup(ig+1) - 1
            CALL c2g%set_ith_conn(ig,c2g_%tab(i1:i2))
        END DO

        ! Deallocates local copies of connectivity data
        CALL v2f_%free_table()
        CALL v2c_%free_table()
        CALL f2c_%free_table()
        CALL c2g_%free_table()

        WRITE(*,'()')

        !     FORMAT Instructions
050     FORMAT(' ERROR! Failure to open Gambit mesh file.',/,&
        &    ' File expected: ', a)
100     FORMAT(' ERROR! Memory allocation failure in RD_GMSH_MESH')


    END PROCEDURE rd_gmsh_mesh

END SUBMODULE rd_gmsh_implementation
