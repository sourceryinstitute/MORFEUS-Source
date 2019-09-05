!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!

SUBMODULE (tools_mesh) rd_exodus_mesh_implementation
    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE rd_exodus_mesh

        USE ISO_FORTRAN_ENV, ONLY : INT64, REAL64
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_face
        USE class_vertex
        USE tools_math
        USE type_table



        IMPLICIT NONE
        INCLUDE 'exodusII.inc'
        INTEGER, PARAMETER :: nlen = 80


        ! Local copies of table_v2f, table_v2c, F2C_ objects (connectivity class)
        TYPE(table) :: c2g_, table_v2c, f2c_, table_v2f

        INTEGER(kind = int64) :: ierr, titlelen
        INTEGER(kind = int64) :: exodus_file_id, num_dims, num_nodes, num_elems,  &
            num_elem_blks, num_node_sets, num_side_sets

        INTEGER(kind = int64) :: cpu_ws, io_ws, num_props, prop_value

        CHARACTER(len = MXLNLN) :: title
        CHARACTER(len = MXSTLN) :: coord_names(3), qa_record(4,2), var_names(3)
        REAL :: vers
        REAL(kind = real64), ALLOCATABLE :: xx(:), yy(:), zz(:), dist_fact(:,:)
        INTEGER(kind = int64), ALLOCATABLE :: elem_map(:), elem_blk_ids(:), &
            num_elems_in_blk(:), num_nodes_per_elem(:), num_attr_per_elem(:), &
            conn(:,:), node_set_ids(:), num_nodes_per_set(:), &
            num_df_per_node_set(:), num_df_per_side_set(:), side_set_ids(:), &
            num_sides_per_set(:), node_list(:,:), side_list(:,:), elem_list(:,:)
        CHARACTER(len = MXSTLN), ALLOCATABLE :: elem_blk_typ(:), prop_names(:)

        INTEGER :: ncells, nfaces
        INTEGER :: istart, iend, ielem
        INTEGER, ALLOCATABLE :: fnv(:), imaster(:), islave(:), iflag(:)
        INTEGER, ALLOCATABLE :: facenv(:), facemaster(:), faceslave(:), faceflag(:)
        INTEGER, ALLOCATABLE :: cellnv(:), cellnf(:), cellgroup(:)
        CHARACTER(len=3), ALLOCATABLE :: cellgeo(:)
        !
        ! Vertex-related variables
        INTEGER :: nverts
        LOGICAL, ALLOCATABLE :: on_boundary(:)
        !
        ! Group-related variables
        INTEGER :: ngroups
        INTEGER, ALLOCATABLE :: groupnc(:), groupmat(:)
        CHARACTER(len=32), ALLOCATABLE:: groupname(:)
        !
        ! BC-related variables
        INTEGER :: nbcfaces
        TYPE(table) :: f2b
        !
        ! Work variables
        LOGICAL :: found
        INTEGER :: i, i1, i2, info, j, j1, j2, k,k1, k2, l, m, n
        INTEGER :: ib, ic, ig
        INTEGER :: IF, if1, if2
        INTEGER :: iv, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8
        INTEGER, ALLOCATABLE :: perm(:), pinv(:)
        INTEGER, ALLOCATABLE :: aux(:,:), buf(:), work(:)
        INTEGER, ALLOCATABLE :: work1(:), work2(:)
        TYPE(table) :: dmy, f2v

        LOGICAL, PARAMETER :: debug = .FALSE.


        ! open exodusII file for reading

        PRINT *, "Reading mesh from ", TRIM(ADJUSTL(mesh_file))

        cpu_ws = 0
        io_ws = 0

        exodus_file_id = exopen(TRIM(ADJUSTL(mesh_file)), EXREAD, cpu_ws, io_ws, vers, ierr)

        ! read the database parameters from the file
        CALL exgini(exodus_file_id, title, num_dims, num_nodes, num_elems, &
            num_elem_blks, num_node_sets, num_side_sets, ierr)

        ngroups = num_elem_blks
        ncells = num_elems
        nbc = num_side_sets
        nverts = num_nodes
        ncd = num_dims
        ALLOCATE(cellnv(ncells), cellnf(ncells), cellgeo(ncells), cellgroup(ncells))

        ! read node coordinates
        ALLOCATE(xx(num_nodes), yy(num_nodes), zz(num_nodes))
        CALL exgcor(exodus_file_id, xx, yy, zz, ierr)

        ! read coordinate names
        CALL exgcon(exodus_file_id, coord_names, ierr)

        ! read element map
        ALLOCATE(elem_map(num_elems))
        CALL exgmap(exodus_file_id, elem_map, ierr)

        ! read element block parameters
        ALLOCATE(elem_blk_ids(num_elem_blks), elem_blk_typ(num_elem_blks), &
            num_elems_in_blk(num_elem_blks), num_nodes_per_elem(num_elem_blks), &
            num_attr_per_elem(num_elem_blks))
        CALL exgebi(exodus_file_id, elem_blk_ids, ierr)

        DO i = 1, num_elem_blks
            CALL exgelb(exodus_file_id, elem_blk_ids(i), elem_blk_typ(i), &
                num_elems_in_blk(i), num_nodes_per_elem(i), &
                num_attr_per_elem(i), ierr)
        END DO

        ALLOCATE(conn(MAXVAL(num_nodes_per_elem), MAXVAL(num_elems_in_blk)))
        ALLOCATE(work(MAXVAL(num_nodes_per_elem)*num_elems))
        ! Read element block conn

        CALL table_v2c%alloc_table(nel=ncells)
        istart = 1
        i1 = 1; i2 = 0; table_v2c%lookup(1) = 1
        DO i = 1,num_elem_blks
            conn = -1
            CALL exgelc(exodus_file_id, elem_blk_ids(i), conn, ierr)
            iend = istart + num_elems_in_blk(i) - 1
            do ielem = 1, num_elems_in_blk(i)
                ic = istart + ielem - 1
                SELECT CASE(num_nodes_per_elem(i))
                CASE(4)
                    cellnv(ic) = 4
                    cellnf(ic) = 4
                    cellgroup(ic) = i
                    cellgeo(ic) = 'tet'
                CASE(5)
                    cellnv(ic) = 5
                    cellnf(ic) = 5
                    cellgroup(ic) = i
                    cellgeo(ic) = 'pyr'
                CASE(6)
                    cellnv(ic) = 6
                    cellnf(ic) = 5
                    cellgroup(ic) = i
                    cellgeo(ic) = 'pri'
                CASE(8)
                    cellnv(ic) = 8
                    cellnf(ic) = 6
                    cellgroup(ic) = i
                    cellgeo(ic) = 'hex'
                CASE default
                    WRITE(*,*) 'WARNING! Unsupported number of nodes per cell: ', num_nodes_per_elem(i), ' IGNORING'
                    CONTINUE
                END SELECT
                i2 = i1 + cellnv(ic) - 1
                work(i1:i2) = conn(:,ielem)
                i1 = i2 + 1
                table_v2c%lookup(ic+1) = i1
            end do
            istart = iend + 1
        ENDDO

        CALL table_v2c%alloc_table(ntab=i2)
        table_v2c%tab = work(1:i2)

        ! Reads groups definition.
        CALL c2g_%alloc_table(nel=ngroups,ntab=ncells)

        i1 = 1
        DO ig = 1, num_elem_blks
            c2g_%lookup(ig) = i1
            i2 = i1 + num_elems_in_blk(ig) - 1
            DO ic = i1, i2
                c2g_%tab(ic) = ic
            END DO
            i1 = i2 + 1
        END DO

        ! Creates face-cell connectivity
        ! NOTE: from here until reading of bcs the code is equal to import_neu.f90
        !
        ! Estimation of elements number for f2c%tab and v2f%tab
        nfaces = 0; n = 0
        DO ic = 1, ncells
            j = cellnf(ic)
            nfaces = nfaces + j
            ! ---
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
        CALL table_v2f%alloc_table(nel=nfaces, ntab=n)

        ! Allocation of temporary face-related arrays
        ALLOCATE(iflag(nfaces), fnv(nfaces),  &
            &   imaster(nfaces), islave(nfaces), stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        i1 = 1; i2 = 0; f2c_%lookup(1) = 1
        j1 = 1; j2 = 0; table_v2f%lookup(1) = 1
        DO ic = 1, ncells
            i2 = i1 + cellnf(ic) - 1
            !
            iv = table_v2c%lookup(ic) - 1
            iv1 = table_v2c%tab(iv+1)
            iv2 = table_v2c%tab(iv+2)
            iv3 = table_v2c%tab(iv+3)
            IF(cellnv(ic) > 3) iv4 = table_v2c%tab(iv+4) ! quad, tet
            IF(cellnv(ic) > 4) iv5 = table_v2c%tab(iv+5) ! pyr
            IF(cellnv(ic) > 5) iv6 = table_v2c%tab(iv+6) ! pri
            IF(cellnv(ic) > 7) THEN ! hex
                iv7 = table_v2c%tab(iv+7)
                iv8 = table_v2c%tab(iv+8)
            END IF
            DO IF = i1, i2
                f2c_%tab(IF) = IF
                iflag(IF) = 0
                imaster(IF) = ic
                islave(IF) = 0 ! Further defined
            END DO

            ! NOTE: standard connectivity vertices-cells defined in SIDS_guide
            SELECT CASE(cellgeo(ic))
            CASE('qua')
                j2 = j1 + 8 - 1
                table_v2f%tab(j1:j2) = (/ iv1,iv2, iv2,iv3, iv3,iv4, iv4,iv1 /)
                table_v2f%lookup(i1:i2) = (/ j1, (j1+2), (j1+4), (j1+6) /)
                fnv(i1:i2) = (/ 2, 2, 2, 2 /)
            CASE('tri')
                j2 = j1 + 6 - 1
                table_v2f%tab(j1:j2) = (/ iv1,iv2, iv2,iv3, iv3,iv1 /)
                table_v2f%lookup(i1:i2) = (/ j1, (j1+2), (j1+4) /)
                fnv(i1:i2) = (/ 2, 2, 2 /)
            CASE('hex')
                j2 = j1 + 24 - 1
                table_v2f%tab(j1:j2) = (/ iv1,iv2,iv6,iv5, iv2,iv3,iv7,iv6, &
                    &              iv3,iv4,iv8,iv7, iv1,iv5,iv8,iv4, iv1,iv4,iv3,iv2, iv5,iv6,iv7,iv8 /)
                table_v2f%lookup(i1:i2) = (/ j1, (j1+4), (j1+8), (j1+12), (j1+16), (j1+20) /)
                fnv(i1:i2) = (/ 4, 4, 4, 4, 4, 4 /)
            CASE('pri')
                j2 = j1 + 18 - 1
                table_v2f%tab(j1:j2) = (/ iv1,iv2,iv5,iv4, iv2,iv3,iv6,iv5, iv1,iv4,iv6,iv3, &
                    &              iv1,iv3,iv2, iv4,iv5,iv6 /)
                table_v2f%lookup(i1:i2) = (/ j1, (j1+4), (j1+8), (j1+12), (j1+15) /)
                fnv(i1:i2) = (/ 4, 4, 4, 3, 3  /)
            CASE('tet')
                j2 = j1 + 12 - 1
                table_v2f%tab(j1:j2) = (/ iv1,iv2,iv4, iv2,iv3,iv4, iv1,iv4,iv3, iv1,iv3,iv2 /)
                table_v2f%lookup(i1:i2) = (/ j1, (j1+3), (j1+6), (j1+9) /)
                fnv(i1:i2) = (/ 3, 3, 3, 3  /)
            CASE('pyr')
                j2 = j1 + 16 - 1
                table_v2f%tab(j1:j2) = (/ iv1,iv2,iv5, iv2,iv3,iv5, &
                    &              iv3,iv4,iv5, iv1,iv5,iv4, iv1,iv4,iv3,iv2 /)
                table_v2f%lookup(i1:i2) = (/ j1, (j1+3), (j1+6), (j1+9), (j1+12) /)
                fnv(i1:i2) = (/ 3, 3, 3, 3, 4  /)
            END SELECT
            ! ---
            i1 = i2 + 1
            j1 = j2 + 1
            f2c_%lookup(ic+1) = i1
            table_v2f%lookup(i2+1) = j1
        END DO

        ! #############################################################################
        IF(debug) THEN
            DO ic = 1, ncells
                i1 = f2c_%lookup(ic)
                i2 = f2c_%lookup(ic+1) - 1
                j1 = f2c_%tab(i1)
                j2 = f2c_%tab(i2)
                DO IF = j1, j2
                    iv1 = table_v2f%lookup(IF)
                    iv2 = table_v2f%lookup(IF+1) - 1
                    WRITE(*,300) ic, cellgeo(ic), IF, table_v2f%tab(iv1:iv2);
                END DO
            END DO
            WRITE(*,*)
    300     FORMAT(' cell: ',i7,' type: ',a3, ' face: ',i7,' vertices:', 8(1x,i5))
        END IF
        ! ###########################################################################


        ! Builds connection table f2v, dual of v2f

        CALL table_v2f%get_dual_table(f2v)

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
                i = f2v%lookup(iv) - 1 ! Offset
                n = f2v%lookup(iv+1) - f2v%lookup(iv) ! # of elements
                if1_2D: DO j = 1, n
                    if1 = f2v%tab(i+j)
                    IF(iflag(if1) == -1 .OR. islave(if1) /= 0) CYCLE
                    i1 = table_v2f%lookup(if1) - 1
                    iv1 = table_v2f%tab(i1+1) ! 1st iv of if1
                    iv2 = table_v2f%tab(i1+2) ! 2nd iv of if1
                    if2_2D: DO k = j + 1, n
                        if2 = f2v%tab(i+k)
                        i2 = table_v2f%lookup(if2)-1
                        iv3 = table_v2f%tab(i2+1) ! 1st iv of if2
                        iv4 = table_v2f%tab(i2+2) ! 2nd iv of if2
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
                if1_3D: DO j = 1, n - 1
                    if1 = f2v%tab(i+j)
                    IF(iflag(if1) == -1 .OR. islave(if1) /= 0) CYCLE
                    ! iflag(if1) = -1  --> repeated face
                    ! islave(if1) /= 0 --> already examinated face
                    i1 = table_v2f%lookup(if1) - 1
                    iv1 = table_v2f%tab(i1+1)
                    iv2 = table_v2f%tab(i1+2)
                    iv3 = table_v2f%tab(i1+3)
                    found = .FALSE.
                    if2_3D: DO k = j + 1, n
                        if2 = f2v%tab(i+k)
                        i2 = table_v2f%lookup(if2)
                        j2 = table_v2f%lookup(if2+1)-1
                        l = fnv(if2)
                        ! Copies if2 iv-sequence in buf
                        buf(1) = table_v2f%tab(j2)
                        buf(2:(l+1)) = table_v2f%tab(i2:j2)
                        buf(l+2) = table_v2f%tab(i2)
                        IF(l == 3) buf(6)=0
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

        ! Completes definition of permutation array
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
                i1 = table_v2f%lookup(IF)
                k = fnv(IF)
                WRITE(*,320) IF, iflag(IF), perm(IF), (table_v2f%tab(j), j = i1, i1+k-1)
            END DO
            WRITE(*,*)
            DO IF = 1,nfaces
                i1 = table_v2f%lookup(IF)
                k = fnv(IF)
            END DO
    320     FORMAT(' face: ',i5,' flag: ',i2,' perm: ',i5,' vertices',4i5)
        END IF
        ! ###########################################################################

        ! Reordering of face based arrays

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

        ! Reordering and definition of FACE face-object
        DO IF = 1, nfaces
            j = pinv(IF)

            facenv(IF)     = fnv(j)
            faceflag(IF)   = 0
            facemaster(IF) = imaster(j)
            faceslave(IF)  = islave(j)
        END DO
        DEALLOCATE(iflag,imaster,islave,fnv)

        ! Reordering and reallocation of v2f
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
            j1 = table_v2f%lookup(j)
            j2 = j1 + facenv(IF) - 1
            dmy%tab(i1:i2) = table_v2f%tab(j1:j2)
            i1 = i2 + 1
            dmy%lookup(IF+1) = i1
        END DO
        CALL table_v2f%free_table()

        CALL table_v2f%alloc_table(nel=nfaces,ntab=k)
        table_v2f%lookup = dmy%lookup(1:nfaces+1)
        table_v2f%tab = dmy%tab(1:k)

        ! Reordering of f2c%tab
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
                i1 = table_v2f%lookup(IF)
                i2 = table_v2f%lookup(IF+1) - 1
                WRITE(*,330) IF, facemaster(IF), faceslave(IF), table_v2f%tab(i1:i2)
            END DO
            WRITE(*,*)
            DO ic = 1, ncells
                i1 = f2c_%lookup(ic)
                i2 = f2c_%lookup(ic+1) - 1
                WRITE(*,340) ic, f2c_%tab(i1:i2)
            END DO
            WRITE(*,*)
    330     FORMAT(' face: ',i5,' master: ',i5,' slave: ',i5,' vertices: ',4i5)
    340     FORMAT(' cell: ',i5, ' faces:',6(1x,i5))
        END IF
        ! ###########################################################################

        ! Read and assign boundary conditions

        ! Counts boundary faces (faceslave = 0)

        k = COUNT(faceslave == 0)
        PRINT *, "FACESLAVES", k

        IF (num_side_sets > 0) THEN
            ALLOCATE(side_set_ids(num_side_sets), num_sides_per_set(num_side_sets), &
                num_df_per_side_set(num_side_sets))

            CALL exgssi(exodus_file_id, side_set_ids, ierr)

            DO i = 1, num_side_sets
                CALL exgsp(exodus_file_id, side_set_ids(i), num_sides_per_set(i), &
                    num_df_per_side_set(i), ierr)
            END DO

            ALLOCATE(side_list(num_side_sets, MAXVAL(num_sides_per_set)), &
                elem_list(num_side_sets, MAXVAL(num_sides_per_set)))
            DO i = 1, num_side_sets
                CALL exgss(exodus_file_id, side_set_ids(i), elem_list(i,:), side_list(i,:), ierr)
            ENDDO

            ALLOCATE(f2b%lookup(nbc+1),f2b%tab(k))

            i1 = 1
            DO ib = 1, num_side_sets
                print *, ib, i1, num_side_sets, num_sides_per_set(ib)
                f2b%lookup(ib) = i1
                i2 = i1 + num_sides_per_set(ib) - 1
                DO j = 1, num_sides_per_set(ib)
                    ic = elem_list(ib, j)
                    l = side_list(ib, j)
                    m = f2c_%lookup(ic) - 1
                    if = f2c_%tab(m + l)

                    f2b%tab(i1+j-1) = if
                    faceflag(if) = ib

                END DO
                i1 = i2 + 1
            END DO

        ELSE
            WRITE( *, *) "ERROR! No boundary conditions in exodus file."
            CALL abort_psblas

        END IF

        PRINT *, "Finished reading EXODUS mesh file"
        CALL exclos(exodus_file_id, ierr)


        ! Reordering of if-indexing following the flag sequence: 0,1,2,...,nbc
        ! f2b%tab is used for building the inverse permutation array pinv
        ! old = pinv(new)
        ! new = perm(old)

        k=SIZE(table_v2f%tab)
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
            PRINT *, nfaces, k, nfaces-k, size(f2b%tab)
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

        ! Reordering of table_v2f%tab and table_v2f%lookup
        k  = SIZE(table_v2f%tab)
        dmy%lookup  = table_v2f%lookup(1:nfaces+1)
        dmy%tab = table_v2f%tab(1:k)
        i1 = 1; i2 = 0; table_v2f%tab(i1) = 1
        DO IF = 1, nfaces
            i2 = i1 + facenv(IF) - 1

            j = pinv(IF)
            j1 = dmy%lookup(j)
            j2 = j1 + facenv(IF) - 1

            table_v2f%tab(i1:i2) = dmy%tab(j1:j2)

            i1 = i2 + 1
            table_v2f%lookup(IF+1) = i1
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
                i2 = i1 + num_sides_per_set(ib) - 1
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
    400     FORMAT(' Fatal Error in RD_EXODUS_MESH:  inconsistent face, boundary lists.')
        END IF
        ! ###########################################################################

        ! f2b is no longer useful. It can be deallocated
        CALL f2b%free_table()
        CALL dmy%free_table()
        DEALLOCATE(perm,pinv,aux)

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
            i1 = table_v2f%lookup(IF)
            i2 = table_v2f%lookup(IF+1) - 1
            DO i = i1, i2
                iv = table_v2f%tab(i)
                on_boundary(iv) = .TRUE.
            END DO
        END DO

        ! Copies vertex-related data in to the object VERTS of VERTEX class
        CALL alloc_vertex(verts,nverts)
        verts = vertex_(xx,yy,zz,on_boundary)
        DEALLOCATE(xx,yy,zz,on_boundary)

        ! Copies face-related data in to the object FACES of FACE class
        CALL alloc_face(faces,nfaces)
        faces = face_(facenv,facemaster,faceslave,faceflag)
        DEALLOCATE(facenv,facemaster,faceslave,faceflag)

        ! Copies cell-related data in to the object CELLS of CELL class
        CALL alloc_cell(cells,ncells)
        cells = cell_(cellnv,cellnf,cellgroup,cellgeo)
        DEALLOCATE(cellnv,cellnf,cellgroup,cellgeo)

        ! Copies connectivities  V2F, V2C, F2C to dummy arguments table_v2f, table_v2c, F2C_
        ! (objects of CONNECTIVITY class)
        CALL alloc_conn(v2f,nel=nfaces,nconn=SIZE(table_v2f%tab))
        CALL alloc_conn(v2c,nel=ncells,nconn=SIZE(table_v2c%tab))
        CALL alloc_conn(f2c,nel=ncells,nconn=SIZE(f2c_%tab))
        CALL alloc_conn(c2g,nel=ngroups,nconn=SIZE(c2g_%tab))

        DO IF = 1, nfaces
            i1 = table_v2f%lookup(IF)
            i2 = table_v2f%lookup(IF+1) - 1
            CALL v2f%set_ith_conn(IF,table_v2f%tab(i1:i2))
        END DO

        DO ic = 1, ncells
            i1 = table_v2c%lookup(ic)
            i2 = table_v2c%lookup(ic+1) - 1
            CALL v2c%set_ith_conn(ic,table_v2c%tab(i1:i2))
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
        CALL table_v2f%free_table()
        CALL table_v2c%free_table()
        CALL f2c_%free_table()
        CALL c2g_%free_table()

        ! Currently group data GROUPNC, GROUPMAT, GROUPNAME are not exported
        ! to the calling program => deallocation
        ! DEALLOCATE(groupnc,groupmat,groupname)

        WRITE(*,'()')

        !     FORMAT Instructions
050     FORMAT(' ERROR! Failure to open Gambit mesh file.',/,&
        &    ' File expected: ', a)
100     FORMAT(' ERROR! Memory allocation failure in RD_EXODUS_MESH')
500     FORMAT(a80)
510     FORMAT(6(1x,i9))
520     FORMAT(i10,3(e20.11))
530     FORMAT(i8,1x,i2,1x,i2,1x)
540     FORMAT(7i8:/(t16,7i8:))
550     FORMAT(t8,i10,tr11,i10,tr11,i10,tr9,i10)
560     FORMAT(a32)
570     FORMAT((10i8))
580     FORMAT(a32,2i8)
590     FORMAT(i10,2i7)


    END PROCEDURE rd_exodus_mesh

END SUBMODULE rd_exodus_mesh_implementation
