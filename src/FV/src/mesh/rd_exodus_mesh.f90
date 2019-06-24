!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
SUBROUTINE rd_exodus_mesh(mesh_file, mesh_id, nbc, ncd, &
    & verts, faces, cells, v2f_, v2c_, f2c_, c2g_)

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

    CHARACTER(len = nlen), INTENT(IN) :: mesh_file
    CHARACTER(len = nlen), INTENT(OUT) :: mesh_id
    INTEGER, INTENT(OUT) :: nbc
    INTEGER, INTENT(OUT) :: ncd
    TYPE(vertex), ALLOCATABLE :: verts(:)
    TYPE(face), ALLOCATABLE :: faces(:)
    TYPE(cell), ALLOCATABLE :: cells(:)
    TYPE(connectivity), INTENT(OUT) :: v2f_, v2c_, f2c_, c2g_

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

    INTEGER :: nverts, ncells, nfaces
    INTEGER, ALLOCATABLE :: fnv(:), imaster(:), islave(:), iflag(:)
    INTEGER, ALLOCATABLE :: facenv(:), facemaster(:), faceslave(:), faceflag(:)
    INTEGER, ALLOCATABLE :: cellnv(:), cellnf(:)
    CHARACTER(len=3), ALLOCATABLE :: cellgeo(:)

    ! Local copies of V2F_, V2C_, F2C_ objects (connectivity class)
    TYPE(table) :: v2f, v2c, f2c

    ! Work variables
    LOGICAL :: found
    INTEGER :: i, i1, i2, info, j, j1, j2, k,k1, k2, l, m, n
    INTEGER :: ib, ic, ig
    INTEGER :: IF, if1, if2
    INTEGER :: iv, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8
    INTEGER, ALLOCATABLE :: perm(:), pinv(:)
    INTEGER, ALLOCATABLE :: aux(:,:), buf(:)
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

    PRINT *, "Element block data"
    DO i = 1, num_elem_blks
        CALL exgelb(exodus_file_id, elem_blk_ids(i), elem_blk_typ(i), &
            num_elems_in_blk(i), num_nodes_per_elem(i), &
            num_attr_per_elem(i), ierr)
        PRINT *, i, elem_blk_ids(i), elem_blk_typ(i), num_elems_in_blk(i), &
            num_nodes_per_elem(i), num_attr_per_elem(i)
    END DO


    ALLOCATE(conn(MAXVAL(num_nodes_per_elem), MAXVAL(num_elems_in_blk)))
    PRINT *, MAXVAL(num_nodes_per_elem), MAXVAL(num_elems_in_blk), SIZE(conn), num_elems_in_blk

    ! Read element block conn
    DO i = 1,num_elem_blks
        conn = -1
        CALL exgelc(exodus_file_id, elem_blk_ids(i), conn, ierr)
        PRINT *, "****************"
        PRINT *, i, elem_blk_ids(i), SIZE(conn(1:num_nodes_per_elem(i), 1:num_elems_in_blk(i)))
    ENDDO

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
    !call alloc_table(f2c,nel=ncells,ntab=nfaces)
    CALL v2f%alloc_table(nel=nfaces,ntab=n)

    ! Allocation of temporary face-related arrays
    ALLOCATE(iflag(nfaces), fnv(nfaces),  &
        &   imaster(nfaces), islave(nfaces), stat=info)
    IF(info /= 0) THEN
        WRITE(*,100)
        CALL abort_psblas
    END IF

    i1 = 1; i2 = 0; f2c%lookup(1) = 1
    j1 = 1; j2 = 0; v2f%lookup(1) = 1
    DO ic = 1, ncells
        i2 = i1 + cellnf(ic) - 1
        !
        iv = v2c%lookup(ic) - 1
        iv1 = v2c%tab(iv+1)
        iv2 = v2c%tab(iv+2)
        iv3 = v2c%tab(iv+3)
        IF(cellnv(ic) > 3) iv4 = v2c%tab(iv+4) ! quad, tet
        IF(cellnv(ic) > 4) iv5 = v2c%tab(iv+5) ! pyr
        IF(cellnv(ic) > 5) iv6 = v2c%tab(iv+6) ! pri
        IF(cellnv(ic) > 7) THEN ! hex
            iv7 = v2c%tab(iv+7)
            iv8 = v2c%tab(iv+8)
        END IF
        DO IF = i1, i2
            f2c%tab(IF) = IF
            iflag(IF) = 0
            imaster(IF) = ic
            islave(IF) = 0 ! Further defined
        END DO

        ! NOTE: standard connectivity vertices-cells defined in SIDS_guide
        SELECT CASE(cellgeo(ic))
        CASE('qua')
            j2 = j1 + 8 - 1
            v2f%tab(j1:j2) = (/ iv1,iv2, iv2,iv3, iv3,iv4, iv4,iv1 /)
            v2f%lookup(i1:i2) = (/ j1, (j1+2), (j1+4), (j1+6) /)
            fnv(i1:i2) = (/ 2, 2, 2, 2 /)
        CASE('tri')
            j2 = j1 + 6 - 1
            v2f%tab(j1:j2) = (/ iv1,iv2, iv2,iv3, iv3,iv1 /)
            v2f%lookup(i1:i2) = (/ j1, (j1+2), (j1+4) /)
            fnv(i1:i2) = (/ 2, 2, 2 /)
        CASE('hex')
            j2 = j1 + 24 - 1
            v2f%tab(j1:j2) = (/ iv1,iv4,iv3,iv2, iv1,iv2,iv6,iv5, iv2,iv3,iv7,iv6, &
                &              iv3,iv4,iv8,iv7, iv1,iv5,iv8,iv4, iv5,iv6,iv7,iv8 /)
            v2f%lookup(i1:i2) = (/ j1, (j1+4), (j1+8), (j1+12), (j1+16), (j1+20) /)
            fnv(i1:i2) = (/ 4, 4, 4, 4, 4, 4 /)
        CASE('pri')
            j2 = j1 + 18 - 1
            v2f%tab(j1:j2) = (/ iv1,iv2,iv5,iv4, iv2,iv3,iv6,iv5, iv3,iv1,iv4,iv6, &
                &              iv1,iv3,iv2, iv4,iv5,iv6 /)
            v2f%lookup(i1:i2) = (/ j1, (j1+4), (j1+8), (j1+12), (j1+15) /)
            fnv(i1:i2) = (/ 4, 4, 4, 3, 3  /)
        CASE('tet')
            j2 = j1 + 12 - 1
            v2f%tab(j1:j2) = (/ iv1,iv3,iv2, iv1,iv2,iv4, iv2,iv3,iv4, iv3,iv1,iv4 /)
            v2f%lookup(i1:i2) = (/ j1, (j1+3), (j1+6), (j1+9) /)
            fnv(i1:i2) = (/ 3, 3, 3, 3  /)
        CASE('pyr')
            j2 = j1 + 16 - 1
            v2f%tab(j1:j2) = (/ iv1,iv4,iv3,iv2, iv1,iv2,iv5, iv2,iv3,iv5, &
                &              iv3,iv4,iv5, iv4,iv1,iv5 /)
            v2f%lookup(i1:i2) = (/ j1, (j1+4), (j1+7), (j1+10), (j1+13) /)
            fnv(i1:i2) = (/ 4, 3, 3, 3, 3  /)
        END SELECT
        ! ---
        i1 = i2 + 1
        j1 = j2 + 1
        f2c%lookup(ic+1) = i1
        v2f%lookup(i2+1) = j1
    END DO

    ! #############################################################################
    IF(debug) THEN
        DO ic = 1, ncells
            i1 = f2c%lookup(ic)
            i2 = f2c%lookup(ic+1) - 1
            j1 = f2c%tab(i1)
            j2 = f2c%tab(i2)
            DO IF = j1, j2
                iv1 = v2f%lookup(IF)
                iv2 = v2f%lookup(IF+1) - 1
                WRITE(*,320) ic, cellgeo(ic), IF, v2f%tab(iv1:iv2);
            END DO
        END DO
        WRITE(*,*)
320     FORMAT(' cell: ',i7,' type: ',a3, ' face: ',i7,' vertices:', 8(1x,i5))
    END IF
    ! ###########################################################################


    ! Builds connection table f2v, dual of v2f

    CALL v2f%get_dual_table(f2v)

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
                i1 = v2f%lookup(if1) - 1
                iv1 = v2f%tab(i1+1) ! 1st iv of if1
                iv2 = v2f%tab(i1+2) ! 2nd iv of if1
                if2_2D: DO k = j + 1, n
                    if2 = f2v%tab(i+k)
                    i2 = v2f%lookup(if2)-1
                    iv3 = v2f%tab(i2+1) ! 1st iv of if2
                    iv4 = v2f%tab(i2+2) ! 2nd iv of if2
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
                i1 = v2f%lookup(if1) - 1
                iv1 = v2f%tab(i1+1)
                iv2 = v2f%tab(i1+2)
                iv3 = v2f%tab(i1+3)
                found = .FALSE.
                if2_3D: DO k = j + 1, n
                    if2 = f2v%tab(i+k)
                    i2 = v2f%lookup(if2)
                    j2 = v2f%lookup(if2+1)-1
                    l = fnv(if2)
                    ! Copies if2 iv-sequence in buf
                    buf(1) = v2f%tab(j2)
                    buf(2:(l+1)) = v2f%tab(i2:j2)
                    buf(l+2) = v2f%tab(i2)
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
            WRITE(*,330) IF, iflag(IF), imaster(IF), islave(IF)
        END DO
        WRITE(*,*)
330     FORMAT(' face: ',i4,' flag: ',i2,' master: ',i4,' slave: ',i4)
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
            i1 = v2f%lookup(IF)
            k = fnv(IF)
            WRITE(*,340) IF, iflag(IF), perm(IF), (v2f%tab(j), j = i1, i1+k-1)
        END DO
        WRITE(*,*)
        DO IF = 1,nfaces
            i1 = v2f%lookup(IF)
            k = fnv(IF)
        END DO
340     FORMAT(' face: ',i5,' flag: ',i2,' perm: ',i5,' vertices',4i5)
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
        j1 = v2f%lookup(j)
        j2 = j1 + facenv(IF) - 1
        dmy%tab(i1:i2) = v2f%tab(j1:j2)
        i1 = i2 + 1
        dmy%lookup(IF+1) = i1
    END DO
    CALL v2f%free_table()

    CALL v2f%alloc_table(nel=nfaces,ntab=k)
    v2f%lookup = dmy%lookup(1:nfaces+1)
    v2f%tab = dmy%tab(1:k)

    ! Reordering of f2c%tab
    DO ic = 1, ncells
        i1 = f2c%lookup(ic)
        i2 = f2c%lookup(ic+1) - 1
        DO i = i1, i2
            j = f2c%tab(i)
            f2c%tab(i) = perm(j)
        END DO
    END DO

    CALL dmy%free_table()
    DEALLOCATE(perm,pinv)

    ! ###########################################################################
    IF(debug) THEN
        DO IF = 1, nfaces
            i1 = v2f%lookup(IF)
            i2 = v2f%lookup(IF+1) - 1
            WRITE(*,350) IF, facemaster(IF), faceslave(IF), v2f%tab(i1:i2)
        END DO
        WRITE(*,*)
        DO ic = 1, ncells
            i1 = f2c%lookup(ic)
            i2 = f2c%lookup(ic+1) - 1
            WRITE(*,360) ic, f2c%tab(i1:i2)
        END DO
        WRITE(*,*)
350     FORMAT(' face: ',i5,' master: ',i5,' slave: ',i5,' vertices: ',4i5)
360     FORMAT(' cell: ',i5, ' faces:',6(1x,i5))
    END IF
    ! ###########################################################################

    PRINT *, "Node Set Data"

    IF (num_node_sets > 0) THEN
        ALLOCATE(node_set_ids(num_node_sets), num_nodes_per_set(num_node_sets), &
            num_df_per_node_set(num_node_sets))

        CALL exgnsi(exodus_file_id, node_set_ids, ierr)

        DO i = 1, num_node_sets
            CALL exgnp(exodus_file_id, node_set_ids(i), num_nodes_per_set(i), &
                num_df_per_node_set(i), ierr)
        END DO

        ALLOCATE(node_list(num_node_sets, MAXVAL(num_nodes_per_set)))
        DO i = 1, num_node_sets
            CALL exgns(exodus_file_id, node_set_ids(i), node_list(i,:), ierr)
            ! call exgnsd(exodus_file_id, node_set_ids(i), dist_fact(i,:), ierr)
            ! print *, i, node_set_ids(i), num_nodes_per_set(i), node_list(i,1:num_nodes_per_set(i))
        ENDDO

    END IF

    PRINT *, "Side Set Data"

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
            ! print *, i, side_set_ids(i), num_sides_per_set(i), elem_list(i,1:num_sides_per_set(i))
            !print *, i, side_set_ids(i), num_sides_per_set(i), side_list(i,1:num_sides_per_set(i))
        ENDDO

    END IF

    CALL exclos(exodus_file_id, ierr)

100 FORMAT(' ERROR! Memory allocation failure in RD_EXODUS_MESH')

END SUBROUTINE rd_exodus_mesh
