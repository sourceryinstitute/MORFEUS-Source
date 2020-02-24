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
! $Id: rd_cgns_mesh.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE (tools_mesh) rd_cgns_mesh_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE rd_cgns_mesh
        !! @note Compatible only with CGNS v3.3.0 and above. Earlier versions are not supported.
        !! Note: CGNS functionality is not working for Morfeus. Currently set to pre-process entire module procedure
        !!       so that the implementation exists but is empty. IP - 6/6/2019
        USE cgns, ONLY : cg_vertex       => vertex, & ! avoids "vertex" name conflict with class_vertex.f90
            & CG_Unstructured => Unstructured, &
            & CG_BAR_2        => BAR_2, &
            & CG_TRI_3        => TRI_3, &
            & CG_QUAD_4       => QUAD_4, &
            & CG_TETRA_4      => TETRA_4, &
            & CG_PYRA_5       => PYRA_5, &
            & CG_PENTA_6      => PENTA_6, &
            & CG_HEXA_8       => HEXA_8, &
            & CG_MODE_READ, &
            & cgenum_t, &
            & cgsize_t, &
            & CG_OK, &
            & cg_nbases_f, &
            & cg_nzones_f, &
            & cg_nbocos_f, &
            & cg_error_exit_f, &
            & cg_nsections_f, &
            & cg_close_f
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_face
        USE class_vertex
        USE tools_math
        USE type_table, ONLY : table, alloc_table, free_table, get_dual_table
        IMPLICIT NONE
        !
        INTEGER, PARAMETER :: nlen = 80
        !
        ! --- Local variables ---
        !
        ! Parameters
        LOGICAL, PARAMETER :: debug = .FALSE.
        INTEGER, PARAMETER :: mesh = 10
        INTEGER, PARAMETER :: MAX_LEN = 32
        !
        ! CGNS variables
        INTEGER :: index_file, index_base, index_zone, index_sect, ier
        INTEGER :: cell_dim, phys_dim, nbocos, nsections
        INTEGER :: itype, istart, iend, nbndry, iparent_flag, iparentdata
        INTEGER :: bocotype, ptset_type
        INTEGER :: NormalIndex, NormalListFlag, NormalDataType, NormalList, ndataset
        INTEGER :: isect_cell, isect_bc
        CHARACTER(len=32) :: basename, zonename, sectionname
        INTEGER(cgenum_t) :: type_coordinate_data
        INTEGER(cgsize_t) :: isize(3)
        CHARACTER(len=MAX_LEN) :: name_coordinate
        INTEGER :: type_file
        !
        ! Local copies of V2F_, V2C_, F2C_ objects (connectivity class)
        TYPE(table) :: v2f, v2c, f2c
        !
        ! Vertex-related variables
        INTEGER :: nverts
        INTEGER(cgsize_t) :: vertsmax !! number of vertices
        INTEGER(cgsize_t) :: vertsmin=1
        LOGICAL, ALLOCATABLE :: on_boundary(:)
        REAL(psb_dpk_), ALLOCATABLE :: xv(:), yv(:), zv(:)
        !
        ! Face-related variables
        INTEGER :: nfaces
        INTEGER, ALLOCATABLE :: fnv(:), imaster(:), islave(:), iflag(:)
        INTEGER, ALLOCATABLE :: facenv(:), facemaster(:), faceslave(:), faceflag(:)
        !
        ! Cell-related variables
        INTEGER :: ncells
        INTEGER, ALLOCATABLE :: cellnv(:), cellnf(:), cellgroup(:)
        CHARACTER(len=3), ALLOCATABLE :: cellgeo(:)
        !
        ! Group-related variables
        INTEGER :: ngroups
        INTEGER, ALLOCATABLE :: groupnc(:), groupmat(:)
        CHARACTER(len=32), ALLOCATABLE :: groupname(:)
        TYPE(table) :: c2g
        !
        ! BC-related variables
        INTEGER :: nbcfaces
        INTEGER, ALLOCATABLE :: bcnf(:)
        INTEGER(cgsize_t), ALLOCATABLE :: bcnf_cg(:)
        CHARACTER(len=32), ALLOCATABLE :: bcname(:)
        TYPE(table) :: f2b, bcface
        !
        ! Work variables
        LOGICAL :: found
        INTEGER :: i, i1, i2, info, j, j1, j2, k,k1, k2, l, m, n
        INTEGER :: ib, ic, ig
        INTEGER :: IF, if1, if2
        INTEGER :: iv, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8
        INTEGER, ALLOCATABLE :: perm(:), pinv(:)
        INTEGER, ALLOCATABLE :: aux(:,:), buf(:), work1(:), work2(:)
        INTEGER(cgsize_t), ALLOCATABLE :: work1_cg(:), f2btab_cg(:), bcfacetab_cg(:)
        TYPE(table) :: dmy, f2v

        ! ------------------------------------------------------------------

        WRITE(*,*) 'Importing mesh in .CGNS format on P0'

        ! Check file type
        CALL cg_is_cgns_f(TRIM(ADJUSTL(mesh_file)),type_file,ier)
        IF(ier /= CG_OK) THEN
            WRITE (*,*) 'ERROR! Unsupported file format'
            CALL abort_psblas
        END IF

        ! Opens file *.cgns
        CALL cg_open_f(TRIM(ADJUSTL(mesh_file)),CG_MODE_READ,index_file,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f

        ! Finds number of bases and zones
        ! NOTE: Only one base and one zone per base are assumed respectively
        CALL cg_nbases_f(index_file,index_base,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f

        IF (index_base /= 1) THEN
            WRITE (*,*) 'ERROR! Unsupported number of CGNS bases'
            CALL abort_psblas
        END IF
        CALL cg_nzones_f(index_file,index_base,index_zone,ier)
        IF (index_zone /= 1) THEN
            WRITE (*,*) 'ERROR! Unsupported number of CGNS zones'
            CALL abort_psblas
        END IF

        ! mesh ID and global parameters
        CALL cg_base_read_f(index_file,index_base,basename,cell_dim,phys_dim,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f
        CALL cg_zone_read_f(index_file,index_base,index_zone,zonename,isize,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f
        CALL cg_nbocos_f(index_file,index_base,index_zone,nbocos,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f
        mesh_id = zonename
        nverts = isize(1)
        ncells = isize(2)
        vertsmax = nverts

        ! By default isize(3) is equal to 0
        ngroups = 1 ! TO DO: how to manage multiple groups?
        nbc = nbocos
        ncd = cell_dim
        IF(ncd /= 2.AND.ncd /= 3) THEN
            WRITE(*,*) 'ERROR! NCD in CGNS file must be 2 or 3'
            CALL abort_psblas
        END IF

        ! #############################################################################
        IF(debug) THEN
            WRITE(*,*)
            WRITE(*,300) TRIM(mesh_id), nverts, ncells, ngroups, nbc, ncd
            WRITE(*,*)
300         FORMAT(' meshid: ',a,' nverts: ',i6,' ncells: ',i6,' ngroups ',i2,&
                & ' nbc: ',i2,' ncd: ',i2)
        END IF
        ! #############################################################################

        ALLOCATE (xv(nverts), yv(nverts), zv(nverts), &
            & cellnv(ncells), cellnf(ncells), cellgeo(ncells), cellgroup(ncells), &
            & stat=info)
        IF(info /= 0) THEN
            WRITE (*,100)
            CALL abort_psblas
        END IF

        ! Check coordinate data type
        CALL cg_coord_info_f(index_file, index_base, index_zone, 1, type_coordinate_data, name_coordinate, ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f()

        ! Reads vertices coordinates
        CALL cg_coord_read_f(index_file,index_base,index_zone,'CoordinateX',&
            &  type_coordinate_data,vertsmin,vertsmax,xv,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f
        CALL cg_coord_read_f(index_file,index_base,index_zone,'CoordinateY',&
            & type_coordinate_data,vertsmin,vertsmax,yv,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f
        CALL cg_coord_read_f(index_file,index_base,index_zone,'CoordinateZ',&
            & type_coordinate_data,vertsmin,vertsmax,zv,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f

        ! #############################################################################
        IF(debug) THEN
            DO iv = 1, nverts
                WRITE(*,310) iv, xv(iv), yv(iv), zv(iv)
            END DO
            WRITE(*,*)
310         FORMAT(' vertex: ',i5,' | ','coordinates: ',3(1x,d13.6))
        END IF
        ! #############################################################################

        ! Finds number of sections in CGNS file
        CALL cg_nsections_f(index_file,index_base,index_zone,nsections,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f

        ! k = number of sections containing cell elements
        DO index_sect = 1, nsections
            CALL cg_section_read_f(index_file,index_base,index_zone,index_sect,&
                & sectionname,itype,istart,iend,nbndry,iparent_flag,ier)
            IF(ier /= CG_OK) CALL cg_error_exit_f
            IF (iend == ncells) THEN
                isect_cell = index_sect
                EXIT
            END IF
        END DO
        isect_bc=nsections

        ! Reads sections containing vertex-to-cell connectivities
        IF(ncd == 2) THEN
            n = 4 * ncells
        ELSE IF(ncd == 3) THEN
            n = 8 * ncells
        END IF
        ALLOCATE(work1(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ALLOCATE(work1_cg(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL v2c%alloc_table(nel=ncells)

        work1 = 0
        work1_cg = 0
        i1 = 1; i2 = 0; v2c%lookup(i1) = 1
        DO index_sect = 1, isect_cell
            CALL cg_section_read_f(index_file,index_base,index_zone,index_sect,&
                & sectionname,itype,istart,iend,nbndry,iparent_flag,ier)
            IF(ier /= CG_OK) CALL cg_error_exit_f
            DO ic = istart, iend
                SELECT CASE(itype)
                CASE(5)
                    cellnv(ic) = 3
                    cellnf(ic) = 3
                    cellgroup(ic) = -1
                    cellgeo(ic) ='tri'
                CASE(7)
                    cellnv(ic) = 4
                    cellnf(ic) = 4
                    cellgroup(ic) = -1
                    cellgeo(ic) = 'qua'
                CASE(10)
                    cellnv(ic) = 4
                    cellnf(ic) = 4
                    cellgroup(ic) = -1
                    cellgeo(ic) = 'tet'
                CASE(12)
                    cellnv(ic) = 5
                    cellnf(ic) = 5
                    cellgroup(ic) = -1
                    cellgeo(ic) = 'pyr'
                CASE(14)
                    cellnv(ic) = 6
                    cellnf(ic) = 5
                    cellgroup(ic) = -1
                    cellgeo(ic) = 'pri'
                CASE(17)
                    cellnv(ic) = 8
                    cellnf(ic) = 6
                    cellgroup(ic) = -1
                    cellgeo(ic) = 'hex'
                CASE default
                    WRITE(*,*) 'WARNING! Unsupported CGNS type of cell: ', itype, ' IGNORING'
                    CONTINUE
                END SELECT
                i2 = i1 + cellnv(ic) - 1
                i1 = i2 + 1
                v2c%lookup(ic+1) = i1
            END DO
            j1 = v2c%lookup(istart)
            j2 = v2c%lookup(iend+1)-1
            CALL cg_elements_read_f(index_file,index_base,index_zone,index_sect,&
                & work1_cg(j1:j2),iparentdata,ier)
            IF(ier /= CG_OK) CALL cg_error_exit_f
        END DO
        CALL v2c%alloc_table(ntab=i2)

        work1 = work1_cg
        v2c%tab = work1(1:i2)
        DEALLOCATE(work1)
        DEALLOCATE(work1_cg)

        ! Defines group
        ! NOTE: only one group is currently supported
        CALL c2g%alloc_table(nel=ngroups,ntab=ncells)

        ALLOCATE(groupnc(ngroups),groupmat(ngroups),&
            &   groupname(ngroups),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        i1 = 1; i2 = 0; c2g%lookup(1) = 1
        DO ig = 1, ngroups
            groupnc(ig) = ncells
            groupmat(ig) = 2        ! Default of Gambit
            groupname(ig) = 'fluid' ! Default of Gambit
            i2 = i1 + groupnc(ig) - 1
            DO ic = i1, i2
                c2g%tab(ic) = ic
            END DO
            i1 = i2 + 1
            c2g%lookup(ig+1) = i1
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
        CALL f2c%alloc_table(nel=ncells,ntab=nfaces)
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
320         FORMAT(' cell: ',i5,' type: ',a3, ' face: ',i6,' vertices:', 8(1x,i5))
        END IF
        ! ###########################################################################

        ! Builds connection table f2v, dual of v2f
        CALL get_dual_table(v2f,f2v)

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
330         FORMAT(' face: ',i4,' flag: ',i2,' master: ',i4,' slave: ',i4)
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
340         FORMAT(' face: ',i5,' flag: ',i2,' perm: ',i5,' vertices',4i5)
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
350         FORMAT(' face: ',i5,' master: ',i5,' slave: ',i5,' vertices: ',4i5)
360         FORMAT(' cell: ',i5, ' faces:',6(1x,i5))
        END IF
        ! ###########################################################################

        ! --- Reads and assigns boundary conditions. Boundary flags are integer
        ! in the range [1-NBC]. If NBC=0 from import, automatically NBC=1 and where
        ! faceslave=0 a boundary face is asssumed.
        ! --- Builds list of faces, grouping them by their flag value.
        !
        ! iflag(if)=0 --> fluid face
        ! iflag(if)=i --> i-th boundary face
        !

        ! Counts boundary faces (faceslave = 0)
        k = COUNT(faceslave == 0)

        IF(nbc /= 0) THEN
            CALL f2b%alloc_table(nel=nbc,ntab=k)
            ALLOCATE(bcnf(nbc),bcname(nbc),stat=info)
        ELSE
            CALL f2b%alloc_table(nel=1,ntab=k)
            ALLOCATE(bcnf(1),bcname(1),stat=info)
        END IF
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        ALLOCATE(f2btab_cg(k), stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        ALLOCATE(bcnf_cg(SIZE(bcnf)), stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Reads bc section
        IF(nbc /= 0) THEN
            i1 = 1; i2 = 0; f2b%lookup(1) = i1
            DO ib = 1, nbc
                CALL cg_boco_info_f(index_file,index_base,index_zone,ib,&
                    & bcname(ib),bocotype,ptset_type,bcnf_cg(ib),&
                    & NormalIndex,NormalListFlag,NormalDataType,ndataset,ier)
                IF(ier /= CG_OK) CALL cg_error_exit_f
                ! NOTE: NORMAL* and NDATASET are unused
                bcnf(ib) = bcnf_cg(ib)

                i2 = i1 + bcnf(ib) - 1
                CALL cg_boco_read_f(index_file,index_base,index_zone,ib,&
                    & f2btab_cg(i1:i2),NormalList,ier)
                IF(ier /= CG_OK) CALL cg_error_exit_f
                f2b%tab(i1:i2) = f2btab_cg(i1:i2)
                i1 = i2 + 1
                f2b%lookup(ib+1) = i1
            END DO
            ! N = number of BC faces
            n = f2b%lookup(nbc+1) - 1

            ! Check on the consistency of bc faces
            IF(n /= COUNT(faceslave == 0)) THEN
                WRITE(*,*) 'ERROR! Inconsistent number of bc faces in CGNS file'
                CALL abort_psblas
            END IF

            !
            ! Counts records and elements of bc-sections; to use in the successive
            ! allocation of bcface%lookup, bcface%tab.
            ! NOTE: this is necessary because it can happen that the total amount of records in
            ! presumed 'bc'-section is greater than the effective number of bc faces.
            ! --> to be investigated more deeply in ICEM!
            k1 = 0; k2 = 0
            DO index_sect = isect_cell + 1, isect_bc
                CALL cg_section_read_f(index_file,index_base,index_zone,index_sect,&
                    & sectionname,itype,istart,iend,nbndry,iparent_flag,ier)
                IF(ier /= CG_OK) CALL cg_error_exit_f
                n = iend - istart + 1 ! Number of records, i.e. faces, inside the bc-section
                k1 = k1 + n           ! Increases bcface%lookup elements number
                SELECT CASE(itype)
                CASE(3) ! bar
                    l=2
                CASE(5) ! tri
                    l=3
                CASE(7) ! quad
                    l=4
                CASE default
                    WRITE(*,*) 'WARNING! Unsupported type of CGNS bc element: itype =', itype, 'IGNORING'
                    CONTINUE
                END SELECT
                k2=k2+n*l ! Increases bcface%tab elements number
            END DO

            ! Allocates bcface%lookup, bcface%tab
            CALL bcface%alloc_table(nel=k1,ntab=k2)
            ALLOCATE(bcfacetab_cg(k2), stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            ! Reads sections containing boundary elements
            bcface%lookup(1) = 1
            DO index_sect = isect_cell + 1, isect_bc
                CALL cg_section_read_f(index_file,index_base,index_zone,index_sect,&
                    & sectionname,itype,istart,iend,nbndry,iparent_flag,ier)
                IF(ier /= CG_OK) CALL cg_error_exit_f
                SELECT CASE(itype)
                CASE(3) ! bar
                    l = 2
                CASE(5) ! tri
                    l = 3
                CASE(7) ! quad
                    l = 4
                CASE default
                    WRITE(*,*) 'WARNING! Unsupported type of CGNS bc element: itype =', itype, 'IGNORING'
                    CONTINUE
                END SELECT
                i1 = istart - ncells
                i2 = iend - ncells
                DO i = i1, i2
                    bcface%lookup(i+1) = bcface%lookup(i) + l
                END DO
                j1 = bcface%lookup(i1)
                j2 = bcface%lookup(i2+1) - 1
                CALL cg_elements_read_f(index_file,index_base,index_zone,index_sect,&
                    & bcfacetab_cg(j1:j2),iparentdata,ier)
                IF(ier /= CG_OK) CALL cg_error_exit_f
                bcface%tab(j1:j2) = bcfacetab_cg(j1:j2)
            END DO

            DEALLOCATE(bcfacetab_cg)

            ! Builds connection table f2v, dual of v2f
            CALL get_dual_table(v2f,f2v)

            ! Associates bc elements to corresponding faces
            IF(ncd == 2) THEN
                ib_2D: DO ib = 1, nbc
                    i1 = f2b%lookup(ib)
                    i2 = f2b%lookup(ib+1) - 1
                    bcface_2D: DO i = i1, i2
                        if1 = f2b%tab(i) - ncells ! Offset of cell-sections
                        k = bcface%lookup(if1)
                        ! NOTE: ICEM don't respect on bc elements the right hand ordering!!
                        ! => for successive benchmark forces iv1 < iv2
                        iv1 = bcface%tab(k)
                        iv2 = bcface%tab(k+1)
                        IF(iv2 < iv1) THEN
                            m = iv1
                            iv1 = iv2
                            iv2 = m
                        END IF

                        j1 = f2v%lookup(iv1)
                        j2 = f2v%lookup(iv1+1) - 1
                        v2f_2D: DO j = j1, j2
                            if2 = f2v%tab(j)
                            IF(faceslave(if2) /= 0 .OR. &
                                & faceflag(if2) > 0) CYCLE
                            l = v2f%lookup(if2)
                            iv3 = v2f%tab(l)
                            iv4 = v2f%tab(l+1)
                            ! NOTE: ICEM don't respect on bc elements the right hand ordering!!
                            ! => for successive benchmark forces iv3 < iv4
                            IF(iv4 < iv3) THEN
                                m = iv4
                                iv4 = iv3
                                iv3 = m
                            END IF
                            IF(iv1 == iv3 .AND. iv2 == iv4 ) THEN
                                faceflag(if2) = ib
                                f2b%tab(i) = if2
                                EXIT v2f_2D
                            END IF
                        END DO v2f_2D
                    END DO bcface_2D
                END DO ib_2D
            ELSE IF(ncd == 3) THEN
                ALLOCATE(work1(4),work2(4),stat=info)
                IF(info /= 0) THEN
                    WRITE(*,100)
                    CALL abort_psblas
                END IF
                ib_3D: DO ib = 1, nbc
                    i1 = f2b%lookup(ib)
                    i2 = f2b%lookup(ib+1) - 1
                    bcface_3D: DO i = i1, i2
                        if1 = f2b%tab(i) - ncells ! Offset of cell-sections
                        ! NOTE: ICEM don't respect on bc elements the right hand ordering!!
                        ! Vertices sequences of if1 and if2 faces are copied respectively in
                        ! arrrays work1 and work2 then ordered.
                        j1 = bcface%lookup(if1)
                        j2 = bcface%lookup(if1+1) - 1
                        l = 0
                        DO j = j1, j2
                            l = l + 1
                            work1(l) = bcface%tab(j)
                        END DO
                        IF(l == 3) work1(4) = 0
                        CALL sort(work1(1:l))

                        iv1 = work1(1)
                        j1 = f2v%lookup(iv1)
                        j2 = f2v%lookup(iv1+1) - 1
                        ! NOTE: l=number of vertex of if1. DO NOT CHANGE during next inner loop!!

                        v2f_3D: DO j = j1, j2
                            if2 = f2v%tab(j)
                            IF(faceslave(if2) /= 0   .OR.&
                                & faceflag(if2) > 0 .OR.&
                                & facenv(if2) /= l) CYCLE
                            k1 = v2f%lookup(if2)
                            k2 = v2f%lookup(if2+1) - 1
                            m = 0
                            DO k = k1, k2
                                m = m + 1
                                work2(m) = v2f%tab(k)
                            END DO
                            CALL sort(work2(1:m))
                            IF(    work1(1) == work2(1) .AND. &
                                & work1(2) == work2(2) .AND. &
                                & work1(3) == work2(3)) THEN
                                faceflag(if2) = ib
                                f2b%tab(i) = if2
                                EXIT v2f_3D
                            END IF
                        END DO v2f_3D
                    END DO bcface_3D
                END DO ib_3D
                DEALLOCATE(work1,work2)
            END IF
            CALL f2v%free_table()
            CALL bcface%free_table()
        ELSE
            ! NBC = 0
            nbc = 1
            bcname(1) = 'bc1'; bcnf(1) = COUNT(faceslave == 0)
            f2b%lookup(1) = k + 1
            DO IF = 1, nfaces
                IF(faceslave(IF) == 0) THEN
                    k = k + 1
                    f2b%tab(k) = IF
                END IF
            END DO
            f2b%lookup(2) = k + 1
            WHERE (faceslave == 0)
                faceflag = 1
            END WHERE
        END IF
        DEALLOCATE(f2btab_cg)

        ! Reordering of if-indexing following the flag sequence: 0,1,2,...,nbc
        ! f2b%tab is used for building the inverse permutation array pinv
        ! old = pinv(new)
        ! new = perm(old)

        k = SIZE(v2f%tab)
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
            facenv(IF)     = aux(1,j)
            faceflag(IF)   = aux(2,j)
            facemaster(IF) = aux(3,j)
            faceslave(IF)  = aux(4,j)
        END DO

        ! Reordering of v2f%tab and v2f%lookup
        k  = SIZE(v2f%tab)
        dmy%lookup  = v2f%lookup(1:nfaces+1)
        dmy%tab = v2f%tab(1:k)
        i1 = 1; i2 = 0; v2f%tab(i1) = 1
        DO IF = 1, nfaces
            i2 = i1 + facenv(IF) - 1

            j = pinv(IF)
            j1 = dmy%lookup(j)
            j2 = j1 + facenv(IF) - 1

            v2f%tab(i1:i2) = dmy%tab(j1:j2)

            i1 = i2 + 1
            v2f%lookup(IF+1) = i1
        END DO

        ! Reordering of f2c%tab
        k = SIZE(f2c%tab)
        DO i = 1, k
            j = f2c%tab(i)
            f2c%tab(i) = perm(j)
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
            WRITE(*,370) (f2b%tab(i), i = 1, f2b%lookup(1) - 1)
            WRITE(*,*)

            DO ib = 1, nbc
                i1 = f2b%lookup(ib)
                i2 = i1 + bcnf(ib) - 1
                WRITE(*,380) ib, bcname(ib), f2b%tab(i1:i2)
            END DO
            WRITE(*,*)

            DO IF = 1, nfaces
                WRITE(*,390) IF, faceflag(IF), facemaster(IF), faceslave(IF)
            END DO
            WRITE(*,*)

            DO i = 0, nbc
                WRITE(*,400) i, COUNT(faceflag == i)
            END DO

            WRITE(*,*)
370         FORMAT(' Fluid faces: ',5i8:/(14x,5i8))
380         FORMAT(' Bcset: ',i2,' Name: ',a10,' List of faces: ', 5i8:/(43x,5i8:))
390         FORMAT(' faces: ',i5,' Flag: ',i2,' Master: ',i5,' Slave: ',i5)
400         FORMAT(' Total amount of faces with flag',i3,':',i8)
        END IF
        ! ###########################################################################

        ! f2b is no longer useful. It can be deallocated
        CALL f2b%free_table()
        CALL dmy%free_table()
        DEALLOCATE(bcname,bcnf)
        DEALLOCATE(bcnf_cg)
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
            i1 = v2f%lookup(IF)
            i2 = v2f%lookup(IF+1) - 1
            DO i = i1, i2
                iv = v2f%tab(i)
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
        CALL alloc_conn(v2f_,nel=nfaces,nconn=SIZE(v2f%tab))
        CALL alloc_conn(v2c_,nel=ncells,nconn=SIZE(v2c%tab))
        CALL alloc_conn(f2c_,nel=ncells,nconn=SIZE(f2c%tab))
        CALL alloc_conn(c2g_,nel=ngroups,nconn=SIZE(c2g%tab))

        DO IF = 1, nfaces
            i1 = v2f%lookup(IF)
            i2 = v2f%lookup(IF+1) - 1
            CALL v2f_%set_ith_conn(IF,v2f%tab(i1:i2))
        END DO

        DO ic = 1, ncells
            i1 = v2c%lookup(ic)
            i2 = v2c%lookup(ic+1) - 1
            CALL v2c_%set_ith_conn(ic,v2c%tab(i1:i2))
        END DO

        DO ic = 1, ncells
            i1 = f2c%lookup(ic)
            i2 = f2c%lookup(ic+1) - 1
            CALL f2c_%set_ith_conn(ic,f2c%tab(i1:i2))
        END DO

        DO ig = 1, ngroups
            i1 = c2g%lookup(ig)
            i2 = c2g%lookup(ig+1) - 1
            CALL c2g_%set_ith_conn(ig,c2g%tab(i1:i2))
        END DO

        ! Deallocates local copies of connectivity data
        CALL v2f%free_table()
        CALL v2c%free_table()
        CALL f2c%free_table()
        CALL c2g%free_table()

        ! Currently group data GROUPNC, GROUPMAT, GROUPNAME are not exported
        ! to the calling program => deallocation
        DEALLOCATE(groupnc,groupmat,groupname)

        CALL cg_close_f(index_file,ier)
        IF(ier /= CG_OK) CALL cg_error_exit_f

        WRITE(*,'()')

100     FORMAT(' ERROR! Memory allocation failure in RD_CGNS_MESH')

        END PROCEDURE rd_cgns_mesh

END SUBMODULE rd_cgns_mesh_implementation
