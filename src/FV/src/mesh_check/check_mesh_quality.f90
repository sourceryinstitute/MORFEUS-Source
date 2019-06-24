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
! $Id: check_mesh_quality.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    Given a mesh and a quality threshold, returns a scalar field containing the
!    quality of local cells (including halos) and the list of the local indices
!    of the strictly local bad cells (quality(ic) < tol).
!
! Provides:
!    QUALITY     (local + halos)   cell quality. Range [0.d0:1.d0]
!    BAD_CELLS   (strictly local)  local indices of bad cells.
!
SUBMODULE (tools_mesh_check) check_mesh_quality_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE check_mesh_quality

    USE class_psblas
    USE class_cell
    USE class_connectivity
    USE class_mesh
    USE class_scalar_field
    USE class_vector
    USE class_vertex
    USE tools_math
    USE tools_mesh_basics
    USE tools_mesh_check, ONLY : check_tet_quality

    IMPLICIT NONE
    !
    LOGICAL, PARAMETER :: debug =.FALSE.
    CHARACTER(len=32), PARAMETER :: WHERE = 'check_mesh_quality'
    !
    LOGICAL :: quiet_
    INTEGER :: err_act, info, icontxt, mypnum
    INTEGER :: nbad, nbad_glob, ncells, ndmy, nfaces, nloc, ntet, nhex
    INTEGER :: i, j, k, ic, ic_glob, IF, itet, ihex
    INTEGER :: firsthex, lasthex
    INTEGER :: iv1, iv2, iv3, iv4, iface,jface
    INTEGER, ALLOCATABLE :: ictype(:) , nctype(:)
    INTEGER, POINTER :: if2c(:) => NULL()
    INTEGER, POINTER :: iv2f(:) => NULL(),jv2f(:) => NULL()
    INTEGER, POINTER :: iv2c(:) => NULL()
    INTEGER, ALLOCATABLE :: iloc_to_glob(:)
    REAL(psb_dpk_) :: max_angle, min_angle, largest, smallest, min_q
    REAL(psb_dpk_), ALLOCATABLE :: quality(:)
    TYPE(vertex)     :: hex_verts(8)
    TYPE(vector)     :: tet_af(4), hex_af(6)
    INTEGER          :: edge, adjacent(12,2)  ! for noting which two faces in a hex share an edge
    LOGICAL          :: match

    ! Sets error handling for PSBLAS-2 routines
    CALL psb_erractionsave(err_act)

    icontxt = icontxt_()
    mypnum  = mypnum_()


    ! WARNING! 2D meshes are NOT supported: to be extended!
    IF(msh%ncd == 2) THEN
        IF(mypnum == 0) THEN
            WRITE(*,'()')
            WRITE(*,100)
            WRITE(*,'()')
        END IF

        fquality = 1.d0

        RETURN ! EXIT
    END IF


    ! ----- Preliminary Checks and Memory Allocation -----

    ncells = SIZE(msh%cells)
    nfaces = SIZE(msh%faces)

    ! Sets variable for data writeout.
    IF(PRESENT(quiet)) THEN
        quiet_ = quiet
    ELSE
        quiet_ = .FALSE. ! Default: write results out
    END IF


    IF (.NOT.quiet_ .AND. (mypnum == 0) ) WRITE(*,*) 'Checking mesh quality'


    IF ((PRESENT(tol).NEQV.PRESENT(bad_cells))) THEN  ! exclusive or
        WRITE(*,200)
        CALL abort_psblas
    ENDIF

    IF (PRESENT(tol)) THEN
        IF(tol < 0.d0 .OR. tol > 1.d0) THEN
            WRITE(*,300)
            CALL abort_psblas
        ENDIF
    END IF

    ALLOCATE(quality(ncells),stat=info)
    IF(info /= 0) THEN
        WRITE(*,400)
        CALL abort_psblas
    END IF

    ! Gets local to global list for cells
    CALL psb_get_loc_to_glob(msh%desc_c,iloc_to_glob)

    ! Sets intial value to maximum quality
    quality(:) = 1.d0

    CALL get_cells_type(msh%cells,nctype,ictype)
    ntet = nctype(itet_)
    nhex = nctype(ihex_)
    firsthex = ntet + nctype(ipyr_) + nctype(ipri_) + 1
    lasthex  = firsthex + nhex - 1

    !--------------------------------------------------------------------|
    ! WARNING! 3D meshes: only tet and hex cells are currently supported |
    !          To be extended to polyhedral cells!                       |
    !--------------------------------------------------------------------|

    ! Connectivity checks on tetrahedral cells.
    DO itet = 1, ntet
        ic = ictype(itet)

        ! Gets which faces belong to this cell
        CALL get_ith_conn(if2c,msh%f2c,ic)

        ! Checks that each tet has four faces.
        IF (SIZE(if2c) /= 4) THEN
            ic_glob = iloc_to_glob(ic)
            WRITE(*,500) ic_glob, SIZE(if2c), 4
            CALL abort_psblas
        ENDIF

        ! Checks that each face is within the bounds of the number of faces
        DO i = 1, 4
            IF = if2c(i)
            IF( IF < 1 .OR. IF > nfaces) THEN
                ic_glob = iloc_to_glob(ic)
                WRITE(*,600) ic, mypnum, ic_glob, IF, i
                CALL abort_psblas
            END IF
        END DO
    END DO

    ! Connectivity checks on hex cells.
    DO ihex = firsthex, lasthex
        ic = ictype(ihex)

        ! Gets which faces belong to this cell
        CALL get_ith_conn(if2c,msh%f2c,ic)

        ! Checks that each hex has six faces.
        IF (SIZE(if2c) /= 6) THEN
            ic_glob = iloc_to_glob(ic)
            WRITE(*,500) ic_glob, SIZE(if2c), 6
            CALL abort_psblas
        ENDIF

        ! Checks that each face is within the bounds of the number of faces
        DO i = 1, 6
            IF = if2c(i)
            IF( IF < 1 .OR. IF > nfaces) THEN
                ic_glob = iloc_to_glob(ic)
                WRITE(*,600) ic, mypnum, ic_glob, IF, i
                CALL abort_psblas
            END IF
        END DO
    END DO

    ! ----- Computes Mesh Quality For Tets and Hexes -----

    DO itet = 1, ntet
        ic = ictype(itet)

        CALL get_ith_conn(iv2c,msh%v2c,ic)

        iv1 = iv2c(1)
        iv2 = iv2c(2)
        iv3 = iv2c(3)
        iv4 = iv2c(4)

        quality(ic) = geom_tet_quality( &
            & msh%verts(iv1),msh%verts(iv2),msh%verts(iv3),msh%verts(iv4),&
            & msh%vol(ic))
    END DO

    DO ihex = firsthex, lasthex
        ic = ictype(ihex)

        CALL get_ith_conn(iv2c,msh%v2c,ic)

        DO i = 1, 8
            hex_verts(i) = msh%verts(iv2c(i))
        ENDDO

        !assumes verts are in CGNS order
        quality(ic) = geom_hex_quality( hex_verts, msh%vol(ic))
    END DO

    min_q = MINVAL(quality)
    CALL psb_amn(icontxt,min_q,info)
    CALL psb_check_error(info,TRIM(WHERE),'psb_amn',icontxt)


    ! ----- Counts local bad cells -----

    IF (PRESENT(bad_cells).AND. PRESENT(tol)) THEN
        nloc = psb_cd_get_local_rows(msh%desc_c)  ! number of strictly local cells
        nbad = COUNT(quality(1:nloc) <= tol)      ! number of strictly local bad cells
        ndmy = nbad                               ! makes a temporary copy of NBAD
        CALL psb_sum(icontxt,nbad,info)           ! sum NBAD over all processes
        CALL psb_check_error(info,TRIM(WHERE),'psb_sum',icontxt)
        nbad_glob = nbad
        nbad = ndmy

        IF ( ALLOCATED(bad_cells) ) DEALLOCATE(bad_cells)

        IF(nbad /= 0) THEN

            ALLOCATE(bad_cells(nbad),stat=info)
            IF(info /= 0) THEN
                WRITE(*,400)
                CALL abort_psblas
            END IF

            ! Groups bad cells indices
            ! WARNING! Only tets and hexes are currently supported.
            !          To be extended to polyhedral cells.
            k = 0
            DO itet = 1, ntet
                ic = ictype(itet)
                IF(    ic <= nloc .AND. &       ! checks whether a cell is strictly local
                    & quality(ic) <= tol) THEN ! checks the cell quality
                    k = k + 1
                    bad_cells(k) = ic
                END IF
            END DO

            IF(debug) THEN
                DO itet = 1, nbad
                    ic = bad_cells(itet)
                    CALL check_tet_quality(msh,ic,2)
                ENDDO
            END IF

            DO ihex = firsthex, lasthex
                ic = ictype(ihex)
                IF(    ic <= nloc .AND. &       ! checks whether a cell is strictly local
                    & quality(ic) <= tol) THEN ! checks the cell quality
                    k = k + 1
                    bad_cells(k) = ic
                END IF
            END DO

        END IF

    ENDIF


    quiet_check: IF (.NOT.quiet_) THEN

        ! ----- Dihedral Angle Check (Only 3D) -----

        ! Computes the largest and smallest dihedral angles
        largest = 0.d0
        smallest = pi

        ! tetrahedral cells
        DO itet = 1, ntet
            ic = ictype(itet)

            CALL get_ith_conn(if2c,msh%f2c,ic)

            DO i = 1,4
                tet_af(i) = msh%af(if2c(i))
            ENDDO

            CALL geom_tet_dihedral_angle( tet_af, max_angle, min_angle)

            largest  = MAX(largest,max_angle)
            smallest = MIN(smallest,min_angle)
        END DO

        ! hex cells
        DO ihex = firsthex, lasthex
            ic = ictype(ihex)

            CALL get_ith_conn(if2c,msh%f2c,ic)

            ! get face normals
            DO i = 1,6
                iface = if2c(i)
                hex_af(i) = msh%af(iface)
            ENDDO


            edge = 0
            adjacent(:,:) = 0
            ! find two faces who share no vertices;  these two faces are opposites

            DO i = 1,5       ! loop over a face
                iface = if2c(i)
                CALL get_ith_conn(iv2f, msh%v2f,iface)

                DO j = i+1,6  ! and all possible pairings
                    jface = if2c(j)
                    CALL get_ith_conn(jv2f, msh%v2f,jface)

                    match = .FALSE.  ! assume there are no shared vertices

                    DO k = 1,SIZE(jv2f) ! check each vertex from the 2nd face against
                        IF (  ANY(iv2f == jv2f(k)) ) THEN ! all the verts from the 1st face
                            match = .TRUE.                 ! and note matches
                        ENDIF

                        IF (match) THEN
                            edge = edge + 1
                            adjacent(edge,1) = i ! note that j is always > i
                            adjacent(edge,2) = j
                            EXIT  ! if there is a match, quit the k loop
                        ENDIF
                    ENDDO

                ENDDO

            ENDDO

            CALL geom_hex_dihedral_angle(hex_af, adjacent, max_angle, min_angle)

            largest  = MAX(largest,max_angle)
            smallest = MIN(smallest,min_angle)
        END DO


        ! Computes the global max. and min. values
        CALL psb_amx(icontxt,largest,info)
        CALL psb_check_error(info,TRIM(WHERE),'psb_amx',icontxt)

        CALL psb_amn(icontxt,smallest,info)
        CALL psb_check_error(info,TRIM(WHERE),'psb_amn',icontxt)

        IF(mypnum == 0) THEN
            WRITE(*,700) ' - The largest dihedral angle is:',180.d0 * largest  / pi
            WRITE(*,700) ' - The smallest dihedral angle is:',180.d0 * smallest / pi
            WRITE(*,700) ' - The lowest cell quality is:', min_q
            IF (PRESENT(bad_cells).AND. PRESENT(tol))  THEN
                WRITE(*,900) ' - Number of bad cells (quality < ',tol,'):', nbad_glob
            ENDIF
            WRITE(*,'()')
        END IF


    ENDIF quiet_check


    ! Sets scalar field
    fquality = quality


    ! Frees storage
    NULLIFY(iv2f,jv2f)
    NULLIFY(if2c,iv2c)
    DEALLOCATE(nctype,ictype)
    DEALLOCATE(iloc_to_glob)
    DEALLOCATE(quality)

    ! ----- Normal termination -----
    CALL psb_erractionrestore(err_act)


100 FORMAT(' WARNING! 2D meshes are not yet supported by CHECK_MESH_QUALITY')
200 FORMAT(' ERROR! Both TOL and BAD_CELLS or neither shall be present in', &
        & ' CHECK_MESH_QUALITY.')
300 FORMAT(' ERROR! Illegal value of TOL argument in CHECK_MESH_QUALITY')
400 FORMAT(' ERROR! Memory allocation failure in CHECK_MESH_QUALITY')
500 FORMAT(' ERROR! Cell ',i5,' has',i2,' faces where ',i2,' are expected.')
600 FORMAT(' ERROR! Illegal face indexing detected in CHECK_MESH_QUALITY.',/,  &
        &        '        Cell ',i5,' on proc ',i2,' (ic_glob=',i5,') references',/, &
        &        '        a face number ',i5,' as its ',i1,' face.')
700 FORMAT(a,f8.3)
900 FORMAT(a,f5.2,a,i5)

        END PROCEDURE check_mesh_quality

END SUBMODULE check_mesh_quality_implementation