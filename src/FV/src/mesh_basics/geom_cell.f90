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
! $Id: geom_cell.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    Evaluates coordinates of centers of mass and volumes of a cell set.
!    The computation is performed by means of a decomposition of each
!    cell in elemental triangles (2D) or tetraheadra (3D) for which exact
!    formulas are available.
!
MODULE PROCEDURE geom_cell
    USE class_psblas
    USE class_cell
    USE class_connectivity
    USE class_face
    USE class_vector
    USE class_vertex
    USE tools_mesh_basics, ONLY : geom_tet_center, geom_tet_volume

    IMPLICIT NONE
    !
    LOGICAL, PARAMETER :: debug = .FALSE.
    INTEGER :: i, ic, IF, info, is, k
    INTEGER :: iv, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8, iv_1, iv_i, iv_im1
    INTEGER :: nfaces, ncells
    LOGICAL :: quiet_
    !
    ! connectivity of I-th element
    INTEGER :: nv2c, nf2c
    INTEGER, POINTER :: iv2c(:) => NULL()
    INTEGER, POINTER :: if2c(:) => NULL()
    INTEGER, POINTER :: iv2f(:) => NULL()
    !
    INTEGER, ALLOCATABLE :: islave(:), nvc(:), nvf(:)
    REAL(psb_dpk_) :: r, w
    REAL(psb_dpk_), PARAMETER :: t3 = 1.d0 / 3.d0
    CHARACTER(len=3), ALLOCATABLE :: geoc(:)
    TYPE(vector) :: p, v, v1, v2

    IF (PRESENT(quiet) ) THEN
        quiet_ = quiet
    ELSE
        quiet_ = .FALSE.
    ENDIF

    ! Preliminary checks
    IF(    .NOT.ALLOCATED(verts) .OR. &
        & .NOT.ALLOCATED(faces) .OR. &
        & .NOT.ALLOCATED(cells)) THEN
        WRITE(*,100)
        CALL abort_psblas
    END IF
    ! Can do without since they are INTENT(OUT).
!!$  if(    allocated(cell_cntr) .or. &
!!$       & allocated(vol)) then
!!$    call abort_psblas
!!$  end if

    ncells = SIZE(cells)
    nfaces = SIZE(faces)

    CALL alloc_vector(cell_cntr,ncells)
    ALLOCATE(islave(nfaces),nvf(nfaces), &
        &   geoc(ncells),nvc(ncells),   &
        &   vol(ncells),stat=info)
    IF(info /= 0) THEN
        WRITE(*,200)
        CALL abort_psblas
    END IF

    ! Extracts slave cell index of every face
    islave(:) = slave_(faces)

    ! Extracts number of vertices of every face
    nvf(:) = nv_(faces)

    ! Extracts geometry type of every cell
    DO ic = 1, ncells
        geoc(ic) = geo_(cells(ic))
    END DO

    ! Extracts number of vertices of every cell
    nvc(:) = nv_(cells)

    ! Loop over cells
    vol(:) = 0.d0
    cell_cntr(:) = vector_(0.d0,0.d0,0.d0)

    IF(ncd == 2) THEN
        DO ic = 1, ncells

            ! Retrieves V2C connectivity of the cell IC.
            CALL get_ith_conn(iv2c,v2c,ic)

            iv_1 = iv2c(1)        ! 1st IV in V2C connectivity of the cell IC
            DO i = 3, nvc(ic)
                iv_im1 = iv2c(i-1) ! (I-1)-th IV in V2C connectivity of the cell IC
                v1 = verts(iv_im1) - verts(iv_1)

                iv_i = iv2c(i)     ! I-th IV in V2C connectivity of the cell IC
                v2 = verts(iv_i) - verts(iv_1)

                p = v1 .cross. v2
                p = 0.5d0 * p         ! Normal vector of the I-th triangle
                w = mag(p)            ! `Volume' of the I-th triangle
                vol(ic) = vol(ic) + w ! `Volume' of the generic polyside 2D cell

                p = verts(iv_1) + verts(iv_im1) + verts(iv_i)
                p = t3 * p                            ! Barycenter of the I-th triangle
                cell_cntr(ic) = cell_cntr(ic) + w * p ! Barycenter of the 2D polygonal cell
            END DO
            r = 1.d0 / vol(ic)
            cell_cntr(ic) = r * cell_cntr(ic)
        END DO

    ELSEIF(ncd == 3) THEN
        DO ic = 1, ncells

            ! Retrieves V2C connectivity of the cell IC.
            CALL get_ith_conn(iv2c,v2c,ic)

            iv1 = iv2c(1)
            iv2 = iv2c(2)
            iv3 = iv2c(3)
            iv4 = iv2c(4)
            SELECT CASE(geoc(ic))
            CASE('tet')
                cell_cntr(ic) = geom_tet_center(verts(iv1),verts(iv2),verts(iv3),verts(iv4))
                vol(ic)       = geom_tet_volume(verts(iv1),verts(iv2),verts(iv3),verts(iv4))
            CASE('pyr')
                iv5 = iv2c(5)

                ! 1st tetrahedron: 1, 2, 3, 5
                v = geom_tet_center(verts(iv1),verts(iv2),verts(iv3),verts(iv5))
                r = geom_tet_volume(verts(iv1),verts(iv2),verts(iv3),verts(iv5))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                ! 2nd tetrahedron: 1, 3, 4, 5
                v = geom_tet_center(verts(iv1),verts(iv3),verts(iv4),verts(iv5))
                r = geom_tet_volume(verts(iv1),verts(iv3),verts(iv4),verts(iv5))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) =cell_cntr(ic) + r * v

                r = 1.d0 / vol(ic)
                cell_cntr(ic) = r * cell_cntr(ic)
            CASE('pri')
                iv5 = iv2c(5)
                iv6 = iv2c(6)

                ! 1st tetrahedron: 1, 2, 3 ,4
                v = geom_tet_center(verts(iv1),verts(iv2),verts(iv3),verts(iv4))
                r = geom_tet_volume(verts(iv1),verts(iv2),verts(iv3),verts(iv4))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                ! 2nd tetrahedron: 4, 6, 5, 2
                v = geom_tet_center(verts(iv4),verts(iv6),verts(iv5),verts(iv2))
                r = geom_tet_volume(verts(iv4),verts(iv6),verts(iv5),verts(iv2))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                ! 3rd tetrahedron: 2, 3, 4, 6
                v = geom_tet_center(verts(iv2),verts(iv3),verts(iv4),verts(iv6))
                r = geom_tet_volume(verts(iv2),verts(iv3),verts(iv4),verts(iv6))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                r = 1.d0 / vol(ic)
                cell_cntr(ic) = r * cell_cntr(ic)
            CASE('hex')
                iv5 = iv2c(5)
                iv6 = iv2c(6)
                iv7 = iv2c(7)
                iv8 = iv2c(8)

                ! 1st tetrahedron: 1, 2, 4, 5
                v = geom_tet_center(verts(iv1),verts(iv2),verts(iv4),verts(iv5))
                r = geom_tet_volume(verts(iv1),verts(iv2),verts(iv4),verts(iv5))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                ! 2nd tetrahedron: 2, 3, 4, 7
                v = geom_tet_center(verts(iv2),verts(iv3),verts(iv4),verts(iv7))
                r = geom_tet_volume(verts(iv2),verts(iv3),verts(iv4),verts(iv7))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                ! 3rd tetrahedron: 7, 6, 5, 2
                v = geom_tet_center(verts(iv7),verts(iv6),verts(iv5),verts(iv2))
                r = geom_tet_volume(verts(iv7),verts(iv6),verts(iv5),verts(iv2))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                ! 4th tetrahedron: 8, 7, 5, 4
                v = geom_tet_center(verts(iv8),verts(iv7),verts(iv5),verts(iv4))
                r = geom_tet_volume(verts(iv8),verts(iv7),verts(iv5),verts(iv4))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                ! 5th tetrahedron: 4, 5, 2, 7
                v = geom_tet_center(verts(iv4),verts(iv5),verts(iv2),verts(iv7))
                r = geom_tet_volume(verts(iv4),verts(iv5),verts(iv2),verts(iv7))
                vol(ic) = vol(ic) + r
                cell_cntr(ic) = cell_cntr(ic) + r * v

                r = 1.d0 / vol(ic)
                cell_cntr(ic) = r * cell_cntr(ic)
            CASE default
                ! Generic Polyhedron

                ! Evaluation of tetrahedra "4th" point P, set equal to the
                ! average of polyhedron vertices coordinates

                ! Retrieves V2C connectivity of the cell IC.
                CALL get_ith_conn(iv2c,v2c,ic)
                nv2c = SIZE(iv2c)

                p =  vector_(0.d0,0.d0,0.d0)

                DO i = 1, nv2c
                    iv = iv2c(i)
                    p = p + verts(iv)
                END DO
                r = 1.d0 / nvc(ic)
                p = r * p

                ! Loop over cell faces

                ! Retrieves F2C connectivity of the cell IC.
                CALL get_ith_conn(if2c,f2c,ic)
                nf2c = SIZE(if2c)

                DO k = 1, nf2c
                    IF = if2c(k)
                    is = islave(IF)
                    !
                    ! Loop over face sub-triangles
                    CALL get_ith_conn(iv2f,v2f,IF)
                    iv_1 = iv2f(1)        ! 1st IV in V2F connectivity of the face IF
                    DO i = 3, nvf(IF)
                        iv_im1 = iv2f(i-1) ! (I-1)-th IV in V2F connectivity of the face IF
                        iv_i = iv2f(i)     ! I-th IV in V2F connectivity of the face IF
                        v = geom_tet_center(verts(iv_1),verts(iv_i),verts(iv_im1),vertex_(p))
                        r = geom_tet_volume(verts(iv_1),verts(iv_i),verts(iv_im1),vertex_(p))
                        IF(is == ic) THEN
                            ! IC is slave fo IF --> invert order of tetrahedron base vertices
                            v = geom_tet_center(verts(iv_1),verts(iv_im1),verts(iv_i),vertex_(p))
                            r = geom_tet_volume(verts(iv_1),verts(iv_im1),verts(iv_i),vertex_(p))
                        END IF
                        vol(ic) = vol(ic) + r
                        cell_cntr(ic) = cell_cntr(ic) + r * v
                    END DO
                END DO
                r = 1.d0 / vol(ic)
                cell_cntr(ic) = r * cell_cntr(ic)
            END SELECT
        END DO
    END IF

    IF (.NOT. quiet_) THEN
        ! Check for non-positive volume cell
        k = COUNT(vol <= 0.d0)
        IF(k > 0) THEN
            WRITE(*,300) k
            CALL abort_psblas
        END IF
    ENDIF

    !#############################################################################
    IF(debug) THEN
        DO ic = 1, ncells
            WRITE(*,500) ic, vol(ic)
        END DO
        WRITE(*,*)
500     FORMAT('cell: ',i5,' || VOL:',d13.6)
    END IF
    !#############################################################################

    !#############################################################################
    IF(debug) THEN
        DO ic = 1, ncells
            WRITE(*,510) ic, x_(cell_cntr(ic)), y_(cell_cntr(ic)), z_(cell_cntr(ic))
        END DO
        WRITE(*,*)
510     FORMAT('cell:',1x,i6,' | Coordinates:',3(1x,d13.6))
    END IF
    !############################################################################

    DEALLOCATE(islave,nvf,geoc,nvc)
    NULLIFY(iv2c,if2c,iv2f)

100 FORMAT(' ERROR! Input array not allocated in GEOM_CELL')
200 FORMAT(' ERROR! Memory allocation failure in GEOM_CELL')
300 FORMAT(' ERROR!',i8,' cells with NON-POSITIVE VOLUME in GEOM_CELL')

END PROCEDURE geom_cell
