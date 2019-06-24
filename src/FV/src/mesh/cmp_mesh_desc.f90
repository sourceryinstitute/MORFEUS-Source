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
SUBMODULE (tools_mesh) cmp_mesh_desc_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE cmp_mesh_desc
        USE class_connectivity
        USE class_psblas
        USE tools_part
        IMPLICIT NONE
        !!
        !! $Id: cmp_mesh_desc.f90 8157 2014-10-09 13:02:44Z sfilippo $
        !!
        !! Description:
        !!
        INTEGER :: err_act, info, icontxt,mypnum
        INTEGER :: i, ic, IF, iv, k
        INTEGER :: ncells, nfaces, nverts, nconn, nelst
        INTEGER :: ilv, ilf, ilc, nlv, nlf, nlc
        INTEGER, ALLOCATABLE :: ia(:), ja(:)
        INTEGER, ALLOCATABLE :: loc_to_glob(:)
        INTEGER, POINTER :: ic2c(:) => NULL()
        INTEGER, POINTER :: if2f(:) => NULL()
        INTEGER, POINTER :: iv2v(:) => NULL()

        CALL tic(sw_dsc)

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()
        mypnum  = mypnum_()

        ! CELL broadcasting
        CALL bcast_conn(c2c)

        ncells = nel_(c2c)

        ! Checks status and size of PART_CELLS vector
        IF (mypnum == 0) THEN
            IF(.NOT.ALLOCATED(part_cells)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
            IF(SIZE(part_cells) /= ncells) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF
        ELSE
            IF(ALLOCATED(part_cells)) DEALLOCATE(part_cells)
            ALLOCATE(part_cells(ncells),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        END IF

        CALL psb_bcast(icontxt,part_cells,0)


        ! ----- CELL Descriptor Allocation -----
        CALL psb_cdall(icontxt,desc_c,info,vg=part_cells)
        CALL psb_check_error(info,'cmp_mesh_desc%desc_c','psb_cdall',icontxt)


        ! ----- CELL Descriptor Inserting -----

        ! Maximum number of elements in the CELL connectivity stencil
        nelst = max_conn(c2c) + 1

        ALLOCATE(ia(nelst),ja(nelst),stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

        ! Insertion loop on DESC_C descriptor
        nlc = psb_cd_get_local_rows(desc_c)
        CALL psb_get_loc_to_glob(desc_c,loc_to_glob)
        DO ilc = 1, nlc
            ic = loc_to_glob(ilc)
            k = 1

            ! First main diagonal elements...
            ia(k) = ic
            ja(k) = ic

            ! ... then off-diagonal elements...
            CALL get_ith_conn(ic2c,c2c,ic)
            nconn = SIZE(ic2c)
            DO i = 1, nconn
                k = k + 1
                ia(k) = ic
                ja(k) = ic2c(i)
            END DO

            !... and last, inserts the ic-th row in the PSBLAS descriptor.
            CALL psb_cdins(k,ia,ja,desc_c,info)
            IF(info /= 0) EXIT
        END DO
        CALL psb_check_error(info,'cmp_mesh_desc%desc_c','psb_cdins',icontxt)

        ! Frees storage
        NULLIFY(ic2c)
        DEALLOCATE(ia,ja,loc_to_glob,stat=info)
        IF(info /= 0) THEN
            WRITE(*,400)
            CALL abort_psblas
        END IF


        ! CELL Descriptors Assembling
        CALL psb_cdasb(desc_c,info) ! DESC_C
        CALL psb_check_error(info,'cmp_mesh_desc%desc_c','psb_cdasb',icontxt)

        ! CELL global to local reallocation
        CALL g2l_conn(c2c,desc_c,is_a2a=.TRUE.)

        !! FACES brodcasting

        CALL bcast_conn(f2f)
        CALL bcast_conn(f2c)

        nfaces = nel_(f2f)

        ! Creates C2F connectivities, stored in the TOOLS_PART module
        CALL get_dual_conn(f2c,c2f)

        ! ----- FACE  Descriptor Allocation -----
        CALL psb_cdall(icontxt,desc_f,info,mg=nfaces,parts=part_faces)
        CALL psb_check_error(info,'cmp_mesh_desc%desc_f','psb_cdall',icontxt)

        ! ----- FACE Descriptor Inserting -----

        ! Maximum number of elements in the FACE connectivity stencil
        nelst = max_conn(f2f) + 1

        ALLOCATE(ia(nelst),ja(nelst),stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

        ! Insertion loop on DESC_F descriptor
        nlf = psb_cd_get_local_rows(desc_f)
        CALL psb_get_loc_to_glob(desc_f,loc_to_glob)
        DO ilf = 1, nlf
            IF = loc_to_glob(ilf)
            k = 1

            ! First main diagonal elements...
            ia(k) = IF
            ja(k) = IF

            ! ... then off-diagonal elements...
            CALL get_ith_conn(if2f,f2f,IF)
            nconn = SIZE(if2f)
            DO i = 1, nconn
                k = k + 1
                ia(k) = IF
                ja(k) = if2f(i)
            END DO

            !... and last, inserts the if-th row in the PSBLAS descriptor.
            CALL psb_cdins(k,ia,ja,desc_f,info)
            IF(info /= 0) EXIT
        END DO
        CALL psb_check_error(info,'cmp_mesh_desc%desc_f','psb_cdins',icontxt)

        ! Frees storage
        NULLIFY(if2f)
        DEALLOCATE(ia,ja,loc_to_glob,stat=info)
        IF(info /= 0) THEN
            WRITE(*,400)
            CALL abort_psblas
        END IF

        ! FACES Descriptors Assembling
        CALL psb_cdasb(desc_f,info) ! DESC_F
        CALL psb_check_error(info,'cmp_mesh_desc%desc_f','psb_cdasb',icontxt)

        ! Frees memory oheld by shared variables in TOOLS_PART module
        CALL free_conn(c2f)

        ! FACES global to ocal reallocation
        CALL g2l_conn(f2c,desc_f,desc_c)
        CALL g2l_conn(f2f,desc_f,is_a2a=.TRUE.)

        ! VERTEX broadcasting

        CALL bcast_conn(v2v)
        CALL bcast_conn(v2c)

        nverts = nel_(v2v)

        ! Creates C2V , stored in the TOOLS_PART module
        CALL get_dual_conn(v2c,c2v)

        ! ----- VERTEX Descriptor Allocation -----
        CALL psb_cdall(icontxt,desc_v,info,mg=nverts,parts=part_verts)
        CALL psb_check_error(info,'cmp_mesh_desc%desc_v','psb_cdall',icontxt)

        ! ----- VERTEX Descriptor Inserting -----

        ! Maximum number of elements in the VERTEX connectivity stencil
        nelst = max_conn(v2v) + 1

        ALLOCATE(ia(nelst),ja(nelst),stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

        ! Insertion loop on DESC_V descriptor
        nlv = psb_cd_get_local_rows(desc_v)
        CALL psb_get_loc_to_glob(desc_v,loc_to_glob)
        DO ilv = 1, nlv
            iv = loc_to_glob(ilv)
            k = 1

            ! First main diagonal elements...
            ia(k) = iv
            ja(k) = iv

            ! ... then off-diagonal elements...
            CALL get_ith_conn(iv2v,v2v,iv)
            nconn = SIZE(iv2v)
            DO i = 1, nconn
                k = k + 1
                ia(k) = iv
                ja(k) = iv2v(i)
            END DO

            !... and last, inserts the iv-th row in the PSBLAS descriptor.
            CALL psb_cdins(k,ia,ja,desc_v,info)
            IF(info /= 0) EXIT
        END DO
        CALL psb_check_error(info,'cmp_mesh_desc%desc_v','psb_cdins',icontxt)

        ! Frees storage
        NULLIFY(iv2v)
        DEALLOCATE(ia,ja,loc_to_glob,stat=info)
        IF(info /= 0) THEN
            WRITE(*,400)
            CALL abort_psblas
        END IF

        ! VERTEX Descriptors Assembling
        CALL psb_cdasb(desc_v,info) ! DESC_V
        CALL psb_check_error(info,'cmp_mesh_desc%desc_v','psb_cdasb',icontxt)

        ! Frees memory oheld by shared variables in TOOLS_PART module
        CALL free_conn(c2v)

        !VERTEX global to local reallocation
        CALL g2l_conn(v2v,desc_v,is_a2a=.TRUE.)

        DEALLOCATE(part_cells,stat=info)
        IF(info /= 0) THEN
            WRITE(*,400)
            CALL abort_psblas
        END IF

        ! ----- Normal termination -----
        CALL psb_erractionrestore(err_act)

        CALL toc(sw_dsc)

100     FORMAT(' ERROR! PART_CELLS pointer in CMP_MESH_DESC is not associated')
200     FORMAT(' ERROR! Size mismatch of PART_CELLS pointer in CMP_MESH_DESC')
300     FORMAT(' ERROR! Memory allocation failure in CMP_MESH_DESC')
400     FORMAT(' ERROR! Memory deallocation failure in CMP_MESH_DESC')

        END PROCEDURE cmp_mesh_desc

END SUBMODULE cmp_mesh_desc_implementation
