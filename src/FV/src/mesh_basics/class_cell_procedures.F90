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
! $Id: class_cell.F90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    Provides cell class functionality.
!
! Includes:
!    CELL class              (all private)
!    CELL_                   basic elemental constructor
!    ALLOC_CELL              construct array pointer for numerous cells
!    FREE_CELL               free array pointers
!    BCAST_CELL              broadcast cells from node 0 to other processes
!    G2L_CELL                global to local reallocation of cells
!    GET_CELL_NV             get number of vertices of cell (elemental)
!    GET_CELL_GEO            get character code indicating type of cell
!    GET_CELLS_TYPE          concatenate lists of cells by type
!

SUBMODULE(class_cell) class_cell_procedures

    USE class_psblas

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_cell_sizeof
        USE psb_base_mod

        nemo_cell_sizeof = 2 * nemo_sizeof_int + LEN(cll%geo)

    END PROCEDURE nemo_cell_sizeof

    ! ----- Constructors -----

    MODULE PROCEDURE cell_
        !! Public default constructor

        !cell_ = cell(nv,nf,group,geo)
        !! Workaround for Intel 18.0.5 error #6053: Structure constructor may not have fields with a PRIVATE attribute
        cell_%nv    = nv
        cell_%nf    = nf
        cell_%group = group
        cell_%geo   = geo

    END PROCEDURE cell_


    MODULE PROCEDURE alloc_cell
        !! Array constructor
        !
        INTEGER :: info

        IF (ALLOCATED(cells))THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ALLOCATE(cells(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Illegal allocation: CELL array already associated')
200     FORMAT(' ERROR! Memory allocation failure in ALLOC_CELL')

    END PROCEDURE alloc_cell


    ! ----- Destructor -----

    MODULE PROCEDURE free_cell
        !
        INTEGER :: info

        IF (ALLOCATED(cells)) THEN
            DEALLOCATE(cells,stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        END IF
100     FORMAT(' ERROR! Memory deallocation failure in FREE_CELL')

    END PROCEDURE free_cell


    ! ----- Parallel Ops. -----

    MODULE PROCEDURE bcast_cell
        INTEGER, ALLOCATABLE :: intbuf(:,:)
        CHARACTER(len=nlen), ALLOCATABLE :: charbuf(:)
        INTEGER :: i, info, n
        INTEGER :: icontxt, mypnum

        mypnum = mypnum_()
        icontxt = icontxt_()

        ! CELLS(:) is supposed to be associated only in P0
        ! => Check and allocation on P > 0
        IF(mypnum == 0) THEN
            IF(.NOT.ALLOCATED(cells)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            n = SIZE(cells)

            CALL psb_bcast(icontxt,n)
        ELSE
            IF (ALLOCATED(cells)) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            CALL psb_bcast(icontxt,n)

            CALL alloc_cell(cells,n)
        END IF

        ALLOCATE(intbuf(3,n),charbuf(n),stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
        END IF

        ! Broadcasting
        IF(mypnum == 0) THEN

            DO i = 1, n
                intbuf(1,i) = cells(i)%nv
                intbuf(2,i) = cells(i)%nf
                intbuf(3,i) = cells(i)%group
                charbuf(i)  = cells(i)%geo
            END DO

            CALL psb_bcast(icontxt,intbuf)
            CALL psb_bcast(icontxt,charbuf)
        ELSE

            CALL psb_bcast(icontxt,intbuf)
            CALL psb_bcast(icontxt,charbuf)

            DO i = 1, n
                cells(i)%nv  = intbuf(1,i)
                cells(i)%nf  = intbuf(2,i)
                cells(i)%group = intbuf(3,i)
                cells(i)%geo = charbuf(i)
            END DO

        END IF

        DEALLOCATE(intbuf,charbuf)

100     FORMAT(' ERROR! Illegal Send: CELL array not associated')
200     FORMAT(' ERROR! Illegal Recv: CELL array already associated')
300     FORMAT(' ERROR! Memory allocation failure in BCAST_CELL')

    END PROCEDURE bcast_cell


    MODULE PROCEDURE g2l_cell
        USE psb_base_mod
        !
        INTEGER :: ic_glob, ic_loc
        INTEGER :: ncells_loc
        INTEGER, ALLOCATABLE :: iloc_to_glob(:)
        TYPE(cell), ALLOCATABLE :: work(:)


        ncells_loc  = psb_cd_get_local_cols(desc_c)

        CALL alloc_cell(work,ncells_loc)

        CALL psb_get_loc_to_glob(desc_c,iloc_to_glob)

        DO ic_loc = 1, ncells_loc
            ic_glob = iloc_to_glob(ic_loc)
            work(ic_loc) = cells(ic_glob)
        END DO
#if defined(HAVE_MOVE_ALLOC)
        CALL move_ALLOC(work,cells)
#else
        IF (SIZE(cells) /= SIZE(work)) THEN
            DEALLOCATE(cells)
            ALLOCATE(cells(SIZE(work)))
        ENDIF
        cells = work
        DEALLOCATE(work)
#endif
        DEALLOCATE(iloc_to_glob)
    END PROCEDURE g2l_cell


    MODULE PROCEDURE l2g_cell
        USE psb_base_mod
        ! WARNING! The global results is allocated only on P0. After its usage
        ! it must be deallocated in the calling unit by means of the statement:
        ! "if(associated(glob_res)) deallocate(glob_res)"

        INTEGER :: icontxt, info, err_act
        INTEGER :: ic, n_glob, n_loc
        INTEGER, ALLOCATABLE  :: ibuf_glob(:,:)
        INTEGER, ALLOCATABLE  :: ibuf_loc(:,:)
        CHARACTER(len=nlen), ALLOCATABLE :: cbuf_glob(:)
        CHARACTER(len=nlen), ALLOCATABLE :: cbuf_loc(:)


        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

        n_loc  = SIZE(cells_loc)            ! # local cells
        n_glob = psb_cd_get_global_cols(desc_c) ! # global cells

        CALL psb_geall(ibuf_loc,desc_c,info,3)
        CALL psb_check_error(info,'l2g_cell','psb_geall',icontxt)

        ALLOCATE(ibuf_glob(n_glob,3),&
            &   cbuf_glob(n_glob),cbuf_loc(n_loc),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF


        DO ic = 1, n_loc
            ibuf_loc(ic,1) = cells_loc(ic)%nv
            ibuf_loc(ic,2) = cells_loc(ic)%nf
            ibuf_loc(ic,3) = cells_loc(ic)%group
            cbuf_loc(ic)   = cells_loc(ic)%geo
        END DO


        CALL psb_gather(ibuf_glob,ibuf_loc,desc_c,info,root=0)
        CALL psb_check_error(info,'l2g_cell','psb_igatherm',icontxt)

        CALL psb_gather(cbuf_glob,cbuf_loc,desc_c,info,root=0)
        CALL psb_check_error(info,'l2g_cell','psb_hgatherv',icontxt)


        IF(mypnum_() == 0) THEN
            IF (ALLOCATED(cells_glob)) CALL free_cell(cells_glob)
            CALL alloc_cell(cells_glob,n_glob)

            DO ic = 1, n_glob
                cells_glob(ic)%nv  = ibuf_glob(ic,1)
                cells_glob(ic)%nf  = ibuf_glob(ic,2)
                cells_glob(ic)%group = ibuf_glob(ic,3)
                cells_glob(ic)%geo = cbuf_glob(ic)
            END DO
        END IF


        ! Frees memory
        DEALLOCATE(cbuf_glob,cbuf_loc,ibuf_glob)

        CALL psb_gefree(ibuf_loc,desc_c,info)
        CALL psb_check_error(info,'l2g_cell','psb_gefree',icontxt)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in L2G_CELL')

    END PROCEDURE l2g_cell


    ! ----- Getters -----

    MODULE PROCEDURE get_cell_nv

        get_cell_nv = c%nv

    END PROCEDURE get_cell_nv


    MODULE PROCEDURE get_cell_geo

        get_cell_geo = c%geo

    END PROCEDURE get_cell_geo

    MODULE PROCEDURE get_cell_group

        get_cell_group = c%group

    END PROCEDURE get_cell_group

    ! This getter function can't be defined as elemental because
    ! of a pending bug in GFortran in handling elemental functions
    ! of character type.

    MODULE PROCEDURE get_cells_type
        USE psb_base_mod
        !
        INTEGER :: err_act, icontxt, info
        INTEGER :: ic, n, ncells_glob, ncells_loc
        REAL(psb_dpk_), ALLOCATABLE :: glob_vect(:)
        REAL(psb_dpk_), ALLOCATABLE :: loc_vect(:)
        !
        ! Polygons and polyhedra numbering
        INTEGER :: ktri, kqua
        INTEGER :: ktet, kpyr, kpri, khex


        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

!!$    No longer necessary because of INTENT(OUT)
!!$    ! Checks that NCTYPE (output) is not associated
!!$    if(associated(nctype)) then
!!$       deallocate(nctype)
!!$    end if
!!$
!!$    ! Checks that ICTYPE (output) is not associated
!!$    if(associated(ictype)) then
!!$       deallocate(ictype)
!!$    end if

        ! Sets total number of supported element types
        n = 6

        IF (PRESENT(desc)) THEN
            ! Local
            ncells_glob = psb_cd_get_global_cols(desc)
            ncells_loc  = psb_cd_get_local_cols(desc)
        ELSE
            ! Global
            ncells_glob = SIZE(cells)
            ncells_loc = ncells_glob
        END IF

        ALLOCATE(nctype(n),ictype(ncells_glob), &
            &   glob_vect(ncells_glob),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        glob_vect = 0.d0

        IF(PRESENT(desc)) THEN
            CALL psb_geall(loc_vect,desc,info)
            CALL psb_check_error(info,'get_cells_type','psb_geall',icontxt)
            loc_vect = 0.d0
        ELSE
            ALLOCATE(loc_vect(SIZE(glob_vect)),stat=info)
        END IF


        WHERE(    cells%geo == 'tri')
            loc_vect = itri_
        ELSEWHERE(cells%geo == 'qua')
            loc_vect = iqua_
        ELSEWHERE(cells%geo == 'tet')
            loc_vect = itet_
        ELSEWHERE(cells%geo == 'pyr')
            loc_vect = ipyr_
        ELSEWHERE(cells%geo == 'pri')
            loc_vect = ipri_
        ELSEWHERE(cells%geo == 'hex')
            loc_vect = ihex_
        END WHERE

        IF(PRESENT(desc)) THEN
!!$       call  psb_gather(glob_vect,loc_vect,desc,info,root=-1) ! BUG!!
            CALL psb_gather(glob_vect,loc_vect,desc,info,root=0)
            CALL psb_check_error(info,'get_cells_type','psb_gather',icontxt)

            IF(mypnum_() == 0) THEN
                CALL psb_bcast(icontxt_(),glob_vect)
            ELSE
                CALL psb_bcast(icontxt_(),glob_vect)
            END IF
            ! WARNING! This could be avoided if there wouldn't be a bug in
            ! PSBLAS such as ROOT=-1 in psb_gather doesn't provide for a
            ! copy of glob_vect on every process.
        ELSE
            glob_vect=loc_vect
        END IF


        ! Computes values for NCTYPE
        nctype(itri_) = COUNT(glob_vect == itri_)
        nctype(iqua_) = COUNT(glob_vect == iqua_)
        nctype(itet_) = COUNT(glob_vect == itet_)
        nctype(ipyr_) = COUNT(glob_vect == ipyr_)
        nctype(ipri_) = COUNT(glob_vect == ipri_)
        nctype(ihex_) = COUNT(glob_vect == ihex_)


        ! Sorts the cells indices according to the cells types
        ktri = 0
        kqua = ktri + nctype(itri_)
        ktet = kqua + nctype(iqua_)
        kpyr = ktet + nctype(itet_)
        kpri = kpyr + nctype(ipyr_)
        khex = kpri + nctype(ipri_)
        !
        DO ic = 1, ncells_glob
            SELECT CASE(INT(glob_vect(ic)))
            CASE(itri_)
                ktri = ktri + 1
                ictype(ktri) = ic
            CASE(iqua_)
                kqua = kqua + 1
                ictype(kqua) = ic
            CASE(itet_)
                ktet = ktet + 1
                ictype(ktet) = ic
            CASE(ipyr_)
                kpyr = kpyr + 1
                ictype(kpyr) = ic
            CASE(ipri_)
                kpri = kpri + 1
                ictype(kpri) = ic
            CASE(ihex_)
                khex = khex + 1
                ictype(khex) = ic
            CASE default
                WRITE(*,200)
                CALL abort_psblas
            END SELECT
        END DO

        IF(PRESENT(desc)) THEN
            CALL psb_gefree(loc_vect,desc,info)
            CALL psb_check_error(info,'get_cells_type','psb_gefree',icontxt)
        ELSE
            DEALLOCATE(loc_vect,stat=info)
        END IF

        IF (info==0) DEALLOCATE(glob_vect,stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
        END IF

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in GET_CELLS_TYPE')
200     FORMAT(' ERROR! Unsupported type of element in GET_CELLS_TYPE')
300     FORMAT(' ERROR! Memory deallocation failure in GET_CELLS_TYPE')

    END PROCEDURE get_cells_type

END SUBMODULE class_cell_procedures
