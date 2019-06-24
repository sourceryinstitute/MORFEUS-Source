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
! $Id$
!
! Description:
!    Contains wrappers and type extensions of PSBLAS routines.
!
! Provides:
!
! PSB_GET_LOC_TO_GLOB: gets the local to global conversion list.
! PSB_GET_GLOB_TO_LOC: gets the globla to local conversion list
! GLOB_TO_LOC_:        single element global to local conversion
! LOC_TO_GLOB_:        single element local to global conversion
! PSB_HBCASTV:         broadcasts an array of strings.
! PSB_HGATHERV:        gathers an array of strings
! PSB_CHECK_ERROR:     checks the error code returned by a PSBLAS call.
!
SUBMODULE(tools_psblas) tools_psblas_procedures 

    USE psb_base_mod
    USE psb_prec_mod
    USE psb_krylov_mod
    IMPLICIT NONE

CONTAINS

    ! ----- Global to Local -----

    MODULE PROCEDURE psb_get_glob_to_loc
        !
        INTEGER :: info, err_act
        INTEGER :: i, icontxt, nglob

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        icontxt = psb_cd_get_context(desc)
        nglob   = psb_cd_get_global_rows(desc)

        CALL psb_realloc(nglob,iglob_to_loc,info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL psb_abort(icontxt)
        END IF

        iglob_to_loc = (/ (i, i = 1, nglob) /)

        CALL psb_glob_to_loc(iglob_to_loc,desc,info,iact='I')
        CALL psb_check_error(info,'psb_get_glob_to_loc','psb_glob_to_loc',&
            & icontxt)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in PSB_GET_GLOB_TO_LOC')

    END PROCEDURE psb_get_glob_to_loc


    MODULE PROCEDURE glob_to_loc_
        !
        INTEGER :: dmy(1), err_act, icontxt, info

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        icontxt = psb_cd_get_context(desc)

        dmy = iglob
        CALL psb_glob_to_loc(dmy,desc,info, iact='I')
        CALL psb_check_error(info,'glob_to_loc_','psb_glob_to_loc',&
            & icontxt)

        glob_to_loc_ = dmy(1)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE glob_to_loc_


    ! ----- Local To Global -----

    MODULE PROCEDURE psb_get_loc_to_glob
        !
        INTEGER :: info, err_act
        INTEGER :: i, icontxt, nloc

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        icontxt = psb_cd_get_context(desc)
        nloc    = psb_cd_get_local_cols(desc)

        CALL psb_realloc(nloc,iloc_to_glob,info)

        IF(info /= 0) THEN
            WRITE(*,100)
            CALL psb_abort(icontxt)
        END IF

        iloc_to_glob = (/ (i, i = 1, nloc) /)

        CALL psb_loc_to_glob(iloc_to_glob,desc,info)
        CALL psb_check_error(info,'psb_get_loc_to_glob','psb_loc_to_glob',&
            icontxt)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in PSB_GET_LOC_TO_GLOB')

    END PROCEDURE psb_get_loc_to_glob


    MODULE PROCEDURE loc_to_glob_
        !
        INTEGER :: dmy(1), err_act, icontxt, info

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        icontxt = psb_cd_get_context(desc)

        dmy = iloc
        CALL psb_loc_to_glob(dmy,desc,info)
        CALL psb_check_error(info,'loc_to_glob_','psb_loc_to_glob',&
            & icontxt)

        loc_to_glob_ = dmy(1)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE loc_to_glob_


!!$  ! ----- Extra PSBLAS Broadcast Routines -----
!!$
!!$  subroutine psb_hbcastv(ictxt,dat,root,length)
!!$    use mpi
!!$    integer, intent(in)             :: ictxt
!!$    character(len=*), intent(inout) :: dat(:)
!!$    integer, intent(in), optional   :: root, length
!!$    !
!!$    integer  :: root_, icomm, length_, info
!!$
!!$    if (present(root)) then
!!$      root_ = root
!!$    else
!!$      root_ = 0
!!$    endif
!!$    if (present(length)) then
!!$      length_ = length
!!$    else
!!$      length_ = len(dat)
!!$    endif
!!$
!!$    call psb_get_mpicomm(ictxt,icomm)
!!$
!!$    call mpi_bcast(dat,length_*size(dat),MPI_CHARACTER,root_,icomm,info)
!!$
!!$  end subroutine psb_hbcastv


    ! ----- Extra PSBLAS Gather Routines -----

    ! The integer versions are already in PSBLAS.

    MODULE PROCEDURE psb_hgatherv
        !
        INTEGER :: icontxt, root_
        INTEGER :: i, j, len_s, n_loc, n_glob
        INTEGER, ALLOCATABLE  :: iglobx(:,:)
        INTEGER, ALLOCATABLE  :: ilocx(:,:)
!!$    integer :: iglobx(size(hglobx),len(hglobx(1)))
!!$    integer :: ilocx(size(hlocx),len(hlocx(1)))

        icontxt = psb_cd_get_context(desc_a)

        IF(PRESENT(root)) THEN
            root_ = root
        ELSE
            root_ = -1
        END IF

        IF(LEN(hlocx(1)) /= LEN(hglobx(1))) THEN
            WRITE(*,100)
            CALL psb_abort(icontxt)
        END IF

        n_loc  = SIZE(hlocx)
        n_glob = SIZE(hglobx)
        len_s = LEN(hglobx(1))
        ALLOCATE(iglobx(n_glob,len_s),ilocx(n_loc,len_s))
        DO j = 1, len_s
            DO i = 1, n_loc
                ilocx(i,j) = IACHAR(hlocx(i)(j:j))
            END DO
        END DO

        CALL psb_gather(iglobx,ilocx,desc_a,info,root=root_)
        DO j = 1, len_s
            DO i = 1, n_glob
                hglobx(i)(j:j) = ACHAR(iglobx(i,j))
            END DO
        END DO

100     FORMAT(' ERROR! String lenght mismatch in PSB_HGATHERV')

    END PROCEDURE psb_hgatherv


    ! ----- Error Handling -----

    MODULE PROCEDURE psb_check_error
        ! Checks in WHERE procedure the error code returned by FROM.
        ! If it's positive prints error stack and forces PSBLAS abortion.
        !

        ! Check value of INFO error code returned by subroutine FROM
        IF(info /= 0) THEN
            info = 4010
            CALL psb_errpush(4013,WHERE,a_err=from,i_err=(/info, 0, 0, 0, 0/))
            CALL psb_error(icontxt)
        END IF

    END PROCEDURE psb_check_error


END SUBMODULE tools_psblas_procedures
