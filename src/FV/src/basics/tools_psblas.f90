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
MODULE tools_psblas

    USE psb_base_mod
    USE psb_prec_mod
    USE psb_krylov_mod
    IMPLICIT NONE
    ! It provides:
    ! - PSBLAS interfaces.
    ! - PSBLAS derived data types, such as PSB_DESCRIPTOR_TYPE,
    !   PSB_PREC_TYPE and PSB_SPMAT_TYPE.
!!$  integer, parameter  :: longndig=12
!!$  integer, parameter  :: nemo_int_long_ = selected_int_kind(longndig)
    INTEGER, PARAMETER  :: nemo_int_long_ = psb_long_int_k_
    INTEGER, PARAMETER  :: nemo_dpk_ = psb_dpk_


    ! ----- Generic Interfaces -----

    INTERFACE psb_gather
        MODULE SUBROUTINE psb_hgatherv(hglobx,hlocx,desc_a,info,root)
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(OUT) :: hglobx(:)
            CHARACTER(len=*), INTENT(IN) :: hlocx(:)
            TYPE(psb_desc_type), INTENT(IN) :: desc_a
            INTEGER, INTENT(OUT) :: info
            INTEGER, INTENT(IN), OPTIONAL :: root
        END SUBROUTINE psb_hgatherv
    END INTERFACE psb_gather


    ! ----- Explicit Interfaces -----


    ! ----- Global to Local -----
    INTERFACE
        MODULE SUBROUTINE psb_get_glob_to_loc(desc,iglob_to_loc)
            IMPLICIT NONE
            TYPE(psb_desc_type), INTENT(IN) :: desc
            INTEGER, ALLOCATABLE  :: iglob_to_loc(:)
        END SUBROUTINE psb_get_glob_to_loc

        MODULE FUNCTION glob_to_loc_(desc,iglob)
            IMPLICIT NONE
            INTEGER :: glob_to_loc_
            TYPE(psb_desc_type), INTENT(IN) :: desc
            INTEGER, INTENT(IN) :: iglob
        END FUNCTION glob_to_loc_
    END INTERFACE

    ! ----- Local To Global -----

    INTERFACE
        MODULE SUBROUTINE psb_get_loc_to_glob(desc,iloc_to_glob)
            IMPLICIT NONE
            TYPE(psb_desc_type), INTENT(IN) :: desc
            INTEGER, ALLOCATABLE  :: iloc_to_glob(:)
        END SUBROUTINE psb_get_loc_to_glob

        MODULE FUNCTION loc_to_glob_(desc,iloc)
            IMPLICIT NONE
            INTEGER :: loc_to_glob_
            TYPE(psb_desc_type), INTENT(IN) :: desc
            INTEGER, INTENT(IN) :: iloc
        END FUNCTION loc_to_glob_
    END INTERFACE

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



    ! ----- Error Handling -----
    INTERFACE
        MODULE SUBROUTINE psb_check_error(info,WHERE,from,icontxt)
        ! Checks in WHERE procedure the error code returned by FROM.
        ! If it's positive prints error stack and forces PSBLAS abortion.
        !
            IMPLICIT NONE
            INTEGER, INTENT(INOUT) :: info
            CHARACTER(len=*), INTENT(IN) :: WHERE, from
            INTEGER, INTENT(IN) :: icontxt
        END SUBROUTINE psb_check_error
    END INTERFACE

END MODULE tools_psblas
