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
! $Id: psb_build_prec.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
!SUBMODULE (tools_math) psb_build_implementation
!    IMPLICIT NONE
!
!    CONTAINS

!        MODULE PROCEDURE psb_build_prec
    SUBROUTINE psb_build_prec(cprec,nlev,cmethod,A,desc_a,prec)
        USE class_psblas
        USE psb_base_mod
        USE psb_prec_mod
        IMPLICIT NONE
        CHARACTER(len=*),      INTENT(IN) :: cprec
        INTEGER,               INTENT(IN) :: nlev
        CHARACTER(len=*),      INTENT(IN) :: cmethod
        TYPE(psb_dspmat_type), INTENT(IN) :: A
        TYPE(psb_desc_type),   INTENT(INOUT) :: desc_a
        TYPE(psb_dprec_type),  INTENT(INOUT) :: prec
        !
        INTEGER(psb_ipk_) :: info
        INTEGER :: err_act
        INTEGER :: icontxt, mypnum
        INTEGER :: nlev_


        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        ! Temporary dummy argument
        nlev_ = nlev

        icontxt = icontxt_()
        mypnum  = mypnum_()

        CALL sw_pre%tic()

        ! Prepares the preconditioner
        SELECT CASE(psb_toupper(TRIM(cprec)))
        CASE('NOPREC','NONE','DIAG','BJAC')
            ! These are OK.
            CALL prec%init(TRIM(cprec),info)

    !!$  case('2LDI4')
    !!$    call prec%init('ml',info,nlev=2)
    !!$    call prec%set(sub_restr_,psb_halo_,info)
    !!$    call prec%set(coarse_solve_,ilu_n_,info)
    !!$    call prec%set(coarse_sweeps_,4,info)
    !!$  case('2LDU4')
    !!$    call prec%init('ml',info,nlev=2)
    !!$    call prec%set(sub_restr_,psb_halo_,info)
    !!$    call prec%set(coarse_solve_,umf_,info)
    !!$    call prec%set(coarse_sweeps_,4,info)
    !!$  case('3LDI4')
    !!$    call prec%init('ml',info,nlev=3)
    !!$    call prec%set(sub_restr_,psb_halo_,info)
    !!$    call prec%set(coarse_solve_,ilu_n_,info)
    !!$    call prec%set(coarse_sweeps_,4,info)
    !!$  case('3LDU4')
    !!$    call prec%init('ml',info,nlev=3)
    !!$    call prec%set(sub_restr_,psb_halo_,info)
    !!$    call prec%set(coarse_solve_,umf_,info)
    !!$    call prec%set(coarse_sweeps_,4,info)
    !!$  case('2LI4S')
    !!$    call prec%init('ml',info,nlev=2)
    !!$    call prec%set(coarse_solve_,ilu_n_,info)
    !!$    call prec%set(coarse_sweeps_,4,info)
    !!$  case('2LU4S')
    !!$    call prec%init('ml',info,nlev=2)
    !!$    call prec%set(coarse_solve_,umf_,info)
    !!$    call prec%set(coarse_sweeps_,4,info)
    !!$  case('3LI4S')
    !!$    call prec%init('ml',info,nlev=3)
    !!$    call prec%set(coarse_solve_,ilu_n_,info)
    !!$    call prec%set(coarse_sweeps_,4,info)
    !!$  case('3LU4S')
    !!$    call prec%init('ml',info,nlev=3)
    !!$    call prec%set(coarse_solve_,umf_,info)
    !!$    call prec%set(coarse_sweeps_,4,info)
        CASE DEFAULT
            WRITE(*,100)
            CALL abort_psblas
        END SELECT

        CALL psb_check_error(info,'psb_build_prec','psb_precinit',icontxt)

        ! Checks preconditioner compatibility with CG method.
        IF(psb_toupper(cmethod) == 'CG') THEN
            SELECT CASE(psb_toupper(TRIM(cprec)))
            CASE ('NOPREC','DIAG','BJAC','2LI4S','2LU4S','3LI4S','3LU4S')
                ! Things are OK
            CASE DEFAULT
                WRITE(*,200)
                CALL abort_psblas
            END SELECT
        END IF

    !!$  ! RAS and Multilevel preconditioners requires BICGSTAB.
    !!$  if(    cprec == 'RASI'  .or. &
    !!$       & cprec == 'NLDI' .or. &
    !!$       & cprec == 'NLDU') then
    !!$    if(cmethod /= 'BICGSTAB') then
    !!$      write(*,300)
    !!$      call abort_psblas
    !!$    end if
    !!$  end if

        ! Builds preconditioner
        IF(mypnum == 0) THEN
            WRITE(*,*) '  + Setting preconditioner to: ',TRIM(cprec)
        END IF
        CALL prec%build(A,desc_a,info)
        CALL psb_check_error(info,'psb_build_prec','psb_precbld',icontxt)


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

        CALL sw_pre%toc()

100     FORMAT(' ERROR! Unsupported preconditioner in PSB_BUILD_PREC')
200     FORMAT(' ERROR! Unsymmetric preconditioner cannot be used with CG')
300     FORMAT(' ERROR! RAS(ML) preconditioners need BICGSTAB')
    END SUBROUTINE psb_build_prec
!        END PROCEDURE psb_build_prec

!END SUBMODULE psb_build_implementation
