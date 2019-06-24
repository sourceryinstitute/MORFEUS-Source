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
! $Id: psb_solve_sys.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE (tools_math) psb_solve_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE psb_solve_sys
        USE class_psblas
        USE psb_base_mod
        USE psb_prec_mod
        USE psb_krylov_mod
        IMPLICIT NONE
        !
        CHARACTER(len=20), PARAMETER  :: name_err = 'psb_solve_sys'
        INTEGER, PARAMETER :: itrace = -1, ml = 2
        TYPE(psb_d_vect_type) :: vb, vx
        REAL(psb_dpk_), ALLOCATABLE :: ax(:)
        !
        INTEGER :: info, err_act
        INTEGER :: icontxt, mypnum
        REAL(psb_dpk_) :: t_sol


        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        mypnum  = mypnum_()
        icontxt = icontxt_()

        IF(mypnum == 0) THEN
            WRITE(*,*) '  + Solving linear system with: ', TRIM(cmethod)
        END IF

        CALL sw_sol%tic()
        CALL vx%bld(x)
        CALL vb%bld(b)

        CALL psb_krylov(cmethod,A,prec,vb,vx,eps,desc,info,&
            & itmax=itmax,iter=iter,err=err,&
            & itrace=itrace,irst=ml)

        CALL psb_check_error(info,name_err,'psb_krylov',icontxt)
        IF (info == 0) THEN
            ax=vx%get_vect()
            x(:) = ax(:)
        END IF
        CALL sw_sol%toc()

        t_sol = sw_sol%partial_()

        ! System solving log message
        IF(mypnum == 0) THEN
            WRITE(*,100) '    - elapsed time:   ', t_sol, ' s'
            WRITE(*,100) '    - time/iteration : ', t_sol/MAX(iter,1), ' s'
            WRITE(*,200) '    - # of iterations: ', iter
            WRITE(*,300) '    - error on exit:  ', err
            WRITE(*,*)
        END IF


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(1x,a,es10.3,a)
200     FORMAT(1x,a,i5)
300     FORMAT(1x,a,es10.3)

        END PROCEDURE psb_solve_sys

END SUBMODULE psb_solve_implementation
