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
! $Id: class_pde.f90 9102 2015-04-24 16:06:49Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_pde) class_pde_procedures
    IMPLICIT NONE

    ! ----- Private Named Constants -----

    INTEGER, PARAMETER :: asb_ = 1
    INTEGER, PARAMETER :: bld_ = 2

    CONTAINS

    MODULE PROCEDURE nemo_sizeof
        USE class_psblas
        IMPLICIT NONE

        INTEGER(kind=nemo_int_long_)   :: val

        val = LEN(eqn%name)&
            & + LEN(eqn%cmethod) + LEN(eqn%cprec)
        val = val + 4 * nemo_sizeof_int + nemo_sizeof_dp
        val = val + eqn%dim%nemo_sizeof()
        val = val + psb_sizeof(eqn%a)
        IF (ALLOCATED(eqn%diag))&
            & val = val + nemo_sizeof_dp * SIZE(eqn%diag)
        val = val + psb_sizeof(eqn%prec)
        ! eqn%msh is an independent object

        nemo_sizeof = val

    END PROCEDURE nemo_sizeof


    MODULE PROCEDURE create_pde
        USE class_connectivity
        USE tools_input
        IMPLICIT NONE

        !
        LOGICAL, PARAMETER :: debug = .FALSE.
        !
        INTEGER :: info, err_act, nnz

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        ! Reads input parameters
        eqn%cmethod = TRIM(read_par(input_file,sec,'cmethod','BICGSTAB'))
        eqn%cprec   = TRIM(read_par(input_file,sec,'cprec','BJAC'))
        IF(eqn%cprec == 'NLDI' .OR. eqn%cprec == 'NLDU') THEN
            eqn%nlev = read_par(input_file,sec,'nlev',mandatory_i_)
        ELSE
            eqn%nlev    = 1
        END IF
        eqn%eps_solv   = read_par(input_file,sec,'eps_solv',1.d-05)
        eqn%itmax_solv = read_par(input_file,sec,'itmax_solv',100)
        eqn%mtx_sys    = read_par(input_file,sec,'mtx_sys',.FALSE.)

        ! Name
        eqn%name = TRIM(sec)

        ! Mesh pointer dereferencing
        eqn%msh => msh

        ! Dimensions
        eqn%dim = dim

        ! Computes the amount of non-zero elements in the local part of
        !  sparse matrix A (this is what SPALL is asking for!)
        nnz = msh%c2c%nconn_(gl='l')
        ! Adds diagonal elements
        nnz = nnz + msh%c2c%nel_(gl='l')

        ! Allocates sparse matrix A
        IF (debug) WRITE(0,*) mypnum_(),' Allocating with :',nnz
        CALL psb_spall(eqn%A,msh%desc_c,info,nnz=nnz)
        CALL psb_check_error(info,'create_pde','psb_spall',icontxt_())

        eqn%status = bld_

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

        IF(debug) THEN
            WRITE(*,*)
            WRITE(*,100) mypnum_()
            WRITE(*,200) 'Pde%name       = ', sec
            WRITE(*,200) 'Pde%cmethod    = ', eqn%cmethod
            WRITE(*,200) 'Pde%cprec      = ', eqn%cprec
            WRITE(*,300) 'Pde%nlev       = ', eqn%nlev
            WRITE(*,400) 'Pde%eps_solv   = ', eqn%eps_solv
            WRITE(*,300) 'Pde%itmax_solv = ', eqn%itmax_solv
            WRITE(*,500) 'Pde%mtx_sys    = ', eqn%mtx_sys
            WRITE(*,*)
        END IF

100     FORMAT(' ----- Process ID = ',i2,' -----')
200     FORMAT(1x,a,a)
300     FORMAT(1x,a,i5)
400     FORMAT(1x,a,es10.3)
500     FORMAT(1x,a,l5)

    END PROCEDURE create_pde


    ! ----- Destructor -----

    MODULE PROCEDURE free_pde
        IMPLICIT NONE
        !
        INTEGER :: err_act, info

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        CALL psb_spfree(eqn%A,eqn%msh%desc_c,info)
        CALL psb_check_error(info,'free_pde','psb_spfree',icontxt_())

        NULLIFY(eqn%msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE free_pde


    ! ----- Getters -----

    MODULE PROCEDURE get_pde_name
        IMPLICIT NONE

        get_pde_name = eqn%name

    END PROCEDURE get_pde_name


    MODULE PROCEDURE get_pde_dim
        IMPLICIT NONE

        get_pde_dim = eqn%dim

    END PROCEDURE get_pde_dim

    !----------------------------------------

    MODULE PROCEDURE get_pde_A
        IMPLICIT NONE

        B = eqn%A

    END PROCEDURE get_pde_A

    !----------------------------------------

    MODULE PROCEDURE get_pde_msh_fun
        IMPLICIT NONE

        get_pde_msh_fun => eqn%msh

    END PROCEDURE get_pde_msh_fun


    MODULE PROCEDURE get_pde_msh_sub
        IMPLICIT NONE

        msh => eqn%msh

    END PROCEDURE get_pde_msh_sub

    MODULE PROCEDURE get_pde_diag
        IMPLICIT NONE

        IF (.NOT. ALLOCATED(eqn%diag)) &
            & CALL eqn%update_diag()
        IF (.NOT. ALLOCATED(d)) &
            & ALLOCATE(d(SIZE(eqn%diag)))
        d = eqn%diag
    END PROCEDURE get_pde_diag

    MODULE PROCEDURE update_pde_diag
        IMPLICIT NONE
        INTEGER :: info, n

        !    n = psb_cd_get_local_rows(eqn%msh%desc_c)
        n = eqn%msh%desc_c%get_local_cols()

        IF (.NOT.ALLOCATED(eqn%diag)) THEN
            eqn%diag = eqn%a%get_diag(info)
        ENDIF

    END PROCEDURE update_pde_diag


    ! ----- Status Inquirer -----

    MODULE PROCEDURE is_pde_bld
        IMPLICIT NONE

        is_pde_bld = (eqn%status == bld_)

    END PROCEDURE is_pde_bld


    MODULE PROCEDURE is_pde_asb
        IMPLICIT NONE

        is_pde_asb = (eqn%status == asb_)

    END PROCEDURE is_pde_asb


    ! ----- Linear System Solving -----

    MODULE PROCEDURE spins_pde
        IMPLICIT NONE
        !! Inserts a ``cloud'' of coefficients into eqn%A
        INTEGER :: info, err_act

        CALL sw_ins%tic()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        ! Inserts CLOUD into sparse matrix eqn%A
        CALL psb_spins(n,ia,ja,cloud,eqn%A,eqn%msh%desc_c,info)
        CALL psb_check_error(info,'spins_pde','psb_spins',icontxt_())

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

        CALL sw_ins%toc()

    END PROCEDURE spins_pde


    MODULE PROCEDURE asb_pde_
        IMPLICIT NONE
        !
        INTEGER :: err_act, info
        LOGICAL, PARAMETER :: debug=.FALSE.

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        IF (eqn%is_pde_asb()) THEN
            ! What if this is called twice? Make it a no-op for now.
            RETURN
        END IF

        IF(mypnum_() == 0) THEN
            WRITE(*,*)
            WRITE(*,*) '* ', TRIM(eqn%name), ':'
            WRITE(*,*) '  - assembling sparse matrix'
        END IF

        CALL sw_asb%tic()

        ! Assemblies sparse matrix associated to eqn
        IF (.FALSE.) THEN
            CALL psb_spasb(eqn%A,eqn%msh%desc_c,info,afmt='CSR  ', &
                & upd=psb_upd_perm_,dupl=psb_dupl_add_)
        ELSE
            CALL psb_spasb(eqn%A,eqn%msh%desc_c,info,afmt='CSR  ', &
                & upd=psb_upd_srch_,dupl=psb_dupl_add_)
        END IF
        IF (debug) THEN
            WRITE(0,*) mypnum_(),' After asb: ',eqn%a%get_nrows(),&
                & eqn%a%get_ncols(), eqn%a%get_nzeros(), eqn%a%sizeof()
        END IF
        CALL psb_check_error(info,'asb_pde','psb_spasb',icontxt_())

        ! Changes eqn%status
        eqn%status = asb_

        IF (ALLOCATED(eqn%diag)) THEN
            DEALLOCATE(eqn%diag)
        ENDIF
        CALL sw_asb%toc()


        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE asb_pde_


    MODULE PROCEDURE build_pde_prec
        USE class_psblas, ONLY : abort_psblas
        USE tools_math,   ONLY : build_prec
        IMPLICIT NONE

        ! Make sure we are in the assembled state
        IF (eqn%is_pde_bld()) THEN
            CALL eqn%asb_pde_()
        END IF
        ! And if not, raise an error.
        IF (.NOT.eqn%is_pde_asb()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        IF(mypnum_() == 0) THEN
            WRITE(*,*)
            WRITE(*,*) '* Solving linear system in ', TRIM(eqn%name)
        END IF

        ! Builds preconditioner
        CALL build_prec(eqn%cprec,eqn%nlev,eqn%cmethod,&
            & eqn%A,eqn%msh%desc_c,eqn%prec)

100     FORMAT(' ERROR! PDE is in an invalid state on BUILD_PDE_PREC')

    END PROCEDURE build_pde_prec


    MODULE PROCEDURE free_pde_prec
        IMPLICIT NONE
        !
        INTEGER :: info, err_act

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        CALL psb_precfree(eqn%prec,info)
        CALL psb_check_error(info,'free_pde_prec','psb_precfree',icontxt_())

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE free_pde_prec


    MODULE PROCEDURE solve_pde_sys
        USE tools_math
        IMPLICIT NONE
        INTEGER :: i
        IF (debug_mat_bld) THEN
            WRITE(0,*) 'Scalar PDE A'
            CALL psb_csprt(0,eqn%a)
            WRITE(0,*) 'Scalar PDE B'
            DO i=1, eqn%a%get_nrows()
                WRITE(0,*) i,b(i)
            END DO
            WRITE(0,*) 'Initial guess'
            DO i=1, eqn%a%get_nrows()
                WRITE(0,*) i,x(i)
            END DO
        END IF
        CALL solve_sys(eqn%A,eqn%prec,b,x,eqn%msh%desc_c,&
            & eqn%cmethod,eqn%eps_solv,eqn%itmax_solv,iter,err)
        IF (debug_mat_bld) THEN
            WRITE(0,*) 'Solution'
            DO i=1, eqn%a%get_nrows()
                WRITE(0,*) i,x(i)
            END DO
        END IF
    END PROCEDURE solve_pde_sys


    MODULE PROCEDURE reinit_pde
        IMPLICIT NONE
        !
        INTEGER :: err_act, info

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        ! If EQN%STATUS is ASB_, then reinit the object
        IF(eqn%status == asb_) THEN

            CALL psb_sprn(eqn%A,eqn%msh%desc_c,info)
            CALL psb_check_error(info,'reinit_pde','psb_sprn',icontxt_())

            eqn%status = bld_
        END IF

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE reinit_pde


    ! ----- Output -----

    MODULE PROCEDURE write_pde
        USE tools_output_basics, ONLY : wr_mtx_matrix
        IMPLICIT NONE

        mtx_rhs = .FALSE.
        IF(.NOT.eqn%mtx_sys) RETURN

        CALL wr_mtx_matrix(eqn%A,eqn%msh%desc_c,mat)
        mtx_rhs = .TRUE.

    END PROCEDURE write_pde

END SUBMODULE class_pde_procedures
