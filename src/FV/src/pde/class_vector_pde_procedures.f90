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
! $Id: class_vector_pde.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_vector_pde) class_vector_pde_procedures
    USE class_pde
    USE class_stopwatch
    USE class_psblas, ONLY : icontxt_, mypnum_, nemo_sizeof_dp, psb_dupl_add_, psb_genrm2, sw_asb, sw_ins, &
        &                    abort_psblas, psb_geall, psb_geasb, psb_gefree, psb_geins,          &
        &                    psb_erractionsave, psb_erractionrestore, psb_check_error
    USE class_mesh
    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_sizeof
        USE psb_base_mod

        INTEGER(kind=nemo_int_long_) :: val

        val = eqn%pde%nemo_sizeof()
        IF (ALLOCATED(eqn%b)) &
            & val = val + nemo_sizeof_dp * SIZE(eqn%b)

        nemo_sizeof = val

    END PROCEDURE nemo_sizeof

    ! ----- Constructor -----

    MODULE PROCEDURE create_pde
        USE class_dimensions
        USE class_mesh
        !
        INTEGER :: err_act, info

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        ! Create PDE%BASE member
        CALL eqn%pde%create_pde(input_file,sec,msh,dim)

        ! Allocates RHS member
        CALL psb_geall(eqn%b,msh%desc_c,info,3)
        CALL psb_check_error(info,'create_vector_pde','psb_geall',icontxt_())

        eqn%b = 0.d0

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE create_pde


    ! ----- Destructor -----

    MODULE PROCEDURE free_pde
        USE class_mesh
        !
        INTEGER :: err_act, info
        TYPE(mesh), POINTER :: msh => NULL()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        CALL eqn%pde%get_mesh(msh)

        ! Frees storage of RHS member
        CALL psb_gefree(eqn%b,msh%desc_c,info)
        CALL psb_check_error(info,'free_pde','psb_gefree',icontxt_())

        NULLIFY(msh)

        ! Frees storage of BASE member
        CALL eqn%pde%free_pde()

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE free_pde


    ! ----- Getters -----

    !-----------------------------------------------

    MODULE PROCEDURE geins_vector_pde_r
        !! Inserts a ``cloud'' of RHS terms into pde%b
        USE class_mesh

        !
        INTEGER :: info, err_act
        TYPE(mesh), POINTER :: msh => NULL()

        CALL sw_ins%tic()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        IF(SIZE(cloud,2) /= 3) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pde%pde%get_mesh(msh)

        ! Inserts CLOUD into RHS pde%b
        CALL psb_geins(n,ia,cloud,pde%b,msh%desc_c,info,psb_dupl_add_)
        CALL psb_check_error(info,'geins_vector_pde','psb_geins',icontxt_())

        NULLIFY(msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

        CALL sw_ins%toc()

100     FORMAT(' ERROR! Size mismatch in GEINS_VECTOR_PDE')

    END PROCEDURE geins_vector_pde_r


    MODULE PROCEDURE geins_vector_pde_v
        !! Wrapper for ``clouds'' of VECTOR type
        USE class_mesh
        USE class_vector

        REAL(psb_dpk_) :: b(SIZE(cloud),3)

        CALL sw_ins%tic()

        b(:,1) = cloud%x_()
        b(:,2) = cloud%y_()
        b(:,3) = cloud%z_()

        CALL geins_vector_pde_r(n,ia,b,pde)

        CALL sw_ins%toc()

    END PROCEDURE geins_vector_pde_v


    MODULE PROCEDURE asb_pde_
        USE class_mesh

        INTEGER :: err_act, info
        TYPE(mesh), POINTER :: msh => NULL()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        ! Assemblies sparse matrix in PDE%BASE member
        CALL eqn%pde%asb_pde_()

        CALL sw_asb%tic()

        CALL eqn%pde%get_mesh(msh)

        ! Assemblies RHS member
        IF(mypnum_() == 0) WRITE(*,*) '  - assembling RHS'
        CALL psb_geasb(eqn%b,msh%desc_c,info)
        CALL psb_check_error(info,'abs_vector_pde','psb_geasb',icontxt_())

        CALL sw_asb%toc()

        NULLIFY(msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE asb_pde_


    MODULE PROCEDURE solve_vector_pde
        USE class_mesh
        USE class_vector
        USE class_vector_field

        !
        INTEGER :: err_act, icontxt, info, iter
        REAL(psb_dpk_), ALLOCATABLE :: x(:,:)
        REAL(psb_dpk_), ALLOCATABLE :: x_old(:,:)
        REAL(psb_dpk_) :: err
        TYPE(mesh), POINTER :: msh => NULL()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

        ! Assemblies PSBLAS sparse matrix and RHS
        CALL pde%asb_pde_()

        ! Build preconditioner in PDE%BASE member
        CALL pde%pde%build_pde_prec()

        ! Is PHI cell-centered?
        IF(phi%on_faces_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pde%get_mesh(msh)
        CALL phi%get_x(x)

        IF (PRESENT(var)) THEN
            ! Allocates X_OLD
            ALLOCATE(x_old(SIZE(x,1),3),stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            ! Copy X_OLD <- X
            x_old = x
        END IF

        ! Solves linear system associated to PDE%BASE
        CALL pde%pde%solve_pde_sys(pde%b(:,1),x(:,1),iter,err) ! Solve for X
        CALL pde%pde%solve_pde_sys(pde%b(:,2),x(:,2),iter,err) ! Solve for Y
        CALL pde%pde%solve_pde_sys(pde%b(:,3),x(:,3),iter,err) ! Solve for Z

        ! Assigns the solution to the vector field
        phi = vector_(x(:,1),x(:,2),x(:,3))

        ! Free preconditioner storage in base member
        CALL pde%pde%free_pde_prec()

        ! Updates phi
        CALL phi%update_field()

        IF(PRESENT(var)) THEN
            ! Computes the norm of the relative variation
            ! variation = || x_new - x_old || / ||x_old||
            var = psb_genrm2(x_old,msh%desc_c,info)
            CALL psb_check_error(info,'solve_vector_pde','psb_genrm2',icontxt)

            x_old = x - x_old

            var = psb_genrm2(x_old,msh%desc_c,info) / var
            CALL psb_check_error(info,'solve_vector_pde','psb_genrm2',icontxt)

            DEALLOCATE(x_old)
        END IF

        DEALLOCATE(x)
        NULLIFY(msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Unknown PHI in SOLVE_VECTOR_PDE is face-centered')
200     FORMAT(' ERROR! Memory allocation failure in SOLVE_VECTOR_PDE')

    END PROCEDURE solve_vector_pde


    MODULE PROCEDURE reinit_pde

!!$    write(0,*) 'Reinit Vector PDE do:',is_pde_asb(pde%base)
        IF(eqn%pde%is_pde_asb()) THEN
            eqn%b = 0.d0
!!$      write(0,*) 'Reinit Vector PDE out ',pde%b
            CALL eqn%pde%reinit_pde()
        END IF

    END PROCEDURE reinit_pde


    ! ----- Output -----

    MODULE PROCEDURE write_vector_pde
        USE class_mesh
        USE tools_output_basics

        !
        LOGICAL :: mtx_rhs
        TYPE(mesh), POINTER :: msh => NULL()

        CALL eqn%pde%write_pde(mat,mtx_rhs)

        IF(.NOT.mtx_rhs) RETURN

        CALL eqn%get_mesh(msh)
        CALL wr_mtx_vector(eqn%b(:,1),msh%desc_c,rhs//'_x')
        CALL wr_mtx_vector(eqn%b(:,2),msh%desc_c,rhs//'_y')
        CALL wr_mtx_vector(eqn%b(:,3),msh%desc_c,rhs//'_z')
        IF(mypnum_() == 0) WRITE(*,*)

        NULLIFY(msh)

    END PROCEDURE write_vector_pde

END SUBMODULE class_vector_pde_procedures
