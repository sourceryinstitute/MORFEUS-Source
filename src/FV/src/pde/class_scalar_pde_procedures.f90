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
! $Id: class_scalar_pde.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_scalar_pde) class_scalar_pde_procedures

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_scalar_pde_sizeof
        USE psb_base_mod

        INTEGER(kind=nemo_int_long_)   :: val

        val = eqn%base%nemo_sizeof()
        IF (ALLOCATED(eqn%b)) &
            & val = val + nemo_sizeof_dp * SIZE(eqn%b)

        nemo_scalar_pde_sizeof = val

    END PROCEDURE nemo_scalar_pde_sizeof

    ! ----- Constructor -----

    MODULE PROCEDURE create_scalar_pde
        USE class_dimensions
        USE class_mesh

        INTEGER :: err_act, info

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        ! Create PDE%BASE member
        CALL create_pde(pde%base,input_file,sec,msh,dim)

        ! Allocates RHS member
        CALL psb_geall(pde%b,msh%desc_c,info)
        CALL psb_check_error(info,'create_scalar_pde','psb_geall',icontxt_())

        pde%b = 0.d0

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE create_scalar_pde


    ! ----- Destructor -----



    MODULE PROCEDURE free_scalar_pde
        USE class_mesh

        INTEGER :: err_act, info
        TYPE(mesh), POINTER :: msh => NULL()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        CALL pde%base%get_mesh(msh)

        ! Frees storage of RHS member
        CALL psb_gefree(pde%b,msh%desc_c,info)
        CALL psb_check_error(info,'free_scalar_pde','psb_gefree',icontxt_())

        NULLIFY(msh)

        ! Frees storage of BASE member
        CALL free_pde(pde%base)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE free_scalar_pde


    ! ----- Getters -----

    MODULE PROCEDURE get_scalar_pde_name

        get_scalar_pde_name = name_(pde%base)

    END PROCEDURE get_scalar_pde_name


    MODULE PROCEDURE get_scalar_pde_dim
        USE class_dimensions

        get_scalar_pde_dim = dim_(pde%base)

    END PROCEDURE get_scalar_pde_dim


    MODULE PROCEDURE get_scalar_pde_msh_fun
        USE class_mesh

        get_scalar_pde_msh_fun => msh_(pde%base)

    END PROCEDURE get_scalar_pde_msh_fun

    MODULE PROCEDURE get_scalar_pde_msh_sub
        USE class_mesh

        CALL pde%base%get_mesh(msh)

    END PROCEDURE get_scalar_pde_msh_sub


    ! ----- Linear System Solving -----

    MODULE PROCEDURE spins_scalar_pde
        !! Inserts a ``cloud'' of coefficients into pde%A

        CALL spins_pde(n,ia,ja,cloud,pde%base)

    END PROCEDURE spins_scalar_pde


    MODULE PROCEDURE geins_scalar_pde
        !! Inserts a ``cloud'' of RHS terms into pde%b
        USE class_mesh


        INTEGER :: info, err_act
        TYPE(mesh), POINTER :: msh => NULL()

        CALL tic(sw_ins)

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        CALL pde%base%get_mesh(msh)

        ! Inserts CLOUD into RHS pde%b
        CALL psb_geins(n,ia,cloud,pde%b,msh%desc_c,info,psb_dupl_add_)
        CALL psb_check_error(info,'geins_scalar_pde','psb_geins',icontxt_())

        NULLIFY(msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

        CALL toc(sw_ins)

    END PROCEDURE geins_scalar_pde

    MODULE PROCEDURE asb_scalar_pde
        USE class_mesh

        INTEGER :: err_act, info
        TYPE(mesh), POINTER :: msh => NULL()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        ! Assemblies sparse matrix in PDE%BASE member
        CALL pde%base%asb_pde()

        CALL tic(sw_asb)

        CALL pde%base%get_mesh(msh)

        ! Assemblies RHS member
        IF(mypnum_() == 0) WRITE(*,*) '  - assembling RHS'
        CALL psb_geasb(pde%b,msh%desc_c,info)
        CALL psb_check_error(info,'abs_scalar_pde','psb_geasb',icontxt_())

        CALL toc(sw_asb)

        NULLIFY(msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE asb_scalar_pde

    MODULE PROCEDURE solve_scalar_pde
        USE class_mesh
        USE class_scalar_field

        INTEGER :: err_act, icontxt, info, iter
        REAL(psb_dpk_), ALLOCATABLE :: x(:)
        REAL(psb_dpk_), ALLOCATABLE :: x_old(:)
        REAL(psb_dpk_) :: err
        TYPE(mesh), POINTER :: msh => NULL()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

        ! Assemblies PSBLAS sparse matrix and RHS
        CALL pde%asb_pde()

        ! Build preconditioner in PDE%BASE member
        CALL build_pde_prec(pde%base)

        ! Is PHI cell-centered?
        IF(on_faces_(phi)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        CALL pde%get_mesh(msh)
        CALL get_x(phi,x)

        IF (PRESENT(var)) THEN
            ! Allocates X_OLD
            ALLOCATE(x_old(SIZE(x)),stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            ! Copy X_OLD <- X
            x_old = x
        END IF

        ! Solves linear system associated to PDE%BASE
        CALL solve_pde_sys(pde%base,pde%b,x,iter,err)

        ! Assigns solution to scalar field
        phi = x

        ! Free preconditioner storage in base member
        CALL free_pde_prec(pde%base)

        ! Updates phi
        CALL update_field(phi, mats)

        IF(PRESENT(var)) THEN
            ! Computes the norm of the relative variation
            ! variation = || x_new - x_old || / ||x_old||
            var = psb_genrm2(x_old,msh%desc_c,info)
            CALL psb_check_error(info,'solve_scalar_pde','psb_genrm2',icontxt)

            x_old = x - x_old

            var = psb_genrm2(x_old,msh%desc_c,info) / var
            CALL psb_check_error(info,'solve_scalar_pde','psb_genrm2',icontxt)

            DEALLOCATE(x_old)
        END IF

        DEALLOCATE(x)
        NULLIFY(msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Unknown PHI in SOLVE_SCALAR_PDE is face-centered')
200     FORMAT(' ERROR! Memory allocation failure in SOLVE_SCALAR_PDE')

    END PROCEDURE solve_scalar_pde


    MODULE PROCEDURE reinit_scalar_pde

!!$    write(0,*) 'Scalar_pde_reinit :',is_pde_asb(pde%base)
        IF(is_pde_asb(pde%base)) THEN
            pde%b = 0.d0
        END IF
        CALL reinit_pde(pde%base)

    END PROCEDURE reinit_scalar_pde


    ! ----- Output -----

    MODULE PROCEDURE write_scalar_pde
        USE class_mesh
        USE tools_output_basics

        LOGICAL :: mtx_rhs
        TYPE(mesh), POINTER :: msh => NULL()

        CALL write_pde(pde%base,mat,mtx_rhs)

        IF(.NOT.mtx_rhs) RETURN

        CALL pde%get_mesh(msh)

        CALL wr_mtx_vector(pde%b,msh%desc_c,rhs)
        IF(mypnum_() == 0) WRITE(*,*)

        NULLIFY(msh)

    END PROCEDURE write_scalar_pde

END SUBMODULE class_scalar_pde_procedures
