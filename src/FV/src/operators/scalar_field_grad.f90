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
!    Computes the gradient of a ScalarField object.
!    The result is a VECTOR object with rank 1.
!
!    REMARK: remember to allocate the LHS before assigning the function result
!
SUBMODULE(op_grad) scalar_field_grad_implementation

    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE scalar_field_grad
            USE class_psblas
            USE class_connectivity
            USE class_face
            USE class_least_squares
            USE class_mesh
            USE class_scalar_field
            USE class_vector

            IMPLICIT NONE
            !
            INTEGER :: i, ib, ibf, ic, IF, im, info, ib_offset
            INTEGER :: n, nb, nbc, ncd, ncells
            INTEGER, POINTER :: ic2c(:) => NULL(), if2b(:) => NULL()
            REAL(psb_dpk_), ALLOCATABLE :: phi_x(:)
            REAL(psb_dpk_), ALLOCATABLE :: phi_b(:)
            REAL(psb_dpk_), ALLOCATABLE :: rhs(:,:)
            TYPE(mesh), POINTER :: msh => NULL()

            ! Gets internal and boundary values
            CALL phi%get_x(phi_x)
            CALL phi%get_bx(phi_b)

            ! Face-centered fields are not supported
            IF(phi%on_faces_()) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            ! Ncells = local + halos
            ncells = SIZE(phi_x)

            IF(.NOT.ALLOCATED(grad)) CALL alloc_vector(grad,ncells)
            ! REMARK: compiler-dependent behaviour
            ! Intel:   allocated(grad) = F
            ! Gfortran: allocated(grad) = T

            ! Gets PHI mesh
            CALL phi%get_mesh(msh)
            ncd = msh%ncd

            ! Allocates RHS for least squares regression
            ALLOCATE(rhs(ncd+1,ncells),stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            ! Initialization
            grad(:) = vector_(0.d0,0.d0,0.d0)
            rhs = 0.d0

            ! Next loops sweep only strictly local cells => NCELLS redefinition.
            ncells = psb_cd_get_local_rows(msh%desc_c)
            nbc = msh%nbc

            ! Builds RHS
            DO ic = 1, ncells
                ! IC
                rhs(1,ic) = phi_x(ic)
                rhs(2,ic) = phi_x(ic) * msh%cell_cntr(ic)%x_()
                rhs(3,ic) = phi_x(ic) * msh%cell_cntr(ic)%y_()

                ! IC's neighbors
                CALL msh%c2c%get_ith_conn(ic2c,ic)
                n = SIZE(ic2c)
                DO i = 1, n
                    nb = ic2c(i)
                    rhs(1,ic) = rhs(1,ic) + phi_x(nb)
                    rhs(2,ic) = rhs(2,ic) + phi_x(nb) * msh%cell_cntr(nb)%x_()
                    rhs(3,ic) = rhs(3,ic) + phi_x(nb) * msh%cell_cntr(nb)%y_()
                END DO
            END DO

            ! Center(s) of possible boundary faces
            DO ib = 1, nbc
                CALL msh%f2b%get_ith_conn(if2b,ib)
                n = SIZE(if2b)
                ib_offset = COUNT(msh%faces%flag_() > 0 .AND. msh%faces%flag_() < ib)

                DO i = 1, n
                    IF = if2b(i)
                    ibf = ib_offset + i
                    im = msh%faces(IF)%master_()
                    rhs(1,im) = rhs(1,im) + phi_b(ibf)
                    rhs(2,im) = rhs(2,im) + phi_b(ibf) * msh%face_cntr(IF)%x_()
                    rhs(3,im) = rhs(3,im) + phi_b(ibf) * msh%face_cntr(IF)%y_()
                END DO
            END DO

            ! Extra elements for 3D case only
            IF(ncd == 3) THEN
                DO ic = 1, ncells
                    ! IC
                    rhs(4,ic) = phi_x(ic) * msh%cell_cntr(ic)%z_()

                    ! IC's neighbors
                    CALL msh%c2c%get_ith_conn(ic2c,ic)
                    n = SIZE(ic2c)
                    DO i = 1, n
                        nb = ic2c(i)
                        rhs(4,ic) = rhs(4,ic) + phi_x(nb) * msh%cell_cntr(nb)%z_()
                    END DO
                END DO

                ! Center(s) of possible boundary faces
                DO ib = 1, nbc
                    CALL msh%f2b%get_ith_conn(if2b,ib)
                    n = SIZE(if2b)
                    ib_offset = COUNT(msh%faces%flag_() > 0 .AND. msh%faces%flag_() < ib)

                    DO i = 1, n
                        IF = if2b(i)
                        ibf = ib_offset + i
                        im = msh%faces(IF)%master_()
                        rhs(4,im) = rhs(4,im) + phi_b(ibf) * msh%face_cntr(IF)%z_()
                    END DO
                END DO
            END IF

            ! Solves least squares problem
            DO ic = 1, ncells
                CALL msh%lsr(ic)%solve_least_squares(rhs(:,ic))
            END DO

            ! Sets grad values
            SELECT CASE(ncd)
            CASE(2)
                DO ic = 1, ncells
                    grad(ic) = vector_(rhs(2,ic),rhs(3,ic),0.d0)
                END DO
            CASE(3)
                DO ic = 1, ncells
                    grad(ic) = vector_(rhs(2,ic),rhs(3,ic),rhs(4,ic))
                END DO
            END SELECT

            ! Updates halo elements
            CALL update_vector_halo(grad,msh%desc_c)

            NULLIFY(if2b)
            NULLIFY(ic2c)
            DEALLOCATE(rhs)
            NULLIFY(msh)
            DEALLOCATE(phi_b,phi_x)

        100 FORMAT(' ERROR! Face-centered fields not supported by SCALAR_FIELD_GRAD')
        200 FORMAT(' ERROR! Memory allocation failure in SCALAR_FIELD_GRAD')

        END PROCEDURE  scalar_field_grad

END SUBMODULE scalar_field_grad_implementation
