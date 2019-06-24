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
! $Id: vector_field_grad.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Computes the gradient of a VectorField object.
!    The result is a VECTOR object with rank 1.
!
!    REMARK: remember to allocate the LHS before assigning the function result
!
SUBMODULE(op_grad) vector_field_grad_implementation

    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE vector_field_grad
            USE class_psblas
            USE class_connectivity
            USE class_face
            USE class_least_squares
            USE class_mesh
            USE class_vector_field
            USE class_vector

            IMPLICIT NONE
            !
            INTEGER :: i, ib, ibf, ic, icoo, IF, im, info, ib_offset
            INTEGER :: n, nb, nbc, ncd, ncells
            INTEGER, POINTER :: ic2c(:) => NULL(), if2b(:) => NULL()
            REAL(psb_dpk_) :: xc, yc, zc, xf, yf, zf, xnb, ynb, znb
            REAL(psb_dpk_), ALLOCATABLE :: phi_x(:,:)
            REAL(psb_dpk_), ALLOCATABLE :: phi_b(:,:)
            REAL(psb_dpk_), ALLOCATABLE :: rhs(:,:,:)
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
            ncells = SIZE(phi_x,1)

            IF(.NOT.ALLOCATED(grad)) THEN
                ALLOCATE(grad(3,ncells),stat=info)
                IF(info /= 0) THEN
                    WRITE(*,200)
                    CALL abort_psblas
                END IF
            END IF
            ! REMARK: compiler-dependent behaviour in SCALAR_FIELD_GRAD
            ! Intel:   allocated(grad) = F
            ! Gfortran: allocated(grad) = T

            ! Gets PHI mesh
            CALL phi%get_mesh(msh)
            ncd = msh%ncd

            ! Allocates RHS for least squares regression
            ALLOCATE(rhs(ncd+1,3,ncells),stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            ! Initialization
            grad(:,:) = vector_(0.d0,0.d0,0.d0)

            ! Next loops sweep only strictly local cells => NCELLS redefinition.
            ncells = psb_cd_get_local_rows(msh%desc_c)
            nbc = msh%nbc

            rhs = 0.d0


            ! ----- Builds RHS -----

            DO ic = 1, ncells

                ! --- IC ---
                xc = msh%cell_cntr(ic)%x_()
                yc = msh%cell_cntr(ic)%y_()

                DO icoo = 1, 3
                    rhs(1,icoo,ic) = phi_x(ic,icoo)
                    rhs(2,icoo,ic) = phi_x(ic,icoo) * xc
                    rhs(3,icoo,ic) = phi_x(ic,icoo) * yc
                END DO

                ! --- IC's neighbors ---
                CALL msh%c2c%get_ith_conn(ic2c,ic)
                n = SIZE(ic2c)

                DO i = 1, n
                    nb = ic2c(i)
                    xnb = msh%cell_cntr(nb)%x_()
                    ynb = msh%cell_cntr(nb)%y_()

                    DO icoo = 1, 3
                        rhs(1,icoo,ic) = rhs(1,icoo,ic) + phi_x(nb,icoo)
                        rhs(2,icoo,ic) = rhs(2,icoo,ic) + phi_x(nb,icoo) * xnb
                        rhs(3,icoo,ic) = rhs(3,icoo,ic) + phi_x(nb,icoo) * ynb
                    END DO
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

                    xf = msh%face_cntr(IF)%x_()
                    yf = msh%face_cntr(IF)%y_()

                    DO icoo = 1, 3
                        rhs(1,icoo,im) = rhs(1,icoo,im) + phi_b(ibf,icoo)
                        rhs(2,icoo,im) = rhs(2,icoo,im) + phi_b(ibf,icoo) * xf
                        rhs(3,icoo,im) = rhs(3,icoo,im) + phi_b(ibf,icoo) * yf
                    END DO
                END DO
            END DO


            ! Extra elements for 3D case only
            IF(ncd == 3) THEN

                DO ic = 1, ncells

                    ! --- IC ---
                    zc = msh%cell_cntr(ic)%z_()

                    DO icoo = 1, 3
                        rhs(4,icoo,ic) = phi_x(ic,icoo) * zc
                    END DO

                    ! --- IC's neighbors ---
                    CALL msh%c2c%get_ith_conn(ic2c,ic)
                    n = SIZE(ic2c)

                    DO i = 1, n
                        nb = ic2c(i)
                        znb = msh%cell_cntr(nb)%z_()

                        DO icoo = 1, 3
                            rhs(4,icoo,ic) = rhs(4,icoo,ic) + phi_x(nb,icoo) * znb
                        END DO
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

                        zf = msh%face_cntr(IF)%z_()

                        DO icoo = 1, 3
                            rhs(4,icoo,im) = rhs(4,icoo,im) + phi_b(ibf,icoo) * zf
                        END DO
                    END DO
                END DO

            END IF

            ! Solves least squares problem
            DO ic = 1, ncells
                DO icoo = 1, 3
                    CALL msh%lsr(ic)%solve_least_squares(rhs(:,icoo,ic))
                END DO
            END DO

            ! Sets grad values
            SELECT CASE(ncd)
            CASE(2)
                DO ic = 1, ncells
                    DO icoo = 1, 3
                        grad(icoo,ic) = vector_(rhs(2,icoo,ic),rhs(3,icoo,ic),0.d0)
                    END DO
                END DO
            CASE(3)
                DO ic = 1, ncells
                    DO icoo = 1, 3
                        grad(icoo,ic) = vector_(rhs(2,icoo,ic),rhs(3,icoo,ic),rhs(4,icoo,ic))
                    END DO
                END DO
            END SELECT


            ! Updates halo elements
            DO icoo = 1, 3
                CALL update_vector_halo(grad(icoo,:),msh%desc_c)
            END DO


            NULLIFY(if2b)
            NULLIFY(ic2c)
            DEALLOCATE(rhs)
            NULLIFY(msh)
            DEALLOCATE(phi_b,phi_x)

        100 FORMAT(' ERROR! Face-centered fields not supported by VECTOR_FIELD_GRAD')
        200 FORMAT(' ERROR! Memory allocation failure in VECTOR_FIELD_GRAD')

        END PROCEDURE vector_field_grad

END SUBMODULE vector_field_grad_implementation
