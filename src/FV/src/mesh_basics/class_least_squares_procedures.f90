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
! $Id: class_least_squares.f90 3175 2008-06-13 12:59:07Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_least_squares) class_least_squares_procedures
    USE class_psblas
    USE class_vector
    USE class_face

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_least_squares_sizeof
        USE psb_base_mod

        nemo_least_squares_sizeof = nemo_sizeof_dp * SIZE(lsq%A)

    END PROCEDURE nemo_least_squares_sizeof


    ! ----- Constructor -----

    MODULE PROCEDURE alloc_least_squares
        !
        INTEGER :: i, ierr, info(n)

        info = 0

        ALLOCATE(lsr(n),stat=ierr)
        IF(ierr /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        SELECT CASE(ncd)
        CASE(2)
            DO i = 1, n
                ALLOCATE(lsr(i)%A(6),stat=ierr)
                info(i) = ierr
            END DO
        CASE(3)
            DO i = 1, n
                ALLOCATE(lsr(i)%A(10),stat=ierr)
                info(i) = ierr
            END DO
        END SELECT
        IF(ANY(info /= 0)) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF


        DO i = 1, n
            lsr(i)%A(:) = 0.d0
        ENDDO

100     FORMAT(' ERROR! Memory allocation failure in ALLOC_LEAST_SQUARES')

    END PROCEDURE alloc_least_squares


    ! ----- Destructor -----

    MODULE PROCEDURE free_least_squares
        !
        INTEGER :: i, ierr, info(SIZE(lsr)), n

        n = SIZE(lsr)
        info = 0

        IF (ALLOCATED(lsr)) THEN
            DO i = 1, n
                DEALLOCATE(lsr(i)%A,stat=info(i)) ! Not supported by GFortran!
!!$       deallocate(lsr(i)%A,stat=ierr)
!!$       info(i) = ierr
            END DO
            IF(ANY(info /= 0)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            DEALLOCATE(lsr,stat=ierr)
            IF(ierr /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
        END IF

100     FORMAT(' ERROR! Memory deallocation failure in FREE_LEAST_SQUARES')

    END PROCEDURE free_least_squares


    ! ----- Setter -----

    MODULE PROCEDURE set_least_squares
        USE class_connectivity
        USE tools_math
        USE psb_base_mod

        !
        INTEGER :: i, ib, ic, IF, im, n, nb, nbc, ncells
        INTEGER :: i11, i21, i31, i41
        INTEGER :: i22, i32, i42
        INTEGER :: i33, i43
        INTEGER :: i44
        INTEGER, POINTER :: ic2c(:) => NULL(), if2b(:) => NULL()

        CALL sw_lsr%tic()

        IF(mypnum_() == 0) THEN
            WRITE(*,*) 'Setting metrics for cell-centered least squares regression'
            WRITE(*,*)
        END IF

        ! LSR is sized according to the number of local + halo cells
        ncells = c2c%nel_()
        nbc    = f2b%nel_() - 2

        CALL alloc_least_squares(lsr,ncells,ncd)

        ! Next loops sweep only strictly local cells => NCELLS redefinition.
        ncells = psb_cd_get_local_rows(desc)

        ! Builds the lower part of the symmetric matrix. The points involved
        ! linked together in the basic stencil:
        ! - IC
        ! - IC's neighbots which share a fluid face
        ! - Center(s) of possible boundary face(s)

        ! Set indices for compact matrix storage
        i11 = 1; i21 = 2; i31 = 3
        SELECT CASE(ncd)
        CASE(2)
            i22 = 4; i32 = 5
            i33 = 6
        CASE(3)
            i41 = 4
            i22 = 5; i32 = 6; i42 = 7
            i33 = 8; i43 = 9
            i44 = 10
        END SELECT

        DO ic = 1, ncells
            ! IC
            lsr(ic)%A(i11) = 1.d0
            lsr(ic)%A(i21) = cell_cntr(ic)%x_()
            lsr(ic)%A(i31) = cell_cntr(ic)%y_()
            lsr(ic)%A(i22) = cell_cntr(ic)%x_() ** 2
            lsr(ic)%A(i32) = cell_cntr(ic)%x_() * cell_cntr(ic)%y_()
            lsr(ic)%A(i33) = cell_cntr(ic)%y_() ** 2

            ! IC's neighbors
            CALL c2c%get_ith_conn(ic2c,ic)
            n = SIZE(ic2c)
            lsr(ic)%A(i11) = lsr(ic)%A(i11) + n
            DO i = 1, n
                nb = ic2c(i)
                lsr(ic)%A(i21) = lsr(ic)%A(i21) + cell_cntr(nb)%x_()
                lsr(ic)%A(i31) = lsr(ic)%A(i31) + cell_cntr(nb)%y_()
                lsr(ic)%A(i22) = lsr(ic)%A(i22) + cell_cntr(nb)%x_() ** 2
                lsr(ic)%A(i32) = lsr(ic)%A(i32) + cell_cntr(nb)%x_() * cell_cntr(nb)%y_()
                lsr(ic)%A(i33) = lsr(ic)%A(i33) + cell_cntr(nb)%y_() ** 2
            END DO
        END DO

        ! Center(s) of possible boundary face(s)
        DO ib = 1, nbc
            CALL f2b%get_ith_conn(if2b,ib)
            n = SIZE(if2b)
            DO i = 1, n
                IF = if2b(i)
                im = faces(IF)%master_()
                lsr(im)%A(i11) = lsr(im)%A(i11) + 1.d0
                lsr(im)%A(i21) = lsr(im)%A(i21) + face_cntr(IF)%x_()
                lsr(im)%A(i31) = lsr(im)%A(i31) + face_cntr(IF)%y_()
                lsr(im)%A(i22) = lsr(im)%A(i22) + face_cntr(IF)%x_() ** 2
                lsr(im)%A(i32) = lsr(im)%A(i32) + face_cntr(IF)%x_() * face_cntr(IF)%y_()
                lsr(im)%A(i33) = lsr(im)%A(i33) + face_cntr(IF)%y_() ** 2
            END DO
        END DO


        ! Extra elements for 3D case only
        IF(ncd == 3) THEN
            DO ic = 1, ncells
                ! IC
                lsr(ic)%A(i41) = cell_cntr(ic)%z_()
                lsr(ic)%A(i42) = cell_cntr(ic)%x_() * cell_cntr(ic)%z_()
                lsr(ic)%A(i43) = cell_cntr(ic)%y_() * cell_cntr(ic)%z_()
                lsr(ic)%A(i44) = cell_cntr(ic)%z_() ** 2

                ! IC's neighbors
                CALL c2c%get_ith_conn(ic2c,ic)
                n = SIZE(ic2c)
                DO i = 1, n
                    nb = ic2c(i)
                    lsr(ic)%A(i41) = lsr(ic)%A(i41) + cell_cntr(nb)%z_()
                    lsr(ic)%A(i42) = lsr(ic)%A(i42) + cell_cntr(nb)%x_() * cell_cntr(nb)%z_()
                    lsr(ic)%A(i43) = lsr(ic)%A(i43) + cell_cntr(nb)%y_() * cell_cntr(nb)%z_()
                    lsr(ic)%A(i44) = lsr(ic)%A(i44) + cell_cntr(nb)%z_() ** 2
                END DO
            END DO

            ! Center(s) of possible boundary face(s)
            DO ib = 1, nbc
                CALL f2b%get_ith_conn(if2b,ib)
                n = SIZE(if2b)
                DO i = 1, n
                    IF = if2b(i)
                    im = faces(IF)%master_()
                    lsr(im)%A(i41) = lsr(im)%A(i41) + face_cntr(IF)%z_()
                    lsr(im)%A(i42) = lsr(im)%A(i42) + face_cntr(IF)%x_() * face_cntr(IF)%z_()
                    lsr(im)%A(i43) = lsr(im)%A(i43) + face_cntr(IF)%y_() * face_cntr(IF)%z_()
                    lsr(im)%A(i44) = lsr(im)%A(i44) + face_cntr(IF)%z_() ** 2
                END DO
            END DO
        END IF

        ! Factorize lsr%A with Cholesky decomposition
        DO ic = 1, ncells
            CALL factorize(lsr(ic)%A)
        END DO

        NULLIFY(ic2c)

        CALL sw_lsr%toc()

    END PROCEDURE set_least_squares


    ! ----- Computational Routines -----

    MODULE PROCEDURE solve_least_squares
        USE tools_math

        CALL sw_lsr%tic()

        ! WARNING!
        ! 2D -> size(RHS) = 3
        ! 3D -> size(RHS) = 4
        ! Check removed for efficiency reason

        CALL solve_sys(lsr%A,rhs)

        CALL sw_lsr%toc()

    END PROCEDURE solve_least_squares


END SUBMODULE class_least_squares_procedures
