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
! $Id: vector_pde_laplacian.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Adds to PDE the contribution of the Laplacian of PHI: div(GAMMA*grad(PHI)).
!    Remark: the diffusivity GAMMA (optional) is face-centered.
!
SUBMODULE (op_laplacian) vector_pde_laplacian_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE vector_pde_laplacian
        USE class_psblas
        USE class_bc
        USE class_connectivity
        USE class_dimensions
        USE class_face
        USE class_scalar_field
        USE class_vector_field
        USE class_mesh
        USE class_vector_pde
        USE class_vector
        USE op_grad
        USE tools_bc
        USE tools_operators

        IMPLICIT NONE
        !
        CHARACTER(len=30) :: op_name = 'VECTOR_PDE_LAPLACIAN'
        INTEGER :: i, ib, ibf, id_math, IF, ifirst, info, ib_offset, k
        INTEGER ::  im, im_glob, is, is_glob
        INTEGER :: n, nbc, ncoeff, nel, nfaces, nf_boundary, nf_fluid, nmax
        INTEGER, POINTER :: if2b(:) => NULL()
        INTEGER, ALLOCATABLE :: ia(:), ja(:)
        INTEGER, ALLOCATABLE :: iloc_to_glob(:)
        REAL(psb_dpk_) :: a_nb, fact, fsign, side_, r, w
        REAL(psb_dpk_), ALLOCATABLE :: A(:)
        REAL(psb_dpk_), ALLOCATABLE :: bc_a(:), bc_b(:)
        REAL(psb_dpk_), ALLOCATABLE :: sp_bc(:,:)
        REAL(psb_dpk_), ALLOCATABLE :: gamma_x(:), gamma_bx(:)
        TYPE(bc_poly), POINTER :: bc(:) => NULL()
        TYPE(dimensions) :: dim, dim_temp
        TYPE(mesh), POINTER :: msh => NULL()
        TYPE(mesh), POINTER :: msh_phi => NULL(), msh_gamma => NULL()
        TYPE(vector) :: corr, delta_rm, delta_rs, nf, proj_m, proj_s, proj_temp
        TYPE(vector) :: phi_proj_m, phi_proj_s, s_no, s_no1, s_no2
        TYPE(vector), ALLOCATABLE :: bc_c(:), grad(:,:), sc_bc(:,:), b(:)
        TYPE(vector), ALLOCATABLE ::  phi_x(:)

        CALL sw_pde%tic()

        IF(mypnum_() == 0) THEN
            WRITE(*,*) '* ', TRIM(pde%name_()), ': applying the Laplacian ',&
                & 'operator to the ', TRIM(phi%name_()), ' field'
        END IF

        ! Possible reinit of PDE
        CALL pde%reinit_pde()

        ! Is PHI cell-centered?
        IF(phi%on_faces_()) THEN
            WRITE(*,100) TRIM(op_name)
            CALL abort_psblas
        END IF

        ! Is GAMMA face-centered?
        IF(PRESENT(gamma)) THEN
            IF(.NOT.gamma%on_faces_()) THEN
                WRITE(*,200) TRIM(op_name)
                CALL abort_psblas
            END IF
        END IF

        ! Checks mesh consistency PDE vs. PHI
        CALL pde%get_mesh(msh)
        CALL phi%get_mesh(msh_phi)
        BLOCK
            !! This is in a block statement due to gfortran-8.3 issue w/ not being able to find this interface
            USE class_mesh, ONLY : check_mesh_consistency
            CALL check_mesh_consistency(msh,msh_phi,op_name)
        END BLOCK

        ! Checks mesh consistency PHI vs. GAMMA
        IF(PRESENT(gamma)) THEN
            CALL gamma%get_mesh(msh_gamma)
            BLOCK
                !! This is in a block statement due to gfortran-8.3 issue w/ not being able to find this interface
                USE class_mesh, ONLY : check_mesh_consistency
                CALL check_mesh_consistency(msh_phi,msh_gamma,op_name)
            END BLOCK
        END IF

        NULLIFY(msh_gamma)
        NULLIFY(msh_phi)

        ! Equation dimensional check
        dim = phi%dim_() * volume_ / (length_ * length_)
        IF(PRESENT(gamma)) dim = gamma%dim_() * dim
        IF(dim /= pde%dim_()) THEN
            CALL dim%debug_dim()
            dim_temp=pde%dim_()
            !CALL pde%dim_()%debug_dim()
            CALL dim_temp%debug_dim()
            WRITE(*,300)  TRIM(op_name)
            CALL abort_psblas
        END IF

        ! Computes sign factor
        IF(PRESENT(side)) THEN
            side_ = side
        ELSE
            side_ = lhs_ ! Default = LHS
        END IF
        fsign = - pde_sign(sign,side_) ! WARNING it's MINUS!


        ! Gets PHI boundary conditions
        bc => phi%bc_()

        ! Gets PHI "x" internal values
        CALL phi%get_x(phi_x)

        IF(PRESENT(gamma)) THEN
            CALL gamma%get_x(gamma_x)
            CALL gamma%get_bx(gamma_bx)
        ELSE
            nf_fluid    = COUNT(msh%faces%flag_() <= 0)
            nf_boundary = COUNT(msh%faces%flag_() > 0)
            ALLOCATE(gamma_x(nf_fluid),gamma_bx(nf_boundary),stat=info)
            IF(info /= 0) THEN
                WRITE(*,400)  TRIM(op_name)
                CALL abort_psblas
            END IF
            gamma_x(:)  = 1.d0
            gamma_bx(:) = 1.d0
        END IF

        ! Computes gradient for non-orthogonality correction
        ! First allocates the result ...
        ALLOCATE(grad(3,SIZE(phi_x)),stat=info)
        IF(info /= 0) THEN
            WRITE(*,400) TRIM(op_name)
            CALL abort_psblas
        END IF

        ! ... then computes it.
        grad = fld_grad(phi)
        ! QUESTION: is it possible to avoid the explicit allocation and,
        ! instead, allocate automatically the result by means of the
        ! assignment? What about move_alloc and realloc?
        ! REMARK: svn diff -r 356:357 scalar_pde_laplacian.f90
        ! - Intel does not support assignment w/o preliminary allocation.
        ! - Gfortran does it.

        ! Total number of matrix coefficients associated to the fluid faces
        ncoeff = 4 * COUNT(msh%faces%flag_() == 0)

        ! Computes maximum size of blocks to be inserted
        nmax = size_blk(1,ncoeff)

        ALLOCATE(A(nmax),b(nmax),ia(nmax),ja(nmax),stat=info)
        IF(info /= 0) THEN
            WRITE(*,400) TRIM(op_name)
            CALL abort_psblas
        END IF

        ! Gets local to global conversion list
        CALL psb_get_loc_to_glob(msh%desc_c,iloc_to_glob)

        ! ----- Insert contributions for fluid faces -----
        A = 0.d0
        ia = 0
        ja = 0

        CALL msh%f2b%get_ith_conn(if2b,0)

        ifirst = 1; i = 0
        insert_fluid: DO
            IF(ifirst > ncoeff) EXIT insert_fluid
            nel = size_blk(ifirst,ncoeff)

            block_fluid: DO k = 1, nel, 4
                ! Local indices
                i = i + 1
                IF = if2b(i)
                im = msh%faces(IF)%master_()
                is = msh%faces(IF)%slave_()
                ! Global indices
                im_glob = iloc_to_glob(im)
                is_glob = iloc_to_glob(is)

                w = msh%interp(IF)

                ! A possible non-conjunctionality correction for gamma is
                ! applied in the update_field routine.

                a_nb = fsign * gamma_x(i) * msh%area(IF) / msh%dist(IF)

                ! WARNING!
                ! GAMMA_X first stores field values for faces with flag == 0
                ! The same subscript I is used for accessing IF2B and GAMMA_X.

                ! ----- Coefficients -----

                ! Diagonal coefficients
                ! Master
                A(k)  = a_nb
                ia(k) = im_glob
                ja(k) = im_glob
                !
                ! Slave
                A(k+1)  = a_nb
                ia(k+1) = is_glob
                ja(k+1) = is_glob

                ! Off-diagonal coefficients
                ! Master row
                A(k+2) = -a_nb
                ia(k+2) = im_glob
                ja(k+2) = is_glob
                !
                ! Slave row
                A(k+3) = -a_nb
                ia(k+3) = is_glob
                ja(k+3) = im_glob
            END DO block_fluid

            CALL pde%spins_pde(nel,ia,ja,A)
            IF (debug_mat_bld) THEN
                WRITE(0,*) 'From vector_pde_laplacian FLUID: SPINS ',nel
                DO k=1,nel
                    WRITE(0,*) ia(k),ja(k),a(k)
                END DO
            END IF

            ifirst = ifirst + nel
        END DO insert_fluid


        ! ----- Insert source terms for non-orthogonality correction -----
        b = vector_(0.d0,0.d0,0.d0)
        ia = 0

        ! Total number of RHS coefficients associated to the correction
        ! of the non-orthogonality.
        ncoeff = 2 * COUNT(msh%faces%flag_() == 0)

        ifirst = 1; i = 0
        insert_nonortho: DO
            IF(ifirst > ncoeff) EXIT insert_nonortho
            nel = size_blk(ifirst,ncoeff)

            block_nonortho: DO k = 1, nel, 2
                ! Local indices
                i = i + 1
                IF = if2b(i)
                im = msh%faces(IF)%master_()
                is = msh%faces(IF)%slave_()
                ! Global indices
                im_glob = iloc_to_glob(im)
                is_glob = iloc_to_glob(is)

                r= 1.d0 / msh%area(IF)
                nf = r * msh%af(IF)

                ! Computes S_NO1 term
                ! Finds PROJ_P and PROJ_NB, projections of P and NB on the
                ! straight line normal to the face and passing by its center
                r = (msh%cell_cntr(im) - msh%face_cntr(IF)) .dot. nf
                proj_m = msh%face_cntr(IF) + r * nf

                r = (msh%cell_cntr(is) - msh%face_cntr(IF)) .dot. nf
                proj_s = msh%face_cntr(IF) + r * nf

                delta_rm = proj_m - msh%cell_cntr(im)
                delta_rs = proj_s - msh%cell_cntr(is)

                ! Interpolates temp value from cell-center to projection point
                phi_proj_m = phi_x(im) + (grad(:,im) .dot. delta_rm)
                phi_proj_s = phi_x(is) + (grad(:,is) .dot. delta_rs)

                ! Distance between projection points
                proj_temp = proj_s - proj_m
                !r = 1.d0 / (proj_s - proj_m)%mag()
                r = 1.d0 / proj_temp%mag()
                s_no1 = r * (phi_proj_s - phi_proj_m)

                ! Computes S_NO2 term
                r = 1.d0 / msh%dist(IF)
                s_no2 = r * (phi_x(is) - phi_x(im))

                corr = s_no1 - s_no2
                s_no = fsign * gamma_x(i) * msh%area(IF) * corr

                ! WARNING!
                ! GAMMA_X first stores field values for faces with flag == 0
                ! The same subscript I is used for accessing IF2B and GAMMA_X.


                ! ----- RHS Coefficients -----

                ! Master
                b(k) = s_no
                ia(k) = im_glob

                ! Slave
                b(k+1) = (-1.d0) * s_no ! WARNING "-"
                ia(k+1) = is_glob

            END DO block_nonortho

            CALL pde%geins_pde(nel,ia,b)
            IF (debug_mat_bld) THEN
                WRITE(0,*) 'From vector_pde_laplacian NON_ORTH_C: geins ',nel
                DO k=1,nel
                    WRITE(0,*) ia(k),b(k)%x_(),b(k)%y_(),b(k)%z_()
                END DO
            END IF
            ifirst = ifirst + nel

        END DO insert_nonortho


        ! ----- Insert source terms for boundary conditions -----

        A = 0.d0
        b = vector_(0.d0,0.d0,0.d0)
        ia = 0
        ja = 0
        nbc = msh%nbc
        nf_fluid = COUNT(msh%faces%flag_() <= 0)
        nfaces   = SIZE(msh%faces)

        nmax = msh%f2b%max_conn(lb=1)
        ALLOCATE(bc_a(nmax),bc_b(nmax),bc_c(nmax),&
            &   sc_bc(nmax,nbc),sp_bc(nmax,nbc),stat=info)
        IF(info /= 0) THEN
            WRITE(*,400)
            CALL abort_psblas
        END IF
        ! WARNING: BC_A, BC_B, BC_C are oversized.

        ! Initialization
        bc_a = 0.d0
        bc_b = 0.d0
        bc_c = vector_(0.d0,0.d0,0.d0)
        sc_bc = vector_(0.d0,0.d0,0.d0)
        sp_bc = 0.d0

        bc_loop: DO ib = 1, nbc
            CALL msh%f2b%get_ith_conn(if2b,ib)
            n = SIZE(if2b)
            ! Gets analytical boundary conditions and builds source terms
            IF (n == 0) CYCLE bc_loop

            CALL bc(ib)%get_abc(phi%dim_(),id_math,bc_a(1:n),bc_b(1:n),bc_c(1:n))

            ! Boundary faces with flag < IB
            ib_offset = COUNT(msh%faces%flag_() > 0 .AND. msh%faces%flag_() < ib)

            SELECT CASE(id_math)
            CASE(bc_dirichlet_, bc_dirichlet_map_)
                ! Dirichlet
                DO i = 1, n
                    IF = if2b(i)
                    ibf = ib_offset + i
                    fact = fsign * gamma_bx(ibf) * msh%area(IF) / msh%dist(IF)
                    sp_bc(i,ib) = fact
                    sc_bc(i,ib) = fact * bc_c(i)
    !!$        write(0,*) 'Vector pde laplacian',i, if, ibf,&
    !!$             & gamma_bx(ibf), msh%area(if), msh%dist(if)
                END DO

            CASE(bc_neumann_)
                ! Neumann, gradient
                DO i = 1, n
                    IF = if2b(i)
                    ibf = ib_offset + i
                    fact = fsign * gamma_bx(ibf) * msh%area(IF)
                    sp_bc(i,ib) = 0.d0
                    sc_bc(i,ib) = fact * bc_c(i)
                END DO

            CASE(bc_neumann_flux_)
                ! Neumann, flux
                DO i = 1, n
                    IF = if2b(i)
                    fact = fsign * msh%area(IF)
                    sp_bc(i,ib) = 0.d0
                    sc_bc(i,ib) = fact * bc_c(i)
                END DO

            CASE(bc_robin_, bc_robin_convection_)
                ! Robin
                DO i = 1, n
                    IF = if2b(i)
                    ibf = ib_offset + i
                    fact = bc_a(i) * msh%dist(IF) + bc_b(i)  !gamma_bx(ibf)
                    fact = fsign * msh%area(IF) * gamma_bx(ibf) / fact
                    sp_bc(i,ib) = fact * bc_a(i)
                    sc_bc(i,ib) = fact * bc_c(i)
                END DO

            CASE default
                WRITE(*,500) TRIM(op_name)
                CALL abort_psblas
            END SELECT
        END DO bc_loop

        DEALLOCATE(bc_a,bc_b,bc_c)

        DO ib = 1, nbc
            CALL msh%f2b%get_ith_conn(if2b,ib)
            n = SIZE(if2b)

            ifirst = 1; i = 0
            insert_boundary: DO
                IF(ifirst > n) EXIT insert_boundary
                nel = size_blk(ifirst,n)

                block_boundary: DO k = 1, nel
                    ! Local index
                    i = i + 1
                    IF = if2b(i)
                    im = msh%faces(IF)%master_()
                    ! Global index
                    im_glob = iloc_to_glob(im)

                    A(k) = sp_bc(i,ib)
                    b(k) = sc_bc(i,ib)

                    ia(k) = im_glob
                    ja(k) = im_glob
                END DO block_boundary

                CALL pde%spins_pde(nel,ia,ja,A)
                CALL pde%geins_pde(nel,ia,b)
                IF (debug_mat_bld) THEN
                    WRITE(0,*) 'From vector_pde_laplacian BC: SPINS ',nel
                    DO k=1,nel
                        WRITE(0,*) ia(k),ja(k),a(k)
                    END DO
                    WRITE(0,*) 'From vector_pde_laplacian BC: geins ',nel
                    DO k=1,nel
                        WRITE(0,*) ia(k),b(k)%x_(),b(k)%y_(),b(k)%z_()
                    END DO
                END IF

                ifirst = ifirst + nel

            END DO insert_boundary
        END DO



        ! Frees memory storage
        DEALLOCATE(sc_bc,sp_bc)
        NULLIFY(if2b)
        DEALLOCATE(iloc_to_glob)
        DEALLOCATE(A,b,ia,ja)
        DEALLOCATE(grad)
        DEALLOCATE(gamma_bx,gamma_x)
        DEALLOCATE(phi_x)
        NULLIFY(msh)
        NULLIFY(bc)

        CALL sw_pde%toc()

    100 FORMAT(' ERROR! Unknown PHI in ',a,' is face-centered')
    200 FORMAT(' ERROR! Diffusivity GAMMA in ',a,' is cell-centered')
    300 FORMAT(' ERROR! Dimensional check failure in ',a)
    400 FORMAT(' ERROR! Memory allocation failure in ',a)
    500 FORMAT(' ERROR! Unsupported BC type in ',a)

    END PROCEDURE vector_pde_laplacian

END SUBMODULE vector_pde_laplacian_implementation
