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
! $Id: $
!
! Description:
!    Adds to PDE the contribution of the Divergence of PHI: div(FLUX*PHI)).
!    Remark: FLUX us face-centered
!
SUBMODULE (op_div) vector_pde_div_implementation
    IMPLICIT NONE

    CONTAINS

    MODULE PROCEDURE vector_pde_div
    USE class_psblas
    USE class_bc
    USE class_connectivity
    USE class_discretization
    USE class_dimensions
    USE class_face
    USE class_vector
    USE class_scalar_field
    USE class_vector_field
    USE class_mesh
    USE class_vector_pde
    USE tools_bc
    USE tools_operators
    USE op_field

    IMPLICIT NONE
    !
    CHARACTER(len=20), PARAMETER :: op_name = 'VECTOR_PDE_DIV'
    !
    INTEGER :: i, ib, ibf, ib_offset, IF, ifirst, info, ka, id_math
    INTEGER :: im, im_glob, is, is_glob
    INTEGER :: n, ncoeff, nel, nmax, nbc
    INTEGER, ALLOCATABLE :: irow_a(:), icol_a(:), iloc_to_glob(:)
    INTEGER, POINTER :: if2b(:) => NULL()
    REAL(psb_dpk_) :: coeff_a(4), fsign, side_, w, fact
    TYPE(bc_poly), POINTER :: bc(:) => NULL()
    TYPE(dimensions) :: dim
    TYPE(discretization) :: ds_
    !
    REAL(psb_dpk_), ALLOCATABLE :: A(:), flux_x(:)
    REAL(psb_dpk_), ALLOCATABLE :: bc_a(:), bc_b(:)
    REAL(psb_dpk_), ALLOCATABLE :: sp_bc(:,:)
    TYPE(mesh), POINTER :: msh => NULL()
    TYPE(mesh), POINTER :: msh_phi => NULL(), msh_flux => NULL()
    TYPE(vector), ALLOCATABLE :: bc_c(:), sc_bc(:,:), b(:)



    CALL tic(sw_pde)

    IF(mypnum_() == 0) THEN
        WRITE(*,*) '* ', TRIM(name_(pde)), ': applying the Divergence ',&
            & 'operator to the ', TRIM(name_(phi)), ' field'
    END IF

    ! Possible reinit of PDE
    CALL reinit_pde(pde)

    ! Is PHI cell-centered?
    IF(on_faces_(phi)) THEN
        WRITE(*,100) TRIM(op_name)
        CALL abort_psblas
    END IF

    ! Is FLUX face-centered?
    IF(.NOT.on_faces_(flux)) THEN
        WRITE(*,200) TRIM(op_name)
        CALL abort_psblas
    END IF

    ! Checks mesh consistency PDE vs. PHI
    CALL pde%get_mesh(msh)
    CALL phi%get_mesh(msh_phi)
    BLOCK
        !! This is in a block statement due to gfortran-8.3 issue w/ not being able to find this interface
        USE class_mesh, ONLY : check_mesh_consistency
        CALL check_mesh_consistency(msh,msh_phi,op_name)
    END BLOCK

    ! Checks mesh consistency PHI vs. FLUX
    CALL flux%get_mesh(msh_flux)
    BLOCK
        !! This is in a block statement due to gfortran-8.3 issue w/ not being able to find this interface
        USE class_mesh, ONLY : check_mesh_consistency
        CALL check_mesh_consistency(msh_phi,msh_flux,op_name)
    END BLOCK

    NULLIFY(msh_flux)
    NULLIFY(msh_phi)

    ! Equation dimensional check
    dim = dim_(flux) * dim_(phi)
    IF(dim /= dim_(pde)) THEN
        WRITE(*,300) TRIM(op_name)
        CALL debug_dim(dim_(pde))
        CALL debug_dim(dim)
        CALL debug_dim(dim_(flux))
        CALL debug_dim(dim_(phi))
        CALL abort_psblas
    END IF

    ! Computes sign factor
    IF(PRESENT(side)) THEN
        side_ = side
    ELSE
        side_ = lhs_ ! Default = LHS
    END IF
    fsign = pde_sign(sign,side_)

    ! Set discretization scheme
    IF(PRESENT(ds)) THEN
        ds_ = ds
    ELSE
        ds_ = cd_ ! Default central differencing
    END IF

    ! Gets FLUX "x" internal values
    CALL get_x(flux,flux_x)

    ! Total number of matrix coefficients associated to the fluid faces
    ncoeff = 4 * COUNT(flag_(msh%faces) == 0)

    ! Computes maximum size of blocks to be inserted
    nmax = size_blk(1,ncoeff)

    ALLOCATE(A(nmax),b(nmax),irow_a(nmax),icol_a(nmax),&
        & stat=info)
    IF(info /= 0) THEN
        WRITE(*,400) TRIM(op_name)
        CALL abort_psblas
    END IF

    ! Gets local to global conversion list
    CALL psb_get_loc_to_glob(msh%desc_c,iloc_to_glob)

    A = 0.d0
    b = vector_(0.d0,0.d0,0.d0)
    irow_a = 0
    icol_a = 0

    CALL get_ith_conn(if2b,msh%f2b,0)

    ifirst = 1; i = 0
    insert_fluid: DO
        IF(ifirst > ncoeff) EXIT insert_fluid
        nel = size_blk(ifirst,ncoeff)

        IF (id_(ds_) == id_(up_)) THEN
            block_fluid_upwind: DO ka = 1, nel, 4
                ! Local indices
                i = i + 1
                IF = if2b(i)
                im = master_(msh%faces(IF))
                is = slave_(msh%faces(IF))
                ! Global indices
                im_glob = iloc_to_glob(im)
                is_glob = iloc_to_glob(is)

                ! --- Sparse Matrix A ---

                ! Gets A's coefficients
                coeff_a(1) = fsign * MAX(flux_x(i),0.d0)
                coeff_a(2) = fsign * MIN(flux_x(i),0.d0)
                coeff_a(3) = -coeff_a(2)
                coeff_a(4) = -coeff_a(1)

                ! Diagonal coefficients
                ! Master
                A(ka)  = coeff_a(1)
                irow_a(ka) = im_glob
                icol_a(ka) = im_glob
                !
                ! Slave
                A(ka+1)  = coeff_a(2)
                irow_a(ka+1) = im_glob
                icol_a(ka+1) = is_glob

                ! Off-diagonal coefficients
                ! Master row
                A(ka+2)  = coeff_a(3)
                irow_a(ka+2) = is_glob
                icol_a(ka+2) = is_glob
                !
                ! Slave row
                A(ka+3)  = coeff_a(4)
                irow_a(ka+3) = is_glob
                icol_a(ka+3) = im_glob

            END DO block_fluid_upwind

        ELSE IF (id_(ds_) == id_(cd_)) THEN
            block_fluid_cd: DO ka = 1, nel, 4
                ! Local indices
                i = i + 1
                IF = if2b(i)
                im = master_(msh%faces(IF))
                is = slave_(msh%faces(IF))
                ! Global indices
                im_glob = iloc_to_glob(im)
                is_glob = iloc_to_glob(is)

                ! --- Sparse Matrix A ---

                ! Gets A's coefficients

                w  = msh%interp(IF)
                fact = fsign * flux_x(i)

                coeff_a(1) = fact * (1-w)
                coeff_a(2) = fact * w
                coeff_a(3) = -coeff_a(1)
                coeff_a(4) = -coeff_a(2)
!!$        write(0,*) 'From vector_pde_div: scale for A',i,im_glob,is_glob,&
!!$             & coeff_a(1:4),fsign, flux_x(i),msh%interp(if)

                ! Diagonal coefficients
                ! Master
                A(ka)  = coeff_a(1)
                irow_a(ka) = im_glob
                icol_a(ka) = im_glob
                !
                ! Slave
                A(ka+1)  = coeff_a(2)
                irow_a(ka+1) = im_glob
                icol_a(ka+1) = is_glob

                ! Off-diagonal coefficients
                ! Master row
                A(ka+2)  = coeff_a(3)
                irow_a(ka+2) = is_glob
                icol_a(ka+2) = is_glob
                !
                ! Slave row
                A(ka+3)  = coeff_a(4)
                irow_a(ka+3) = is_glob
                icol_a(ka+3) = im_glob

            END DO block_fluid_cd
        ELSE
            WRITE(*,401) TRIM(op_name)
            CALL abort_psblas

        ENDIF
        CALL spins_pde(nel,irow_a,icol_a,A,pde)
        IF (debug_mat_bld) THEN
            WRITE(0,*) 'From vector_pde_div FLUID: SPINS ',nel
            DO ka=1,nel
                WRITE(0,*) irow_a(ka),icol_a(ka),a(ka)
            END DO
        END IF
        ifirst = ifirst + nel
    END DO insert_fluid

    ! Free storage
    DEALLOCATE(flux_x)

    ! ----- Insert source terms for boundary conditions -----

    ! Gets PHI boundary conditions
    bc => bc_(phi)
    ! Gets boundary values of FLUX fields
    CALL get_bx(flux,flux_x)

    A = 0.d0
    b = vector_(0.d0,0.d0,0.d0)
    irow_a = 0
    icol_a = 0

    nbc = msh%nbc

    nmax = max_conn(msh%f2b,lb=1)
    ALLOCATE(bc_a(nmax),bc_b(nmax),bc_c(nmax),&
        & sc_bc(nmax,nbc), sp_bc(nmax,nbc),stat=info)
    IF(info /= 0) THEN
        WRITE(*,400) TRIM(op_name)
        CALL abort_psblas
    END IF
    ! WARNING: BC_A, BC_B, BC_C are oversized.

    ! Initialization
    bc_a = 0.d0
    bc_b = 0.d0
    bc_c = vector_(0.d0,0.d0,0.d0)
    sc_bc = vector_(0.d0,0.d0,0.d0)
    sp_bc = 0.d0

    ! Boundary faces with flag < IB
    ib_offset = 0

    bc_loop: DO ib = 1, nbc
        CALL get_ith_conn(if2b,msh%f2b,ib)
        n = SIZE(if2b)
        ! Gets analytical boundary conditions and builds source terms
        IF(n == 0) CYCLE bc_loop

        CALL get_abc(bc(ib),dim_(phi),id_math,bc_a(1:n),bc_b(1:n),bc_c(1:n))

        SELECT CASE(id_math)
        CASE(bc_dirichlet_, bc_dirichlet_map_)
            ! Dirichlet
            DO i = 1, n
                IF = if2b(i)
                ibf = ib_offset + i
                fact = fsign * flux_x(ibf)
                sc_bc(i,ib) = fact * bc_c(i)
                sp_bc(i,ib) = 0.d0
            END DO
        CASE(bc_neumann_)
            ! Neumann, gradient
            DO i = 1, n
                IF = if2b(i)
                ibf = ib_offset + i
                IF (id_(ds_) == id_(up_)) THEN
                    sc_bc(i,ib) = fsign * MIN(flux_x(ibf),0.d0) * msh%dist(IF) * bc_c(i)
                    sp_bc(i,ib) = fsign * flux_x(ibf)
                ELSE
                    fact = fsign * flux_x(ibf)
                    sc_bc(i,ib) = fact *  msh%dist(IF) * bc_c(i)
                    sp_bc(i,ib) = fact
                ENDIF
            END DO
        CASE default
            WRITE(*,402)
            CALL abort_psblas
        END SELECT
        ib_offset = ib_offset + n
    END DO bc_loop

    DEALLOCATE(bc_a,bc_b,bc_c)

    DO ib = 1, nbc
        CALL get_ith_conn(if2b,msh%f2b,ib)
        n = SIZE(if2b)

        ifirst = 1; i = 0
        insert_boundary: DO
            IF(ifirst > n) EXIT insert_boundary
            nel = size_blk(ifirst,n)

            block_boundary: DO ka = 1, nel
                ! Local index
                i = i + 1
                IF = if2b(i)
                im = master_(msh%faces(IF))
                ! Global index
                im_glob = iloc_to_glob(im)

                A(ka) = sp_bc(i,ib)
                b(ka) = -sc_bc(i,ib)
                irow_a(ka) = im_glob
                icol_a(ka) = im_glob
            END DO block_boundary


            CALL spins_pde(nel,irow_a,icol_a,A,pde)
            CALL geins_pde(nel,irow_a,b,pde)
            IF (debug_mat_bld) THEN
                WRITE(0,*) 'From vector_pde_div BC: SPINS ',nel
                DO ka=1,nel
                    WRITE(0,*) irow_a(ka),a(ka)
                END DO
                WRITE(0,*) 'From vector_pde_div BC: geins ',nel
                DO ka=1,nel
                    WRITE(0,*) irow_a(ka),x_(b(ka)),y_(b(ka)),z_(b(ka))
                END DO
            END IF
            ifirst = ifirst + nel

        END DO insert_boundary
    END DO

    ! Free storage
    DEALLOCATE(sc_bc,sp_bc)
    NULLIFY(if2b)
    DEALLOCATE(flux_x)
    DEALLOCATE(iloc_to_glob)
    DEALLOCATE(A,b,irow_a,icol_a)
    NULLIFY(msh)
    NULLIFY(bc)

    CALL toc(sw_pde)

100 FORMAT(' ERROR! PHI field in ',a,' is not cell centered')
200 FORMAT(' ERROR! FLUX field in ',a,' is cell centered')
300 FORMAT(' ERROR! Dimensional check failure in ',a)
400 FORMAT(' ERROR! Memory allocation failure in ',a)
401 FORMAT(' ERROR! Unimplemented discretization scheme in ',a)
402 FORMAT(' ERROR! Unimplemented BC  in ',a)

    END PROCEDURE vector_pde_div

END SUBMODULE vector_pde_div_implementation
