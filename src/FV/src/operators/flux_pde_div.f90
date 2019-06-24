!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
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
! $Id: flux_pde_div.f90 3878 2010-01-04 12:26:01Z sfilippo $
!
! Description:
!    Adds to PDE the contribution of the Divergence of FLUX as a source term.
SUBMODULE (op_div) flux_pde_div_implementation
    IMPLICIT NONE

    CONTAINS

    MODULE PROCEDURE flux_pde_div
    USE class_psblas
    USE class_connectivity
    USE class_discretization
    USE class_dimensions
    USE class_face
    USE class_vector
    USE class_scalar_field
    USE class_vector_field
    USE class_mesh, ONLY : mesh
    USE class_scalar_pde
    USE tools_operators
    USE op_field

    IMPLICIT NONE
    !
    CHARACTER(len=*), PARAMETER :: op_name = 'PDE_DIV_AS_SOURCE_TERM'
    !
    INTEGER :: i, ib, ibf, ib_offset, IF, ifirst, info, ka
    INTEGER :: im, im_glob, is, is_glob
    INTEGER :: n, ncoeff, nel, nmax, nbc
    INTEGER, ALLOCATABLE :: irow_a(:), iloc_to_glob(:)
    INTEGER, POINTER :: if2b(:) => NULL()
    REAL(psb_dpk_) :: fsign, side_, fact
    REAL(psb_dpk_), ALLOCATABLE :: b(:), flux_x(:)
    TYPE(dimensions) :: dim
    TYPE(discretization) :: ds_
    TYPE(mesh), POINTER :: msh => NULL()
    TYPE(mesh), POINTER :: msh_flux => NULL()
    TYPE(vector) :: vf
    CALL sw_pde%tic()

    IF(mypnum_() == 0) THEN
        WRITE(*,*) '* ', TRIM(pde%name_()), ': applying the Divergence ',&
            & 'operator to the ', TRIM(flux%name_()), ' field'
    END IF

    ! Possible reinit of PDE
    CALL pde%reinit_pde()

    ! Is FLUX face-centered?
    IF(.NOT. flux%on_faces_()) THEN
        WRITE(*,200) TRIM(op_name)
        CALL abort_psblas
    END IF

    ! Checks mesh consistency PDE vs. FLUX
    CALL pde%get_mesh(msh)
    CALL flux%get_mesh(msh_flux)
    BLOCK
        !! This is in a block statement due to gfortran-8.3 issue w/ not being able to find this interface
        USE class_mesh, ONLY : check_mesh_consistency
        CALL check_mesh_consistency(msh,msh_flux,op_name)
    END BLOCK

    NULLIFY(msh_flux)

    ! Equation dimensional check
    dim = flux%dim_()
    IF(dim /= pde%dim_()) THEN
        WRITE(*,300) TRIM(op_name)
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
    CALL flux%get_x(flux_x)

    ! Total number of matrix coefficients associated to the fluid faces
    ncoeff = 2 * COUNT(msh%faces%flag_() == 0)

    ! Computes maximum size of blocks to be inserted
    nmax = size_blk(1,ncoeff)

    ALLOCATE(b(nmax),irow_a(nmax), stat=info)
    IF(info /= 0) THEN
        WRITE(*,400) TRIM(op_name)
        CALL abort_psblas
    END IF

    ! Gets local to global conversion list
    CALL psb_get_loc_to_glob(msh%desc_c,iloc_to_glob)


    ! ----- Insert contributions for fluid faces -----

    b = 0.d0
    irow_a = 0

    CALL msh%f2b%get_ith_conn(if2b,0)

    ifirst = 1; i = 0
    insert_fluid: DO
        IF(ifirst > ncoeff) EXIT insert_fluid
        nel = size_blk(ifirst,ncoeff)

        block_fluid: DO ka = 1, nel, 2
            ! Local indices
            i = i + 1
            IF = if2b(i)
            im = msh%faces(IF)%master_()
            is = msh%faces(IF)%slave_()
            ! Global indices
            im_glob = iloc_to_glob(im)
            is_glob = iloc_to_glob(is)

            ! --- Sparse Matrix A ---

            ! Gets A's coefficients

            fact = (-1.d0 * fsign) * flux_x(i)


            ! Master
            b(ka)  = fact
            irow_a(ka) = im_glob
            !
            ! Slave
            b(ka+1)  = (-1.d0) * fact
            irow_a(ka+1) = is_glob

        END DO block_fluid

        CALL geins_pde(nel,irow_a,b,pde)
        IF (debug_mat_bld) THEN
            WRITE(0,*) 'From scalar_vector_pde_div FLUID: geins ',nel
            DO ka=1,nel
                WRITE(0,*) irow_a(ka),b(ka)
            END DO
        END IF

        ifirst = ifirst + nel
    END DO insert_fluid

    ! Free storage
    DEALLOCATE(flux_x)


    ! ----- Insert source terms for boundary conditions -----

    ! Gets boundary values of FLUX fields
    CALL flux%get_bx(flux_x)

    b = 0.d0
    irow_a = 0

    ! Boundary faces with flag < IB
    ib_offset = 0

    nbc = msh%nbc

    bc_loop: DO ib = 1, nbc
        CALL msh%f2b%get_ith_conn(if2b,ib)
        n = SIZE(if2b)
        IF(n == 0) CYCLE bc_loop

        ifirst = 1; i = 0
        insert_boundary: DO
            IF(ifirst > n) EXIT insert_boundary
            nel = size_blk(ifirst,n)

            block_boundary: DO ka = 1, nel
                ! Local index
                i = i + 1
                IF = if2b(i)
                ibf = ib_offset + i
                im = msh%faces(IF)%master_()
                ! Global index
                im_glob = iloc_to_glob(im)

                b(ka) =  (-1.d0 * fsign) * flux_x(ibf)
                irow_a(ka) = im_glob
            END DO block_boundary

            CALL geins_pde(nel,irow_a,b,pde)
            IF (debug_mat_bld) THEN
                WRITE(0,*) 'From scalar_vector_pde_div BC: geins ',nel
                DO ka=1,nel
                    WRITE(0,*) irow_a(ka),b(ka)
                END DO
            END IF
            ifirst = ifirst + nel

        END DO insert_boundary

        ib_offset = ib_offset + n

    END DO bc_loop


    ! Free storage
    NULLIFY(if2b)
    DEALLOCATE(flux_x)
    DEALLOCATE(iloc_to_glob)
    DEALLOCATE(b,irow_a)
    NULLIFY(msh)

    CALL sw_pde%toc()

100 FORMAT(' ERROR! PHI field in ',a,' is not cell centered')
200 FORMAT(' ERROR! FLUX field in ',a,' is not cell centered')
300 FORMAT(' ERROR! Dimensional check failure in ',a)
400 FORMAT(' ERROR! Memory allocation failure in ',a)
401 FORMAT(' ERROR! Unimplemented discretization scheme in ',a)
402 FORMAT(' ERROR! Unimplemented BC  in ',a)

    END PROCEDURE flux_pde_div

END SUBMODULE flux_pde_div_implementation
