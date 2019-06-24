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
!    Gradient operator as a source term in a vector PDE.
!
SUBMODULE(op_grad) vector_pde_grad_implementation

    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE vector_pde_grad
            USE class_psblas
            USE class_connectivity
            USE class_dimensions
            USE class_face
            USE class_scalar_field
            USE class_mesh
            USE class_vector_pde
            USE class_vector
            USE tools_math
            USE tools_operators

            IMPLICIT NONE
            !
            CHARACTER(len=*), PARAMETER :: op_name = 'VECTOR_PDE_GRAD'
            INTEGER :: i, ib, ibf, ib_offset, IF, ifirst, info, ic
            INTEGER :: ic_glob
            INTEGER :: n, nbc, nmax, nel, ncells
            INTEGER, ALLOCATABLE :: ia(:), iloc_to_glob(:)
            INTEGER, POINTER :: if2b(:) => NULL()
            REAL(psb_dpk_) :: fsign, side_, w
            REAL(psb_dpk_), ALLOCATABLE :: x(:), bx(:)
            TYPE(dimensions)   :: dim
            TYPE(mesh), POINTER :: msh => NULL(), msh_phi => NULL()
            TYPE(vector) :: fact
            TYPE(vector), ALLOCATABLE :: b(:), grad(:)

            CALL sw_pde%tic()

            IF(mypnum_() == 0) THEN
                WRITE(*,*) '* ', TRIM(pde%name_()), ': applying the Gradient ',&
                    & 'operator to the ', TRIM(phi%name_()), ' field'
            END IF

            ! Possible reinit of PDE
            CALL pde%reinit_pde()

            ! Is PHI cell-centered?
            IF(phi%on_faces_()) THEN
                WRITE(*,100) TRIM(op_name)
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
            NULLIFY(msh_phi)

            ! Equation dimensional check
            dim = phi%dim_() / length_ * volume_
            IF(dim /= pde%dim_()) THEN
                WRITE(*,200) TRIM(op_name)
                CALL abort_psblas
            END IF

            ! Computes sign factor
            IF(PRESENT(side)) THEN
                side_ = side
            ELSE
                side_ = lhs_ ! Default = LHS
            END IF
            fsign = pde_sign(sign,side_)

            ! Gets PHI internal ("x") and boundary ("bx") values
            CALL phi%get_x(x)
            CALL phi%get_bx(bx)

            ! Computes gradient for non-orthogonality correction
            ! First allocates the result ...
            ALLOCATE(grad(SIZE(x)),stat=info)
            IF(info /= 0) THEN
                WRITE(*,300) TRIM(op_name)
                CALL abort_psblas
            END IF

            ! ... then computes it.
            grad = fld_grad(phi)

            nbc = msh%nbc

            ! Number of strictly local cells
            ncells = psb_cd_get_local_rows(msh%desc_c)

            ! Gets local to global list for cell indices
            CALL psb_get_loc_to_glob(msh%desc_c,iloc_to_glob)

            ! Computes maximum size of blocks to be inserted
            nmax = size_blk(1,ncells)

            ALLOCATE(b(nmax),ia(nmax),stat=info)
            IF(info /= 0) THEN
                WRITE(*,300) TRIM(op_name)
                CALL abort_psblas
            END IF




            ifirst = 1; ic = 0
            insert : DO
                IF(ifirst > ncells) EXIT insert
                nel = size_blk(ifirst,ncells)

                BLOCK: DO i = 1, nel
                    ! Local indices
                    ic = ic + 1
                    fact = fsign * msh%vol(ic) * grad(ic)
                    ! Global indices
                    ic_glob = iloc_to_glob(ic)

                    ! WARNING! If the source term is applied to the RHS then
                    ! pde_sign returns FSIGN < 0. Hence, in order to have a
                    ! positive contribution one has to use - FSIGN

                    b(i) = -fact
                    ia(i) = ic_glob

                END DO BLOCK

                CALL geins_pde(nel,ia,b,pde)

                ifirst = ifirst + nel
            END DO insert



            NULLIFY(if2b)
            DEALLOCATE(iloc_to_glob)
            DEALLOCATE(b,ia)
            NULLIFY(msh)
            DEALLOCATE(x,bx)


            CALL sw_pde%toc()

        100 FORMAT(' ERROR! Operand in ',a,' is not cell centered')
        200 FORMAT(' ERROR! Dimensional check failure in ',a)
        300 FORMAT(' ERROR! Memory allocation failure in ',a)

        END PROCEDURE vector_pde_grad

END SUBMODULE vector_pde_grad_implementation
