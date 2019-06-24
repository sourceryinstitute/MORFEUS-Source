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
!    Adds to PDE the source term SRC = sc + sp * phi_P
!
SUBMODULE(op_source) scalar_pde_source_implementation

    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE scalar_pde_source
            USE class_psblas
            USE class_dimensions
            USE class_mesh
            USE class_material
            USE class_scalar_field
            USE class_scalar_source
            USE class_scalar_pde
            USE tools_operators

            IMPLICIT NONE
            INTEGER :: i, ic, im, ic_glob, ifirst, info, ncells, nel, nmax
            INTEGER, ALLOCATABLE :: ia(:), ja(:)
            INTEGER, ALLOCATABLE  :: iloc_to_glob(:)
            REAL(psb_dpk_) :: fact, fsign, sc, sp, side_
            REAL(psb_dpk_), ALLOCATABLE :: A(:), b(:)
            TYPE(dimensions) :: dim
            TYPE(mesh), POINTER :: msh => NULL()

            CALL tic(sw_pde)

            IF(mypnum_() == 0) THEN
                WRITE(*,*) '* ', TRIM(name_(pde)), ': applying the Source term'
            END IF

            ! Possible reinit of pde
            CALL reinit_pde(pde)


            ! Dimensional check
            dim = dim_(src) * volume_
            IF(dim /= dim_(pde)) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            ! Computes sign factor
            IF(PRESENT(side)) THEN
                side_ = side
            ELSE
                side_ = lhs_ ! Default = LHS
            END IF
            fsign = pde_sign(sign,side_)


            ! Gets PDE mesh
            CALL pde%get_mesh(msh)

            ! Gets components of source term once for all
            sc = sc_(src)
            sp = sp_(src)

            ! Number of strictly local cells
            ncells = psb_cd_get_local_rows(msh%desc_c)

            ! Gets local to global list for cell indices
            CALL psb_get_loc_to_glob(msh%desc_c,iloc_to_glob)

            ! Computes maximum size of blocks to be inserted
            nmax = size_blk(1,ncells)

            ALLOCATE(A(nmax),b(nmax),ia(nmax),ja(nmax),stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            ifirst = 1; ic = 0
            insert: DO
                IF(ifirst > ncells) EXIT insert
                nel = size_blk(ifirst,ncells)

                BLOCK: DO i = 1, nel
                    ! Local indices
                    ic = ic + 1

                    im = get_scalar_field_mat_id(phi,i+ifirst-1)
                    IF (mat_id_(mats(im)%mat) < 100) THEN
                        fact = fsign * msh%vol(ic)
                    ELSE
                        fact = 0.0
                    END IF

                    ! WARNING! If one assumes that SRC is applied to RHS then
                    ! pde_sign returns FSIGN < 0.
                    A(i) = fact * sp
                    b(i) = -fact * sc

                    ! Global indices in COO format
                    ic_glob = iloc_to_glob(ic)
                    ia(i)= ic_glob
                    ja(i)= ic_glob
                END DO BLOCK

                CALL spins_pde(nel,ia,ja,A,pde)
                CALL geins_pde(nel,ia,b,pde)

                ifirst = ifirst + nel

            END DO insert
            DEALLOCATE(A,b,ia,ja)
            DEALLOCATE(iloc_to_glob)
            NULLIFY(msh)

            CALL toc(sw_pde)

        100 FORMAT(' ERROR! Dimensional check failure in SCALAR_PDE_SOURCE')
        200 FORMAT(' ERROR! Memory allocation failure in SCALAR_PDE_SOURCE')

        END PROCEDURE scalar_pde_source

END SUBMODULE scalar_pde_source_implementation
