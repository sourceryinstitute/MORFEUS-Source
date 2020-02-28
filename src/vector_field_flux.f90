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
! $Id: vector_field_flux.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(op_field) vector_field_flux_implementation

    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE vector_field_flux
            USE class_psblas
            USE class_connectivity
            USE class_dimensions
            USE class_field
            USE class_scalar_field
            USE class_vector_field
            USE class_vector
            USE class_mesh

            IMPLICIT NONE
            !
            INTEGER :: i, ib, ibf, ib_offset, IF, info, n
            INTEGER, POINTER :: if2b(:) => NULL()
            REAL(psb_dpk_), ALLOCATABLE :: r_x(:), r_bx(:)
            TYPE(dimensions) :: dim
            TYPE(field) :: base
            TYPE(vector_field) :: fld_f
            TYPE(mesh), POINTER :: msh => NULL()
            TYPE(vector), ALLOCATABLE :: fld_x(:), fld_bx(:)


            ! Gets FLD mesh
            CALL fld%get_mesh(msh)

            ! Gets BASE, X and BX members of the face-centered vector field
            IF(fld%on_faces_()) THEN
                CALL fld%get_x(fld_x)
                CALL fld%get_bx(fld_bx)
                CALL fld%get_base(base)

            ELSE
                ! If FLD is cell-centered, first interpolate it
                fld_f = fld%interp_on_faces()

                CALL fld_f%get_x(fld_x)
                CALL fld_f%get_bx(fld_bx)
                CALL fld_f%get_base(base)
                CALL fld_f%free_field()
            END IF

            ! Computes result dimensions
            dim = fld%dim_() * surface_

            ! Sets DIM member in the base field object
            CALL base%set_field_dim(dim)

            ! Allocates arrays for storing result values
            ALLOCATE(r_x(SIZE(fld_x)),r_bx(SIZE(fld_bx)),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            ! First computes fluxes on faces with flag = 0 and flag = -1...
            ib_offset = 0
            DO ib = 0, 0 !-1, -1
                CALL msh%f2b%get_ith_conn(if2b,ib)
                n = SIZE(if2b)

                DO i = 1, n
                    IF = if2b(i)
                    ibf = ib_offset + i
                    r_x(ibf) = fld_x(ibf) .dot. msh%af(IF)
                END DO
                ib_offset = ib_offset + n
            END DO

            ! ... then computes fluxes on boundary faces.
            ib_offset = 0
            DO ib = 1, msh%nbc
                CALL msh%f2b%get_ith_conn(if2b,ib)
                n = SIZE(if2b)

                DO i = 1, n
                    IF = if2b(i)
                    ibf = ib_offset + i
                    r_bx(ibf) = fld_bx(ibf) .dot. msh%af(IF)
                END DO

                ib_offset = ib_offset + n
            END DO

            ! Eventually construct the result
            res = scalar_field(base,x=r_x,bx=r_bx)

            NULLIFY(if2b)
            DEALLOCATE(r_x,r_bx)
            CALL base%free_field()
            CALL free_vector(fld_x)
            CALL free_vector(fld_bx)
            NULLIFY(msh)

        100 FORMAT(' ERROR! Memory allocation failure in VECTOR_FIELD_FLUX')

        END PROCEDURE vector_field_flux

END SUBMODULE vector_field_flux_implementation
