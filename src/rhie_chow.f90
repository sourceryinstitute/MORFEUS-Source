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
!
! Description:
!    To be added...
!
SUBMODULE(op_field) rhie_chow_implementation

    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE rhie_chow
            USE class_psblas
            USE class_connectivity
            USE class_dimensions
            USE class_field
            USE class_scalar_field
            USE class_vector_field
            USE class_vector
            USE class_mesh
            USE class_face
            USE class_vector

            IMPLICIT NONE
            !
            INTEGER :: i, ib, ibf, ib_offset, IF, info, n, im, is
            INTEGER, POINTER :: if2b(:) => NULL()
            REAL(psb_dpk_), ALLOCATABLE :: r_x(:), r_bx(:)
            TYPE(dimensions) :: dim
            TYPE(field) :: base
            TYPE(scalar_field) :: fld_f, phi_f
            TYPE(mesh), POINTER :: msh => NULL()
            REAL(psb_dpk_), ALLOCATABLE :: fld_x(:), fld_bx(:), phi_x(:), phi_bx(:)
            REAL(psb_dpk_), ALLOCATABLE :: phi_fx(:), phi_fbx(:)
            TYPE(vector) :: fact
            REAL(psb_dpk_) :: factu,factv,factz

            ! Gets FLD mesh
            CALL phi%get_mesh(msh)

            ! Gets BASE, X and BX members of the face-centered scalar field
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

            CALL phi%get_x(phi_x)
            CALL phi%get_bx(phi_bx)

            phi_f = phi%interp_on_faces()
            CALL phi_f%get_x(phi_fx)
            CALL phi_f%get_bx(phi_fbx)


            ! Computes result dimensions
            dim = phi%dim_() * fld%dim_() * surface_ * surface_ /mass_

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
            DO ib = 0, -1, -1
                CALL msh%f2b%get_ith_conn(if2b,ib)
                n = SIZE(if2b)

                DO i = 1, n
                    IF = if2b(i)
                    im = msh%faces(IF)%master_()
                    is = msh%faces(IF)%slave_()
                    ibf = ib_offset + i
                    r_x(ibf) =  msh%area(IF)*((msh%df(IF) .dot. msh%af(IF))*fld_x(ibf)*&
                        & ((phi_x(is)-phi_x(im))/msh%dist(IF)-phi_fx(ibf)))

                END DO
                ib_offset = ib_offset + n
            END DO

            r_bx=0.d0

            ! ... then computes fluxes on boundary faces.
            !  ib_offset = 0
            !  do ib = 1, msh%nbc
            !     call get_ith_conn(if2b,msh%f2b,ib)
            !     n = size(if2b)
            !     do i = 1, n
            !        if = if2b(i)
            !        im = master_(msh%faces(if))
            !        ibf = ib_offset + i
            !        r_bx(ibf) =  ((msh%df(if) .dot. msh%af(if))*fld_bx(ibf)*&
            !                  & ((phi_bx(ibf)-phi_x(im))/msh%dist(if)-phi_fbx(ibf)))&
            !                  &  *msh%area(if)
            !     end do
            !     ib_offset = ib_offset + n
            !  end do

            ! Eventually construct the result
            res = scalar_field(base,x=r_x,bx=r_bx)

            NULLIFY(if2b)
            DEALLOCATE(r_x,r_bx, fld_x, fld_bx, phi_x, phi_bx, phi_fx,phi_fbx)
            CALL base%free_field()
            NULLIFY(msh)

        100 FORMAT(' ERROR! Memory allocation failure in VECTOR_FIELD_FLUX')

        END PROCEDURE rhie_chow

END SUBMODULE rhie_chow_implementation
