!! label: Morfeus-FV
!!
!! Define and manipulate scalar-field values associated with a field-parent mesh

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
MODULE class_scalar_field
  !! label: Morfeus-FV
  !!
  !! Define and manipulate scalar-field values associated with a field-parent mesh
    USE class_psblas, ONLY : psb_dpk_, nemo_int_long_, nemo_sizeof_dp, nemo_sizeof_int,&
        & icontxt_, psb_erractionsave, abort_psblas, psb_check_error, psb_erractionrestore,&
        & psb_halo
    USE class_field, ONLY : field
    USE class_mesh, ONLY : mesh
    USE class_dimensions, ONLY : dimensions
    USE class_bc, ONLY : bc_poly
    USE class_material, ONLY : material, matptr

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: scalar_field           !! Class

    TYPE, EXTENDS(field) :: scalar_field
      !! label: Morfeus-FV
      !!
      !! Encapsulate scalar-field state and operations, storing a mesh association in the parent field
      !!
        PRIVATE
        REAL(psb_dpk_), ALLOCATABLE :: x(:)  !! field nodal values
        REAL(psb_dpk_), ALLOCATABLE :: xp(:) !! previous nodal values
        REAL(psb_dpk_), ALLOCATABLE :: bx(:) !! face values
        INTEGER, ALLOCATABLE :: mat(:)
        INTEGER, ALLOCATABLE :: bmat(:)
    CONTAINS
        PROCEDURE, PUBLIC :: create_scalar_field               !! Constructor
        PROCEDURE, PUBLIC ::  free_field                       !! Destructor
        PROCEDURE, PUBLIC :: get_base
        PROCEDURE, PRIVATE :: get_scalar_field_x, get_scalar_field_element
        GENERIC, PUBLIC :: get_x => get_scalar_field_x, get_scalar_field_element
        PROCEDURE, PRIVATE :: get_scalar_field_xp, get_scalar_field_element_prev
        GENERIC, PUBLIC :: get_xp => get_scalar_field_xp, get_scalar_field_element_prev
        PROCEDURE, PRIVATE :: get_scalar_field_bx
        GENERIC, PUBLIC :: get_bx => get_scalar_field_bx
        PROCEDURE :: get_scalar_field_mat_id
        PROCEDURE, PRIVATE :: update_scalar_field
        GENERIC, PUBLIC :: update_field => update_scalar_field
        PROCEDURE, PRIVATE :: set_scalar_field_element, set_scalar_field_group
        GENERIC, PUBLIC :: set_field_element => set_scalar_field_element
        GENERIC, PUBLIC :: set_field_group => set_scalar_field_group
        PROCEDURE,PRIVATE :: interp_on_faces_s
        GENERIC, PUBLIC :: interp_on_faces => interp_on_faces_s
        PROCEDURE, PRIVATE :: nemo_scalar_field_normi, nemo_scalar_field_norm1
        GENERIC, PUBLIC :: field_normi => nemo_scalar_field_normi
        GENERIC, PUBLIC :: field_norm1 => nemo_scalar_field_norm1
        PROCEDURE, PRIVATE :: scalar_field_sum, scalar_field_dif, scalar_field_dif_s
        PROCEDURE, PRIVATE :: scalar_field_div
        PROCEDURE, PASS(f), PRIVATE :: scalar_field_scal
        PROCEDURE, PRIVATE :: scalar_field_mul
        PROCEDURE, PRIVATE :: assign_scalar_field_s, assign_scalar_field_v
        GENERIC :: ASSIGNMENT(=) => assign_scalar_field_s, assign_scalar_field_v  !! User-defined assignment
        GENERIC :: OPERATOR(*) => scalar_field_scal, scalar_field_mul             !! Algebra operations
        GENERIC :: OPERATOR(+) => scalar_field_sum
        GENERIC :: OPERATOR(-) => scalar_field_dif, scalar_field_dif_s
        GENERIC :: OPERATOR(/) => scalar_field_div
        PROCEDURE, PUBLIC :: nemo_sizeof
        PROCEDURE, PRIVATE :: check_mesh_consistency_sf
        GENERIC, PUBLIC :: check_mesh_consistency => check_mesh_consistency_sf
    END TYPE scalar_field

    INTERFACE scalar_field

        MODULE FUNCTION scalar_field_(base,x,bx)
          !! label: Morfeus-FV
          !!
          !! Result is new scalar_field from parent field
            IMPLICIT NONE
            TYPE(scalar_field) :: scalar_field_
            TYPE(field),      INTENT(IN) :: base
            REAL(psb_dpk_), INTENT(IN) :: x(:)
            REAL(psb_dpk_), INTENT(IN) :: bx(:)
        END FUNCTION scalar_field_

    END INTERFACE scalar_field

    INTERFACE

        MODULE FUNCTION nemo_sizeof(fld)
          !! label: Morfeus-FV
          !!
          !! Result is the scalar_field size
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN) :: fld
            INTEGER(kind=nemo_int_long_)   :: nemo_sizeof
        END FUNCTION nemo_sizeof

        MODULE SUBROUTINE create_scalar_field(fld,msh,dim,bc,mats,on_faces,x0)
          !! label: Morfeus-FV
          !!
          !! Define scalar_field from mesh, dimensions, boundary conditions, materials, face-centeredness, and values
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(OUT)        :: fld
            TYPE(mesh),          INTENT(IN), TARGET :: msh
            TYPE(dimensions), INTENT(IN), OPTIONAL         :: dim
            TYPE(bc_poly),    INTENT(IN), OPTIONAL, TARGET :: bc(:)
            TYPE(matptr),     INTENT(IN), OPTIONAL, TARGET :: mats(:)
            LOGICAL,          INTENT(IN), OPTIONAL         :: on_faces
            REAL(psb_dpk_),   INTENT(IN), OPTIONAL         :: x0
        END SUBROUTINE create_scalar_field

        MODULE SUBROUTINE free_field(fld)
          !! label: Morfeus-FV
          !!
          !! Deallocate scalar_field components: field parent and component arrays
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(INOUT) :: fld
        END SUBROUTINE free_field

        INTEGER MODULE FUNCTION get_scalar_field_mat_id(fld,i)
          !! label: Morfeus-FV
          !!
          !! Result is scalar_field material identifier
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN) :: fld
            INTEGER, INTENT(IN) :: i
        END FUNCTION get_scalar_field_mat_id

        MODULE SUBROUTINE get_base(fld,base)
          !! label: Morfeus-FV
          !!
          !! Get parent field
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN)  :: fld
            TYPE(field),         INTENT(OUT) :: base
        END SUBROUTINE get_base

        MODULE SUBROUTINE get_scalar_field_x(fld,x)
          !! label: Morfeus-FV
          !!
          !! get all scalar_field values
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN) :: fld
            REAL(psb_dpk_),      INTENT(OUT), ALLOCATABLE :: x(:)
        END SUBROUTINE get_scalar_field_x

        MODULE SUBROUTINE get_scalar_field_element(fld,x,i)
          !! label: Morfeus-FV
          !!
          !! get one scalar_field value
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN)  :: fld
            REAL(psb_dpk_),      INTENT(OUT) :: x
            INTEGER,             INTENT(IN)  :: i
        END SUBROUTINE get_scalar_field_element

        MODULE SUBROUTINE get_scalar_field_xp(fld,xp)
          !! label: Morfeus-FV
          !!
          !! Get all scalar_field previous values
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN) :: fld
            REAL(psb_dpk_),      INTENT(OUT), ALLOCATABLE :: xp(:)
        END SUBROUTINE get_scalar_field_xp

        MODULE SUBROUTINE get_scalar_field_element_prev(fld,xp,i)
          !! label: Morfeus-FV
          !!
          !! get one scalar_field previous value
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN)  :: fld
            REAL(psb_dpk_),      INTENT(OUT) :: xp
            INTEGER,             INTENT(IN)  :: i
        END SUBROUTINE get_scalar_field_element_prev

        MODULE SUBROUTINE get_scalar_field_bx(fld,bx)
          !! label: Morfeus-FV
          !!
          !! get scalar_field boundary values
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN) :: fld
            REAL(psb_dpk_),      INTENT(OUT), ALLOCATABLE :: bx(:)
        END SUBROUTINE get_scalar_field_bx

        MODULE SUBROUTINE update_scalar_field(fld,mats,temp)
          !! label: Morfeus-FV
          !!
          !! Set all scalar_field values
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(INOUT) :: fld
            TYPE(scalar_field),  INTENT(IN), OPTIONAL :: temp     !! temperature
            TYPE(matptr),        INTENT(IN), POINTER  :: mats(:)
        END SUBROUTINE update_scalar_field

        MODULE SUBROUTINE set_scalar_field_element(f,i,x)
          !! label: Morfeus-FV
          !!
          !! Set one scalar_field value
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(INOUT) :: f
            INTEGER,             INTENT(IN)    :: i !! element number to set
            REAL(psb_dpk_),      INTENT(IN)    :: x !! value
        END SUBROUTINE set_scalar_field_element

        MODULE SUBROUTINE set_scalar_field_group(f,ig,x)
          !! label: Morfeus-FV
          !!
          !! Set a group of scalar field values
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(INOUT) :: f
            INTEGER,             INTENT(IN)    :: ig
            REAL(psb_dpk_),      INTENT(IN)    :: x
        END SUBROUTINE set_scalar_field_group

        !! ----- Algebra Operations -----

        MODULE FUNCTION nemo_scalar_field_normi(fld) RESULT(norm)
          !! label: Morfeus-FV
          !!
          !! Result is infinity-norm of scalar_field
            IMPLICIT NONE
            ClASS(scalar_field), INTENT(IN) :: fld
            REAL(psb_dpk_)                  :: norm
        END FUNCTION nemo_scalar_field_normi

        MODULE FUNCTION nemo_scalar_field_norm1(fld) RESULT(norm)
          !! label: Morfeus-FV
          !!
          !! Result is L1 norm of scalar_field
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN) :: fld
            REAL(psb_dpk_)                  :: norm
        END FUNCTION nemo_scalar_field_norm1

        MODULE FUNCTION scalar_field_sum(f1,f2)RESULT(r)
          !! label: Morfeus-FV
          !!
          !! Result is the arithmetic sum of the two scalar_field operands
            IMPLICIT NONE
            TYPE(scalar_field) :: r
            CLASS(scalar_field), INTENT(IN) :: f1
            TYPE(scalar_field),  INTENT(IN) :: f2
        END FUNCTION scalar_field_sum

        MODULE FUNCTION scalar_field_dif(f1,f2)RESULT(r)
          !! label: Morfeus-FV
          !!
          !! Result is the arithmetic difference of the two scalar_field operands
            IMPLICIT NONE
            TYPE(scalar_field) :: r
            CLASS(scalar_field), INTENT(IN) :: f1
            TYPE(scalar_field),  INTENT(IN) :: f2
        END FUNCTION scalar_field_dif

        MODULE FUNCTION scalar_field_dif_s(f1,f2)RESULT(r)
          !! label: Morfeus-FV
          !!
          !! Result is the arithmetic difference of the scalar_field operand and the constant operand
            IMPLICIT NONE
            TYPE(scalar_field) :: r
            CLASS(scalar_field), INTENT(IN) :: f1
            REAL(psb_dpk_),      INTENT(IN) :: f2
        END FUNCTION scalar_field_dif_s

        MODULE FUNCTION scalar_field_div(f1,f2)RESULT(r)
          !! label: Morfeus-FV
          !!
          !! Result is the ratio of the left-hand-side numerator and right-hand-side denominator scalar_field objects
            IMPLICIT NONE
            TYPE(scalar_field) :: r
            CLASS(scalar_field), INTENT(IN) :: f1
            TYPE(scalar_field),  INTENT(IN) :: f2
        END FUNCTION scalar_field_div

        MODULE FUNCTION interp_on_faces_s(fld)RESULT(r)
          !! label: Morfeus-FV
          !!
          !! Result contains scalar_field cell values interpolated onto cell faces
            IMPLICIT NONE
            TYPE(scalar_field) :: r
            CLASS(scalar_field), INTENT(IN) :: fld
        END FUNCTION interp_on_faces_s

        MODULE SUBROUTINE check_mesh_consistency_sf(f1,f2,WHERE)
          !! label: Morfeus-FV
          !!
          !! Abort iff the scalar_field arguments are not associated with the same mesh
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(IN) :: f1
            TYPE(scalar_field),  INTENT(IN) :: f2
            CHARACTER(len=*),    INTENT(IN) :: WHERE
        END SUBROUTINE check_mesh_consistency_sf

        MODULE FUNCTION scalar_field_scal(a,f)RESULT(r)
          !! label: Morfeus-FV
          !!
          !! Rescale the scalar_field argument by the provided multiplicative constant
            IMPLICIT NONE
            TYPE(scalar_field) :: r
            REAL(psb_dpk_),      INTENT(IN) :: a
            CLASS(scalar_field), INTENT(IN) :: f
        END FUNCTION scalar_field_scal

        MODULE FUNCTION scalar_field_mul(f1,f2)RESULT(r)
          !! label: Morfeus-FV
          !!
          !! REsult is the arithmetic product of two scalar_field objects
            IMPLICIT NONE
            TYPE(scalar_field) :: r
            CLASS(scalar_field), INTENT(IN) :: f1
            TYPE(scalar_field),  INTENT(IN) :: f2
        END FUNCTION scalar_field_mul

        MODULE SUBROUTINE assign_scalar_field_s(f,x)
          !! label: Morfeus-FV
          !!
          !! Set the left-hand-side scalar_field values to the provided constant
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(INOUT) :: f
            REAL(psb_dpk_),      INTENT(IN)    :: x
        END SUBROUTINE assign_scalar_field_s

        MODULE SUBROUTINE assign_scalar_field_v(f,x)
          !! label: Morfeus-FV
          !!
          !! Set the left-hand-side scalar_field values to the provided values
            IMPLICIT NONE
            CLASS(scalar_field), INTENT(INOUT) :: f
            REAL(psb_dpk_),      INTENT(IN)    :: x(:)
        END SUBROUTINE assign_scalar_field_v

    END INTERFACE

END MODULE class_scalar_field
