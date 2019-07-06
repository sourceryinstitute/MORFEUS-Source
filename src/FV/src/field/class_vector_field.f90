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
! $Id: class_vector_field.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_vector_field

    USE class_psblas, ONLY : psb_dpk_, nemo_int_long_, nemo_sizeof_dp, nemo_sizeof_int,&
        & icontxt_, psb_erractionsave, abort_psblas, psb_check_error, psb_erractionrestore,&
        & psb_halo
    USE class_field, ONLY : field
    USE class_vector, ONLY : vector
    USE class_mesh, ONLY : mesh
    USE class_dimensions, ONLY : dimensions
    USE class_bc, ONLY : bc_poly
    USE class_material, ONLY : matptr, material

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: vector_field      !! Class

    TYPE, EXTENDS(field) :: vector_field
        PRIVATE
        TYPE(vector), ALLOCATABLE :: x(:)
        TYPE(vector), ALLOCATABLE :: xp(:)
        TYPE(vector), ALLOCATABLE :: bx(:)
        INTEGER, ALLOCATABLE :: mat(:)
    CONTAINS
        PROCEDURE, PUBLIC :: create_vector_field
        PROCEDURE, PUBLIC :: free_field
        PROCEDURE, PRIVATE :: get_vector_field_base
        GENERIC, PUBLIC :: get_base => get_vector_field_base
        PROCEDURE, PRIVATE :: get_vector_field_x_r, get_vector_field_x_v
        GENERIC, PUBLIC :: get_x => get_vector_field_x_r, get_vector_field_x_v
        PROCEDURE, PRIVATE :: get_vector_field_xp_r, get_vector_field_xp_v
        GENERIC, PUBLIC :: get_xp => get_vector_field_xp_r, get_vector_field_xp_v
        PROCEDURE, PRIVATE :: get_vector_field_bx_r, get_vector_field_bx_v
        GENERIC, PUBLIC :: get_bx => get_vector_field_bx_r, get_vector_field_bx_v
        PROCEDURE, PRIVATE :: update_vector_field
        GENERIC, PUBLIC :: update_field => update_vector_field
        PROCEDURE, PRIVATE :: set_vector_field_element, set_vector_field_group
        GENERIC, PUBLIC :: set_field_element => set_vector_field_element
        GENERIC, PUBLIC :: set_field_group => set_vector_field_group
        PROCEDURE, PRIVATE :: set_vector_field_bound_element
        GENERIC, PUBLIC :: set_field_bound_element => set_vector_field_bound_element
        PROCEDURE, PRIVATE :: set_vector_field_x
        GENERIC, PUBLIC :: set_x => set_vector_field_x
        PROCEDURE, PRIVATE :: interp_on_faces_v
        GENERIC, PUBLIC :: interp_on_faces => interp_on_faces_v
        PROCEDURE, PRIVATE :: vector_field_sum, vector_field_dif
        PROCEDURE, PASS(f2), PRIVATE :: vector_field_scal     !! Algebra operations
        PROCEDURE, PRIVATE :: assign_vector_field_s, assign_vector_field_v
        GENERIC :: ASSIGNMENT(=) => assign_vector_field_s, assign_vector_field_v !! User-defined assignemnts
        GENERIC :: OPERATOR(+) => vector_field_sum
        GENERIC :: OPERATOR(-) => vector_field_dif
        GENERIC :: OPERATOR(*) => vector_field_scal     !! Algebra operations
        PROCEDURE, PUBLIC:: nemo_sizeof
        PROCEDURE, PRIVATE :: check_mesh_consistency_vf
        GENERIC, PUBLIC :: check_mesh_consistency => check_mesh_consistency_vf
    END TYPE vector_field

    ! ----- Generic Interfaces -----

    INTERFACE vector_field

      MODULE FUNCTION vector_field_(base,x,bx)
        !! Default public constructor, necessary with ifort
          IMPLICIT NONE
          TYPE(vector_field) :: vector_field_
          TYPE(field),      INTENT(IN) :: base
          TYPE(vector), INTENT(IN) :: x(:)
          TYPE(vector), INTENT(IN) :: bx(:)
      END FUNCTION vector_field_

    END INTERFACE vector_field

    INTERFACE

      MODULE FUNCTION nemo_sizeof(fld)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(IN) :: fld
          INTEGER(kind=nemo_int_long_)   :: nemo_sizeof
      END FUNCTION nemo_sizeof

      MODULE SUBROUTINE create_vector_field(fld,msh,dim,bc,mats,on_faces,x0)
          !! Constructor
          IMPLICIT NONE
          ! Mandatory arguments
          CLASS(vector_field), INTENT(OUT)        :: fld
          TYPE(mesh),         INTENT(IN), TARGET :: msh
          !
          ! Optional arguments
          TYPE(dimensions), INTENT(IN), OPTIONAL         :: dim
          TYPE(bc_poly),    INTENT(IN), OPTIONAL, TARGET :: bc(:)
          TYPE(matptr),   INTENT(IN), OPTIONAL, TARGET :: mats(:)
          LOGICAL,          INTENT(IN), OPTIONAL         :: on_faces
          TYPE(vector),     INTENT(IN), OPTIONAL         :: x0
      END SUBROUTINE create_vector_field

    ! ----- Destructor -----

      MODULE SUBROUTINE free_field(fld)
          !! Destructor
          IMPLICIT NONE
          CLASS(vector_field), INTENT(INOUT) :: fld
      END SUBROUTINE free_field

    ! ----- Getters for Inherited Members -----

      MODULE SUBROUTINE get_vector_field_base(fld,base)
          CLASS(vector_field), INTENT(IN)  :: fld
          TYPE(field),        INTENT(OUT) :: base
      END SUBROUTINE get_vector_field_base

    ! ----- Getters for Additional Members -----

      !! Getters for Additional Members

      MODULE SUBROUTINE get_vector_field_x_r(fld,x)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(IN) :: fld
          REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: x(:,:)
      END SUBROUTINE get_vector_field_x_r

      MODULE SUBROUTINE get_vector_field_x_v(fld,x)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(IN) :: fld
          TYPE(vector),       INTENT(OUT), ALLOCATABLE :: x(:)
      END SUBROUTINE get_vector_field_x_v

      MODULE SUBROUTINE get_vector_field_xp_r(fld,xp)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(IN) :: fld
          REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: xp(:,:)
      END SUBROUTINE get_vector_field_xp_r

      MODULE SUBROUTINE get_vector_field_xp_v(fld,xp)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(IN) :: fld
          TYPE(vector),       INTENT(OUT), ALLOCATABLE :: xp(:)
      END SUBROUTINE get_vector_field_xp_v


      MODULE SUBROUTINE get_vector_field_bx_r(fld,bx)
          CLASS(vector_field), INTENT(IN) :: fld
          REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: bx(:,:)
      END SUBROUTINE get_vector_field_bx_r

      MODULE SUBROUTINE get_vector_field_bx_v(fld,bx)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(IN) :: fld
          TYPE(vector),       INTENT(OUT), ALLOCATABLE :: bx(:)
      END SUBROUTINE get_vector_field_bx_v

      ! ----- Setters -----

      MODULE SUBROUTINE set_vector_field_x(fld,x)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(INOUT) :: fld
          REAL(psb_dpk_),   INTENT(IN), ALLOCATABLE :: x(:,:)
      END SUBROUTINE set_vector_field_x

      MODULE SUBROUTINE update_vector_field(fld)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(INOUT) :: fld
      END SUBROUTINE update_vector_field

      MODULE SUBROUTINE set_vector_field_element(f,i,x)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(INOUT) :: f
          INTEGER, INTENT(IN) :: i
          TYPE(vector), INTENT(IN) :: x
      END SUBROUTINE set_vector_field_element

      MODULE SUBROUTINE set_vector_field_bound_element(f,i,x)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(INOUT) :: f
          INTEGER, INTENT(IN) :: i
          TYPE(vector), INTENT(IN) :: x
      END SUBROUTINE set_vector_field_bound_element

      MODULE SUBROUTINE set_vector_field_group(f,ig,x)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(INOUT) :: f
          INTEGER, INTENT(IN) :: ig
          TYPE(vector), INTENT(IN) :: x
      END SUBROUTINE set_vector_field_group

      ! ----- Algebra Operations -----

      MODULE FUNCTION vector_field_sum(f1,f2)RESULT(r)
          IMPLICIT NONE
          TYPE(vector_field) :: r
          CLASS(vector_field), INTENT(IN) :: f1
          TYPE(vector_field), INTENT(IN) :: f2
      END FUNCTION vector_field_sum

      MODULE FUNCTION vector_field_dif(f1,f2) RESULT(r)
          IMPLICIT NONE
          TYPE(vector_field) :: r
          CLASS(vector_field), INTENT(IN) :: f1
          TYPE(vector_field), INTENT(IN) :: f2
      END FUNCTION vector_field_dif

      MODULE FUNCTION interp_on_faces_v(fld)RESULT(r)
          IMPLICIT NONE
          TYPE(vector_field) :: r
          CLASS(vector_field), INTENT(IN) :: fld
      END FUNCTION interp_on_faces_v

      ! ----- Check Procedures -----

      MODULE SUBROUTINE check_mesh_consistency_vf(f1,f2,WHERE)
          IMPLICIT NONE
          CLASS(vector_field), INTENT(IN) :: f1
          TYPE(vector_field), INTENT(IN) :: f2
          CHARACTER(len=*), INTENT(IN) :: WHERE
      END SUBROUTINE check_mesh_consistency_vf

      MODULE FUNCTION vector_field_scal(a,f2)RESULT(r)
          USE class_dimensions
          IMPLICIT NONE
          TYPE(vector_field) :: r
          REAL(psb_dpk_), INTENT(IN) :: a
          CLASS(vector_field), INTENT(IN) :: f2
      END FUNCTION vector_field_scal

      MODULE SUBROUTINE assign_vector_field_s(f,x)
          USE class_vector
          IMPLICIT NONE
          CLASS(vector_field), INTENT(INOUT) :: f
          TYPE(vector),       INTENT(IN)    :: x
      END SUBROUTINE assign_vector_field_s

      MODULE SUBROUTINE assign_vector_field_v(f,x)
          USE class_vector
          IMPLICIT NONE
          CLASS(vector_field), INTENT(INOUT) :: f
          TYPE(vector),       INTENT(IN)    :: x(:)
      END SUBROUTINE assign_vector_field_v

    END INTERFACE

END MODULE class_vector_field
