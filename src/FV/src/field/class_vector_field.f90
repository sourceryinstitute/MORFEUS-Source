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
    PUBLIC :: vector_field_     !! Constructor
    PUBLIC :: ASSIGNMENT(=)     !! Setters
    PUBLIC :: OPERATOR(*)       !! Algebra operations

    TYPE vector_field
        PRIVATE
        TYPE(field) :: base
        TYPE(vector), ALLOCATABLE :: x(:)
        TYPE(vector), ALLOCATABLE :: xp(:)
        TYPE(vector), ALLOCATABLE :: bx(:)
        INTEGER, ALLOCATABLE :: mat(:)
    CONTAINS
        PROCEDURE, PRIVATE :: create_vector_field, free_vector_field ! Constructor/destructor
        GENERIC, PUBLIC :: create_field => create_vector_field
        GENERIC, PUBLIC :: free_field => free_vector_field
        PROCEDURE, PRIVATE :: get_vector_field_on_faces, get_vector_field_bc
        GENERIC, PUBLIC :: on_faces_ => get_vector_field_on_faces
        GENERIC, PUBLIC :: bc_ => get_vector_field_bc
        PROCEDURE, PRIVATE :: get_vector_field_mat, get_vector_field_mat_sub
        GENERIC, PUBLIC :: mat_ => get_vector_field_mat
        GENERIC, PUBLIC :: get_material => get_vector_field_mat_sub
        PROCEDURE, PRIVATE :: get_vector_field_dim, get_vector_field_msh_fun ! Getters
        GENERIC, PUBLIC :: dim_ => get_vector_field_dim
        GENERIC, PUBLIC :: msh_ => get_vector_field_msh_fun
        PROCEDURE, PRIVATE :: get_vector_field_name
        GENERIC, PUBLIC :: name_ => get_vector_field_name
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
        PROCEDURE :: vector_field_sum, vector_field_dif
        GENERIC :: OPERATOR(+) => vector_field_sum
        GENERIC :: OPERATOR(-) => vector_field_dif
        PROCEDURE, PRIVATE :: nemo_vector_field_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_vector_field_sizeof
        PROCEDURE, PRIVATE :: get_vector_field_msh_sub
        GENERIC, PUBLIC :: get_mesh => get_vector_field_msh_sub
        PROCEDURE, PRIVATE :: check_mesh_consistency_vf
        GENERIC, PUBLIC :: check_mesh_consistency => check_mesh_consistency_vf
    END TYPE vector_field

    ! ----- Generic Interfaces -----

    INTERFACE

        MODULE FUNCTION nemo_vector_field_sizeof(fld)
            IMPLICIT NONE
            CLASS(vector_field), INTENT(IN) :: fld
            INTEGER(kind=nemo_int_long_)   :: nemo_vector_field_sizeof
        END FUNCTION nemo_vector_field_sizeof

        ! ----- Constructor -----

        ! Default public constructor, necessary with ifort
        MODULE FUNCTION vector_field_(base,x,bx)
            IMPLICIT NONE
            TYPE(vector_field) :: vector_field_
            TYPE(field),      INTENT(IN) :: base
            TYPE(vector), INTENT(IN) :: x(:)
            TYPE(vector), INTENT(IN) :: bx(:)
        END FUNCTION vector_field_

        ! Constructor
        MODULE SUBROUTINE create_vector_field(fld,msh,dim,bc,mats,on_faces,x0)
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

        !! Destructor
        MODULE SUBROUTINE free_vector_field(fld)
            IMPLICIT NONE
            CLASS(vector_field), INTENT(INOUT) :: fld
        END SUBROUTINE free_vector_field


        ! ----- Getters for Inherited Members -----

        MODULE FUNCTION get_vector_field_name(fld)
            IMPLICIT NONE
            CHARACTER(len=32) :: get_vector_field_name
            CLASS(vector_field), INTENT(IN) :: fld
        END FUNCTION get_vector_field_name

        MODULE FUNCTION get_vector_field_dim(fld)
            IMPLICIT NONE
            TYPE(dimensions) :: get_vector_field_dim
            CLASS(vector_field), INTENT(IN) :: fld
        END FUNCTION get_vector_field_dim

        MODULE FUNCTION get_vector_field_msh_fun(fld)
            IMPLICIT NONE
            TYPE(mesh), POINTER :: get_vector_field_msh_fun
            CLASS(vector_field), INTENT(IN), TARGET  :: fld
        END FUNCTION get_vector_field_msh_fun

        ! ----- Temporary up to Gfortran patch -----
        MODULE SUBROUTINE get_vector_field_msh_sub(fld,msh)
            IMPLICIT NONE
            CLASS(vector_field), INTENT(IN) :: fld
            TYPE(mesh), POINTER :: msh
        END SUBROUTINE get_vector_field_msh_sub
        ! ------------------------------------------


        MODULE FUNCTION get_vector_field_on_faces(fld)
            IMPLICIT NONE
            LOGICAL :: get_vector_field_on_faces
            CLASS(vector_field), INTENT(IN) :: fld
        END FUNCTION get_vector_field_on_faces

        MODULE FUNCTION get_vector_field_bc(fld)
            IMPLICIT NONE
            TYPE(bc_poly), POINTER :: get_vector_field_bc(:)
            CLASS(vector_field), INTENT(IN) :: fld
        END FUNCTION get_vector_field_bc

        MODULE FUNCTION get_vector_field_mat(fld)
            IMPLICIT NONE
            TYPE(material), POINTER :: get_vector_field_mat
            CLASS(vector_field), INTENT(IN) :: fld
        END FUNCTION get_vector_field_mat

        ! ----- Temporary up to Gfortran patch -----
        MODULE SUBROUTINE get_vector_field_mat_sub(fld,i,mat)
            IMPLICIT NONE
            CLASS(vector_field), INTENT(IN) :: fld
            INTEGER, INTENT(IN), OPTIONAL :: i
            TYPE(material), POINTER :: mat
        END SUBROUTINE get_vector_field_mat_sub
        ! ------------------------------------------

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

        !! Setters
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

    END INTERFACE

    INTERFACE ASSIGNMENT(=)
        !! User-defined assignemnts

        MODULE SUBROUTINE assign_vector_field_s(f,x)
            USE class_vector
            IMPLICIT NONE
            TYPE(vector_field), INTENT(INOUT) :: f
            TYPE(vector),       INTENT(IN)    :: x
        END SUBROUTINE assign_vector_field_s

        MODULE SUBROUTINE assign_vector_field_v(f,x)
            USE class_vector
            IMPLICIT NONE
            TYPE(vector_field), INTENT(INOUT) :: f
            TYPE(vector),       INTENT(IN)    :: x(:)
        END SUBROUTINE assign_vector_field_v

    END INTERFACE ASSIGNMENT(=)

    INTERFACE OPERATOR(*)
        MODULE FUNCTION vector_field_scal(a,f2)RESULT(r)
            USE class_dimensions
            !use class_vector
            IMPLICIT NONE
            TYPE(vector_field) :: r
            REAL(psb_dpk_), INTENT(IN) :: a
            TYPE(vector_field), INTENT(IN) :: f2
        END FUNCTION vector_field_scal
    END INTERFACE OPERATOR(*)

END MODULE class_vector_field
