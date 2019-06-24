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
    USE class_field
    USE class_vector

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: vector_field                               ! Class
    PUBLIC :: vector_field_, create_field, free_field    ! Constructor/destructor
    PUBLIC :: name_, dim_, msh_, on_faces_, bc_, &       ! Getters
        &    get_material, get_base, &        !   "
        &    get_x, get_xp, get_bx                      !   "
    PUBLIC :: update_field, ASSIGNMENT(=), &             ! Setters
        &    set_field_element, set_field_group,&       !   "
        &    set_field_bound_element                    !   "
    PUBLIC :: set_x
    PUBLIC :: OPERATOR(*)   ! Algebra operations
    PUBLIC :: interp_on_faces                            !   "

    TYPE vector_field
        PRIVATE
        TYPE(field) :: base
        TYPE(vector), ALLOCATABLE :: x(:)
        TYPE(vector), ALLOCATABLE :: xp(:)
        TYPE(vector), ALLOCATABLE :: bx(:)
        INTEGER, ALLOCATABLE :: mat(:)
    CONTAINS
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
        !use psb_base_mod
        IMPLICIT NONE
        CLASS(vector_field), INTENT(IN) :: fld
        INTEGER(kind=nemo_int_long_)   :: nemo_vector_field_sizeof
    END FUNCTION nemo_vector_field_sizeof

    ! ----- Constructor -----

    ! Default public constructor, necessary with ifort
    MODULE FUNCTION vector_field_(base,x,bx)
        TYPE(vector_field) :: vector_field_
        TYPE(field),      INTENT(IN) :: base
        TYPE(vector), INTENT(IN) :: x(:)
        TYPE(vector), INTENT(IN) :: bx(:)
    END FUNCTION vector_field_
  END INTERFACE


  ! Constructor
  INTERFACE create_field
    MODULE SUBROUTINE create_vector_field(fld,msh,dim,bc,mats,on_faces,x0)
        USE class_bc
        USE class_connectivity
        USE class_dimensions
        USE class_material
        USE class_mesh
        USE class_vector
        USE tools_material
        IMPLICIT NONE
        ! Mandatory arguments
        TYPE(vector_field), INTENT(OUT)        :: fld
        TYPE(mesh),         INTENT(IN), TARGET :: msh
        !
        ! Optional arguments
        TYPE(dimensions), INTENT(IN), OPTIONAL         :: dim
        TYPE(bc_poly),    INTENT(IN), OPTIONAL, TARGET :: bc(:)
        TYPE(matptr),   INTENT(IN), OPTIONAL, TARGET :: mats(:)
        LOGICAL,          INTENT(IN), OPTIONAL         :: on_faces
        TYPE(vector),     INTENT(IN), OPTIONAL         :: x0
    END SUBROUTINE create_vector_field
  END INTERFACE create_field


  ! ----- Destructor -----

  INTERFACE free_field
    !! Destructor
    MODULE SUBROUTINE free_vector_field(fld)
        IMPLICIT NONE
        TYPE(vector_field), INTENT(INOUT) :: fld
    END SUBROUTINE free_vector_field
  END INTERFACE free_field


  ! ----- Getters for Inherited Members -----

  INTERFACE name_
    MODULE FUNCTION get_vector_field_name(fld)
        IMPLICIT NONE
        CHARACTER(len=32) :: get_vector_field_name
        TYPE(vector_field), INTENT(IN) :: fld
    END FUNCTION get_vector_field_name
  END INTERFACE name_

  INTERFACE dim_
    MODULE FUNCTION get_vector_field_dim(fld)
        USE class_dimensions
        IMPLICIT NONE
        TYPE(dimensions) :: get_vector_field_dim
        TYPE(vector_field), INTENT(IN) :: fld
    END FUNCTION get_vector_field_dim
  END INTERFACE dim_


  INTERFACE msh_
    MODULE FUNCTION get_vector_field_msh_fun(fld)
        USE class_mesh
        IMPLICIT NONE
        TYPE(mesh), POINTER :: get_vector_field_msh_fun
        TYPE(vector_field), INTENT(IN), TARGET  :: fld
    END FUNCTION get_vector_field_msh_fun
  END INTERFACE msh_




  ! ----- Temporary up to Gfortran patch -----
  INTERFACE
    MODULE SUBROUTINE get_vector_field_msh_sub(fld,msh)
        USE class_mesh
        IMPLICIT NONE
        CLASS(vector_field), INTENT(IN) :: fld
        TYPE(mesh), POINTER :: msh
    END SUBROUTINE get_vector_field_msh_sub
  END INTERFACE
  ! ------------------------------------------


  INTERFACE on_faces_
    MODULE FUNCTION get_vector_field_on_faces(fld)
        IMPLICIT NONE
        LOGICAL :: get_vector_field_on_faces
        TYPE(vector_field), INTENT(IN) :: fld
    END FUNCTION get_vector_field_on_faces
  END INTERFACE on_faces_

  INTERFACE bc_
    MODULE FUNCTION get_vector_field_bc(fld)
        USE class_bc
        IMPLICIT NONE
        TYPE(bc_poly), POINTER :: get_vector_field_bc(:)
        TYPE(vector_field), INTENT(IN) :: fld
    END FUNCTION get_vector_field_bc
  END INTERFACE bc_


  INTERFACE mat_
    MODULE FUNCTION get_vector_field_mat(fld)
        USE class_material
        IMPLICIT NONE
        TYPE(material), POINTER :: get_vector_field_mat
        TYPE(vector_field), INTENT(IN) :: fld
    END FUNCTION get_vector_field_mat
  END INTERFACE mat_


  ! ----- Temporary up to Gfortran patch -----
  INTERFACE get_material
    MODULE SUBROUTINE get_vector_field_mat_sub(fld,i,mat)
        USE class_material
        IMPLICIT NONE
        TYPE(vector_field), INTENT(IN) :: fld
        INTEGER, INTENT(IN), OPTIONAL :: i
        TYPE(material), POINTER :: mat
    END SUBROUTINE get_vector_field_mat_sub
  END INTERFACE get_material
  ! ------------------------------------------


  INTERFACE get_base
    MODULE SUBROUTINE get_vector_field_base(fld,base)
        TYPE(vector_field), INTENT(IN)  :: fld
        TYPE(field),        INTENT(OUT) :: base
    END SUBROUTINE get_vector_field_base
  END INTERFACE get_base


  ! ----- Getters for Additional Members -----

  INTERFACE get_x
    !! Getters for Additional Members

    MODULE SUBROUTINE get_vector_field_x_r(fld,x)
        IMPLICIT NONE
        TYPE(vector_field), INTENT(IN) :: fld
        REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: x(:,:)
    END SUBROUTINE get_vector_field_x_r

    MODULE SUBROUTINE get_vector_field_x_v(fld,x)
        IMPLICIT NONE
        TYPE(vector_field), INTENT(IN) :: fld
        TYPE(vector),       INTENT(OUT), ALLOCATABLE :: x(:)
    END SUBROUTINE get_vector_field_x_v

  END INTERFACE get_x


  INTERFACE get_xp

    MODULE SUBROUTINE get_vector_field_xp_r(fld,xp)
        IMPLICIT NONE
        TYPE(vector_field), INTENT(IN) :: fld
        REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: xp(:,:)
    END SUBROUTINE get_vector_field_xp_r

    MODULE SUBROUTINE get_vector_field_xp_v(fld,xp)
        IMPLICIT NONE
        TYPE(vector_field), INTENT(IN) :: fld
        TYPE(vector),       INTENT(OUT), ALLOCATABLE :: xp(:)
    END SUBROUTINE get_vector_field_xp_v

  END INTERFACE get_xp


  INTERFACE get_bx

    MODULE SUBROUTINE get_vector_field_bx_r(fld,bx)
        TYPE(vector_field), INTENT(IN) :: fld
        REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: bx(:,:)
    END SUBROUTINE get_vector_field_bx_r

    MODULE SUBROUTINE get_vector_field_bx_v(fld,bx)
        IMPLICIT NONE
        TYPE(vector_field), INTENT(IN) :: fld
        TYPE(vector),       INTENT(OUT), ALLOCATABLE :: bx(:)
    END SUBROUTINE get_vector_field_bx_v

  END INTERFACE get_bx


  ! ----- Setters -----

  INTERFACE set_x
    !! Setters
    MODULE SUBROUTINE set_vector_field_x(fld,x)
        IMPLICIT NONE
        TYPE(vector_field), INTENT(INOUT) :: fld
        REAL(psb_dpk_),   INTENT(IN), ALLOCATABLE :: x(:,:)
    END SUBROUTINE set_vector_field_x
  END INTERFACE set_x


  INTERFACE update_field
    MODULE SUBROUTINE update_vector_field(fld)
        USE class_bc
        USE class_dimensions
        USE class_material
        USE class_mesh
        IMPLICIT NONE
        TYPE(vector_field), INTENT(INOUT) :: fld
    END SUBROUTINE update_vector_field
  END INTERFACE update_field


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


  INTERFACE set_field_element
    MODULE SUBROUTINE set_vector_field_element(f,i,x)
        USE class_vector
        IMPLICIT NONE
        TYPE(vector_field), INTENT(INOUT) :: f
        INTEGER, INTENT(IN) :: i
        TYPE(vector), INTENT(IN) :: x
    END SUBROUTINE set_vector_field_element
  END INTERFACE set_field_element

  INTERFACE set_field_bound_element
    MODULE SUBROUTINE set_vector_field_bound_element(f,i,x)
        USE class_vector
        IMPLICIT NONE
        TYPE(vector_field), INTENT(INOUT) :: f
        INTEGER, INTENT(IN) :: i
        TYPE(vector), INTENT(IN) :: x
    END SUBROUTINE set_vector_field_bound_element
  END INTERFACE set_field_bound_element


  INTERFACE set_field_group
    MODULE SUBROUTINE set_vector_field_group(f,ig,x)
        USE class_connectivity
        USE class_mesh
        USE class_vector
        IMPLICIT NONE
        TYPE(vector_field), INTENT(INOUT) :: f
        INTEGER, INTENT(IN) :: ig
        TYPE(vector), INTENT(IN) :: x
    END SUBROUTINE set_vector_field_group
  END INTERFACE set_field_group


  ! ----- Algebra Operations -----

  INTERFACE
    MODULE FUNCTION vector_field_sum(f1,f2)RESULT(r)
        USE class_dimensions
        !use class_vector
        IMPLICIT NONE
        TYPE(vector_field) :: r
        CLASS(vector_field), INTENT(IN) :: f1
        TYPE(vector_field), INTENT(IN) :: f2
    END FUNCTION vector_field_sum
  END INTERFACE


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

  INTERFACE
    MODULE FUNCTION vector_field_dif(f1,f2)RESULT(r)
        USE class_dimensions
        !use class_vector
        IMPLICIT NONE
        TYPE(vector_field) :: r
        CLASS(vector_field), INTENT(IN) :: f1
        TYPE(vector_field), INTENT(IN) :: f2
    END FUNCTION vector_field_dif
  END INTERFACE


  INTERFACE interp_on_faces
    MODULE FUNCTION interp_on_faces_v(fld)RESULT(r)
        USE class_connectivity
        USE class_face
        USE class_mesh
        USE tools_math
        IMPLICIT NONE
        TYPE(vector_field) :: r
        TYPE(vector_field), INTENT(IN) :: fld
    END FUNCTION interp_on_faces_v
  END INTERFACE interp_on_faces


    ! ----- Check Procedures -----

  INTERFACE
    MODULE SUBROUTINE check_mesh_consistency_vf(f1,f2,WHERE)
        IMPLICIT NONE
        CLASS(vector_field), INTENT(IN) :: f1
        TYPE(vector_field), INTENT(IN) :: f2
        CHARACTER(len=*), INTENT(IN) :: WHERE
    END SUBROUTINE check_mesh_consistency_vf
  END INTERFACE

END MODULE class_vector_field
