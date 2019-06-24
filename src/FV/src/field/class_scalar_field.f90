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
! $Id: class_scalar_field.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_scalar_field

    USE class_psblas, ONLY : psb_dpk_, nemo_int_long_, nemo_sizeof_dp, nemo_sizeof_int,&
        & icontxt_, psb_erractionsave, abort_psblas, psb_check_error, psb_erractionrestore,&
        & psb_halo
    USE class_field!, ONLY : field

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: scalar_field                              ! Class
    PUBLIC :: scalar_field_, create_field, free_field   ! Constructor/destructor
    PUBLIC :: name_, dim_, msh_, on_faces_, bc_, &      ! Getters
        &    get_material, get_base,&        !   "
        &    get_x, get_xp, get_bx, get_scalar_field_mat_id                             !   "
    PUBLIC :: update_field, ASSIGNMENT(=), &            ! Setters
        &    set_field_element, set_field_group        !   "
    PUBLIC :: interp_on_faces, OPERATOR(*)                          ! Algebra operations
    PUBLIC :: field_normi, field_norm1

    TYPE scalar_field
        PRIVATE
        TYPE(field) :: base
        REAL(psb_dpk_), ALLOCATABLE :: x(:)
        REAL(psb_dpk_), ALLOCATABLE :: xp(:)
        REAL(psb_dpk_), ALLOCATABLE :: bx(:)
        INTEGER, ALLOCATABLE :: mat(:)
        INTEGER, ALLOCATABLE :: bmat(:)
    CONTAINS
        PROCEDURE :: scalar_field_sum, scalar_field_dif, scalar_field_dif_s
        PROCEDURE :: scalar_field_div
        GENERIC :: OPERATOR(+) => scalar_field_sum
        GENERIC :: OPERATOR(-) => scalar_field_dif, scalar_field_dif_s
        GENERIC :: OPERATOR(/) => scalar_field_div
        PROCEDURE, PRIVATE :: nemo_scalar_field_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_scalar_field_sizeof
        PROCEDURE, PRIVATE :: get_scalar_field_msh_sub
        GENERIC, PUBLIC :: get_mesh => get_scalar_field_msh_sub
        PROCEDURE, PRIVATE :: check_mesh_consistency_sf
        GENERIC, PUBLIC :: check_mesh_consistency => check_mesh_consistency_sf
    END TYPE scalar_field

    INTERFACE

    MODULE FUNCTION nemo_scalar_field_sizeof(fld)
        !use psb_base_mod
        IMPLICIT NONE
        CLASS(scalar_field), INTENT(IN) :: fld
        INTEGER(kind=nemo_int_long_)   :: nemo_scalar_field_sizeof
    END FUNCTION nemo_scalar_field_sizeof

    ! ----- Constructor -----

    MODULE FUNCTION scalar_field_(base,x,bx)
      !! Default public constructor, necessary with ifort
        IMPLICIT NONE
        TYPE(scalar_field) :: scalar_field_
        TYPE(field),      INTENT(IN) :: base
        REAL(psb_dpk_), INTENT(IN) :: x(:)
        REAL(psb_dpk_), INTENT(IN) :: bx(:)
    END FUNCTION scalar_field_
  END INTERFACE


  ! ----- Generic Interfaces -----

  ! Constructor
  INTERFACE create_field
    MODULE SUBROUTINE create_scalar_field(fld,msh,dim,bc,mats,on_faces,x0)
        USE class_bc
        USE class_connectivity
        USE class_dimensions
        USE class_material
        USE class_mesh
        USE class_cell
        USE class_face
        USE tools_material
        IMPLICIT NONE
        ! Mandatory arguments
        TYPE(scalar_field), INTENT(OUT)        :: fld
        TYPE(mesh),         INTENT(IN), TARGET :: msh
        !
        ! Optional arguments
        TYPE(dimensions), INTENT(IN), OPTIONAL         :: dim
        TYPE(bc_poly),    INTENT(IN), OPTIONAL, TARGET :: bc(:)
        TYPE(matptr),   INTENT(IN), OPTIONAL, TARGET :: mats(:)
        LOGICAL,          INTENT(IN), OPTIONAL         :: on_faces
        REAL(psb_dpk_), INTENT(IN), OPTIONAL         :: x0
    END SUBROUTINE create_scalar_field
  END INTERFACE create_field

  ! ----- Destructor -----

  INTERFACE free_field
    !! Destructor
    MODULE SUBROUTINE free_scalar_field(fld)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(INOUT) :: fld
    END SUBROUTINE free_scalar_field
  END INTERFACE free_field

  ! ----- Getters for Inherited Members -----

  INTERFACE name_
    MODULE FUNCTION get_scalar_field_name(fld)
        IMPLICIT NONE
        CHARACTER(len=32) :: get_scalar_field_name
        TYPE(scalar_field), INTENT(IN) :: fld
    END FUNCTION get_scalar_field_name
  END INTERFACE name_

  INTERFACE dim_
    MODULE FUNCTION get_scalar_field_dim(fld)
        USE class_dimensions
        IMPLICIT NONE
        TYPE(dimensions) :: get_scalar_field_dim
        TYPE(scalar_field), INTENT(IN) :: fld
    END FUNCTION get_scalar_field_dim
  END INTERFACE dim_

  INTERFACE msh_
    MODULE FUNCTION get_scalar_field_msh_fun(fld)
        USE class_mesh
        IMPLICIT NONE
        TYPE(mesh), POINTER :: get_scalar_field_msh_fun
        TYPE(scalar_field), INTENT(IN), TARGET  :: fld
    END FUNCTION get_scalar_field_msh_fun
  END INTERFACE msh_

  ! ----- Temporary up to Gfortran patch -----
  INTERFACE
    MODULE SUBROUTINE get_scalar_field_msh_sub(fld,msh)
        USE class_mesh
        IMPLICIT NONE
        CLASS(scalar_field), INTENT(IN) :: fld
        TYPE(mesh), POINTER :: msh
    END SUBROUTINE get_scalar_field_msh_sub
  END INTERFACE
  ! ------------------------------------------

  INTERFACE on_faces_
    MODULE FUNCTION get_scalar_field_on_faces(fld)
        IMPLICIT NONE
        LOGICAL :: get_scalar_field_on_faces
        TYPE(scalar_field), INTENT(IN) :: fld
    END FUNCTION get_scalar_field_on_faces
  END INTERFACE on_faces_

  INTERFACE bc_
    MODULE FUNCTION get_scalar_field_bc(fld)
        USE class_bc
        IMPLICIT NONE
        TYPE(bc_poly), POINTER :: get_scalar_field_bc(:)
        TYPE(scalar_field), INTENT(IN) :: fld
    END FUNCTION get_scalar_field_bc
  END INTERFACE bc_


  INTERFACE mat_
    MODULE FUNCTION get_scalar_field_mat(fld, i)
        USE class_material
        IMPLICIT NONE
        TYPE(material), POINTER :: get_scalar_field_mat
        INTEGER, INTENT(IN), OPTIONAL :: i
        TYPE(scalar_field), INTENT(IN) :: fld
    END FUNCTION get_scalar_field_mat
  END INTERFACE mat_


  ! ----- Temporary up to Gfortran patch -----
  INTERFACE get_material
    MODULE SUBROUTINE get_scalar_field_mat_sub(fld,i,mat)
      USE class_material
      IMPLICIT NONE
      TYPE(scalar_field), INTENT(IN) :: fld
      INTEGER, INTENT(IN), OPTIONAL :: i
      TYPE(material), POINTER :: mat
    END SUBROUTINE get_scalar_field_mat_sub
  END INTERFACE get_material
  ! ------------------------------------------

  INTERFACE
    INTEGER MODULE FUNCTION get_scalar_field_mat_id(fld,i)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(IN) :: fld
        INTEGER, INTENT(IN) :: i
    END FUNCTION get_scalar_field_mat_id
  END INTERFACE

  INTERFACE get_base
    MODULE SUBROUTINE get_scalar_field_base(fld,base)
        TYPE(scalar_field), INTENT(IN)  :: fld
        TYPE(field),        INTENT(OUT) :: base
    END SUBROUTINE get_scalar_field_base
  END INTERFACE get_base


  ! ----- Getters for Additional Members -----

  INTERFACE get_x

    MODULE SUBROUTINE get_scalar_field_x(fld,x)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(IN) :: fld
        REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: x(:)
    END SUBROUTINE get_scalar_field_x

    MODULE SUBROUTINE get_scalar_field_element(fld,x,i)
        TYPE(scalar_field), INTENT(IN) :: fld
        REAL(psb_dpk_),   INTENT(OUT) :: x
        INTEGER, INTENT(IN) :: i
    END SUBROUTINE get_scalar_field_element

  END INTERFACE get_x


  INTERFACE get_xp

    MODULE SUBROUTINE get_scalar_field_xp(fld,xp)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(IN) :: fld
        REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: xp(:)
    END SUBROUTINE get_scalar_field_xp

    MODULE SUBROUTINE get_scalar_field_element_prev(fld,xp,i)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(IN) :: fld
        REAL(psb_dpk_),   INTENT(OUT) :: xp
        INTEGER, INTENT(IN) :: i
    END SUBROUTINE get_scalar_field_element_prev

  END INTERFACE get_xp


  INTERFACE get_bx
    MODULE SUBROUTINE get_scalar_field_bx(fld,bx)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(IN) :: fld
        REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: bx(:)
    END SUBROUTINE get_scalar_field_bx
  END INTERFACE get_bx

    ! ----- Setters -----

  INTERFACE update_field
    MODULE SUBROUTINE update_scalar_field(fld,mats,temp)
        USE class_bc
        USE class_connectivity
        USE class_dimensions
        USE class_face
        USE class_material
        USE class_mesh
        USE tools_math
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(INOUT) :: fld
        TYPE(scalar_field), INTENT(IN), OPTIONAL :: temp
        TYPE(matptr), INTENT(IN), POINTER :: mats(:)
    END SUBROUTINE update_scalar_field
  END INTERFACE update_field

  !
  ! REMARK: thanks to TR15581 extensions there is no need of defining
  ! the corresponding routine for the assignment between two fields.
  ! f1 = f2 is legal and it does not require an explicit creation of f1

  INTERFACE ASSIGNMENT(=)
    !! User-defined assignment

    MODULE SUBROUTINE assign_scalar_field_s(f,x)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(INOUT) :: f
        REAL(psb_dpk_),   INTENT(IN)    :: x
    END SUBROUTINE assign_scalar_field_s

    MODULE SUBROUTINE assign_scalar_field_v(f,x)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(INOUT) :: f
        REAL(psb_dpk_),   INTENT(IN)    :: x(:)
    END SUBROUTINE assign_scalar_field_v

  END INTERFACE ASSIGNMENT(=)


  ! Setters
  INTERFACE set_field_element
    MODULE SUBROUTINE set_scalar_field_element(f,i,x)
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(INOUT) :: f
        INTEGER, INTENT(IN) :: i
        REAL(psb_dpk_), INTENT(IN) :: x
    END SUBROUTINE set_scalar_field_element
  END INTERFACE set_field_element


  INTERFACE set_field_group
    MODULE SUBROUTINE set_scalar_field_group(f,ig,x)
        USE class_connectivity
        USE class_mesh
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(INOUT) :: f
        INTEGER, INTENT(IN) :: ig
        REAL(psb_dpk_), INTENT(IN) :: x
    END SUBROUTINE set_scalar_field_group
  END INTERFACE set_field_group


  ! ----- Algebra Operations -----

  INTERFACE field_normi
    MODULE FUNCTION nemo_scalar_field_normi(fld) RESULT(norm)
        USE class_mesh
        USE tools_psblas
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(IN) :: fld
        REAL(psb_dpk_)                 :: norm
    END FUNCTION nemo_scalar_field_normi
  END INTERFACE field_normi


  INTERFACE field_norm1
    MODULE FUNCTION nemo_scalar_field_norm1(fld) RESULT(norm)
        USE class_mesh
        USE tools_psblas
        IMPLICIT NONE
        TYPE(scalar_field), INTENT(IN) :: fld
        REAL(psb_dpk_)                 :: norm
    END FUNCTION nemo_scalar_field_norm1
  END INTERFACE field_norm1


  ! Algebra Operations

  INTERFACE
    MODULE FUNCTION scalar_field_sum(f1,f2)RESULT(r)
        USE class_dimensions
        IMPLICIT NONE
        TYPE(scalar_field) :: r
        CLASS(scalar_field), INTENT(IN) :: f1
        TYPE(scalar_field), INTENT(IN) :: f2
    END FUNCTION scalar_field_sum
  END INTERFACE


  INTERFACE

    MODULE FUNCTION scalar_field_dif(f1,f2)RESULT(r)
        USE class_dimensions
        IMPLICIT NONE
        TYPE(scalar_field) :: r
        CLASS(scalar_field), INTENT(IN) :: f1
        TYPE(scalar_field), INTENT(IN) :: f2
    END FUNCTION scalar_field_dif

    MODULE FUNCTION scalar_field_dif_s(f1,f2)RESULT(r)
        IMPLICIT NONE
        TYPE(scalar_field) :: r
        CLASS(scalar_field), INTENT(IN) :: f1
        REAL(psb_dpk_), INTENT(IN) :: f2
    END FUNCTION scalar_field_dif_s

  END INTERFACE

  INTERFACE OPERATOR(*)

    MODULE FUNCTION scalar_field_scal(a,f)RESULT(r)
        IMPLICIT NONE
        TYPE(scalar_field) :: r
        REAL(psb_dpk_),   INTENT(IN) :: a
        CLASS(scalar_field), INTENT(IN) :: f
    END FUNCTION scalar_field_scal

    MODULE FUNCTION scalar_field_mul(f1,f2)RESULT(r)
        USE class_dimensions
        IMPLICIT NONE
        TYPE(scalar_field) :: r
        CLASS(scalar_field), INTENT(IN) :: f1
        TYPE(scalar_field), INTENT(IN) :: f2
    END FUNCTION scalar_field_mul

  END INTERFACE OPERATOR(*)


  INTERFACE
    MODULE FUNCTION scalar_field_div(f1,f2)RESULT(r)
        USE class_dimensions

        TYPE(scalar_field) :: r
        CLASS(scalar_field), INTENT(IN) :: f1
        TYPE(scalar_field), INTENT(IN) :: f2
    END FUNCTION scalar_field_div
  END INTERFACE


  INTERFACE interp_on_faces
    MODULE FUNCTION interp_on_faces_s(fld)RESULT(r)
        USE class_connectivity
        USE class_face
        USE class_mesh
        USE tools_math
        IMPLICIT NONE
        TYPE(scalar_field) :: r
        TYPE(scalar_field), INTENT(IN) :: fld
    END FUNCTION interp_on_faces_s
  END INTERFACE interp_on_faces


  ! ----- Check Procedures -----

  INTERFACE
    MODULE SUBROUTINE check_mesh_consistency_sf(f1,f2,WHERE)
        CLASS(scalar_field), INTENT(IN) :: f1
        TYPE(scalar_field), INTENT(IN) :: f2
        CHARACTER(len=*), INTENT(IN) :: WHERE
    END SUBROUTINE check_mesh_consistency_sf
  END INTERFACE

END MODULE class_scalar_field
