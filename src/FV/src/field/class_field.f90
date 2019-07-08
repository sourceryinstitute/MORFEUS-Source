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
! $Id: class_field.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Base class for scalar, vector and tensor field
!
MODULE class_field
    USE class_psblas, ONLY : nemo_int_long_, psb_dpk_
    USE class_bc, ONLY : bc_poly
    USE class_dimensions, ONLY : dimensions
    USE class_material, ONLY : matptr, material
    USE class_mesh, ONLY : mesh

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: field                           !! Class
    PUBLIC :: fld_internal_, fld_boundary_    !! Named constants

    TYPE field
        PRIVATE
        CHARACTER(len=32)       :: name
        TYPE(dimensions)        :: dim
        TYPE(mesh),     POINTER :: msh   => NULL()
        LOGICAL                 :: on_faces
        TYPE(bc_poly),  POINTER :: bc(:) => NULL()
        TYPE(matptr),   POINTER :: mats(:) => NULL()
    CONTAINS
        PROCEDURE :: create_field, free_field          !! Constructor/destructor
        PROCEDURE, PUBLIC :: on_faces_                 !! Getters
        PROCEDURE, PRIVATE :: get_field_size, get_field_mat_sub
        PROCEDURE, PUBLIC :: mat_
        PROCEDURE, PUBLIC :: bc_
        GENERIC, PUBLIC :: fld_size => get_field_size
        GENERIC, PUBLIC :: get_material => get_field_mat_sub
        PROCEDURE, PRIVATE :: get_field_dim            !! Getters
        PROCEDURE, PUBLIC :: msh_
        GENERIC, PUBLIC :: dim_ => get_field_dim
        PROCEDURE, PUBLIC :: name_
        PROCEDURE :: set_field_dim, set_field_on_faces !! Setters
        PROCEDURE :: check_field_operands
        PROCEDURE, PUBLIC :: nemo_sizeof
        PROCEDURE, PUBLIC :: get_mesh
        PROCEDURE, PRIVATE :: check_mesh_consistency_bf
        GENERIC, PUBLIC :: check_mesh_consistency => check_mesh_consistency_bf
    END TYPE field

    ! Default FIELD%ON_FACES = .false. => cell-centered

    ! ----- Named Constants -----

    INTEGER, PARAMETER :: fld_internal_ = 1
    INTEGER, PARAMETER :: fld_boundary_ = 2

    INTERFACE

        MODULE FUNCTION nemo_sizeof(fld)
            IMPLICIT NONE
            CLASS(field), INTENT(IN)     :: fld
            INTEGER(kind=nemo_int_long_) :: nemo_sizeof
        END FUNCTION nemo_sizeof

        MODULE SUBROUTINE create_field(fld,msh,dim,bc,mats,on_faces)
            !! Constructor
            IMPLICIT NONE
            !! Mandatory arguments
            CLASS(field),     INTENT(OUT)          :: fld
            TYPE(mesh),       INTENT(IN), TARGET   :: msh
            !! Optional arguments
            TYPE(dimensions), INTENT(IN), OPTIONAL :: dim
            TYPE(bc_poly),    INTENT(IN), OPTIONAL, TARGET :: bc(:)
            TYPE(matptr),     INTENT(IN), OPTIONAL, TARGET :: mats(:)
            LOGICAL,          INTENT(IN), OPTIONAL :: on_faces
        END SUBROUTINE create_field

        !! ----- Destructor -----

        MODULE SUBROUTINE free_field(fld)
            !! Destructor
            IMPLICIT NONE
            CLASS(field), INTENT(INOUT) :: fld
        END SUBROUTINE free_field

        !! ----- Getters -----

        MODULE FUNCTION name_(fld)
            IMPLICIT NONE
            CLASS(field), INTENT(IN) :: fld
            CHARACTER(len=32) :: name_
        END FUNCTION name_

        MODULE FUNCTION get_field_dim(fld)
            IMPLICIT NONE
            CLASS(field), INTENT(IN) :: fld
            TYPE(dimensions) :: get_field_dim
        END FUNCTION get_field_dim

        MODULE FUNCTION msh_(fld)
            IMPLICIT NONE
            CLASS(field), INTENT(IN), TARGET :: fld
            TYPE(mesh), POINTER :: msh_
        END FUNCTION msh_

        MODULE FUNCTION on_faces_(fld)
            IMPLICIT NONE
            CLASS(field), INTENT(IN) :: fld
            LOGICAL :: on_faces_
        END FUNCTION on_faces_

        MODULE FUNCTION bc_(fld)
            IMPLICIT NONE
            CLASS(field), INTENT(IN), TARGET  :: fld
            TYPE(bc_poly), POINTER :: bc_(:)
        END FUNCTION bc_

        MODULE FUNCTION mat_(fld, i)
            IMPLICIT NONE
            CLASS(field), INTENT(IN) :: fld
            INTEGER, INTENT(IN), OPTIONAL :: i
            TYPE(material), POINTER :: mat_
        END FUNCTION mat_

        MODULE FUNCTION get_field_size(fld) RESULT(isize)
            IMPLICIT NONE
            CLASS(field), INTENT(IN) :: fld
            INTEGER :: isize(2)
        END FUNCTION get_field_size

        !! ----- Temporary up to Gfortran patch -----
        MODULE SUBROUTINE get_mesh(fld,msh)
            IMPLICIT NONE
            CLASS(field), INTENT(IN) :: fld
            TYPE(mesh),   POINTER    :: msh
        END SUBROUTINE get_mesh

        MODULE SUBROUTINE get_field_mat_sub(fld,i,mat)
            IMPLICIT NONE
            CLASS(field), INTENT(IN) :: fld
            INTEGER,      INTENT(IN), OPTIONAL :: i
            TYPE(material), POINTER  :: mat
        END SUBROUTINE get_field_mat_sub

        !! Check operations
        MODULE SUBROUTINE check_mesh_consistency_bf(f1,f2,WHERE)
            IMPLICIT NONE
            CLASS(field),     INTENT(IN) :: f1
            TYPE(field),      INTENT(IN) :: f2
            CHARACTER(len=*), INTENT(IN) :: WHERE
        END SUBROUTINE check_mesh_consistency_bf

        !! ----- Setters -----

        MODULE SUBROUTINE set_field_dim(fld,dim)
            IMPLICIT NONE
            CLASS(field),     INTENT(INOUT) :: fld
            TYPE(dimensions), INTENT(IN)    :: dim
        END SUBROUTINE set_field_dim

        MODULE SUBROUTINE set_field_on_faces(fld,on_faces)
            IMPLICIT NONE
            CLASS(field), INTENT(INOUT) :: fld
            LOGICAL,      INTENT(IN)    :: on_faces
        END SUBROUTINE set_field_on_faces

        !! ----- Auxiliary Routines -----

        MODULE SUBROUTINE check_field_operands(f1,f2,WHERE)
            IMPLICIT NONE
            CLASS(field),     INTENT(IN) :: f1
            TYPE(field),      INTENT(IN) :: f2
            CHARACTER(len=*), INTENT(IN) :: WHERE
        END SUBROUTINE check_field_operands

    END INTERFACE

END MODULE class_field
