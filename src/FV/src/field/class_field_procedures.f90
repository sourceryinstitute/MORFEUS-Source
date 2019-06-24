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
SUBMODULE(class_field) class_field_procedures

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_field_sizeof
        ! use psb_base_mod
        USE class_psblas

        !
        ! msh, bc and mat are independent objects.
        !
        nemo_field_sizeof = nemo_sizeof_int + LEN(fld%name)

    END PROCEDURE nemo_field_sizeof

    ! ----- Constructor -----

    MODULE PROCEDURE create_field
        USE tools_material

        ! Assigns mandatory arguments
        fld%msh => msh

        ! Assigns optional arguments
        IF(PRESENT(dim)) THEN
            fld%dim = dim
        ELSE
            fld%dim = null_dim_
        END IF

        IF(PRESENT(bc)) THEN
            fld%bc => bc
        ELSE
            fld%bc => NULL()
        END IF

        IF(PRESENT(mats)) THEN
            fld%mats => mats
        ELSE
            fld%mats => NULL()
        END IF

        IF(PRESENT(on_faces)) THEN
            fld%on_faces = on_faces
        ELSE
            fld%on_faces = .FALSE. ! Default is cell-centered
        END IF

        ! Field name
        fld%name = quantity(fld%dim)

    END PROCEDURE create_field


    ! ----- Destructor -----

    MODULE PROCEDURE free_field
        NULLIFY(fld%msh)
        NULLIFY(fld%bc)
        NULLIFY(fld%mats)
    END PROCEDURE free_field


    ! ----- Getters -----

    MODULE PROCEDURE get_field_name

        get_field_name = fld%name

    END PROCEDURE get_field_name


    MODULE PROCEDURE get_field_dim

        get_field_dim = fld%dim

    END PROCEDURE get_field_dim


    MODULE PROCEDURE get_field_msh_sub

        msh => fld%msh

    END PROCEDURE get_field_msh_sub

    MODULE PROCEDURE get_field_msh_fun

        get_field_msh_fun => fld%msh

    END PROCEDURE get_field_msh_fun


    MODULE PROCEDURE get_field_on_faces

        get_field_on_faces = fld%on_faces

    END PROCEDURE get_field_on_faces


    MODULE PROCEDURE get_field_bc

        get_field_bc => fld%bc

    END PROCEDURE get_field_bc


    MODULE PROCEDURE get_field_mat_sub

        IF (PRESENT(i) .AND. i < SIZE(fld%mats)) THEN
            mat => fld%mats(i)%mat
        ELSE
            mat => fld%mats(1)%mat
        END IF

    END PROCEDURE get_field_mat_sub
    !
    MODULE PROCEDURE get_field_mat_fun

        IF (PRESENT(i) .AND. i < SIZE(fld%mats)) THEN
            get_field_mat_fun => fld%mats(i)%mat
        ELSE
            get_field_mat_fun => fld%mats(1)%mat
        END IF

    END PROCEDURE get_field_mat_fun


    MODULE PROCEDURE get_field_size
        USE class_face

        ! Number of internal elements
        IF(fld%on_faces) THEN
            ! Face-centered
            isize(fld_internal_) = COUNT(flag_(fld%msh%faces) <= 0)
        ELSE
            ! Cell-centered
            isize(fld_internal_) = SIZE(fld%msh%cells)
        END IF

        ! Number of boundary faces
        isize(fld_boundary_) = COUNT(flag_(fld%msh%faces) > 0)

    END PROCEDURE get_field_size

    ! ----- Setter -----

    MODULE PROCEDURE set_field_dim

        fld%dim = dim

    END PROCEDURE set_field_dim

    MODULE PROCEDURE set_field_on_faces

        fld%on_faces = on_faces

    END PROCEDURE set_field_on_faces


    ! ----- Auxiliary Routines -----

    MODULE PROCEDURE check_mesh_consistency_bf
        USE class_mesh, ONLY : check_mesh_consistency

        CALL check_mesh_consistency(f1%msh,f2%msh,WHERE)

        ! bf = Base Field
    END PROCEDURE check_mesh_consistency_bf


    MODULE PROCEDURE check_field_operands
        USE class_mesh, ONLY : check_mesh_consistency
        INTEGER :: i

        CALL check_mesh_consistency(f1%msh,f2%msh,'CHECk_FIELD_OPERANDS')

        IF (SIZE(f1%mats) /= SIZE(f2%mats)) THEN
            WRITE(*,100) TRIM(WHERE)
            CALL abort_psblas
        END IF

        CALL check_material_consistency(f1%mats,f2%mats,'CHECK_FIELD_OPERANDS')

        IF(on_faces_(f1) .AND. .NOT.(on_faces_(f2)) &
            & .OR. (on_faces_(f2) .AND. .NOT.(on_faces_(f1)))) THEN
            WRITE(*,100) TRIM(WHERE)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Face-centered field operand in ',a)

    END PROCEDURE check_field_operands


END SUBMODULE class_field_procedures
