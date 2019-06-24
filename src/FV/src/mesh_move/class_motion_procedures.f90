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
! $Id: class_motion.f90 3323 2008-08-28 15:44:18Z sfilippo $
!
! Description:
!    To be added...
!
SUBMODULE(class_motion) class_motion_procedures

    USE class_vector
    USE class_psblas

    IMPLICIT NONE

CONTAINS


    MODULE PROCEDURE nemo_motion_sizeof
        USE psb_base_mod
        USE class_psblas
        INTEGER(kind=nemo_int_long_)   :: val

        val = 3 * nemo_sizeof_int
        IF (ALLOCATED(mot%law_x)) &
            & val = val + nemo_sizeof_dp * SIZE(mot%law_x)
        IF (ALLOCATED(mot%law_y)) &
            & val = val + SUM(mot%law_y%nemo_sizeof())
        nemo_motion_sizeof = val

    END PROCEDURE nemo_motion_sizeof



    ! ----- Constructor -----

    MODULE PROCEDURE create_motion
        USE class_psblas
        USE tools_mesh_move

        mot%surface_motion = surface_motion
        mot%vertex_motion = vertex_motion

        IF(surface_motion == moving_) THEN
            IF(TRIM(ml_file) == '') THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF
            CALL rd_inp_motion_law(ml_file,mot%iml,mot%law_x,mot%law_y)
        END IF

100     FORMAT(' ERROR! Motion-law file not specified for moving boundary')

    END PROCEDURE create_motion


    ! ----- Destructor -----

    MODULE PROCEDURE free_motion

        IF(ALLOCATED(mot%law_x)) DEALLOCATE(mot%law_x)
        IF(ALLOCATED(mot%law_y)) CALL free_vector(mot%law_y)

    END PROCEDURE free_motion


    ! ----- Getters -----

    MODULE PROCEDURE surface_motion_

        surface_motion_ = mot%surface_motion

    END PROCEDURE surface_motion_


    MODULE PROCEDURE vertex_motion_

        vertex_motion_ = mot%vertex_motion

    END PROCEDURE vertex_motion_


    MODULE PROCEDURE get_motion_displacement
        USE tools_math
        USE tools_mesh_move

        INTEGER :: i1, i2, i3, i4
        TYPE(vector) :: p1, p2, v_A, v_B, v_C, v_D
        TYPE(vector) :: i_trap_1, i_trap_2, i_trap_3

        SELECT CASE(mot%surface_motion)
        CASE(stationary_)
            get_motion_displacement = vector_(0.d0,0.d0,0.d0)
        CASE(moving_)
            SELECT CASE(mot%iml)
            CASE(ml_position_)

                CALL pwl_interp(p1,x1,mot%law_y,mot%law_x)
                CALL pwl_interp(p2,x2,mot%law_y,mot%law_x)
                get_motion_displacement = p2 - p1

            CASE(ml_velocity_)
                ! Integration: trapezoidal rule

                CALL pwl_nearest(x1,mot%law_x,i1,i2)
                CALL pwl_nearest(x2,mot%law_x,i3,i4)

                CALL pwl_interp(v_A,x1,mot%law_y,mot%law_x)
                v_B = mot%law_y(i2)
                v_C = mot%law_y(i3)
                CALL pwl_interp(v_D,x2,mot%law_y,mot%law_x)

                i_trap_1 = 0.5d0 * (mot%law_x(i2) - x1) * (v_A + v_B)
                i_trap_2 = 0.5d0 * (mot%law_x(i3) - mot%law_x(i2)) * (v_B + v_C)
                i_trap_3 = 0.5d0 * (x2 - mot%law_x(i3)) * (v_C + v_D)

                ! If i1 = i3 and i2 = i4: i_trap_2 < 0
                ! If i2 = i3            : i_trap_3 = 0

                get_motion_displacement = i_trap_1 + i_trap_2 + i_trap_3
            END SELECT
        END SELECT

    END PROCEDURE get_motion_displacement


    MODULE PROCEDURE get_motion_velocity
        USE tools_math
        USE tools_mesh_move
        IMPLICIT NONE
        !
        TYPE(vector) :: v

        SELECT CASE(mot%surface_motion)
        CASE(stationary_)
            get_motion_velocity = vector_(0.d0,0.d0,0.d0)
        CASE(moving_)
            SELECT CASE(mot%iml)
            CASE(ml_position_)
                CALL pwl_deriv(v,x,mot%law_y,mot%law_x)
            CASE(ml_velocity_)
                CALL pwl_interp(v,x,mot%law_y,mot%law_x)
            END SELECT
            get_motion_velocity = v
        END SELECT

    END PROCEDURE get_motion_velocity


    ! Moves the boundary vertices and the associated surface
    MODULE PROCEDURE move_boundary

        USE class_psblas
        USE class_connectivity
        USE class_mesh
        USE class_surface
        USE class_vertex
        USE tools_mesh_move, ONLY: stationary_, moving_, sticky_, sliding_

        ! Local variables

        INTEGER, POINTER :: iv2b(:) => NULL()
        TYPE(vertex)  :: bndry_vert      ! vertex attached to this boundary
        TYPE(vector)  :: vert_displ
        TYPE(vector)  :: location
        TYPE(vector)             :: normal
        INTEGER                  :: i,iv,nverts

        IF ( this_motion%surface_motion == stationary_ ) RETURN ! nothing to do !

        ! trap error:  we can't move unknown surface types with the sliding flag,
        ! since the normal is ill-defined

        IF (  ( type_(msh%surf(ib)) == iunknown_ ) .AND. &
            &( this_motion%vertex_motion == sliding_ ) ) THEN
            WRITE(6,100)
            CALL abort_psblas
        ENDIF

        CALL translate_surface(msh%surf(ib),displacement)

        ! Get list of vertices
        CALL get_ith_conn(iv2b,msh%v2b,ib)

        nverts = SIZE(iv2b)

        !loop over vertices on boundary
        DO i = 1,nverts
            iv = iv2b(i)

            ! get vertices on this boundary
            bndry_vert = msh%verts(iv)

            location  = position_(bndry_vert)

            IF ( this_motion%vertex_motion == sliding_ ) THEN

                ! have to put this in loop, since normal may vary
                normal = normal_( msh%surf(ib),location )
                vert_displ = (displacement .dot. normal) * normal

                ! sliding verts move with only the normal displacement
                bndry_vert = location + vert_displ

            ELSE

                !sticky verts move with total displacement
                bndry_vert = location + displacement

            ENDIF

            ! now clean up vertex locations...this would be necessary for sliding
            ! motion where a curved surface moves in a general direction
            ! if the motion is not orthogonal to the surface normal, the points
            ! move off the surface
            IF ( type_(msh%surf(ib)) /= iunknown_ ) &
                &call reform_vertex(msh%surf(ib), bndry_vert)

            ! save vertex location in the mesh
            msh%verts(iv) = bndry_vert

        ENDDO

        NULLIFY(iv2b)


100     FORMAT('Error! Cannot move boundary of unknown type with sliding vertices.')
    END PROCEDURE move_boundary

END SUBMODULE class_motion_procedures
