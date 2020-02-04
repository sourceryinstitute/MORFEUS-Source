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
! $Id: class_bc.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Run-time polymorphism for BC_* classes
!
SUBMODULE(class_bc) class_bc_procedures
    USE class_face
    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_bc_poly_sizeof
        USE psb_base_mod

        INTEGER(kind=nemo_int_long_)   :: val

        val = nemo_sizeof_int
        val = val + bc%mot%nemo_sizeof()
        IF (ASSOCIATED(bc%math)) val = val + bc%math%nemo_sizeof()
        IF (ASSOCIATED(bc%wall)) val = val + bc%wall%nemo_sizeof()
        nemo_bc_poly_sizeof = val

    END PROCEDURE nemo_bc_poly_sizeof


    ! ----- Constuctors -----

    MODULE PROCEDURE create_bc
    !! Global constructor
        USE class_mesh
        USE tools_bc
        USE tools_input
        !
        CHARACTER(LEN=2) :: ib_string
        CHARACTER(len=32) :: bc_sec
        INTEGER :: info, mypnum
        INTEGER :: ib, nbf
        INTEGER :: digits
        CHARACTER(len=8) :: aformat

        mypnum = mypnum_()

        ! Checks whether MSH has been already set.
        IF(.NOT.msh%set) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Allocates and sizes BC(:) pointer
        IF(ALLOCATED(bc)) THEN
            WRITE(*,200)
            CALL abort_psblas
        ELSE
            ALLOCATE(bc(msh%nbc),stat=info)
            IF(info /= 0) THEN
                WRITE(*,300)
                CALL abort_psblas
            END IF
        END IF

        IF(mypnum == 0) THEN
            WRITE(*,*) 'Reading boundary conditions from ', TRIM(input_file)
        END IF

        ! Reads BC(:)%ID
        CALL rd_inp_bc(input_file,sec,msh%nbc,bc(:)%id,bc(:)%mot)


        ! ----- POLYMORPHISM -----
        ! Allocates and reads proper BC according to BC(:)%ID
        DO ib = 1, msh%nbc

            ! Composes BC section name
            WRITE(ib_string, '(i0)') ib
            bc_sec = "MORFEUS_FV.BCS("//trim(ib_string)//")"

            ! Counts number of boundary faces
            nbf = COUNT(msh%faces%flag_() == ib)
            SELECT CASE(bc(ib)%id)
            CASE(bc_math_)
                CALL create_bc_math(bc(ib)%math,input_file,bc_sec,nbf)
            CASE(bc_wall_)
                CALL create_bc_wall(bc(ib)%wall,input_file,bc_sec,nbf)
            CASE DEFAULT
                WRITE(*,400)
                CALL abort_psblas
            END SELECT
        END DO
        ! ------------------------

        IF(mypnum == 0) THEN
            WRITE(*,*)
        END IF

100     FORMAT(' ERROR! MESH object must be set before calling CREATE_BC')
200     FORMAT(' ERROR! Illegal call to CREATE_BC: BC pointer is already associated')
300     FORMAT(' ERROR! Memory allocation failure in CREATE_BC')
400     FORMAT(' ERROR! Unsupported type of boundary condition in CREATE_BC')

    END PROCEDURE create_bc


    ! ----- Destructor -----

    MODULE PROCEDURE free_bc
        USE tools_bc

        INTEGER :: info
        INTEGER :: ib, nbc

        nbc = SIZE(bc)

        DO ib = 1, nbc

            ! ----- POLYMORPHISM -----
            SELECT CASE(bc(ib)%id)
            CASE(bc_math_)
                CALL free_bc_math(bc(ib)%math)
            CASE(bc_wall_)
                CALL free_bc_wall(bc(ib)%wall)
            END SELECT
            ! ------------------------

            CALL bc(ib)%mot%free_motion()
        END DO

        DEALLOCATE(bc,stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! Memory deallocation failure in FREE_BC')

    END PROCEDURE free_bc


    ! ----- Getter -----

    MODULE PROCEDURE get_abc_s
        USE class_dimensions
        USE tools_bc

        ! ----- POLYMORPHISM -----
        SELECT CASE(bc%id)
        CASE(bc_math_)
            CALL get_abc_math(bc%math,id,a,b,c)
        CASE(bc_wall_)
            CALL bc%wall%get_abc_wall(dim,id,a,b,c)
        END SELECT
        ! ------------------------

    END PROCEDURE get_abc_s


    MODULE PROCEDURE get_abc_v
        USE class_dimensions
        USE class_vector
        USE tools_bc

        ! ----- POLYMORPHISM -----
        SELECT CASE(bc%id)
        CASE(bc_wall_)
            CALL bc%wall%get_abc_wall(dim,id,a,b,c)
        END SELECT
        ! ------------------------

    END PROCEDURE get_abc_v


    MODULE PROCEDURE get_bc_surface_motion

        get_bc_surface_motion = bc%mot%surface_motion_()

    END PROCEDURE get_bc_surface_motion


    MODULE PROCEDURE get_bc_vertex_motion

        get_bc_vertex_motion = bc%mot%vertex_motion_()

    END PROCEDURE get_bc_vertex_motion


    MODULE PROCEDURE get_bc_motion_displacement
        USE class_vector

        res = bc%mot%get_displacement(x1,x2)

    END PROCEDURE get_bc_motion_displacement


    MODULE PROCEDURE get_bc_motion_velocity
        USE class_vector

        res = bc%mot%get_velocity(x)

    END PROCEDURE get_bc_motion_velocity


    ! ----- Setters -----

    MODULE PROCEDURE set_bc_poly_map_s
        USE tools_bc

        ! ----- POLYMORPHISM -----
        SELECT CASE(bc%id)
        CASE(bc_math_)
            CALL bc%math%set_bc_math_map(i,a,b,c)
        CASE(bc_wall_)
            CALL bc%wall%set_bc_wall_map(i,a,b,c)
        CASE DEFAULT
            WRITE(*,100)
            CALL abort_psblas
        END SELECT
        ! ------------------------

100     FORMAT(' ERROR! Unsupported type of BC in SET_BC_POLY_MAP')

    END PROCEDURE set_bc_poly_map_s

    MODULE PROCEDURE set_bc_poly_map_v
        USE class_vector
        USE tools_bc

        ! ----- POLYMORPHISM -----
        SELECT CASE(bc%id)
        CASE(bc_wall_)
            CALL bc%wall%set_bc_wall_map(i,a,b,c)
        CASE DEFAULT
            WRITE(*,100)
            CALL abort_psblas
        END SELECT
        ! ------------------------

100     FORMAT(' ERROR! Unsupported type of BC in SET_BC_POLY_MAP')

    END PROCEDURE set_bc_poly_map_v


    ! ----- Boundary Values Updater -----

    MODULE PROCEDURE update_boundary_s
        USE class_dimensions
        USE class_material
        USE class_mesh
        USE tools_bc

        ! ----- POLYMORPHISM -----
        SELECT CASE(bc%id)
        CASE(bc_math_)
            CALL update_boundary_math(ib,bc%math,msh,x,bx)
        CASE(bc_wall_)
            CALL update_boundary_wall(ib,bc%wall,dim,msh,mats,im,x,bx)
        CASE DEFAULT
            WRITE(*,100)
            CALL abort_psblas
        END SELECT
        ! ------------------------

100     FORMAT(' ERROR! Unsupported BC type in UPDATE_BOUNDARY_S')

    END PROCEDURE update_boundary_s


    MODULE PROCEDURE update_boundary_v
        USE class_dimensions
        USE class_mesh
        USE class_vector
        USE tools_bc

        ! ----- POLYMORPHISM -----
        SELECT CASE(bc%id)
        CASE(bc_wall_)
            CALL update_boundary_wall(ib,bc%wall,dim,msh,x,bx)
        CASE DEFAULT
            WRITE(*,100)
            CALL abort_psblas
        END SELECT
        ! ------------------------

100     FORMAT(' ERROR! Unsupported BC type in UPDATE_BOUNDARY_V')

    END PROCEDURE update_boundary_v


    MODULE PROCEDURE move_boundaries
    !! loop over all boundaries, moving the vertices and conceptual surfaces
    !! from the time interval t1 to t2
        USE class_mesh
        USE class_vector

        ! Local variables
        INTEGER :: ib
        TYPE(motion) :: this_motion
        TYPE(vector) :: offset

        DO ib = 1 , msh%nbc

            this_motion = bc(ib)%mot
            offset = this_motion%get_displacement(t1,t2)

            ! move each single boundary
            CALL move_boundary(ib,this_motion,offset,msh)
        ENDDO

    END PROCEDURE move_boundaries

END SUBMODULE class_bc_procedures
