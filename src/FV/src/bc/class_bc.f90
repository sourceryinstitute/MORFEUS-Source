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
MODULE class_bc
    !! Emulate runtime polymorphism for boundary condition classes.

    USE class_psblas
    USE class_bc_math
    USE class_bc_wall
    USE class_motion

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: bc_poly                           ! Class
    PUBLIC :: create_bc, free_bc                ! Constructor/Destructor
    PUBLIC :: move_boundaries           ! Setter
    PUBLIC :: update_boundary                   ! Updater

    TYPE bc_poly
      !! Implement a runtime polymorphism pattern [1] to emulate an abstract parent of the bc_math/bc_wall classes.
      !! TODO: Modernize by 1. Refactoring math/wall so they extend bc_poly, 2. Replace all TYPE(bc_poly) with CLASS(bc_poly), and
      !! 3. Make bc_poly ABSTRACT.
      !!
      !! [1] Akin, E. (2003) Object-oriented programming via Fortran 90/95. Cambridge University Press.
        PRIVATE
        INTEGER :: id
        TYPE(motion) :: mot
        TYPE(bc_math), POINTER :: math => NULL()
        TYPE(bc_wall), POINTER :: wall => NULL()
    CONTAINS
        PROCEDURE, PRIVATE :: get_abc_s, get_abc_v
        GENERIC, PUBLIC :: get_abc => get_abc_s, get_abc_v
        PROCEDURE, PRIVATE :: set_bc_poly_map_s, set_bc_poly_map_v
        GENERIC, PUBLIC :: set_bc => set_bc_poly_map_s, set_bc_poly_map_v
        PROCEDURE, PRIVATE :: get_bc_surface_motion, get_bc_vertex_motion
        GENERIC, PUBLIC :: surface_motion_ => get_bc_surface_motion
        GENERIC, PUBLIC :: vertex_motion_ => get_bc_vertex_motion
        PROCEDURE, PRIVATE :: get_bc_motion_displacement
        GENERIC, PUBLIC :: get_displacement => get_bc_motion_displacement
        PROCEDURE, PRIVATE :: get_bc_motion_velocity
        GENERIC, PUBLIC :: get_velocity => get_bc_motion_velocity
        PROCEDURE, PRIVATE :: nemo_bc_poly_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_bc_poly_sizeof
    END TYPE bc_poly

    ! ----- Generic Interfaces -----

    INTERFACE

        ELEMENTAL MODULE FUNCTION nemo_bc_poly_sizeof(bc)
            USE class_psblas, ONLY : nemo_int_long_
            IMPLICIT NONE
            CLASS(bc_poly), INTENT(IN) :: bc
            INTEGER(kind=nemo_int_long_)   :: nemo_bc_poly_sizeof
        END FUNCTION nemo_bc_poly_sizeof

        MODULE SUBROUTINE create_bc(bc,input_file,sec,msh)
        !! Global constructor
            USE class_face
            USE class_mesh
            USE tools_bc
            USE tools_input
            TYPE(bc_poly), ALLOCATABLE :: bc(:)
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            TYPE(mesh), INTENT(IN) :: msh
        END SUBROUTINE create_bc

        MODULE SUBROUTINE free_bc(bc)
        !! ----- Destructor -----
            USE tools_bc
            TYPE(bc_poly), ALLOCATABLE  :: bc(:)
        END SUBROUTINE free_bc

        !! ----- Getters -----

        MODULE SUBROUTINE get_abc_s(bc,dim,id,a,b,c)
            USE class_dimensions
            USE tools_bc
            CLASS(bc_poly), INTENT(IN) :: bc
            TYPE(dimensions), INTENT(IN) :: dim
            INTEGER, INTENT(OUT) :: id
            REAL(psb_dpk_), INTENT(INOUT) :: a(:)
            REAL(psb_dpk_), INTENT(INOUT) :: b(:)
            REAL(psb_dpk_), INTENT(INOUT) :: c(:)
        END SUBROUTINE get_abc_s

        MODULE SUBROUTINE get_abc_v(bc,dim,id,a,b,c)
            USE class_dimensions
            USE class_vector
            USE tools_bc
            CLASS(bc_poly), INTENT(IN) :: bc
            TYPE(dimensions), INTENT(IN) :: dim
            INTEGER, INTENT(OUT) :: id
            REAL(psb_dpk_), INTENT(INOUT) :: a(:)
            REAL(psb_dpk_), INTENT(INOUT) :: b(:)
            TYPE(vector),     INTENT(INOUT) :: c(:)
        END SUBROUTINE get_abc_v

        MODULE FUNCTION get_bc_surface_motion(bc)
            INTEGER :: get_bc_surface_motion
            CLASS(bc_poly), INTENT(IN) :: bc
        END FUNCTION get_bc_surface_motion

        MODULE FUNCTION get_bc_vertex_motion(bc)
            IMPLICIT NONE
            CLASS(bc_poly), INTENT(IN) :: bc
            INTEGER :: get_bc_vertex_motion
        END FUNCTION get_bc_vertex_motion

        MODULE FUNCTION get_bc_motion_displacement(bc,x1,x2)RESULT(res)
            USE class_vector
            TYPE(vector) :: res
            CLASS(bc_poly), INTENT(IN) :: bc
            REAL(psb_dpk_), INTENT(IN) :: x1, x2
        END FUNCTION get_bc_motion_displacement

        MODULE FUNCTION get_bc_motion_velocity(bc,x)RESULT(res)
            USE class_vector
            TYPE(vector) :: res
            CLASS(bc_poly), INTENT(IN) :: bc
            REAL(psb_dpk_), INTENT(IN) :: x
        END FUNCTION get_bc_motion_velocity

        !! ----- Setters -----

        MODULE SUBROUTINE set_bc_poly_map_s(bc,i,a,b,c)
            USE tools_bc
            CLASS(bc_poly), INTENT(INOUT) :: bc
            INTEGER, INTENT(IN) :: i
            REAL(psb_dpk_), INTENT(IN) :: a, b, c
        END SUBROUTINE set_bc_poly_map_s

        MODULE SUBROUTINE set_bc_poly_map_v(bc,i,a,b,c)
            USE class_vector
            USE tools_bc
            CLASS(bc_poly), INTENT(INOUT) :: bc
            INTEGER, INTENT(IN) :: i
            REAL(psb_dpk_), INTENT(IN) :: a, b
            TYPE(vector), INTENT(IN) :: c
        END SUBROUTINE set_bc_poly_map_v

        MODULE SUBROUTINE move_boundaries(msh,bc,t1,t2)
        !! loop over all boundaries, moving the vertices and conceptual surfaces
        !! from the time interval t1 to t2
            USE class_mesh
            USE class_vector
            TYPE(mesh), INTENT(INOUT) :: msh
            TYPE(bc_poly), INTENT(IN) :: bc(:)
            REAL(psb_dpk_)         :: t1
            REAL(psb_dpk_)         :: t2
        END SUBROUTINE move_boundaries

    END INTERFACE


    INTERFACE update_boundary
        !! ----- Boundary Values Updater -----

        MODULE SUBROUTINE update_boundary_s(ib,bc,dim,msh,mats,im,x,bx)
            USE class_dimensions
            USE class_material
            USE class_mesh
            USE tools_bc
            INTEGER,          INTENT(IN) :: ib
            TYPE(bc_poly),    INTENT(IN) :: bc
            TYPE(dimensions), INTENT(IN) :: dim
            TYPE(mesh),       INTENT(IN) :: msh
            TYPE(matptr),   INTENT(IN),POINTER :: mats(:)
            INTEGER, INTENT(IN) :: im(:)
            REAL(psb_dpk_), INTENT(IN) :: x(:)
            REAL(psb_dpk_), INTENT(INOUT) :: bx(:)
        END SUBROUTINE update_boundary_s


        MODULE SUBROUTINE update_boundary_v(ib,bc,dim,msh,x,bx)
            USE class_dimensions
            USE class_mesh
            USE class_vector
            USE tools_bc
            INTEGER,          INTENT(IN) :: ib
            TYPE(bc_poly),    INTENT(IN) :: bc
            TYPE(dimensions), INTENT(IN) :: dim
            TYPE(mesh),       INTENT(IN) :: msh
            TYPE(vector),     INTENT(IN) :: x(:)
            TYPE(vector),     INTENT(INOUT) :: bx(:)
        END SUBROUTINE update_boundary_v

    END INTERFACE update_boundary

END MODULE class_bc
