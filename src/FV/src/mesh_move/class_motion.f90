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
! $Id: class_motion.f90 3323 2008-08-28 15:44:18Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_motion

    USE class_vector
    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: motion ! Class
    PUBLIC :: move_boundary                      ! Setters


    TYPE motion
        PRIVATE
        INTEGER :: surface_motion
        INTEGER :: vertex_motion
        !
        INTEGER :: iml
        REAL(psb_dpk_), ALLOCATABLE :: law_x(:)
        TYPE(vector), ALLOCATABLE :: law_y(:)
    CONTAINS
        PROCEDURE, PUBLIC :: create_motion ! Constructor
        PROCEDURE :: free_motion ! Destructor
        PROCEDURE :: surface_motion_
        PROCEDURE, PRIVATE :: get_motion_displacement, get_motion_velocity
        GENERIC, PUBLIC :: get_displacement => get_motion_displacement
        GENERIC, PUBLIC :: get_velocity => get_motion_velocity
        PROCEDURE, PRIVATE :: nemo_motion_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_motion_sizeof
        PROCEDURE, PUBLIC :: vertex_motion_
    END TYPE motion


    INTERFACE

    ELEMENTAL MODULE FUNCTION nemo_motion_sizeof(mot)
        USE psb_base_mod
        USE class_psblas
        IMPLICIT NONE
        CLASS(motion), INTENT(IN) :: mot
        INTEGER(kind=nemo_int_long_)   :: nemo_motion_sizeof
    END FUNCTION nemo_motion_sizeof

    ! ----- Constructor -----

    MODULE SUBROUTINE create_motion(mot,surface_motion,vertex_motion,ml_file)
        IMPLICIT NONE
        CLASS(motion), INTENT(INOUT) :: mot
        INTEGER, INTENT(IN) :: surface_motion
        INTEGER, INTENT(IN) :: vertex_motion
        CHARACTER(len=*), INTENT(IN) :: ml_file
    END SUBROUTINE create_motion


    ! ----- Destructor -----

    MODULE SUBROUTINE free_motion(mot)
        IMPLICIT NONE
        CLASS(motion), INTENT(INOUT) :: mot
    END SUBROUTINE free_motion


    ! ----- Getters -----

    MODULE FUNCTION surface_motion_(mot)
        IMPLICIT NONE
        INTEGER :: surface_motion_
        CLASS(motion), INTENT(IN) :: mot
    END FUNCTION surface_motion_


    MODULE FUNCTION vertex_motion_(mot)
        IMPLICIT NONE
        INTEGER :: vertex_motion_
        CLASS(motion), INTENT(IN) :: mot
    END FUNCTION vertex_motion_

    MODULE FUNCTION get_motion_displacement(mot,x1,x2)
        USE tools_math
        USE tools_mesh_move
        IMPLICIT NONE
        TYPE(vector) :: get_motion_displacement
        CLASS(motion), INTENT(IN) :: mot
        REAL(psb_dpk_), INTENT(IN) :: x1, x2
    END FUNCTION get_motion_displacement

    MODULE FUNCTION get_motion_velocity(mot,x)
        USE tools_math
        USE tools_mesh_move
        IMPLICIT NONE
        TYPE(vector) :: get_motion_velocity
        CLASS(motion), INTENT(IN) :: mot
        REAL(psb_dpk_), INTENT(IN) :: x
    END FUNCTION get_motion_velocity

    MODULE SUBROUTINE move_boundary(ib,this_motion,displacement,msh)
      !! Moves the boundary vertices and the associated surface
        USE class_psblas
        !USE class_connectivity
        USE class_mesh
        !USE class_surface
        !USE class_vertex
        USE tools_mesh_move, ONLY: stationary_, moving_, sticky_, sliding_
        IMPLICIT NONE
        INTEGER, INTENT(IN)        :: ib
        TYPE(motion), INTENT(IN)   :: this_motion
        TYPE(vector),INTENT(IN)    :: displacement
        TYPE(mesh), INTENT(INOUT)  :: msh
    END SUBROUTINE move_boundary

  END INTERFACE

END MODULE class_motion
