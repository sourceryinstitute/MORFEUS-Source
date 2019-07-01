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
! $Id: class_cylinder.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    A class that handles the geometric concept of a cylinder.
!
! Provides:
!    CYLINDER              class.
!    ALLOC_CYLINDER        constructor for CYLINDER class.
!    FREE_CYLINDER         destructor for CYLINDER class.
!    GET_CYLINDER_NORMAL   returns the normal of the cylinder at a certain point
!    GET_CYLINDER_R2       returns how good a fit the cylinder is
!    GET_PT_CYLINDER       returns the nearest point that lies on the cylinder
!    TRANSLATE_CYLINDER    moves a cylinder in 3D space, translation only

! To be done :: a move cylinder function

MODULE class_cylinder
    USE class_vector
    USE class_vertex
    USE class_psblas, ONLY : psb_dpk_

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: cylinder ! Class
    PUBLIC :: alloc_cylinder, free_cylinder          ! Constructor/Destructor

    ! The cylinder is arbitrarily rotated and sized, with infinite length.
    TYPE cylinder
        PRIVATE

        TYPE(vector)      :: center       ! the location of the surface
        TYPE(vector)      :: axis         ! a unit vector pointing along the axis
        REAL(psb_dpk_) :: radius       ! the radius of the cylinder
        REAL(psb_dpk_) :: r2           ! The correlation parameter for the fit
    CONTAINS
        PROCEDURE :: get_cylinder_normal, get_cylinder_r2   ! Getters
        PROCEDURE :: get_pt_cylinder                        ! Getters, cont.
        PROCEDURE :: translate_cylinder                     ! Setters
        PROCEDURE, PRIVATE :: nemo_cylinder_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_cylinder_sizeof
    END TYPE cylinder

    TYPE(vertex),ALLOCATABLE      :: my_vertices(:)
    !! Local copy of vertices, available to function to be optimized.

    ! ----- Generic Interfaces -----

    INTERFACE
        ELEMENTAL MODULE FUNCTION nemo_cylinder_sizeof(cyl)
            USE psb_base_mod
            USE class_psblas
            IMPLICIT NONE
            CLASS(cylinder), INTENT(IN) :: cyl
            INTEGER(kind=nemo_int_long_)   :: nemo_cylinder_sizeof
        END FUNCTION nemo_cylinder_sizeof

        ! ----- Constructor -----

        MODULE SUBROUTINE alloc_cylinder(vertices,this_cylinder)
        !! Constructs cylinder by using steepest descents to fit a cylinder to the vertices' locations.
        !! We make our figure of merit f = 1 - error
            USE class_psblas
            USE class_vector
            USE tools_math ! contains interface to Levinburg-Marquardt algorithm
            IMPLICIT NONE
            TYPE(vertex),INTENT(IN)      :: vertices(:)
            TYPE(cylinder),POINTER       :: this_cylinder       !inout, so that we can check if set.
        END SUBROUTINE alloc_cylinder

        ! ----- Destructor -----

        MODULE SUBROUTINE free_cylinder(this_cylinder)
            USE class_psblas
            IMPLICIT NONE
            TYPE( cylinder ), POINTER  :: this_cylinder
        END SUBROUTINE free_cylinder


        ! ----- Getters -----

        MODULE FUNCTION get_cylinder_normal(this_cylinder,this_point)
        !! Returns the cylinder normal
            USE class_psblas
            USE class_vector
            IMPLICIT NONE
            TYPE(vector)              :: get_cylinder_normal
            CLASS(cylinder), INTENT(IN)   :: this_cylinder
            TYPE(vector),   INTENT(IN)   :: this_point
        END FUNCTION get_cylinder_normal


        ! Returns the point on a cylinder that is closest to the given point
        MODULE FUNCTION get_pt_cylinder(this_cylinder,point)
            USE class_vector
            IMPLICIT NONE
            TYPE(vector)              :: get_pt_cylinder   ! function result
            CLASS(cylinder),  INTENT(IN)  :: this_cylinder     ! the cylinder that we are on
            TYPE(vector) ,INTENT(IN)  :: point          ! the point off the cylinder
        END FUNCTION get_pt_cylinder


        MODULE FUNCTION get_cylinder_r2(this_cylinder)
            IMPLICIT NONE
            REAL(psb_dpk_)  :: get_cylinder_r2
            CLASS(cylinder),INTENT(IN)        :: this_cylinder
        END FUNCTION get_cylinder_r2

        ! ----- Setters -----

        MODULE SUBROUTINE translate_cylinder(this_cylinder,offset)
        !! moves the definition of a cylinder in the direction of the offset vector
            USE class_vector
            IMPLICIT NONE
            CLASS(cylinder), INTENT(INOUT) :: this_cylinder       ! the cylinder we are moving
            TYPE(vector),  INTENT(IN)     :: offset              ! translation vector
        END SUBROUTINE translate_cylinder

    END INTERFACE

END MODULE class_cylinder
