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
! $Id: class_plane.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    A class that handles the geometric concept of a plane.
!
! Provides:
!    PLANE              class.
!    ALLOC_PLANE        constructor for PLANE class.
!    FREE_PLANE         destructor for PLANE class.
!    GET_PLANE_NORMAL   returns the normal of the plane at a certain point
!    GET_PLANE_R2       returns how good a fit the plane is
!    GET_PT_PLANE       returns the nearest point that lies on the plane
!    TRANSLATE_PLANE    moves a plane in 3D space, translation only

! To be done :: a move plane function

MODULE class_plane
    USE class_vector
    USE class_psblas, ONLY : psb_dpk_

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: plane ! Class
    PUBLIC :: alloc_plane, free_plane          ! Constructor/Destructor
    PUBLIC :: get_plane_normal, get_plane_r2   ! Getters
    PUBLIC :: get_pt_plane                     ! Getters, cont.
    PUBLIC :: translate_plane                  ! Setters

    TYPE plane
        PRIVATE

        TYPE(vector)      :: normal       ! the surface unit normal
        REAL(psb_dpk_) :: a,b,c,d      ! linear form of a plane:
        !  a * x + b * y + c * z =d
        REAL(psb_dpk_) :: r2           ! The correlation parameter for the fit
    CONTAINS
        PROCEDURE, PRIVATE :: nemo_plane_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_plane_sizeof
    END TYPE plane


  ! ----- Generic Interfaces -----

    INTERFACE

    ELEMENTAL MODULE FUNCTION nemo_plane_sizeof(pl)
        USE class_psblas, ONLY : nemo_int_long_
        IMPLICIT NONE
        CLASS(plane), INTENT(IN) :: pl
        INTEGER(kind=nemo_int_long_)   :: nemo_plane_sizeof
    END FUNCTION nemo_plane_sizeof

    ! ----- Constructor -----

    MODULE SUBROUTINE alloc_plane(vertices,this_plane)
      !! Constructs plane by a least-squares fit
        USE class_psblas
        USE class_vector
        USE class_vertex
        USE tools_math
        IMPLICIT NONE
        TYPE(vertex),INTENT(IN)      :: vertices(:)
        TYPE(plane),POINTER          :: this_plane       !inout, so that we can check SET.
    END SUBROUTINE alloc_plane

    ! ----- Destructor -----

    MODULE SUBROUTINE free_plane(this_plane)
        USE class_psblas
        IMPLICIT NONE
        TYPE( plane ), POINTER  :: this_plane
    END SUBROUTINE free_plane

    ! ----- Getters -----

    MODULE FUNCTION get_plane_normal(this_plane)
      !! Returns the plane normal
      !! (which is trivial for a plane, since the normal is constant)
        USE class_psblas
        USE class_vector
        IMPLICIT NONE
        TYPE(vector)              :: get_plane_normal
        TYPE(plane), INTENT(IN)   :: this_plane
    END FUNCTION get_plane_normal


    MODULE FUNCTION get_plane_r2(this_plane)
      !! Returns an approximation for the goodness of fit, R2 value, from 0 to 1
        USE class_vector
        IMPLICIT NONE
        REAL(psb_dpk_)  :: get_plane_r2
        TYPE(plane)        :: this_plane
    END FUNCTION get_plane_r2

    MODULE FUNCTION get_pt_plane(this_plane,point)
      !! Returns the point on a plane that is closest to the given point
        USE class_vector
        IMPLICIT NONE
        TYPE(vector)              :: get_pt_plane   ! function result
        TYPE(plane),  INTENT(IN)  :: this_plane     ! the plane that we are on
        TYPE(vector) ,INTENT(IN)  :: point          ! the point off the plane
    END FUNCTION get_pt_plane


    ! ----- Setters -----

    MODULE SUBROUTINE translate_plane(this_plane,offset)
      !! moves the definition of a plane in the direction of the offset vector
        USE class_vector
        IMPLICIT NONE
        TYPE(plane), INTENT(INOUT) :: this_plane       ! the plane we are moving
        TYPE(vector),  INTENT(IN)    :: offset           ! translation vector
    END SUBROUTINE translate_plane

  END INTERFACE

END MODULE class_plane
