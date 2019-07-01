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
! $Id: class_surface.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    A base class that contains various kinds of surfaces on which a boundary may or
!    may not lie.  Provides generic access for getting a surfaces local normal, the
!    goodness of fit between the boundary points and the surface, the surface's
!    constructor, destructor.  Also provides services like
!    the closest point on the surface.  Many functions assume that the shape of
!    the surface is a recognized type.  Check if ( type_ == unknown ) before using
!    these functions.
!
! Provides:
!    SURFACE              class.
!    ALLOC_SURFACE        constructor for SURFACE class.
!    FREE_SURFACE         destructor for SURFACE class.
!    TYPE_                indentifies the type of the surface
!    GET_SURFACE_NORMAL   returns the normal of the surface at a certain point
!    GET_SURFACE_R2       returns how good a fit the surface is
!    GET_CLOSEST_POINT    returns the nearest point that lies on the surface
!    REFORM_VERTEX        moves the vertex to the closest point on the surface
!    TRANSLATE_SURFACE    translates the definition of the surface in 3D space
!    SET_                 true if the surface has been set up

MODULE class_surface

    USE class_plane
    USE class_cylinder
    USE class_connectivity

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: surface ! Class
    PUBLIC :: alloc_surface                        ! Constructor
    PUBLIC :: iunknown_, iplane_, icylinder_, isphere_ ! Named constants

    ! To be done:  a routine that moves the surface

    ! Surface is a base class that contains plane, cylinder, and surface

    TYPE surface
        PRIVATE

        ! Surface holds one of these three kinds of surface.
        ! PLANE is-a SURFACE, A_CYLINDER is-a SURFACE, etc.

        TYPE( plane ),    POINTER  :: my_plane =>null()
        TYPE( cylinder ), POINTER  :: my_cylinder => NULL()
        !     type( sphere   ), pointer :: my_sphere => null()    ! not yet written
        LOGICAL :: set = .FALSE.      ! true if this surface has already been created
        INTEGER :: itype = 0          ! indicates the type of surface
    CONTAINS
        PROCEDURE :: free_surface                                  ! Destructor
        PROCEDURE, PRIVATE :: get_surface_type                     ! Returns Int. ID #
        GENERIC, PUBLIC :: type_ => get_surface_type
        PROCEDURE, PRIVATE :: get_surface_normal, get_surface_r2   ! Getters
        GENERIC, PUBLIC :: normal_ => get_surface_normal
        GENERIC, PUBLIC :: r2_ => get_surface_r2
        PROCEDURE :: get_closest_point                    ! Getters, cont.
        PROCEDURE :: reform_vertex, translate_surface     ! Setters.
        PROCEDURE, PRIVATE :: nemo_surface_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_surface_sizeof
    END TYPE surface

    ! ----- Named Constants -----

    INTEGER, PARAMETER :: iunknown_ = 0,  iplane_ = 1     ! unknown, plane
    INTEGER, PARAMETER :: icylinder_ = 2, isphere_ = 3    ! cylinder, sphere


    ! ----- Generic Interfaces -----

    INTERFACE

        ELEMENTAL MODULE FUNCTION nemo_surface_sizeof(surf)
            USE class_psblas, ONLY : nemo_int_long_
            IMPLICIT NONE
            CLASS(surface), INTENT(IN) :: surf
            INTEGER(kind=nemo_int_long_)   :: nemo_surface_sizeof
        END FUNCTION nemo_surface_sizeof

        ! ----- Constructor -----

        ! Constructs surface by instantiating a plane, cylinder, etc.
        ! and checking the goodness of fit.  We keep this surface pointing
        ! to the first good fit.  If none fit, then we return a null pointer.
        ! The idea is that we are auto-detecting the type of surface by trial-
        ! and-error, and that null means "unknown or irregular surface."

        MODULE SUBROUTINE alloc_surface(v2b,ib,vertices,this_surface)
            USE class_psblas
            !USE class_connectivity
            USE class_cylinder
            USE class_plane
            USE class_vertex
            !    use class_sphere
            IMPLICIT NONE
            TYPE(connectivity), INTENT(IN) :: v2b
            INTEGER, INTENT(IN)            :: ib
            TYPE(vertex),INTENT(IN)        :: vertices(:)    ! All mesh vertices
            TYPE(surface),INTENT(INOUT)    :: this_surface   ! inout, so that we can check SET.
        END SUBROUTINE alloc_surface

        ! ----- Destructor -----

        MODULE SUBROUTINE free_surface(this_surface)
            USE class_psblas
            IMPLICIT NONE
            CLASS(surface), INTENT(INOUT) :: this_surface
        END SUBROUTINE free_surface


        ! ----- Getters -----

        ! Returns a named integer constant indicating type:
        ! iunknown, iplane, icylinder, isphere

        ! Getters

        MODULE FUNCTION get_surface_type(this_surface)
            IMPLICIT NONE
            INTEGER                   :: get_surface_type
            CLASS(surface), INTENT(IN) :: this_surface
        END FUNCTION get_surface_type

        MODULE FUNCTION get_surface_normal(this_surface, this_point)
        !! Returns the surface normal at an appropriately close point
        !! We assume that the point is actually on the surface
            USE class_psblas
            USE class_vector
            IMPLICIT NONE
            TYPE(vector)              :: get_surface_normal
            CLASS(surface), INTENT(IN) :: this_surface
            TYPE(vector),  INTENT(IN) :: this_point         ! the point on the surface where
        END FUNCTION get_surface_normal

        MODULE FUNCTION get_surface_r2(this_surface)
        !! Returns the goodness of fit, R2 value, from 0 to 1
            USE class_psblas, ONLY : psb_dpk_
            USE class_vector
            IMPLICIT NONE
            REAL(psb_dpk_)         :: get_surface_r2
            CLASS(surface), INTENT(IN) :: this_surface
        END FUNCTION get_surface_r2

        MODULE SUBROUTINE translate_surface(this_surface,offset)
        !! Move a surface by translating it in 3D space
            USE class_psblas
            USE class_vector
            USE class_plane
            IMPLICIT NONE
            CLASS(surface), INTENT(INOUT) :: this_surface
            TYPE(vector),  INTENT(IN)    :: offset
        END SUBROUTINE translate_surface

        MODULE FUNCTION get_closest_point(this_surface,point)
        !! Returns the point on a surface that is closest to the given point
            USE class_psblas
            USE class_vector
            IMPLICIT NONE
            TYPE(vector)              :: get_closest_point
            CLASS(surface), INTENT(IN) :: this_surface
            TYPE(vector) , INTENT(IN) :: point
        END FUNCTION get_closest_point


        ! ----- Setters -----

        MODULE SUBROUTINE reform_vertex(this_surface, vtx)
        !! moves the given vertex onto the closest point on the surface
            USE class_psblas
            USE class_vector
            USE class_vertex
            IMPLICIT NONE
            CLASS(surface), INTENT(IN)    :: this_surface
            TYPE(vertex) , INTENT(INOUT) :: vtx
        END SUBROUTINE reform_vertex

    END INTERFACE

    !Below is commented out because this procedure is not used

    !INTERFACE set_
    !  MODULE FUNCTION get_surface_set(this_surface)
    !      ! returns true if the surface has been set up already
    !      LOGICAL :: get_surface_set
    !      TYPE(surface), INTENT(IN)    :: this_surface
    !  END FUNCTION get_surface_set
    !END INTERFACE set_

END MODULE class_surface
