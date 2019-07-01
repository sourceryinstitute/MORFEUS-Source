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

SUBMODULE(class_surface) class_surface_procedures

    USE class_plane
    USE class_cylinder

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_surface_sizeof
        USE psb_base_mod
        USE class_psblas
        INTEGER(kind=nemo_int_long_)   :: val

        val = 2 * nemo_sizeof_int
        IF (ASSOCIATED(surf%my_plane)) &
            & val = val + surf%my_plane%nemo_sizeof()
        IF (ASSOCIATED(surf%my_cylinder)) &
            & val = val + surf%my_cylinder%nemo_sizeof()

        nemo_surface_sizeof = val

    END PROCEDURE nemo_surface_sizeof

    ! ----- Constructor -----

    ! Constructs surface by instantiating a plane, cylinder, etc.
    ! and checking the goodness of fit.  We keep this surface pointing
    ! to the first good fit.  If none fit, then we return a null pointer.
    ! The idea is that we are auto-detecting the type of surface by trial-
    ! and-error, and that null means "unknown or irregular surface."

    MODULE PROCEDURE alloc_surface
        USE class_psblas
        USE class_connectivity
        USE class_cylinder
        USE class_plane
        USE class_vertex

        !    use class_sphere

        IMPLICIT NONE

        !
        ! Local variables
        REAL(psb_dpk_), PARAMETER :: acceptable = 0.9999d0 ! min acceptable goodness of fit

        TYPE(vertex),ALLOCATABLE     :: bndry_verts(:)      ! vertices for this boundary surf.
        INTEGER,POINTER              :: iv2b(:) => NULL()
        INTEGER                      :: nverts, info ,i

        ! Check to see if surface already exists
        IF ( this_surface%set ) THEN
            WRITE(6,100)
            CALL abort_psblas
        ENDIF

        ! Get list of vertices
        CALL v2b%get_ith_conn(iv2b,ib)

        nverts = SIZE(iv2b)

        ALLOCATE( bndry_verts(nverts), stat=info)

        IF ( info /= 0 ) THEN
            WRITE(6,200)
            CALL abort_psblas
        ENDIF

        ! save list of vertices on this boundary
        DO i = 1,nverts
            bndry_verts(i) = vertices(iv2b(i))
        ENDDO

        ! try to fit each type of surface to these points
        CALL alloc_plane(bndry_verts,this_surface%my_plane)
        ! CALL alloc_cylinder(bndry_verts,this_surface%my_cylinder)

        ! check the types of surface by goodness of fit. Note that the acceptable tolerance
        ! must be very high to avoid "false positives."
        IF ( this_surface%my_plane%get_plane_r2() >= acceptable ) THEN

            this_surface%itype = iplane_
            IF ( mypnum_() == 0) &
                & WRITE(6,'(a,i3,a,f10.6)')"  BC: ",ib,"  Autodetected plane with certainty    = ", &
                & this_surface%my_plane%get_plane_r2()

            ! CALL free_cylinder(this_surface%my_cylinder)

        ! ELSEIF ( get_cylinder_r2( this_surface%my_cylinder ) >= acceptable ) THEN

        !     this_surface%itype = icylinder_
        !     IF ( mypnum_() == 0) &
        !         & WRITE(6,'(a,i3,a,f10.6)')"  BC: ",ib,"  Autodetected cylinder with certainty = ", &
        !         & get_cylinder_r2( this_surface%my_cylinder )

        !     CALL free_plane(this_surface%my_plane)

        ELSE
            this_surface%itype = iunknown_
            IF ( mypnum_() == 0) &
                & WRITE(6,'(a,i3,a)')"  BC: ",ib,"  Unknown surface detected."
            ! this surface is unrecognized

            CALL free_plane(this_surface%my_plane)
            ! CALL free_cylinder(this_surface%my_cylinder)

            !       nullify(this_surface%sphere)
        ENDIF

        ! I am not yet sure if this is a good idea or not...
        ! but for now it is at least useful for debugging
        ! Clean up the surface positions of the vertices

        IF ( this_surface%itype /= iunknown_ ) THEN
            DO i = 1, SIZE(bndry_verts)
                CALL this_surface%reform_vertex(bndry_verts(i))
            ENDDO
        ENDIF

        ! altnerate exits would create memory leaks--make sure code goes through here

        NULLIFY(iv2b)
        DEALLOCATE(bndry_verts)

100     FORMAT(' ERROR! Surface already exists...cannot allocate memory.')
200     FORMAT(' ERROR! Failure to allocate memory in ALLOC_SURFACE.')

    END PROCEDURE alloc_surface

    ! ----- Destructor -----

    MODULE PROCEDURE free_surface
        USE class_psblas

        IMPLICIT NONE

        IF ( ASSOCIATED(this_surface%my_plane) )    &
            & CALL free_plane( this_surface%my_plane )

        IF ( ASSOCIATED(this_surface%my_cylinder) ) &
            & CALL free_cylinder( this_surface%my_cylinder )

    END PROCEDURE free_surface


    ! ----- Getters -----

    ! Returns a named integer constant indicating type:
    ! iunknown, iplane, icylinder, isphere

    MODULE PROCEDURE get_surface_type

        IMPLICIT NONE

        get_surface_type = this_surface%itype

    END PROCEDURE get_surface_type

    ! Returns the surface normal at an appropriately close point
    ! We assume that the point is actually on the surface
    MODULE PROCEDURE get_surface_normal
        USE class_psblas
        USE class_vector

        ! we want to know the normal

        SELECT CASE ( this_surface%itype )

        CASE (iplane_)
            get_surface_normal = this_surface%my_plane%get_plane_normal()

        CASE (icylinder_)
            get_surface_normal = this_surface%my_cylinder%get_cylinder_normal(this_point)

        CASE default
            WRITE(6,100)
            CALL abort_psblas

        END SELECT

100     FORMAT('Normal not known for surface in GET_SURFACE_NORMAL')

    END PROCEDURE get_surface_normal


    MODULE PROCEDURE get_surface_r2
    !! Returns the goodness of fit, R2 value, from 0 to 1
        USE class_psblas
        USE class_vector

        IMPLICIT NONE

        SELECT CASE ( this_surface%itype)

        CASE (iplane_)
            get_surface_r2 = this_surface%my_plane%get_plane_r2()

        CASE (icylinder_)
            get_surface_r2 = this_surface%my_cylinder%get_cylinder_r2()

        CASE default
            WRITE(6,100)
            CALL abort_psblas

        END SELECT

100     FORMAT('R2 not known for surface in GET_SURFACE_R2')

    END PROCEDURE get_surface_r2

    MODULE PROCEDURE translate_surface
    !! Move a surface by translating it in 3D space
        USE class_psblas
        USE class_vector
        USE class_plane

        IMPLICIT NONE

        SELECT CASE ( this_surface%itype)

        CASE (iplane_)
            CALL this_surface%my_plane%translate_plane(offset)

        CASE (icylinder_)
            CALL this_surface%my_cylinder%translate_cylinder(offset)

        CASE default
            ! do nothing...there are no data describing an irregular surface
        END SELECT

    END PROCEDURE translate_surface

    MODULE PROCEDURE get_closest_point
    !! Returns the point on a surface that is closest to the given point
        USE class_psblas
        USE class_vector

        IMPLICIT NONE

        SELECT CASE ( this_surface%itype )

        CASE (iplane_)
            get_closest_point = this_surface%my_plane%get_pt_plane(point)

        CASE (icylinder_)
            get_closest_point = this_surface%my_cylinder%get_pt_cylinder(point)

        CASE default
            WRITE(6,100)
            CALL abort_psblas

        END SELECT

100     FORMAT('Unrecognized type of surface in GET_CLOSEST_POINT')

    END PROCEDURE get_closest_point


    ! ----- Setters -----

    MODULE PROCEDURE reform_vertex
    !! moves the given vertex onto the closest point on the surface

        USE class_psblas
        USE class_vector
        USE class_vertex

        IMPLICIT NONE
        !
        ! Local variables
        TYPE(vector) :: old_pos,new_pos

        old_pos = vtx%position_()

        new_pos = this_surface%get_closest_point(old_pos)

        vtx = new_pos

    END PROCEDURE reform_vertex


    !Below is commented out because this procedure is not used

    !MODULE PROCEDURE get_surface_set
    !! returns true if the surface has been set up already

    !    get_surface_set = this_surface%set

    !END PROCEDURE get_surface_set

END SUBMODULE class_surface_procedures
