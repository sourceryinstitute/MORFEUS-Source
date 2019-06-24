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

SUBMODULE(class_plane) class_plane_procedures
    USE class_vertex
    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_plane_sizeof
        USE psb_base_mod
        USE class_psblas

        nemo_plane_sizeof = 5 * nemo_sizeof_dp + pl%normal%nemo_sizeof()

    END PROCEDURE nemo_plane_sizeof


    ! ----- Constructor -----

    ! Constructs plane by a least-squares fit

    MODULE PROCEDURE alloc_plane
        USE class_psblas
        USE class_vector
        USE tools_math

        ! Local variables
        INTEGER                    :: nverts, info
        REAL(psb_dpk_),TARGET,ALLOCATABLE   :: x(:),y(:),z(:)
        REAL(psb_dpk_)          :: MAT(3,3),rhs(3)    !Linear system for least-squares fit
        REAL(psb_dpk_)          :: A,B,C,D
        REAL(psb_dpk_)          :: coeff(0:2)         ! coefficiencts of fit
        REAL(psb_dpk_), POINTER :: dep(:),x1(:),x2(:) ! dep. & 2 indep. variables
        REAL(psb_dpk_)          :: xrange, yrange, zrange
        CHARACTER(len = 1)         :: smallest           ! which of x,y,z has the smallest range
        TYPE(vector)               :: unit_n             ! the unit normal of the plane
        REAL(psb_dpk_)          :: s_t,s_r            ! total and random deviations
        REAL(psb_dpk_)          :: ybar               ! average dep. variable
        REAL(psb_dpk_),ALLOCATABLE  :: dev(:)         ! deviation of each point

        IF ( ASSOCIATED(this_plane) ) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        ALLOCATE(this_plane,stat = info)
        IF ( info /= 0 ) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        nverts = SIZE(vertices)

        ALLOCATE( x(nverts), y(nverts), z(nverts), dev(nverts), stat=info)

        IF ( info /= 0 ) THEN
            WRITE(6,100)
            CALL abort_psblas
        ENDIF

        x = vertices%x_()
        y = vertices%y_()
        z = vertices%z_()

        ! Decide which will be our independent variable and which will be the dependent
        ! varible. Make the variable with the minimum range our independent variable.

        xrange = MAXVAL(x) - MINVAL(x)
        yrange = MAXVAL(y) - MINVAL(y)
        zrange = MAXVAL(z) - MINVAL(z)

        IF     ( ( xrange < yrange ) .AND. (xrange < zrange) ) THEN

            smallest = "x"

        ELSEIF ( ( yrange < xrange ) .AND. (yrange < zrange) ) THEN

            smallest = "y"

        ELSE

            smallest = "z"

        ENDIF


        SELECT CASE (smallest)

        CASE ("x")
            dep    => x
            x1     => y
            x2     => z

        CASE ("y")
            dep    => y
            x1     => z
            x2     => x

        CASE ("z")
            dep    => z
            x1     => x
            x2     => y

        END SELECT

        ! now do a least squares fit : indep = c0 + c1 * x1 + c2 * x2

        ! set up matrix
        MAT(1,1) = nverts

        MAT(2,1) = SUM(x1)
        MAT(1,2) = MAT(2,1)

        MAT(1,3) = SUM(x2)
        MAT(3,1) = MAT(1,3)

        MAT(2,2) = SUM(x1 * x1)

        MAT(2,3) = SUM(x1 * x2)
        MAT(3,2) = MAT(2,3)

        MAT(3,3) = SUM (x2 * x2)

        ! set up rhs
        rhs(1) = SUM( dep )
        rhs(2) = SUM (x1 * dep)
        rhs(3) = SUM (x2 * dep)

        ! solve linear system
        CALL factorize(MAT)
        CALL solve_sys(MAT,rhs)  ! returns solution in rhs

        coeff(0:2) = rhs(1:3)    ! change to notation of Chapra & Canalle

        !starting with the form dep = c0 + c1*x1 + c2*x2
        !now put in the form Ax + By +Cz = D

        d = coeff(0)

        SELECT CASE (smallest)

        CASE ("x")
            a  =  1.0d0
            b  = -coeff(1)  !y is 1st indep.
            c  = -coeff(2)  !z is 2nd indep.

        CASE ("y")
            b  =  1.0d0
            c  = -coeff(1)  !z is 1st indep.
            a  = -coeff(2)  !x is 2nd indep.

        CASE ("z")
            c  =  1.0d0
            a  = -coeff(1)  !x is 1st indep.
            b  = -coeff(2)  !y is 1st indep.

        END SELECT

        unit_n = (1.0d0/SQRT(a*a + b*b + c*c)) * vector_(a,b,c)

        ! scale the constant by the same factor
        d = d / SQRT(a*a + b*b + c*c)

        ! set plane's linear constants

        this_plane%a = a; this_plane%b =b; this_plane%c = c; this_plane%d =d

        ! set this_plane%normal
        this_plane%normal = unit_n

        ybar = SUM ( dep ) /nverts

        dev = dep - ybar

        s_t = SUM (dev * dev)

        s_r = SUM( (dep - coeff(0) - coeff(1)*x1(:) - coeff(2)*x2(:) ) ** 2 )

        IF ( s_t <= 1.0D-15 ) THEN ! the dependent variable is constant

            this_plane % r2 = 1.0d0

        ELSE

            this_plane % r2 = (s_t -s_r ) /s_t

        ENDIF

        DEALLOCATE(x,y,z,dev)
        NULLIFY(dep,x1,x2)

100     FORMAT(' ERROR! Failure to allocate memory in ALLOC_PLANE.')
200     FORMAT(' ERROR! Attempt to allocate previously created plane in ALLOC_PLANE.')

    END PROCEDURE alloc_plane

    ! ----- Destructor -----

    MODULE PROCEDURE free_plane
        USE class_psblas

        DEALLOCATE(this_plane)
    END PROCEDURE free_plane


    ! ----- Getters -----

    ! Returns the plane normal
    ! (which is trivial for a plane, since the normal is constant)
    MODULE PROCEDURE get_plane_normal
        USE class_psblas
        USE class_vector

        get_plane_normal = this_plane%normal

    END PROCEDURE get_plane_normal


    ! Returns an approximation for the goodness of fit, R2 value, from 0 to 1
    MODULE PROCEDURE get_plane_r2
        USE class_vector

        get_plane_r2 = this_plane%r2

    END PROCEDURE get_plane_r2

    ! Returns the point on a plane that is closest to the given point
    MODULE PROCEDURE get_pt_plane
        USE class_vector

        ! Local variables
        REAL(psb_dpk_)  ::    distance  ! the distance from the plane to the point
        TYPE(vector)       ::    offset    ! the vector pointing from the plane to the pt

        distance = ( this_plane%normal .dot. point ) - this_plane%d

        offset = distance * this_plane%normal

        get_pt_plane = point - offset

    END PROCEDURE get_pt_plane


    ! ----- Setters -----

    ! moves the definition of a plane in the direction of the offset vector
    MODULE PROCEDURE translate_plane
        USE class_vector

        this_plane%d = this_plane%d + ( this_plane%normal .dot. offset )

    END PROCEDURE translate_plane

END SUBMODULE class_plane_procedures
