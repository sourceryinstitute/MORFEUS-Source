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

SUBMODULE(class_cylinder) class_cylinder_procedures
    USE class_vector
    USE class_vertex
    USE class_psblas, ONLY : psb_dpk_

    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_cylinder_sizeof
        USE psb_base_mod
        USE class_psblas

        nemo_cylinder_sizeof = 2 * nemo_sizeof_dp &
            & + cyl%center%nemo_sizeof() + cyl%axis%nemo_sizeof()

    END PROCEDURE nemo_cylinder_sizeof
    ! ----- Constructor -----

    ! Constructs cylinder by using steepest descents to fit a cylinder to the vertices' locations.
    ! We make our figure of merit f = 1 - error

    MODULE PROCEDURE alloc_cylinder
        USE class_psblas
        USE class_vector
        USE tools_math ! contains interface to Levinburg-Marquardt algorithm

        IMPLICIT NONE

        ! Local variables
        INTEGER                    :: nverts, info
        REAL(psb_dpk_),TARGET,ALLOCATABLE   :: x(:),y(:),z(:)
        REAL(psb_dpk_)          :: xrange, yrange, zrange

        REAL(psb_dpk_),TARGET   :: unknowns(6)            ! list of unknowns to be optimized
        REAL(psb_dpk_),POINTER  :: xc    => NULL()
        REAL(psb_dpk_),POINTER  :: yc    => NULL()
        REAL(psb_dpk_),POINTER  :: zc    => NULL()
        REAL(psb_dpk_),POINTER  :: alpha => NULL()
        REAL(psb_dpk_),POINTER  :: beta  => NULL()
        REAL(psb_dpk_),POINTER  :: radius=> NULL()

        REAL(psb_dpk_)          :: ax,ay,az
        TYPE(vector)               :: axis

        ! For calling Levingburg-Marquart algorithm DNLS1E
        REAL(psb_dpk_),ALLOCATABLE  :: err(:)
        INTEGER :: IOPT,M,N,NPRINT,LWA
        INTEGER, ALLOCATABLE :: IW(:)
        REAL(psb_dpk_) :: TOL
        REAL(psb_dpk_), ALLOCATABLE :: WA(:)
        INTEGER :: tries

        iopt = 1                 ! Numerically approximate the Jacobian
        m = SIZE(vertices)       ! Number of errors to be minimized
        n = 6                    ! Number of free parameters
        nprint = 0               ! Frequency of output during iterations (0 means never print)
        lwa = N*(M+5)+M          ! Size of real work array
        tol = SQRT(EPSILON(tol)) ! Error tolerance, as recommended by documentation for SLATEC

        ALLOCATE(iw(n),wa(lwa),stat = info)
        IF ( info /= 0 ) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        !establish correspondence of unknowns to descriptive names
        xc => unknowns(1); yc => unknowns(2) ; zc => unknowns(3)

        alpha => unknowns(4)
        beta  => unknowns(5)
        radius=> unknowns(6)

        IF ( ASSOCIATED(this_cylinder) ) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        ALLOCATE(this_cylinder,stat = info)
        IF ( info /= 0 ) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        nverts = SIZE(vertices)

        ALLOCATE(my_vertices(nverts), x(nverts), y(nverts), z(nverts),err(nverts), stat=info)

        IF ( info /= 0 ) THEN
            WRITE(6,100)
            CALL abort_psblas
        ENDIF

        my_vertices = vertices ! my_vertices is a copy with module scope for use in FCN

        ! guess no rotation
        alpha = 0.0d0 ; beta = 0.0d0

        ! used for initial guesses
        x = vertices%x_()
        y = vertices%y_()
        z = vertices%z_()

        xrange = MAXVAL(x) - MINVAL(x)
        yrange = MAXVAL(y) - MINVAL(y)
        zrange = MAXVAL(z) - MINVAL(z)

        ! set error flag to indicate non-fatal, quiet response to errors
        ! CALL xsetf(0)

        tries = 0

        ! Try a few different initial guesses for the cylinder orientation...if the guess is exactly
        ! orthogonal to the true orientation, recognition can fail
        try: DO

            tries = tries + 1

            ! Set initial guesses for cylinder

            ! Use the mean of the x,y,z values for the cylinder center
            xc = SUM(x)/SIZE(x)
            yc = SUM(y)/SIZE(y)
            zc = SUM(z)/SIZE(z)

            ! guess radius based on ranges
            radius = (xrange + yrange + zrange)/ 6.0d0

            ! call main least-squares solver
            ! CALL dnls1e(fcn,iopt,m,n,unknowns,err,tol,nprint, &
            !     &            info,iw,wa,lwa)

            IF (info == 0) THEN  ! for other problems, fail quietly.  But info == 0 should be fatal.
                WRITE(6,300)
                CALL abort_psblas
            ENDIF

            this_cylinder%center = vector_(xc,yc,zc)

            ! calculate axis from angles of rotation

            ax =  SIN(beta)
            ay = -SIN(alpha)*COS(beta)
            az =  COS(alpha)*COS(beta)

            axis = vector_(ax,ay,az)

            this_cylinder%axis   = vector_(ax,ay,az)  ! a unit vector axis

            this_cylinder%radius = radius

            ! set cylinder r2.  If this is approximately 1.0 then we found a cylinder
            ! otherwise, it is likely not a cylinder
            this_cylinder%r2 = try_cylinder_r2(this_cylinder%center,axis,radius,vertices)

            IF ( (this_cylinder%r2 > 0.999 ) .OR. (tries >= 3) ) EXIT try

            ! perturb ortation guess

            IF ( tries == 1 ) THEN
                alpha = 0.785d0   ! approx pi/4
                beta  = 0.0d0
            ENDIF


            IF ( tries == 2 ) THEN
                alpha = 0.0d0
                beta  = 0.785d0  ! approx pi/4
            ENDIF

        END DO try

        DEALLOCATE(x,y,z,my_vertices,err)
        DEALLOCATE(iw,wa)

100     FORMAT(' ERROR! Failure to allocate memory in ALLOC_CYLINDER.')
200     FORMAT(' ERROR! Attempt to allocate previously created cylinder in ALLOC_CYLINDER.')
300     FORMAT(' ERROR! Failure to provide correct inputs to optimizer in ALLOC_CYLINDER.')

    END PROCEDURE alloc_cylinder

    ! ----- Destructor -----

    MODULE PROCEDURE free_cylinder
        USE class_psblas

        IMPLICIT NONE

        DEALLOCATE(this_cylinder)

    END PROCEDURE free_cylinder


    ! ----- Getters -----


    ! Returns the cylinder normal

    MODULE PROCEDURE get_cylinder_normal
        USE class_psblas
        USE class_vector

        ! Local variables

        TYPE(vector)       ::    xrad      ! points from the axis to the point
        TYPE(vector)       ::    xrel      ! points from the cylinder center to the point

        ! calculate position relative to the cylinder's center
        xrel = this_point - this_cylinder%center

        ! get radial vector
        xrad = xrel - ( xrel .dot. this_cylinder%axis) * this_cylinder%axis

        ! return unit vector

        get_cylinder_normal = xrad%unit()

    END PROCEDURE get_cylinder_normal


    ! Returns the point on a cylinder that is closest to the given point
    MODULE PROCEDURE get_pt_cylinder
        USE class_vector

        IMPLICIT NONE

        ! Local variables
        TYPE(vector)       ::    offset    ! points from the closest pt. on cyl. surf. to point
        TYPE(vector)       ::    xrad      ! points from the axis to the point
        TYPE(vector)       ::    xrel      ! points from the cylinder center to the point

        ! calculate position relative to the cylinder's center
        xrel = point - this_cylinder%center

        ! get radial vector from the closest spot on the axis to the point
        ! by subracting off the axial component of the relative position

        xrad = xrel - ( xrel .dot. this_cylinder%axis) * this_cylinder%axis

        ! subtract a vector from the closest spot on the axis to the cyl. surface
        offset = xrad - ( this_cylinder%radius / xrad%mag() )  * xrad

        ! offset is now the difference between the point and the cylinder surface
        get_pt_cylinder = point - offset

    END PROCEDURE get_pt_cylinder


    MODULE PROCEDURE get_cylinder_r2

        IMPLICIT NONE

        get_cylinder_r2 = this_cylinder%r2

    END PROCEDURE get_cylinder_r2

    ! ----- Setters -----

    ! moves the definition of a cylinder in the direction of the offset vector
    MODULE PROCEDURE translate_cylinder
        USE class_vector

        IMPLICIT NONE

        this_cylinder%center = this_cylinder%center + offset

    END PROCEDURE translate_cylinder

    ! ------ Private Functions ------


    ! returns L2 norm of error for a vector of unknowns
    ! the assumed order is center vector, axis vector, and radius
!!$  function f(unknowns,vertices)
!!$
!!$    use class_vector
!!$
!!$    real(kind(1.0d0))                 :: f
!!$    real(kind(1.0d0)),intent(in)      :: unknowns(6)
!!$    type(vertex)     ,intent(in)      :: vertices(:)
!!$
!!$    ! Local variables
!!$    type(vector)         :: center
!!$    type(vector)         :: axis
!!$    real(kind(1.0d0))    :: radius,alpha,beta
!!$    real(kind(1.0d0))    :: ax,ay,az
!!$
!!$    center = vector_(unknowns(1),unknowns(2),unknowns(3))
!!$    alpha  = unknowns(4)
!!$    beta   = unknowns(5)
!!$    radius = unknowns(6)
!!$
!!$    ax =  sin(beta)
!!$    ay = -sin(alpha)*cos(beta)
!!$    az =  cos(alpha)*cos(beta)
!!$
!!$    axis = vector_(ax,ay,az)
!!$
!!$    f =  calc_error(center,axis,radius,vertices)
!!$
!!$  end function f

    ! Returns an approximation for the goodness of fit using specified axis, radius, center
    ! Not exactly an R2 value, since it is not limited from 0 to 1
    FUNCTION try_cylinder_r2(center,axis,radius,vertices)
        USE class_vector

        IMPLICIT NONE

        REAL(psb_dpk_)  :: try_cylinder_r2

        TYPE(vector),INTENT(IN)        :: center
        TYPE(vector),INTENT(IN)        :: axis
        REAL(psb_dpk_),INTENT(IN)   :: radius
        TYPE(vertex),INTENT(IN)        :: vertices(:)

        ! Local variables
        REAL(psb_dpk_)              :: s_t                ! total deviations
        INTEGER                        :: nverts             ! number of vertices
        INTEGER                        :: i
        REAL(psb_dpk_)              :: err                ! cumulative error
        TYPE(vector)                   :: centroid, xrad

        ! Calculate R2 for the cylinder

        nverts = SIZE(vertices)

        err = calc_error(center,axis,radius,vertices)
        centroid = vector_(0.0d0, 0.0d0, 0.0d0)

        DO i = 1,nverts

            centroid = centroid + vertices(i)%position_()

        ENDDO

        centroid = ( 1.0d0/float(nverts) )* centroid

        s_t = 0.0

        ! calculate cumulative error
        DO i = 1,nverts

            ! calculate position relative to the centroid
            xrad = centroid - vertices(i)%position_()

            s_t = s_t + xrad%mag()**2

        ENDDO

        s_t = s_t / float(nverts)

        try_cylinder_r2 = ( s_t - err ) / s_t

    END FUNCTION try_cylinder_r2

    ! Returns an approximation for the goodness of fit: the L2 norm of error
    ! A perfect fit is zero.

    FUNCTION calc_error(center,axis,radius,vertices)
        USE class_vector

        IMPLICIT NONE

        REAL(psb_dpk_)  :: calc_error

        TYPE(vector),INTENT(IN)        :: center
        TYPE(vector),INTENT(IN)        :: axis
        REAL(psb_dpk_),INTENT(IN)   :: radius
        TYPE(vertex),INTENT(IN)        :: vertices(:)

        ! Local variables
        INTEGER                        :: nverts             ! number of vertices
        INTEGER                        :: i
        REAL(psb_dpk_)              :: err                ! cumulative error
        TYPE(vector)                   :: x,xrad,xrel


        nverts = SIZE(vertices)

        err = 0

        ! calculate cumulative error
        DO i = 1,nverts

            x = vertices(i)%position_()

            ! calculate position relative to the cylinder's center
            xrel = x - center

            ! get radial vector
            xrad = xrel - ( xrel .dot. axis) * axis

            ! because our ideal cylinder is infinitely long, the error is only
            ! the vector difference in the radial direction

            err = err + ( xrad%mag() - radius)**2

        ENDDO

        calc_error = err/float(nverts)

    END FUNCTION calc_error

    SUBROUTINE FCN(iflag,nverts,nunknowns,unknowns,err,FJAC,LDFJAC)
    !! This is the function of merit to be optimized
        USE class_psblas
        USE class_vector

        IMPLICIT NONE

        INTEGER,INTENT(IN)            :: iflag
        INTEGER,INTENT(IN)            :: nverts             ! number of vertices (M)
        INTEGER,INTENT(IN)            :: nunknowns          ! number of variables (N)
        REAL(psb_dpk_),INTENT(IN)  :: unknowns(6)
        REAL(psb_dpk_),INTENT(OUT) :: err(nverts)
        REAL(psb_dpk_),INTENT(IN)  :: fjac(1,1)     ! ignored
        INTEGER,INTENT(IN)            :: ldfjac        ! ignored

        TYPE(vector)        :: center
        TYPE(vector)        :: axis
        REAL(psb_dpk_)   :: radius

        ! Local variables
        INTEGER                        :: i
        TYPE(vector)                   :: x,xrad,xrel
        REAL(psb_dpk_)              :: ax,ay,az           ! components of the axis vector
        REAL(psb_dpk_)              :: alpha, beta        !angles of cylinder rotation

        INTEGER                        :: idummy             ! this variable exists to make the
        ! compiler happy, but serves no other use

        idummy = SIZE(fjac) ! since we don't control the interface to this subroutines, we must
        idummy = ldfjac     ! include these two unused variables.  These two lines make them
        ! appear used and supressed "unused variable warnings.

        IF (nunknowns /= SIZE(unknowns) ) THEN
            WRITE(6,*)"Error in argument to FCN"
            CALL abort_psblas
        ENDIF

        IF (iflag == 0) THEN
            WRITE(6,'(a)')"========= Center        X           Y           Z          ALPHA       &
                &BETA        RADIUS"
            WRITE(6,'(a,6(2x,e10.4))')"The unknowns are: ",unknowns
            RETURN
        ENDIF

        center = vector_(unknowns(1),unknowns(2),unknowns(3))
        alpha  = unknowns(4)
        beta   = unknowns(5)
        radius = unknowns(6)

        ax =  SIN(beta)
        ay = -SIN(alpha)*COS(beta)
        az =  COS(alpha)*COS(beta)

        axis = vector_(ax,ay,az)

        !    nverts = size(my_vertices)

        err = 0

        ! calculate cumulative error
        DO i = 1,nverts

            x = my_vertices(i)%position_()

            ! calculate position relative to the cylinder's center
            xrel = x - center

            ! get radial vector
            xrad = xrel - ( xrel .dot. axis) * axis

            ! because our ideal cylinder is infinitely long, the error is only
            ! the vector difference in the radial direction

            err(i) =  ( xrad%mag() - radius)**2

        ENDDO

    END SUBROUTINE FCN

END SUBMODULE class_cylinder_procedures
