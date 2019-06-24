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
! $Id: geom_diff.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    Evaluation of diffusion distance, intersection point between face
!    and master/slave segment, interpolation factor (fluid faces only).
!
MODULE PROCEDURE geom_diff
    USE class_psblas
    USE class_connectivity
    USE class_face
    USE class_vector

    IMPLICIT NONE
    !
    LOGICAL, PARAMETER :: debug = .FALSE.
    INTEGER :: icontxt, mypnum
    INTEGER :: i, ib, IF, im, info, is, m, n, nbc, nfaces
    INTEGER, ALLOCATABLE :: imaster(:), islave(:)
    INTEGER, POINTER :: if2b(:) => NULL()
    REAL(psb_dpk_) :: area, dfs, r, s, t
    TYPE(vector) :: norm, v
    TYPE(vector), ALLOCATABLE :: isect(:)

    icontxt = icontxt_()
    mypnum = mypnum_()

    ! Preliminary check
    IF(    .NOT.ALLOCATED(faces)     .OR. &
        & .NOT.ALLOCATED(face_cntr) .OR. &
        & .NOT.ALLOCATED(af)        .OR. &
        & .NOT.ALLOCATED(cell_cntr)) THEN
        WRITE(*,100)
        CALL abort_psblas
    END IF
    ! Not necessary since they are INTENT(OUT)
!!$  if(    allocated(df)   .or. &
!!$       & allocated(dist) .or. &
!!$       & allocated(int_fact)) then
!!$     call abort_psblas
!!$  end if

    nfaces = SIZE(faces)

    ALLOCATE(imaster(nfaces),islave(nfaces),isect(nfaces), &
        & df(nfaces),dist(nfaces),int_fact(nfaces),stat=info)
    IF(info /= 0) THEN
        WRITE(*,200)
        CALL abort_psblas
    END IF

    ! Extract master and slave indices from faces
    imaster(:) = master_(faces)
    islave(:)  = slave_(faces)

    ! Initialization.
    df(:) = vector_(0.d0,0.d0,0.d0)
    dist(:) = 0.0d0
    isect(:) = vector_(0.d0,0.d0,0.d0)
    int_fact(:) = 0.d0


    ! ----- Internal Fluid Faces (flag = 0) -----

    CALL get_ith_conn(if2b,f2b,0)
    n = SIZE(if2b)

    DO i = 1, n
        IF = if2b(i)
        im = imaster(IF)
        is = islave(IF)

        ! Definition of distance vector between the centers of master and slave cells
        df(IF) = cell_cntr(is) - cell_cntr(im)
        dist(IF) = mag(df(IF))

        ! Definition of intersection point between face and master-slave segment
        v = face_cntr(IF) - cell_cntr(im)
        r = v .dot. af(IF)
        s = df(IF) .dot. af(IF)
        t = r / s

        isect(IF) = cell_cntr(im) + t * df(IF)

        ! Definition of interpolation factor
        v = cell_cntr(is) - isect(IF)
        dfs = mag(v)
        int_fact(IF) = dfs / dist(IF)
    END DO


    ! ----- Boundary Faces (flag > 0) -----

    ! DF is a vector starting from the master center, normal to the face and of
    ! magnitude equal to the distance master-face.

    nbc = nel_(f2b) - 2

    DO ib = 1, nbc
        CALL get_ith_conn(if2b,f2b,ib)
        n = SIZE(if2b)

        DO i = 1, n
            IF = if2b(i)

            ! Intersection point between the face and the normal from master center
            im = imaster(IF)
            area = mag(af(IF))
            norm = (1.d0 / area) * af(IF)
            v = face_cntr(IF) - cell_cntr(im)

            ! Definition of diffusion distance vector
            dist(IF) = v .dot. norm
            df(IF) = dist(IF) * norm

            ! Definition of interpolation factor (redundant, not needed for bc faces)
            int_fact(IF) = 0.5d0
        END DO
    END DO


    ! Check for the boundeness of INT_FACT
    m = COUNT(int_fact < 0.d0)
    n = COUNT(int_fact > 1.d0)
    !
    CALL psb_sum(icontxt,m)
    IF( mypnum == 0 .AND. m /= 0 ) THEN
        WRITE(*,300) m
        CALL abort_psblas
    END IF
    !
    CALL psb_sum(icontxt,n)
    IF( mypnum == 0 .AND. n /= 0 ) THEN
        WRITE(*,400) n
        CALL abort_psblas
    END IF

    !#############################################################################
    IF(debug) THEN
        DO IF = 1, nfaces
            WRITE(*,600) IF, flag_(faces(IF)), master_(faces(IF)), slave_(faces(IF)), int_fact(IF)
        END DO
        WRITE(*,*)
600     FORMAT('Face: ',i6,' | Flag:',i3,' | Master:',i5,' | Slave:',i5,' | Int_fct:',3(1x,d13.6))
    END IF
    !#############################################################################

    NULLIFY(if2b)
    DEALLOCATE(imaster,islave,isect)

100 FORMAT(' ERROR! Input array not allocated in GEOM_DIFF')
200 FORMAT(' ERROR! Memory allocation failure in GEOM_DIFF')
300 FORMAT(' WARNING!',i8,' faces with INT_FACT less than 0.d0')
400 FORMAT(' WARNING!',i8,' faces with INT_FACT more than 1.d0')


END PROCEDURE geom_diff
