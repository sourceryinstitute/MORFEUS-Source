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
! $Id: geom_face.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    Calculates face related geometry, such as face-normal area vectors & area centroid
!
MODULE PROCEDURE geom_face
    USE class_psblas
    USE class_connectivity
    USE class_vector
    USE class_vertex

    IMPLICIT NONE
    !
    LOGICAL, PARAMETER :: debug = .FALSE.
    INTEGER :: i, i1, i2, IF, info, nfaces
    INTEGER :: iv_1, iv_i, iv_im1
    INTEGER, ALLOCATABLE :: nvf(:)
    INTEGER, POINTER :: iv2f(:) => NULL()
    REAL(psb_dpk_), PARAMETER :: t3 = 1.d0 / 3.d0
    REAL(psb_dpk_) :: r, w
    TYPE(vector) :: p, v1, v2

    ! Evaluation of face-normal vectors AF(:) and face-centers coordinates
    ! FACE_CNTR(:), referred as centers of mass.
    ! Relation for arbitrary shape CVs reported in Ferziger & Peric.


    ! Preliminary checks

    IF(.NOT.ALLOCATED(verts)) THEN
        WRITE(*,100)
        CALL abort_psblas
    END IF

    ! Can do without since they are INTENT(OUT)
!!$  if(    allocated(face_cntr) .or. &
!!$       & allocated(af)        .or. &
!!$       & allocated(area)) then
!!$     call abort_psblas
!!$  end if

    nfaces = nel_(v2f)

    CALL alloc_vector(face_cntr,nfaces)
    CALL alloc_vector(af,nfaces)
    ALLOCATE(area(nfaces),nvf(nfaces),stat=info)
    IF(info /= 0) THEN
        WRITE(*,200)
        CALL abort_psblas
    END IF

    af(:)        = vector_(0.d0,0.d0,0.d0)
    face_cntr(:) = vector_(0.d0,0.d0,0.d0)

    IF(ncd == 2) THEN
        v2 = vector_(0.d0,0.d0,1.d0)

        loop_2D: DO IF = 1, nfaces
            CALL get_ith_conn(iv2f,v2f,IF)
            i1 = iv2f(1)
            i2 = iv2f(2)
            v1 = verts(i2) - verts(i1)
            af(IF) = v1 .cross. v2
            area(IF) = mag(af(IF))

            face_cntr(IF) = verts(i1) + verts(i2)
            face_cntr(IF) = 0.5d0 * face_cntr(IF)
        END DO loop_2D

    ELSE IF(ncd == 3) THEN

        ! Exctracts number of face vertices from v2f table
        DO IF = 1, nfaces
            CALL get_ith_conn(iv2f,v2f,IF)
            nvf(IF) = SIZE(iv2f)
        END DO

        loop_3D: DO IF = 1, nfaces

            ! Retrieves V2F connectivity of the face IF
            CALL get_ith_conn(iv2f,v2f,IF)

            iv_1 = iv2f(1)        ! 1st IV in V2F connectivity of the face IF

            DO i = 3, nvf(IF)
                iv_im1 = iv2f(i-1) ! (I-1)-th IV in V2F connectivity of the face IF
                v1 = verts(iv_im1) - verts(iv_1)

                iv_i = iv2f(i)     ! I-th IV in V2F connectivity of the face IF
                v2 = verts(iv_i) - verts(iv_1)

                p = v1 .cross. v2
                p = 0.5d0 * p           ! Normal vector of the I-th triangle
                af(IF) = af(IF) + p     ! Normal vector of the generic polyside face

                w = mag(p)              ! Weighting factor of the I-th triangle
                p = verts(iv_1) + verts(iv_im1) + verts(iv_i)
                p = t3 * p                            ! Barycenter of the I-th triangle
                face_cntr(IF) = face_cntr(IF) + w * p ! Barycenter of the generic polygonal face
            END DO

            area(IF) = mag(af(IF))
            r = 1.d0 / area(IF)
            face_cntr(IF) = r * face_cntr(IF)
        END DO loop_3D

    END IF

    !#############################################################################
    IF(debug) THEN
        DO IF = 1, nfaces
            WRITE(*,500) IF, x_(af(IF)), y_(af(IF)), z_(af(IF)), area(IF)
        END DO
        WRITE(*,*)
500     FORMAT('face: ',i6,' | Normal:',3(1x,d13.6),' | Area: ',d13.6)
    END IF
    !#############################################################################

    NULLIFY(iv2f)
    DEALLOCATE(nvf)

100 FORMAT(' ERROR! Input array not associated in GEOM_FACE')
200 FORMAT(' ERROR! Memory allocation failure in GEOM_FACE')

END PROCEDURE geom_face
