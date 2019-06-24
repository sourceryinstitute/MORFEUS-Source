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
! $Id: class_vector.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_vector

    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: vector                                   ! Class
    PUBLIC :: vector_, alloc_vector, free_vector       ! Constructor/destructor
    PUBLIC :: bcast_vector, g2l_vector, l2g_vector     ! Parallel Operations
    PUBLIC :: x_, y_, z_                               ! Getters
    PUBLIC :: set_x_, set_y_, set_z_                   ! Setters
    PUBLIC :: update_vector_halo                       ! Setters
    PUBLIC :: OPERATOR(+), OPERATOR(-), OPERATOR(*), & ! Vector Algebra
        &    OPERATOR(.dot.), OPERATOR(.cross.),    &
        &    mag, unit, OPERATOR(==),  ASSIGNMENT(=)
    PUBLIC :: mandatory_v_

    TYPE vector
        PRIVATE
        REAL(psb_dpk_) :: x
        REAL(psb_dpk_) :: y
        REAL(psb_dpk_) :: z
    CONTAINS
        PROCEDURE, PRIVATE :: nemo_vector_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_vector_sizeof
    END TYPE vector


    ! ----- Generic Interfaces -----

    ! Getters
    INTERFACE x_
        ELEMENTAL MODULE FUNCTION get_vector_x(vert)
            IMPLICIT NONE
            REAL(psb_dpk_) :: get_vector_x
            TYPE(vector), INTENT(IN) :: vert
        END FUNCTION get_vector_x
    END INTERFACE x_

    INTERFACE y_
        ELEMENTAL MODULE FUNCTION get_vector_y(vert)
            IMPLICIT NONE
            REAL(psb_dpk_) :: get_vector_y
            TYPE(vector), INTENT(IN) :: vert
        END FUNCTION get_vector_y
    END INTERFACE y_

    INTERFACE z_
        ELEMENTAL MODULE FUNCTION get_vector_z(vert)
            IMPLICIT NONE
            REAL(psb_dpk_) :: get_vector_z
            TYPE(vector), INTENT(IN) :: vert
        END FUNCTION get_vector_z
    END INTERFACE z_

    ! Setters

    INTERFACE set_x_
        MODULE SUBROUTINE set_vector_x(vect,r)
            IMPLICIT NONE
            REAL(psb_dpk_)::r
            TYPE(vector), INTENT(INOUT) :: vect
        END SUBROUTINE set_vector_x
    END INTERFACE set_x_

    INTERFACE set_y_
        MODULE SUBROUTINE set_vector_y(vect,r)
            IMPLICIT NONE
            REAL(psb_dpk_)::r
            TYPE(vector), INTENT(INOUT) :: vect
        END SUBROUTINE set_vector_y
    END INTERFACE set_y_

    INTERFACE set_z_
        MODULE SUBROUTINE set_vector_z(vect,r)
            IMPLICIT NONE
            REAL(psb_dpk_)::r
            TYPE(vector), INTENT(INOUT) :: vect
        END SUBROUTINE set_vector_z
    END INTERFACE set_z_

    ! ----- Operators Overloading -----

    INTERFACE OPERATOR(+)
        ELEMENTAL MODULE FUNCTION vec_sum(a,b)
            IMPLICIT NONE
            TYPE(vector) :: vec_sum
            TYPE(vector), INTENT(IN) :: a, b
        END FUNCTION vec_sum
    END INTERFACE OPERATOR(+)

    INTERFACE OPERATOR(-)
        ELEMENTAL MODULE FUNCTION vec_diff(a,b)
            IMPLICIT NONE
            TYPE(vector) :: vec_diff
            TYPE(vector), INTENT(IN) :: a, b
        END FUNCTION vec_diff

        ELEMENTAL MODULE FUNCTION vec_minus(v)
            IMPLICIT NONE
            TYPE(vector) :: vec_minus
            TYPE(vector), INTENT(IN) :: v
        END FUNCTION vec_minus
    END INTERFACE OPERATOR(-)

    INTERFACE OPERATOR(*)
        ELEMENTAL MODULE FUNCTION scalar_vector_prod(alpha,v)
            IMPLICIT NONE
            TYPE(vector) :: scalar_vector_prod
            REAL(psb_dpk_), INTENT(IN) :: alpha
            TYPE(vector), INTENT(IN) :: v
        END FUNCTION scalar_vector_prod
    END INTERFACE OPERATOR(*)

    INTERFACE OPERATOR(.dot.)
        PURE MODULE FUNCTION dot_prod_t(a,b)
         ! Used for GRAD .dot. DELTA in VECTOR_PDE_LAPLACIAN.
         ! GRAD is a tensor. It would need a proper class...
            IMPLICIT NONE
            TYPE(vector) :: dot_prod_t
            TYPE(vector), INTENT(IN) :: a(:), b
        END FUNCTION dot_prod_t


        PURE MODULE FUNCTION dot_prod_v(a,b)
            IMPLICIT NONE
            REAL(psb_dpk_) :: dot_prod_v
            TYPE(vector), INTENT(IN) :: a, b
        END FUNCTION dot_prod_v
    END INTERFACE OPERATOR(.dot.)

    INTERFACE OPERATOR(.cross.)
        PURE MODULE FUNCTION cross_prod(a,b)
            IMPLICIT NONE
            TYPE(vector) :: cross_prod
            TYPE(vector), INTENT(IN) :: a, b
        END FUNCTION cross_prod
    END INTERFACE OPERATOR(.cross.)

    INTERFACE mag
        ELEMENTAL MODULE FUNCTION vec_mag(v)
            IMPLICIT NONE
            REAL(psb_dpk_) :: vec_mag
            TYPE(vector), INTENT(IN) :: v
        END FUNCTION vec_mag
    END INTERFACE mag

    INTERFACE unit
        ELEMENTAL MODULE FUNCTION vec_unit(v)
        ! Returns a unit vector in the direction of V
        !
            IMPLICIT NONE
            TYPE(vector)             :: vec_unit
            TYPE(vector), INTENT(IN) :: v
        END FUNCTION vec_unit
    END INTERFACE unit

    INTERFACE OPERATOR(==)
        PURE MODULE FUNCTION vec_eq(a,b)
            IMPLICIT NONE
            LOGICAL :: vec_eq
            TYPE(vector), INTENT(IN) :: a, b
        END FUNCTION vec_eq
    END INTERFACE OPERATOR(==)

    INTERFACE ASSIGNMENT(=)
        PURE MODULE SUBROUTINE vector_to_array(a,b)
            IMPLICIT NONE
            REAL(psb_dpk_), INTENT(OUT) :: a(3)
            TYPE(vector),      INTENT(IN)  :: b
        END SUBROUTINE vector_to_array
    END INTERFACE ASSIGNMENT(=)

    INTERFACE
        ELEMENTAL MODULE FUNCTION nemo_vector_sizeof(vec)
            USE psb_base_mod
            IMPLICIT NONE
            CLASS(vector), INTENT(IN) :: vec
            INTEGER(kind=nemo_int_long_)   :: nemo_vector_sizeof
        END FUNCTION nemo_vector_sizeof

        ! ----- Constructors -----

        ! Public default constructor
        ELEMENTAL MODULE FUNCTION vector_(x,y,z)
            IMPLICIT NONE
            TYPE(vector) :: vector_
            REAL(psb_dpk_), INTENT(IN) :: x, y, z
        END FUNCTION vector_


        ! Array constructor
        MODULE SUBROUTINE alloc_vector(vect,n)
            IMPLICIT NONE
            TYPE(vector), ALLOCATABLE :: vect(:)
            INTEGER, INTENT(IN) :: n
        END SUBROUTINE alloc_vector 
    END INTERFACE

    ! ----- Destructor -----
    INTERFACE
        MODULE SUBROUTINE free_vector(vect)
            IMPLICIT NONE
            TYPE(vector), ALLOCATABLE :: vect(:)
        END SUBROUTINE free_vector
    END INTERFACE

    ! ----- Parallel Operations -----
    INTERFACE
        MODULE SUBROUTINE bcast_vector(vect)
            IMPLICIT NONE
            TYPE(vector), ALLOCATABLE :: vect(:)
        END SUBROUTINE bcast_vector
        
        MODULE SUBROUTINE g2l_vector(verts,desc_v)
            USE psb_base_mod
            IMPLICIT NONE
            TYPE(vector), ALLOCATABLE  :: verts(:)
            TYPE(psb_desc_type), INTENT(IN) :: desc_v
        END SUBROUTINE g2l_vector

        MODULE SUBROUTINE l2g_vector(verts_loc,verts_glob,desc_v)
            USE psb_base_mod
            ! WARNING! The global results is allocated only on P0. After its usage
            ! it must be deallocated in the calling unit by means of the statement:
            ! "if(associated(glob_res)) deallocate(glob_res)"
            IMPLICIT NONE
            TYPE(vector)               :: verts_loc(:)
            TYPE(vector), ALLOCATABLE  :: verts_glob(:)
            TYPE(psb_desc_type), INTENT(IN) :: desc_v
        END SUBROUTINE l2g_vector
    END INTERFACE

    INTERFACE
         MODULE SUBROUTINE update_vector_halo(v,desc)
             USE psb_base_mod
             IMPLICIT NONE
             TYPE(vector), INTENT(INOUT) :: v(:)
             TYPE(psb_desc_type), INTENT(IN) :: desc
         END SUBROUTINE update_vector_halo 
     END INTERFACE

    ! ----- Named Constants -----

    TYPE(vector), PARAMETER :: mandatory_v_ = vector(-999.d0,999.d0,999.d0)

END MODULE class_vector
