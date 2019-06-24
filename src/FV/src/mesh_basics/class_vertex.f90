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
! $Id: class_vertex.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    to be added...
!
MODULE class_vertex

    USE class_psblas
    USE class_vector

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: vertex                                   ! Class
    PUBLIC :: vertex_, alloc_vertex, free_vertex       ! Constructor/Destructor
    PUBLIC :: bcast_vertex, g2l_vertex, l2g_vertex     ! Parallel Operations
    PUBLIC :: update_vertex_halo                       ! Parallel Operations(2)
    PUBLIC :: x_, y_, z_, position_, on_boundary_      ! Getters
    PUBLIC :: ASSIGNMENT(=)                            ! Setters
    PUBLIC :: OPERATOR(+), OPERATOR(-), OPERATOR(*), & ! Vector Algerbra
        &    OPERATOR(.dot.), OPERATOR(.cross.), mag

    TYPE vertex
        PRIVATE
        TYPE(vector) :: position
        LOGICAL :: on_boundary
    CONTAINS
        PROCEDURE, PRIVATE :: nemo_vertex_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_vertex_sizeof
    END TYPE vertex


    ! ----- Generic Interfaces -----

  INTERFACE
    ELEMENTAL MODULE FUNCTION nemo_vertex_sizeof(vtx)
        USE class_psblas, ONLY : nemo_int_long_
        IMPLICIT NONE
        CLASS(vertex), INTENT(IN) :: vtx
        INTEGER(kind=nemo_int_long_)   :: nemo_vertex_sizeof
    END FUNCTION nemo_vertex_sizeof
  END INTERFACE

  ! ----- Constructors -----

  ! Constructor
  INTERFACE vertex_

    ! Public default constructor
    ELEMENTAL MODULE FUNCTION vertex_1_(x,y,z,on_boundary)RESULT(res)
        IMPLICIT NONE
        TYPE(vertex) :: res
        REAL(psb_dpk_), INTENT(IN) :: x, y, z
        LOGICAL, INTENT(IN), OPTIONAL :: on_boundary
    END FUNCTION vertex_1_

    ELEMENTAL MODULE FUNCTION vertex_2_(position,on_boundary)RESULT(res)
        IMPLICIT NONE
        TYPE(vertex) :: res
        TYPE(vector), INTENT(IN) :: position
        LOGICAL, INTENT(IN), OPTIONAL :: on_boundary
    END FUNCTION vertex_2_

  END INTERFACE vertex_


  INTERFACE 

    ! Array constructor
    MODULE SUBROUTINE alloc_vertex(verts,n)
        IMPLICIT NONE
        TYPE(vertex), ALLOCATABLE :: verts(:)
        INTEGER, INTENT(IN) :: n
    END SUBROUTINE alloc_vertex


    ! ----- Destructor -----

    MODULE SUBROUTINE free_vertex(verts)
        IMPLICIT NONE
        TYPE(vertex), ALLOCATABLE  :: verts(:)
    END SUBROUTINE free_vertex


    ! ----- Parallel Operations -----

    MODULE SUBROUTINE bcast_vertex(verts)
        IMPLICIT NONE
        TYPE(vertex), ALLOCATABLE :: verts(:)
    END SUBROUTINE bcast_vertex


    MODULE SUBROUTINE g2l_vertex(verts,desc_v)
        USE psb_base_mod
        IMPLICIT NONE
        TYPE(vertex), ALLOCATABLE :: verts(:)
        TYPE(psb_desc_type), INTENT(IN) :: desc_v
    END SUBROUTINE g2l_vertex


    MODULE SUBROUTINE l2g_vertex(verts_loc,verts_glob,desc_v)
        USE psb_base_mod
        IMPLICIT NONE
        ! WARNING! The global results is allocated only on P0. After its usage
        ! it must be deallocated in the calling unit by means of the statement:
        ! "if(associated(glob_res)) deallocate(glob_res)"

        TYPE(vertex), ALLOCATABLE :: verts_loc(:)
        TYPE(vertex), ALLOCATABLE :: verts_glob(:)
        TYPE(psb_desc_type), INTENT(IN) :: desc_v
    END SUBROUTINE l2g_vertex


    MODULE SUBROUTINE update_vertex_halo(verts,desc)
      IMPLICIT NONE
      ! synchronizes vertex positions
      ! (other data could be synchronized, too, if needed)
      
      TYPE(vertex)                    :: verts(:)
      TYPE(psb_desc_type), INTENT(IN) :: desc
   END SUBROUTINE update_vertex_halo


    ! ----- Getters -----

    ELEMENTAL MODULE FUNCTION position_(vert)
      IMPLICIT NONE
      TYPE(vector) :: position_
      TYPE(vertex), INTENT(IN) :: vert
    END FUNCTION position_


    ELEMENTAL MODULE FUNCTION on_boundary_(vert)
      IMPLICIT NONE
      LOGICAL :: on_boundary_
      TYPE(vertex), INTENT(IN) :: vert
    END FUNCTION on_boundary_

  END INTERFACE

  ! Getters

  INTERFACE x_
    ELEMENTAL MODULE FUNCTION get_vertex_x(vert)
      IMPLICIT NONE
      REAL(psb_dpk_) :: get_vertex_x
      TYPE(vertex), INTENT(IN) :: vert
    END FUNCTION get_vertex_x
  END INTERFACE x_


  INTERFACE y_
    ELEMENTAL MODULE FUNCTION get_vertex_y(vert)
      IMPLICIT NONE
      REAL(psb_dpk_) :: get_vertex_y
      TYPE(vertex), INTENT(IN) :: vert
    END FUNCTION get_vertex_y
  END INTERFACE y_

  INTERFACE z_
    ELEMENTAL MODULE FUNCTION get_vertex_z(vert)
      IMPLICIT NONE
      REAL(psb_dpk_) :: get_vertex_z
      TYPE(vertex), INTENT(IN) :: vert
    END FUNCTION get_vertex_z
  END INTERFACE z_


  ! ----- Vector Algebra Operations -----

  INTERFACE ASSIGNMENT(=)

    ! ----- Setters -----

    MODULE SUBROUTINE set_vertex_position(vert,position)
      IMPLICIT NONE
      TYPE(vertex), INTENT(INOUT) :: vert
      TYPE(vector), INTENT(IN) :: position
    END SUBROUTINE set_vertex_position

  END INTERFACE ASSIGNMENT(=)


    ! ----- User-defined operators -----

  INTERFACE OPERATOR(+)

    PURE MODULE FUNCTION vert_sum_1(a,b)
      IMPLICIT NONE
      TYPE(vector) :: vert_sum_1
      TYPE(vertex), INTENT(IN) :: a, b
    END FUNCTION vert_sum_1

    PURE MODULE FUNCTION vert_sum_2(a,b)
      IMPLICIT NONE
      TYPE(vector) :: vert_sum_2
      TYPE(vector), INTENT(IN) :: a
      TYPE(vertex), INTENT(IN) :: b
    END FUNCTION vert_sum_2

  END INTERFACE OPERATOR(+)


  INTERFACE OPERATOR(-)
    PURE MODULE FUNCTION vert_diff(a,b)
      IMPLICIT NONE
      TYPE(vector) :: vert_diff
      TYPE(vertex), INTENT(IN) :: a, b
    END FUNCTION vert_diff
  END INTERFACE OPERATOR(-)


  INTERFACE OPERATOR(*)
    ELEMENTAL MODULE FUNCTION scalar_vertex_prod(alpha,v)
      IMPLICIT NONE
      TYPE(vertex) :: scalar_vertex_prod
      REAL(psb_dpk_), INTENT(IN) :: alpha
      TYPE(vertex), INTENT(IN) :: v
    END FUNCTION scalar_vertex_prod
  END INTERFACE OPERATOR(*)


  INTERFACE OPERATOR(.dot.)
    PURE MODULE FUNCTION dot_prod(a,b)
      IMPLICIT NONE
      REAL(psb_dpk_) :: dot_prod
      TYPE(vertex), INTENT(IN) :: a, b
    END FUNCTION dot_prod
  END INTERFACE OPERATOR(.dot.)


  INTERFACE OPERATOR(.cross.)
    PURE MODULE FUNCTION cross_prod(a,b)
      IMPLICIT NONE
      TYPE(vector) :: cross_prod
      TYPE(vertex), INTENT(IN) :: a, b
    END FUNCTION cross_prod
  END INTERFACE OPERATOR(.cross.)

  INTERFACE mag
    MODULE FUNCTION vert_mag(v)
      IMPLICIT NONE
      REAL(psb_dpk_) :: vert_mag
      TYPE(vertex), INTENT(IN) :: v
    END FUNCTION vert_mag
  END INTERFACE mag

END MODULE class_vertex
