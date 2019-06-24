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
! $Id: class_face.F90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    face class and face management
!
! Provides:
!    FACE            class      (no public data)
!    FACE_           constructor (elemental)
!    FREE_FACE       destructor
!    NV_             get number of cell faces
!    MASTER_         get the master cell for a face
!    SLAVE_          get the slave cell for a face
!    FLAG_           indicates of a face is interior, on a processor/physical bndry
!    BCAST_FACE      broadcast face structures from process 0 to other processes
!    G2L_FACE        global to local renumbering and reallocation of faces
!    SET_FACE        resets number of vertices, master, slave, or flag for a face
!
MODULE class_face

    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: face                           ! Class
    PUBLIC :: face_, alloc_face, free_face   ! Constructor/destructor
    PUBLIC :: bcast_face, g2l_face, l2g_face ! Parallel ops.

    TYPE face
        PRIVATE
        INTEGER :: nv
        INTEGER :: master
        INTEGER :: slave
        INTEGER :: flag
    CONTAINS
        PROCEDURE, PRIVATE :: get_face_nv, get_face_master, get_face_slave, get_face_flag   ! Getters
        GENERIC, PUBLIC :: nv_ => get_face_nv
        GENERIC, PUBLIC :: master_ => get_face_master
        GENERIC, PUBLIC :: slave_ => get_face_slave
        GENERIC, PUBLIC :: flag_ => get_face_flag
        PROCEDURE :: set_face                                  ! Setters
        PROCEDURE, PRIVATE :: nemo_face_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_face_sizeof
    END TYPE face

    !-------------------------|
    ! Face Flag | Face Type   |
    !-------------------------|
    !      0    | interior    |
    !    < 0    | proc. bndry |
    !    > 0    | phys. bndry |
    !-------------------------|

    ! ----- Generic Interfaces -----

  INTERFACE

    ELEMENTAL MODULE FUNCTION nemo_face_sizeof(fc)
        USE class_psblas, ONLY : nemo_int_long_
        IMPLICIT NONE
        CLASS(face), INTENT(IN) :: fc
        INTEGER(kind=nemo_int_long_)   :: nemo_face_sizeof
    END FUNCTION nemo_face_sizeof

    ! ----- Constructors -----

    ELEMENTAL MODULE FUNCTION face_(nv,master,slave,flag)
      !! Public default constructor
        IMPLICIT NONE
        TYPE(face) :: face_
        INTEGER, INTENT(IN) :: nv, master, slave, flag
    END FUNCTION face_


    MODULE SUBROUTINE alloc_face(faces,n)
      !! Array constructor
        IMPLICIT NONE
        TYPE(face), ALLOCATABLE :: faces(:)
        INTEGER, INTENT(IN) :: n
    END SUBROUTINE alloc_face


    ! ----- Destructor -----

    MODULE SUBROUTINE free_face(faces)
        IMPLICIT NONE
        TYPE(face), ALLOCATABLE :: faces(:)
    END SUBROUTINE free_face


    ! ----- Parallel Ops. -----

    MODULE SUBROUTINE bcast_face(faces)
        IMPLICIT NONE
        TYPE(face), ALLOCATABLE :: faces(:)
    END SUBROUTINE bcast_face

    MODULE SUBROUTINE g2l_face(faces,desc_f,desc_c)
        USE psb_base_mod
        IMPLICIT NONE
        TYPE(face), ALLOCATABLE :: faces(:)
        TYPE(psb_desc_type), INTENT(IN) :: desc_f, desc_c
    END SUBROUTINE g2l_face

    MODULE SUBROUTINE l2g_face(faces_loc,faces_glob,desc_f,desc_c)
        USE psb_base_mod
        IMPLICIT NONE
        !! WARNING! The global results is allocated only on P0. After its usage
        !! it must be deallocated in the calling unit by means of the statement:
        !! "if(associated(glob_res)) deallocate(glob_res)"
        TYPE(face), ALLOCATABLE :: faces_loc(:)
        TYPE(face), ALLOCATABLE :: faces_glob(:)
        TYPE(psb_desc_type), INTENT(IN) :: desc_f
        TYPE(psb_desc_type), INTENT(IN) :: desc_c
    END SUBROUTINE l2g_face

    ! ----- Getters -----

    ELEMENTAL MODULE FUNCTION get_face_nv(f)
        IMPLICIT NONE
        INTEGER :: get_face_nv
        CLASS(face), INTENT(IN) :: f
    END FUNCTION get_face_nv

    ELEMENTAL MODULE FUNCTION get_face_master(f)
        IMPLICIT NONE
        INTEGER :: get_face_master
        CLASS(face), INTENT(IN) :: f
    END FUNCTION get_face_master

    ! Getters

    ELEMENTAL MODULE FUNCTION get_face_slave(f)
        IMPLICIT NONE
        INTEGER :: get_face_slave
        CLASS(face), INTENT(IN) :: f
    END FUNCTION get_face_slave

    ELEMENTAL MODULE FUNCTION get_face_flag(f)
        IMPLICIT NONE
        INTEGER :: get_face_flag
        CLASS(face), INTENT(IN) :: f
    END FUNCTION get_face_flag

    ! ----- Setter -----
    MODULE SUBROUTINE set_face(f,nv,master,slave,flag)
        IMPLICIT NONE
        CLASS(face), INTENT(INOUT) :: f
        INTEGER, INTENT(IN), OPTIONAL :: nv, master, slave, flag
    END SUBROUTINE set_face

  END INTERFACE

END MODULE class_face
