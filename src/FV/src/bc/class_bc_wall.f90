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
! $Id: class_bc_wall.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    WALL boundary condition class.
!
MODULE class_bc_wall
    USE class_psblas,     ONLY : nemo_int_long_, psb_dpk_
    USE class_bc_math,    ONLY : bc_math
    USE class_dimensions, ONLY : dimensions
    USE class_vector,     ONLY : vector
    USE class_mesh,       ONLY : mesh
    USE class_material,   ONLY : matptr

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: bc_wall                      !! Class
    PUBLIC :: create_bc_wall, free_bc_wall !! Constructor/Destructor
    PUBLIC :: update_boundary_wall         !! Updater

    TYPE bc_wall
        PRIVATE
        INTEGER :: id(5)
        TYPE(bc_math) :: temp
        TYPE(bc_math) :: pressure
        TYPE(bc_math) :: conc
        TYPE(bc_math) :: vel(3)
        TYPE(bc_math) :: pos(3)
        TYPE(bc_math) :: stress(3)
    CONTAINS
        PROCEDURE, PRIVATE :: get_abc_wall_s, get_abc_wall_v
        GENERIC, PUBLIC :: get_abc_wall => get_abc_wall_s, get_abc_wall_v ! Getter
        PROCEDURE, PRIVATE :: set_bc_wall_map_s, set_bc_wall_map_v
        GENERIC, PUBLIC :: set_bc_wall_map => set_bc_wall_map_s, set_bc_wall_map_v  ! Setter
        PROCEDURE, PRIVATE :: nemo_bc_wall_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_bc_wall_sizeof
    END TYPE bc_wall


  ! ----- Generic Interfaces -----

  INTERFACE
    !! ----- Getter -----

    MODULE SUBROUTINE get_abc_wall_s(bc,dim,id,a,b,c)
        IMPLICIT NONE
        CLASS(bc_wall), INTENT(IN) :: bc
        TYPE(dimensions), INTENT(IN) :: dim
        INTEGER, INTENT(OUT) :: id
        REAL(psb_dpk_), INTENT(INOUT) :: a(:)
        REAL(psb_dpk_), INTENT(INOUT) :: b(:)
        REAL(psb_dpk_), INTENT(INOUT) :: c(:)
          !! WARNING! Use intent(inout) for ABC.
    END SUBROUTINE get_abc_wall_s


    MODULE SUBROUTINE get_abc_wall_v(bc,dim,id,a,b,c)
        IMPLICIT NONE
        CLASS(bc_wall),    INTENT(IN) :: bc
        TYPE(dimensions), INTENT(IN) :: dim
        INTEGER,          INTENT(OUT) :: id
        REAL(psb_dpk_), INTENT(INOUT) :: a(:)
        REAL(psb_dpk_), INTENT(INOUT) :: b(:)
        TYPE(vector),     INTENT(INOUT) :: c(:)
        ! WARNING! Use intent(inout) for ABC.
    END SUBROUTINE get_abc_wall_v

    !! ----- Setter -----

    MODULE SUBROUTINE set_bc_wall_map_s(bc,i,a,b,c)
!        USE class_vector
!        USE tools_bc
!        USE class_bc_math
        IMPLICIT NONE
        CLASS(bc_wall), INTENT(INOUT) :: bc
        INTEGER, INTENT(IN) :: i
        REAL(psb_dpk_), INTENT(IN) :: a, b, c
    END SUBROUTINE set_bc_wall_map_s

    MODULE SUBROUTINE set_bc_wall_map_v(bc,i,a,b,c)
!        USE class_vector
!        USE tools_bc
!        USE class_bc_math
        IMPLICIT NONE
        CLASS(bc_wall), INTENT(INOUT) :: bc
        INTEGER, INTENT(IN) :: i
        REAL(psb_dpk_), INTENT(IN) :: a, b
        TYPE(vector), INTENT(IN) :: c
    END SUBROUTINE set_bc_wall_map_v

    ELEMENTAL MODULE FUNCTION nemo_bc_wall_sizeof(bc)
      !! REMARK: the implementation of run-time polymorphism requires
      !! specific BC object as POINTERS
        IMPLICIT NONE
        CLASS(bc_wall), INTENT(IN) :: bc
        INTEGER(kind=nemo_int_long_)   :: nemo_bc_wall_sizeof
    END FUNCTION nemo_bc_wall_sizeof

    MODULE SUBROUTINE create_bc_wall(bc,input_file,sec,nbf)
        IMPLICIT NONE
      !! ----- Constructor -----
        TYPE(bc_wall), POINTER :: bc
        CHARACTER(len=*), INTENT(IN) :: input_file
        CHARACTER(len=*), INTENT(IN) :: sec
        INTEGER, INTENT(IN) :: nbf
    END SUBROUTINE create_bc_wall

    MODULE SUBROUTINE free_bc_wall(bc)
        IMPLICIT NONE
      !! ----- Destructor -----
        TYPE(bc_wall), POINTER :: bc
    END SUBROUTINE free_bc_wall

  END INTERFACE

  INTERFACE update_boundary_wall
    !! ----- Boundary Values Updater -----

    MODULE SUBROUTINE update_boundary_wall_s(ib,bc,dim,msh,mats,im,x,bx)
!        USE class_dimensions
!        USE class_face
!        USE class_material
!        USE class_mesh
!        USE tools_bc
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ib
        TYPE(bc_wall), INTENT(IN) :: bc
        TYPE(dimensions), INTENT(IN) :: dim
        TYPE(mesh), INTENT(IN) :: msh
        TYPE(matptr), INTENT(IN), POINTER:: mats(:)
        INTEGER, INTENT(IN) :: im(:)
        REAL(psb_dpk_), INTENT(IN) :: x(:)
        REAL(psb_dpk_), INTENT(INOUT) :: bx(:)
    END SUBROUTINE update_boundary_wall_s

    MODULE SUBROUTINE update_boundary_wall_v(ib,bc,dim,msh,x,bx)
        !! WARNING!
        !! - Use intent(inout) for BX(:).

        !! Number of boundary faces with flag < IB
 !       USE class_dimensions
 !       USE class_face
 !       USE class_mesh
 !       USE class_vector
 !       USE tools_bc
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ib
        TYPE(bc_wall), INTENT(IN) :: bc
        TYPE(dimensions), INTENT(IN) :: dim
        TYPE(mesh), INTENT(IN) :: msh
        TYPE(vector), INTENT(IN) :: x(:)
        TYPE(vector), INTENT(INOUT) :: bx(:)
    END SUBROUTINE update_boundary_wall_v

  END INTERFACE update_boundary_wall

END MODULE class_bc_wall
