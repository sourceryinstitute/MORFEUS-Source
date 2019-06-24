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
! $Id: class_bc_math.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    MATHematical boundary condition class.
!
MODULE class_bc_math

    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: bc_math                           ! Class
    PUBLIC :: create_bc_math, alloc_bc_math     ! Constructors
    PUBLIC :: free_bc_math, dealloc_bc_math     ! Destructors
    PUBLIC :: get_abc_math, is_allocated        ! Getter
    PUBLIC :: set_bc_math_map!, set_bc_math_uni ! Setter
    PUBLIC :: update_boundary_math, &
        &    apply_abc_to_boundary             ! Updater
    PUBLIC :: debug_bc_math                     ! Debug

    TYPE bc_math
        PRIVATE
        INTEGER :: id
        INTEGER :: nbf
        REAL(psb_dpk_), ALLOCATABLE :: a(:)
        REAL(psb_dpk_), ALLOCATABLE :: b(:)
        REAL(psb_dpk_), ALLOCATABLE :: c(:)
    CONTAINS
        PROCEDURE, PRIVATE :: nemo_bc_math_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_bc_math_sizeof
    END TYPE bc_math


    ! ----- Generic Interface -----


    ! ----- Description -----

    ! Id            | a    | b    | c
    !------------------------------------------
    ! bc_dirichlet_ | 1.d0 | 0.d0 | phi_b
    ! bc_neumann_   | 0.d0 | 1.d0 | grad(phi)_b
    ! bc_robin_     | a    | b    | c


  INTERFACE

    ELEMENTAL MODULE FUNCTION nemo_bc_math_sizeof(bc)
        USE class_psblas, ONLY : nemo_int_long_
        IMPLICIT NONE
        CLASS(bc_math), INTENT(IN) :: bc
        INTEGER(kind=nemo_int_long_)   :: nemo_bc_math_sizeof
    END FUNCTION nemo_bc_math_sizeof

    ! REMARK: the implementation of run-time polymorphism requires
    ! specific BC object as POINTERS!

    ! ----- Constructor -----

    MODULE SUBROUTINE create_bc_math(bc,input_file,sec,nbf)
        USE class_mesh
        USE tools_bc

        TYPE(bc_math), POINTER :: bc
        CHARACTER(len=*), INTENT(IN) :: input_file
        CHARACTER(len=*), INTENT(IN) :: sec
        INTEGER, INTENT(IN) :: nbf
    END SUBROUTINE create_bc_math


    MODULE SUBROUTINE alloc_bc_math(bc,id,nbf,a,b,c)
        USE tools_bc

        TYPE(bc_math), INTENT(INOUT) :: bc
        INTEGER, INTENT(IN) :: id
        INTEGER, INTENT(IN) :: nbf
        REAL(psb_dpk_), INTENT(IN) :: a(:), b(:), c(:)
    END SUBROUTINE alloc_bc_math


    ! ----- Destructor -----

    MODULE SUBROUTINE free_bc_math(bc)
        !! To be invoked when BC_MATH is used as high-level BC.
        TYPE(bc_math), POINTER :: bc
    END SUBROUTINE free_bc_math


    MODULE SUBROUTINE dealloc_bc_math(bc)
        !! To be invoked when BC_MATH is a member of another BC class.
        TYPE(bc_math) :: bc
    END SUBROUTINE dealloc_bc_math

  END INTERFACE


  ! ----- Getter -----

  INTERFACE get_abc_math

    MODULE SUBROUTINE get_abc_math_s(bc,id,a,b,c)
        USE class_connectivity
        USE class_face
        USE class_mesh
        USE tools_bc

        TYPE(bc_math), INTENT(IN) :: bc
        INTEGER, INTENT(OUT) :: id
        REAL(psb_dpk_), INTENT(INOUT) :: a(:)
        REAL(psb_dpk_), INTENT(INOUT) :: b(:)
        REAL(psb_dpk_), INTENT(INOUT) :: c(:)
    END SUBROUTINE get_abc_math_s


    MODULE SUBROUTINE get_abc_math_v(bc,id,a,b,c)
        USE class_connectivity
        USE class_face
        USE class_mesh
        USE class_vector
        USE tools_bc

        TYPE(bc_math), INTENT(IN) :: bc(:)
        INTEGER, INTENT(OUT) :: id
        REAL(psb_dpk_), INTENT(INOUT) :: a(:)
        REAL(psb_dpk_), INTENT(INOUT) :: b(:)
        TYPE(vector),     INTENT(INOUT) :: c(:)
    END SUBROUTINE get_abc_math_v

  END INTERFACE get_abc_math


  INTERFACE

    ELEMENTAL MODULE FUNCTION is_allocated(bc)
        LOGICAL :: is_allocated
        TYPE(bc_math), INTENT(IN) :: bc
    END FUNCTION is_allocated

  END INTERFACE


    ! ----- Setter -----

  INTERFACE set_bc_math

    MODULE SUBROUTINE set_bc_math_map(bc,i,a,b,c)
        USE tools_bc

        TYPE(bc_math), INTENT(INOUT) :: bc
        INTEGER, INTENT(IN) :: i
        REAL(psb_dpk_), INTENT(IN) :: a, b, c
    END SUBROUTINE set_bc_math_map

  END INTERFACE set_bc_math

    ! ----- Boundary Values Updating -----

  INTERFACE apply_abc_to_boundary

    MODULE SUBROUTINE apply_abc_to_boundary_s(id,a,b,c,ib,msh,x,bx)
        !! WARNING!
        !! - Use intent(inout) for BX(:)
        !! - Do loop on the faces subset corresponding to IB bc.
        !! - Only this section of BX(:) is going to be modified.
        !! - BX(:) indexing starts from 1: when BX(:) is referenced an offset
        !!   equal to the # of boundary faces with flag > 0 and < IB must be
        !!   added to the I counter.
        USE class_connectivity
        USE class_face
        USE class_mesh
        USE tools_bc

        INTEGER, INTENT(IN) :: id
        REAL(psb_dpk_), INTENT(IN) :: a(:)
        REAL(psb_dpk_), INTENT(IN) :: b(:)
        REAL(psb_dpk_), INTENT(IN) :: c(:)
        INTEGER, INTENT(IN) :: ib
        TYPE(mesh), INTENT(IN) :: msh
        REAL(psb_dpk_), INTENT(IN) :: x(:)
        REAL(psb_dpk_), INTENT(INOUT) :: bx(:)
    END SUBROUTINE apply_abc_to_boundary_s


    MODULE SUBROUTINE apply_abc_to_boundary_v(id,a,b,c,ib,msh,x,bx)
        USE class_connectivity
        USE class_face
        USE class_mesh
        USE class_vector
        USE tools_bc

        INTEGER, INTENT(IN) :: id
        REAL(psb_dpk_), INTENT(IN) :: a(:)
        REAL(psb_dpk_), INTENT(IN) :: b(:)
        TYPE(vector),     INTENT(IN) :: c(:)
        INTEGER,          INTENT(IN) :: ib
        TYPE(mesh),       INTENT(IN) :: msh
        TYPE(vector),     INTENT(IN) :: x(:)
        TYPE(vector),     INTENT(INOUT) :: bx(:)
    END SUBROUTINE apply_abc_to_boundary_v

  END INTERFACE apply_abc_to_boundary


  INTERFACE

    MODULE SUBROUTINE update_boundary_math(ib,bc,msh,x,bx)
        !! WARNING! Use intent(inout) for BX(:)
        USE class_mesh

        INTEGER, INTENT(IN) :: ib
        TYPE(bc_math), INTENT(IN) :: bc
        TYPE(mesh), INTENT(IN) :: msh
        REAL(psb_dpk_), INTENT(IN) :: x(:)
        REAL(psb_dpk_), INTENT(INOUT) :: bx(:)
    END SUBROUTINE update_boundary_math


    ! ----- Debug -----

    MODULE SUBROUTINE debug_bc_math(bc)
        USE tools_bc

        TYPE(bc_math), INTENT(IN) :: bc
    END SUBROUTINE debug_bc_math

  END INTERFACE
END MODULE class_bc_math
