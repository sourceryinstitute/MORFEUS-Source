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
! $Id: class_vector_pde.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_vector_pde

    USE class_psblas,       ONLY : psb_dspmat_type, psb_dpk_, nemo_int_long_
    USE class_pde,          ONLY : pde
    USE class_mesh,         ONLY : mesh
    USE class_vector,       ONLY : vector
    USE class_vector_field, ONLY : vector_field

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: vector_pde                       ! Class
    PUBLIC :: spins_pde, geins_pde             ! Linear System Solving
    PRIVATE :: pde ! Required by INTEL FC!
    ! INTEL Bug!
    ! The Intel compiler for some reason ignores the default PRIVATE
    ! clause and makes available also outside the module the PDE type.
    !

    TYPE vector_pde
        PRIVATE
        TYPE(pde) :: base
        REAL(psb_dpk_), ALLOCATABLE :: b(:,:)
    CONTAINS
        PROCEDURE, PRIVATE :: create_vector_pde, free_vector_pde
        GENERIC, PUBLIC :: create_pde => create_vector_pde   ! Constructor
        GENERIC, PUBLIC :: free_pde => free_vector_pde       ! Destructor
        PROCEDURE, PRIVATE :: get_vector_pde_dim, get_vector_pde_msh_fun  ! Getters
        GENERIC, PUBLIC :: dim_ => get_vector_pde_dim
        GENERIC, PUBLIC :: msh_ => get_vector_pde_msh_fun
        PROCEDURE, PRIVATE :: get_vector_pde_diag
        GENERIC, PUBLIC :: get_diag => get_vector_pde_diag
        PROCEDURE, PRIVATE :: get_vector_pde_A
        GENERIC, PUBLIC :: get_A => get_vector_pde_A
        PROCEDURE, PRIVATE :: reinit_vector_pde
        GENERIC, PUBLIC :: reinit_pde => reinit_vector_pde
        PROCEDURE, PRIVATE :: get_vector_pde_name
        GENERIC, PUBLIC :: name_ => get_vector_pde_name
        PROCEDURE, PRIVATE :: nemo_vector_pde_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_vector_pde_sizeof
        PROCEDURE, PRIVATE :: solve_vector_pde
        GENERIC, PUBLIC :: solve_pde => solve_vector_pde
        PROCEDURE, PRIVATE :: write_vector_pde
        GENERIC, PUBLIC :: write_pde => write_vector_pde
        PROCEDURE, PRIVATE :: get_vector_pde_msh_sub
        GENERIC, PUBLIC :: get_mesh => get_vector_pde_msh_sub
        PROCEDURE, PRIVATE :: asb_vector_pde
        GENERIC, PUBLIC :: asb_pde => asb_vector_pde
    END TYPE vector_pde

    INTERFACE spins_pde
        MODULE PROCEDURE :: spins_vector_pde
    END INTERFACE spins_pde

    INTERFACE geins_pde
        MODULE PROCEDURE :: geins_vector_pde_v
        MODULE PROCEDURE :: geins_vector_pde_r
    END INTERFACE geins_pde

    ! ----- Generic Interfaces -----

  INTERFACE

    MODULE FUNCTION nemo_vector_pde_sizeof(eqn)
        USE class_psblas, ONLY : nemo_int_long_
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(IN) :: eqn
        INTEGER(kind=nemo_int_long_)   :: nemo_vector_pde_sizeof
     END FUNCTION nemo_vector_pde_sizeof

    !! ----- Constructor -----
    MODULE SUBROUTINE create_vector_pde(pde,input_file,sec,msh,dim)
      !! Constructor
        USE class_dimensions, ONLY : dimensions
        IMPLICIT NONE

        CLASS(vector_pde), INTENT(OUT)           :: pde
        CHARACTER(len=*), INTENT(IN)            :: input_file
        CHARACTER(len=*), INTENT(IN)            :: sec
        TYPE(mesh),       INTENT(INOUT), TARGET :: msh
        TYPE(dimensions), INTENT(IN)            :: dim
    END SUBROUTINE create_vector_pde

    !! ----- Destructor -----
    MODULE SUBROUTINE free_vector_pde(pde)
      !! Destructor
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(INOUT) :: pde
    END SUBROUTINE free_vector_pde

  ! ----- Getters -----

    !! Getters
    MODULE FUNCTION get_vector_pde_name(pde)
        IMPLICIT NONE
        CHARACTER(len=32) :: get_vector_pde_name
        CLASS(vector_pde), INTENT(IN) :: pde
    END FUNCTION get_vector_pde_name

    MODULE FUNCTION get_vector_pde_dim(pde)
        USE class_dimensions, ONLY : dimensions
        IMPLICIT NONE
        TYPE(dimensions) :: get_vector_pde_dim
        CLASS(vector_pde), INTENT(IN) :: pde
    END FUNCTION get_vector_pde_dim

  !-----------------------------------------

    MODULE SUBROUTINE get_vector_pde_A(pde,A)
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(INOUT) :: pde
        TYPE(psb_dspmat_type)  :: A
    END SUBROUTINE get_vector_pde_A

  !-----------------------------------------

    MODULE FUNCTION get_vector_pde_msh_fun(pde)
        IMPLICIT NONE
        TYPE(mesh), POINTER :: get_vector_pde_msh_fun
        CLASS(vector_pde), INTENT(IN) :: pde
    END FUNCTION get_vector_pde_msh_fun

    MODULE SUBROUTINE get_vector_pde_diag(pde,d)
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(INOUT) :: pde
        REAL(psb_dpk_), ALLOCATABLE  :: d(:)
    END SUBROUTINE get_vector_pde_diag

  ! ----- Temporary up to Gfortran patch -----
    MODULE SUBROUTINE get_vector_pde_msh_sub(pde,msh)
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(IN) :: pde
        TYPE(mesh), POINTER :: msh
    END SUBROUTINE get_vector_pde_msh_sub
  ! ------------------------------------------

  ! Linear System Solving

    !! ----- Linear System Solving -----
    MODULE SUBROUTINE spins_vector_pde(n,ia,ja,cloud,pde)
        !! Inserts a ``cloud'' of coefficients into pde%A
        IMPLICIT NONE
        INTEGER,          INTENT(IN)    :: n
        INTEGER,          INTENT(IN)    :: ia(:), ja(:)
        REAL(psb_dpk_), INTENT(IN)    :: cloud(:)
        TYPE(vector_pde), INTENT(INOUT) :: pde
    END SUBROUTINE spins_vector_pde

    MODULE SUBROUTINE geins_vector_pde_v(n,ia,cloud,pde)
      !! Wrapper for ``clouds'' of VECTOR type
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN) :: ia(:)
        TYPE(vector),     INTENT(IN)    :: cloud(:)
        TYPE(vector_pde), INTENT(INOUT) :: pde
    END SUBROUTINE geins_vector_pde_v

    MODULE SUBROUTINE geins_vector_pde_r(n,ia,cloud,pde)
      !! Inserts a ``cloud'' of RHS terms into pde%b
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN) :: ia(:)
        REAL(psb_dpk_), INTENT(IN) :: cloud(:,:)
        TYPE(vector_pde), INTENT(INOUT) :: pde
    END SUBROUTINE geins_vector_pde_r

    MODULE SUBROUTINE asb_vector_pde(pde)
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(INOUT) :: pde
    END SUBROUTINE asb_vector_pde

    MODULE SUBROUTINE solve_vector_pde(pde,phi,var)
        !! Assigns the solution to the vector field
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(INOUT) :: pde
        TYPE(vector_field), INTENT(INOUT) :: phi
        REAL(psb_dpk_), INTENT(OUT), OPTIONAL :: var
    END SUBROUTINE solve_vector_pde

    MODULE SUBROUTINE reinit_vector_pde(pde)
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(INOUT) :: pde
    END SUBROUTINE reinit_vector_pde

  ! Output

    MODULE SUBROUTINE write_vector_pde(pde,mat,rhs)
      !! ----- Output -----
        IMPLICIT NONE
        CLASS(vector_pde), INTENT(IN) :: pde
        CHARACTER(len=*), INTENT(IN) :: mat, rhs
    END SUBROUTINE write_vector_pde
  END INTERFACE

END MODULE class_vector_pde
