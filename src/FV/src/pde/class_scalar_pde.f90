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
! $Id: class_scalar_pde.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_scalar_pde
    USE class_psblas
    USE class_pde
    USE class_material
    USE class_mesh
    USE psb_base_mod
    USE class_dimensions!, ONLY : dimensions
    USE class_scalar_field, ONLY : scalar_field
    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: scalar_pde                       ! Class
    PUBLIC :: spins_pde, geins_pde             ! Linear System Solving
    PRIVATE :: pde ! Requuired by INTEL FC!
    ! INTEL Bug!
    ! The Intel compiler for some reason ignores the default PRIVATE
    ! clause and makes available also outside the module the PDE type.
    !

    TYPE scalar_pde
        PRIVATE
        TYPE(pde) :: base
        REAL(psb_dpk_), ALLOCATABLE :: b(:)
    CONTAINS
        PROCEDURE, PRIVATE :: create_scalar_pde, free_scalar_pde
        GENERIC, PUBLIC :: create_pde => create_scalar_pde   ! Constructor
        GENERIC, PUBLIC :: free_pde =>  free_scalar_pde      ! Destructor
        PROCEDURE, PRIVATE :: get_scalar_pde_dim, get_scalar_pde_msh_fun  ! Getters
        GENERIC, PUBLIC :: dim_ => get_scalar_pde_dim
        GENERIC, PUBLIC :: msh_ => get_scalar_pde_msh_fun
        PROCEDURE, PRIVATE :: reinit_scalar_pde
        GENERIC, PUBLIC :: reinit_pde => reinit_scalar_pde
        PROCEDURE, PRIVATE :: get_scalar_pde_name
        GENERIC, PUBLIC :: name_ => get_scalar_pde_name
        PROCEDURE, PRIVATE :: nemo_scalar_pde_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_scalar_pde_sizeof
        PROCEDURE, PRIVATE :: get_scalar_pde_msh_sub
        GENERIC, PUBLIC :: get_mesh => get_scalar_pde_msh_sub
        PROCEDURE, PRIVATE :: asb_scalar_pde
        GENERIC, PUBLIC :: asb_pde => asb_scalar_pde
        PROCEDURE, PRIVATE :: solve_scalar_pde
        GENERIC, PUBLIC :: solve_pde => solve_scalar_pde
        PROCEDURE, PRIVATE :: write_scalar_pde
        GENERIC, PUBLIC :: write_pde => write_scalar_pde     ! Output
    END TYPE scalar_pde

  ! ----- Generic Interfaces -----

    INTERFACE spins_pde
        PROCEDURE :: spins_scalar_pde
    END INTERFACE spins_pde

    INTERFACE geins_pde
        PROCEDURE :: geins_scalar_pde
    END INTERFACE geins_pde

    INTERFACE

        MODULE FUNCTION nemo_scalar_pde_sizeof(eqn)
        IMPLICIT NONE
        CLASS(scalar_pde), INTENT(IN) :: eqn
        INTEGER(kind=nemo_int_long_)   :: nemo_scalar_pde_sizeof
        END FUNCTION nemo_scalar_pde_sizeof

        MODULE SUBROUTINE create_scalar_pde(pde,input_file,sec,msh,dim)
        !! ----- Constructor -----
        IMPLICIT NONE
        CLASS(scalar_pde), INTENT(OUT)           :: pde
        CHARACTER(len=*), INTENT(IN)            :: input_file
        CHARACTER(len=*), INTENT(IN)            :: sec
        TYPE(mesh),       INTENT(INOUT), TARGET :: msh
        TYPE(dimensions), INTENT(IN)            :: dim
        END SUBROUTINE create_scalar_pde

        MODULE SUBROUTINE free_scalar_pde(pde)
        !! ----- Destructor -----
        IMPLICIT NONE
        CLASS(scalar_pde), INTENT(INOUT) :: pde
        END SUBROUTINE free_scalar_pde

        ! Getters

        MODULE FUNCTION get_scalar_pde_name(pde)
        IMPLICIT NONE
        CHARACTER(len=32) :: get_scalar_pde_name
        CLASS(scalar_pde), INTENT(IN) :: pde
        END FUNCTION get_scalar_pde_name

        MODULE FUNCTION get_scalar_pde_dim(pde)
        IMPLICIT NONE
        TYPE(dimensions) :: get_scalar_pde_dim
        CLASS(scalar_pde), INTENT(IN) :: pde
        END FUNCTION get_scalar_pde_dim

        MODULE FUNCTION get_scalar_pde_msh_fun(pde)
        IMPLICIT NONE
        TYPE(mesh), POINTER :: get_scalar_pde_msh_fun
        CLASS(scalar_pde), INTENT(IN) :: pde
        END FUNCTION get_scalar_pde_msh_fun

        ! ----- Temporary up to Gfortran patch -----

        MODULE SUBROUTINE get_scalar_pde_msh_sub(pde,msh)
        IMPLICIT NONE
        CLASS(scalar_pde), INTENT(IN) :: pde
        TYPE(mesh), POINTER :: msh
        END SUBROUTINE get_scalar_pde_msh_sub

        ! Linear System Solving
        MODULE SUBROUTINE spins_scalar_pde(n,ia,ja,cloud,pde)
        !! ----- Linear System Solving -----
        !! Inserts a ``cloud'' of coefficients into pde%A
        IMPLICIT NONE
        INTEGER,          INTENT(IN)    :: n
        INTEGER,          INTENT(IN)    :: ia(:), ja(:)
        REAL(psb_dpk_), INTENT(IN)    :: cloud(:)
        TYPE(scalar_pde), INTENT(INOUT) :: pde
        END SUBROUTINE spins_scalar_pde

        MODULE SUBROUTINE geins_scalar_pde(n,ia,cloud,pde)
        IMPLICIT NONE
        ! Inserts a ``cloud'' of RHS terms into pde%b
        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(IN) :: ia(:)
        REAL(psb_dpk_), INTENT(IN) :: cloud(:)
        TYPE(scalar_pde), INTENT(INOUT) :: pde
        END SUBROUTINE geins_scalar_pde

        MODULE SUBROUTINE asb_scalar_pde(pde)
        IMPLICIT NONE
        CLASS(scalar_pde), INTENT(INOUT) :: pde
        END SUBROUTINE asb_scalar_pde

        MODULE SUBROUTINE solve_scalar_pde(pde, mats, phi,var)
        IMPLICIT NONE
        CLASS(scalar_pde), INTENT(INOUT) :: pde
        TYPE(matptr),   INTENT(IN), OPTIONAL, POINTER :: mats(:)
        TYPE(scalar_field), INTENT(INOUT) :: phi
        REAL(psb_dpk_), INTENT(OUT), OPTIONAL :: var
        END SUBROUTINE solve_scalar_pde

        MODULE SUBROUTINE reinit_scalar_pde(pde)
        IMPLICIT NONE
        CLASS(scalar_pde), INTENT(INOUT) :: pde
        END SUBROUTINE reinit_scalar_pde

        ! Output

        MODULE SUBROUTINE write_scalar_pde(pde,mat,rhs)
        !! ----- Output -----
        IMPLICIT NONE
        CLASS(scalar_pde), INTENT(IN) :: pde
        CHARACTER(len=*), INTENT(IN) :: mat, rhs
        END SUBROUTINE write_scalar_pde

    END INTERFACE

END MODULE class_scalar_pde
