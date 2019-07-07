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
    PUBLIC  :: scalar_pde                       ! Class
    PRIVATE :: pde ! Requuired by INTEL FC!
    ! INTEL Bug!
    ! The Intel compiler for some reason ignores the default PRIVATE
    ! clause and makes available also outside the module the PDE type.
    !

    TYPE, EXTENDS(PDE) :: scalar_pde
        PRIVATE
        REAL(psb_dpk_), ALLOCATABLE :: b(:)
    CONTAINS
        PROCEDURE, PUBLIC :: create_pde                   !! Constructor
        PROCEDURE, PUBLIC :: free_pde                     !! Destructor
        PROCEDURE, PUBLIC :: write_scalar_pde
        PROCEDURE, PUBLIC, PASS(pde) :: geins_scalar_pde
        GENERIC, PUBLIC :: geins_pde => geins_scalar_pde  !! Linear System Solving
        PROCEDURE, PUBLIC :: nemo_sizeof
        PROCEDURE, PUBLIC :: reinit_pde
        PROCEDURE, PUBLIC :: asb_pde_
        PROCEDURE, PRIVATE :: solve_scalar_pde
        GENERIC, PUBLIC :: solve_pde => solve_scalar_pde
    END TYPE scalar_pde

    !! ----- Generic Interfaces -----

    INTERFACE

        MODULE FUNCTION nemo_sizeof(eqn)
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(IN) :: eqn
            INTEGER(kind=nemo_int_long_)  :: nemo_sizeof
        END FUNCTION nemo_sizeof

        MODULE SUBROUTINE create_pde(eqn,input_file,sec,msh,dim)
            !! ----- Constructor -----
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(OUT)           :: eqn
            CHARACTER(len=*),  INTENT(IN)            :: input_file
            CHARACTER(len=*),  INTENT(IN)            :: sec
            TYPE(mesh),        INTENT(INOUT), TARGET :: msh
            TYPE(dimensions),  INTENT(IN)            :: dim
        END SUBROUTINE create_pde

        MODULE SUBROUTINE free_pde(eqn)
            !! ----- Destructor -----
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(INOUT) :: eqn
        END SUBROUTINE free_pde

        !! Getters

        MODULE FUNCTION get_scalar_pde_name(pde)
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(IN) :: pde
            CHARACTER(len=32) :: get_scalar_pde_name
        END FUNCTION get_scalar_pde_name

        MODULE FUNCTION get_scalar_pde_msh_fun(pde)
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(IN) :: pde
            TYPE(mesh), POINTER :: get_scalar_pde_msh_fun
        END FUNCTION get_scalar_pde_msh_fun

        !! ----- Temporary up to Gfortran patch -----

        MODULE SUBROUTINE get_scalar_pde_msh_sub(pde,msh)
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(IN) :: pde
            TYPE(mesh), POINTER :: msh
        END SUBROUTINE get_scalar_pde_msh_sub

        !! Linear System Solving
        MODULE SUBROUTINE geins_scalar_pde(n,ia,cloud,pde)
            IMPLICIT NONE
            !! Inserts a ``cloud'' of RHS terms into pde%b
            INTEGER,           INTENT(IN)    :: n
            INTEGER,           INTENT(IN)    :: ia(:)
            REAL(psb_dpk_),    INTENT(IN)    :: cloud(:)
            CLASS(scalar_pde), INTENT(INOUT) :: pde
        END SUBROUTINE geins_scalar_pde

        MODULE SUBROUTINE asb_pde_(eqn)
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(INOUT) :: eqn
        END SUBROUTINE asb_pde_

        MODULE SUBROUTINE solve_scalar_pde(pde, mats, phi,var)
            IMPLICIT NONE
            CLASS(scalar_pde),  INTENT(INOUT) :: pde
            TYPE(matptr),       INTENT(IN), OPTIONAL, POINTER :: mats(:)
            TYPE(scalar_field), INTENT(INOUT) :: phi
            REAL(psb_dpk_),     INTENT(OUT), OPTIONAL :: var
        END SUBROUTINE solve_scalar_pde

        MODULE SUBROUTINE reinit_pde(eqn)
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(INOUT) :: eqn
        END SUBROUTINE reinit_pde

        !! Output

        MODULE SUBROUTINE write_scalar_pde(eqn,mat,rhs)
            !! ----- Output -----
            IMPLICIT NONE
            CLASS(scalar_pde), INTENT(IN) :: eqn
            CHARACTER(len=*),  INTENT(IN) :: mat
            CHARACTER(len=*),  INTENT(IN) :: rhs
        END SUBROUTINE write_scalar_pde

    END INTERFACE

END MODULE class_scalar_pde
