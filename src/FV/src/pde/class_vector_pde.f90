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
    PUBLIC :: vector_pde            ! Class
    PRIVATE :: pde ! Required by INTEL FC!
    ! INTEL Bug!
    ! The Intel compiler for some reason ignores the default PRIVATE
    ! clause and makes available also outside the module the PDE type.
    !

    TYPE, EXTENDS(PDE) :: vector_pde
        PRIVATE
        REAL(psb_dpk_), ALLOCATABLE :: b(:,:)
    CONTAINS
        PROCEDURE, PUBLIC :: create_pde
        PROCEDURE, PUBLIC :: free_pde
        PROCEDURE, PUBLIC :: write_vector_pde
        PROCEDURE, PUBLIC :: reinit_pde
        PROCEDURE, PUBLIC :: nemo_sizeof
        PROCEDURE, PRIVATE :: solve_vector_pde
        GENERIC, PUBLIC :: solve_pde => solve_vector_pde
        PROCEDURE, PUBLIC :: asb_pde_
        PROCEDURE, PUBLIC, PASS(pde) :: geins_vector_pde_v
        PROCEDURE, PUBLIC, PASS(pde) :: geins_vector_pde_r
        GENERIC, PUBLIC :: geins_pde =>  geins_vector_pde_v, geins_vector_pde_r      ! Linear System Solving
    END TYPE vector_pde


    INTERFACE

        MODULE FUNCTION nemo_sizeof(eqn)
            USE class_psblas, ONLY : nemo_int_long_
            IMPLICIT NONE
            CLASS(vector_pde), INTENT(IN) :: eqn
            INTEGER(kind=nemo_int_long_)  :: nemo_sizeof
        END FUNCTION nemo_sizeof

        !! ----- Constructor -----

        MODULE SUBROUTINE create_pde(eqn,input_file,sec,msh,dim)
            !! Constructor
            USE class_dimensions, ONLY : dimensions
            IMPLICIT NONE
            CLASS(vector_pde), INTENT(OUT)           :: eqn
            CHARACTER(len=*),  INTENT(IN)            :: input_file
            CHARACTER(len=*),  INTENT(IN)            :: sec
            TYPE(mesh),        INTENT(INOUT), TARGET :: msh
            TYPE(dimensions),  INTENT(IN)            :: dim
        END SUBROUTINE create_pde

        !! ----- Destructor -----

        MODULE SUBROUTINE free_pde(eqn)
            !! Destructor
            IMPLICIT NONE
            CLASS(vector_pde), INTENT(INOUT) :: eqn
        END SUBROUTINE free_pde

        !! ----- Getters -----

        MODULE SUBROUTINE get_vector_pde_A(pde,A)
            IMPLICIT NONE
            CLASS(vector_pde), INTENT(INOUT) :: pde
            TYPE(psb_dspmat_type)  :: A
        END SUBROUTINE get_vector_pde_A

        MODULE SUBROUTINE get_vector_pde_diag(pde,d)
            IMPLICIT NONE
            CLASS(vector_pde), INTENT(INOUT) :: pde
            REAL(psb_dpk_),    ALLOCATABLE   :: d(:)
        END SUBROUTINE get_vector_pde_diag

        !! ----- Linear System Solving -----

        MODULE SUBROUTINE geins_vector_pde_v(n,ia,cloud,pde)
            !! Wrapper for ``clouds'' of VECTOR type
            IMPLICIT NONE
            INTEGER,           INTENT(IN)    :: n
            INTEGER,           INTENT(IN)    :: ia(:)
            TYPE(vector),      INTENT(IN)    :: cloud(:)
            CLASS(vector_pde), INTENT(INOUT) :: pde
        END SUBROUTINE geins_vector_pde_v

        MODULE SUBROUTINE geins_vector_pde_r(n,ia,cloud,pde)
            !! Inserts a ``cloud'' of RHS terms into pde%b
            IMPLICIT NONE
            INTEGER,           INTENT(IN)    :: n
            INTEGER,           INTENT(IN)    :: ia(:)
            REAL(psb_dpk_),    INTENT(IN)    :: cloud(:,:)
            CLASS(vector_pde), INTENT(INOUT) :: pde
        END SUBROUTINE geins_vector_pde_r

        MODULE SUBROUTINE asb_pde_(eqn)
            IMPLICIT NONE
            CLASS(vector_pde), INTENT(INOUT) :: eqn
        END SUBROUTINE asb_pde_

        MODULE SUBROUTINE solve_vector_pde(pde,phi,var)
            !! Assigns the solution to the vector field
            IMPLICIT NONE
            CLASS(vector_pde),  INTENT(INOUT) :: pde
            TYPE(vector_field), INTENT(INOUT) :: phi
            REAL(psb_dpk_),     INTENT(OUT), OPTIONAL :: var
        END SUBROUTINE solve_vector_pde

        MODULE SUBROUTINE reinit_pde(eqn)
            IMPLICIT NONE
            CLASS(vector_pde), INTENT(INOUT) :: eqn
        END SUBROUTINE reinit_pde

        !! Output

        MODULE SUBROUTINE write_vector_pde(eqn,mat,rhs)
            !! ----- Output -----
            IMPLICIT NONE
            CLASS(vector_pde), INTENT(IN) :: eqn
            CHARACTER(len=*),  INTENT(IN) :: mat
            CHARACTER(len=*),  INTENT(IN) :: rhs
        END SUBROUTINE write_vector_pde

    END INTERFACE

END MODULE class_vector_pde
