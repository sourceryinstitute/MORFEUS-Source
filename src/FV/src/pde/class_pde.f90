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
! $Id: class_pde.f90 9102 2015-04-24 16:06:49Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_pde

    USE class_psblas
    USE class_dimensions
    USE class_mesh

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: pde                              ! Class
    PUBLIC :: create_pde, free_pde             ! Constructor/destructor
    PUBLIC :: name_, dim_, msh_                ! Getters
    PUBLIC :: get_diag, get_A, update_diag            ! Getters & Setters
    PUBLIC :: is_pde_bld, is_pde_asb           ! Status inquirer
    PUBLIC :: spins_pde, &            ! Linear System Solving
        &    free_pde_prec, build_pde_prec, & !           "
        &    solve_pde_sys, reinit_pde        !           "
    PUBLIC :: write_pde                        ! Output

    TYPE pde
        PRIVATE
        CHARACTER(len=32) :: name            ! Name
        TYPE(dimensions) :: dim              ! Dimensions
        TYPE(mesh), POINTER :: msh => NULL() ! Mesh
        TYPE(psb_dspmat_type) :: A           ! PSBLAS
        REAL(psb_dpk_), ALLOCATABLE :: diag(:)  ! A's diag

        ! Linear System
        TYPE(psb_dprec_type) :: prec
        CHARACTER(len=10) :: cmethod
        CHARACTER(len=10) :: cprec
        INTEGER :: nlev
        REAL(psb_dpk_) :: eps_solv
        INTEGER :: itmax_solv
        LOGICAL :: mtx_sys

        ! Status
        INTEGER :: status
    CONTAINS
        PROCEDURE, PRIVATE :: nemo_pde_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_pde_sizeof
        PROCEDURE, PRIVATE :: get_pde_msh_sub
        GENERIC, PUBLIC :: get_mesh => get_pde_msh_sub
        PROCEDURE, PUBLIC :: asb_pde
    END TYPE pde


    ! ----- Generic Interfaces -----

  INTERFACE
    MODULE FUNCTION nemo_pde_sizeof(eqn)
        USE class_psblas
        IMPLICIT NONE
        CLASS(pde), INTENT(IN) :: eqn
        INTEGER(kind=nemo_int_long_)   :: nemo_pde_sizeof
    END FUNCTION nemo_pde_sizeof
  END INTERFACE

  ! Constructor
  INTERFACE create_pde
    MODULE SUBROUTINE create_pde(eqn,input_file,sec,msh,dim)
        USE class_connectivity
        USE tools_input
        IMPLICIT NONE
        TYPE(pde),        INTENT(OUT)           :: eqn
        CHARACTER(len=*), INTENT(IN)            :: input_file
        CHARACTER(len=*), INTENT(IN)            :: sec
        TYPE(mesh),       INTENT(INOUT), TARGET :: msh
        TYPE(dimensions), INTENT(IN)            :: dim
    END SUBROUTINE create_pde
  END INTERFACE create_pde

  INTERFACE free_pde
    !! ----- Destructor -----
    MODULE SUBROUTINE free_pde(eqn)
      !!  Destructor
        IMPLICIT NONE
        TYPE(pde), INTENT(INOUT) :: eqn
    END SUBROUTINE free_pde
  END INTERFACE free_pde

  ! ----- Getters -----

  INTERFACE name_
    !! Getters
    MODULE FUNCTION get_pde_name(eqn)
        IMPLICIT NONE
        CHARACTER(len=32) :: get_pde_name
        TYPE(pde), INTENT(IN) :: eqn
    END FUNCTION get_pde_name
  END INTERFACE name_


  INTERFACE dim_
    MODULE FUNCTION get_pde_dim(eqn)
        IMPLICIT NONE
        TYPE(dimensions) :: get_pde_dim
        TYPE(pde), INTENT(IN) :: eqn
    END FUNCTION get_pde_dim
  END INTERFACE dim_


  INTERFACE get_A
    MODULE SUBROUTINE get_pde_A(eqn,B)
        IMPLICIT NONE
        TYPE(psb_dspmat_type) :: B
        TYPE(pde), INTENT(INOUT) :: eqn
    END SUBROUTINE get_pde_A
  END INTERFACE get_A


  INTERFACE msh_
    MODULE FUNCTION get_pde_msh_fun(eqn)
        IMPLICIT NONE
        TYPE(mesh), POINTER :: get_pde_msh_fun
        TYPE(pde), INTENT(IN) :: eqn
    END FUNCTION get_pde_msh_fun
  END INTERFACE msh_


  INTERFACE get_diag
    MODULE SUBROUTINE get_pde_diag(eqn,d)
        IMPLICIT NONE
        TYPE(pde), INTENT(INOUT) :: eqn
        REAL(psb_dpk_), ALLOCATABLE  :: d(:)
    END SUBROUTINE get_pde_diag
  END INTERFACE get_diag


  INTERFACE update_diag
    MODULE SUBROUTINE update_pde_diag(eqn)
        IMPLICIT NONE
        TYPE(pde), INTENT(INOUT) :: eqn
    END SUBROUTINE update_pde_diag
  END INTERFACE update_diag

  ! ----- Temporary up to Gfortran patch -----
  INTERFACE
    MODULE SUBROUTINE get_pde_msh_sub(eqn,msh)
        IMPLICIT NONE
        CLASS(pde), INTENT(IN) :: eqn
        TYPE(mesh), POINTER :: msh
    END SUBROUTINE get_pde_msh_sub
  END INTERFACE
  ! ------------------------------------------

  ! ----- Linear System Solving -----
  INTERFACE spins_pde
    !! Linear System Solving
    MODULE SUBROUTINE spins_pde(n,ia,ja,cloud,eqn)
        !! Inserts a ``cloud'' of coefficients into eqn%A
        IMPLICIT NONE
        INTEGER,          INTENT(IN)    :: n
        INTEGER,          INTENT(IN)    :: ia(:), ja(:)
        REAL(psb_dpk_), INTENT(IN)    :: cloud(:)
        TYPE(pde), INTENT(INOUT) :: eqn
    END SUBROUTINE spins_pde
  END INTERFACE spins_pde

  INTERFACE 

    MODULE SUBROUTINE asb_pde(eqn)
        IMPLICIT NONE
        CLASS(pde), INTENT(INOUT) :: eqn
    END SUBROUTINE asb_pde

  END INTERFACE

  INTERFACE reinit_pde
    MODULE SUBROUTINE reinit_pde(eqn)
        IMPLICIT NONE
        TYPE(pde), INTENT(INOUT) :: eqn
    END SUBROUTINE reinit_pde
  END INTERFACE reinit_pde

  !! ----- Output -----

  INTERFACE write_pde
    !! Output
    MODULE SUBROUTINE write_pde(eqn,mat,mtx_rhs)
        USE tools_output_basics
        IMPLICIT NONE
        TYPE(pde),        INTENT(IN) :: eqn
        CHARACTER(len=*), INTENT(IN) :: mat
        LOGICAL,          INTENT(OUT) :: mtx_rhs
    END SUBROUTINE write_pde
  END INTERFACE write_pde



  INTERFACE

    ! ----- Status Inquirer -----

    MODULE FUNCTION is_pde_bld(eqn)
        IMPLICIT NONE
        LOGICAL :: is_pde_bld
        TYPE(pde), INTENT(IN) :: eqn
    END FUNCTION is_pde_bld

    MODULE FUNCTION is_pde_asb(eqn)
        IMPLICIT NONE
        LOGICAL :: is_pde_asb
        TYPE(pde), INTENT(IN) :: eqn
    END FUNCTION is_pde_asb

    MODULE SUBROUTINE build_pde_prec(eqn)
        IMPLICIT NONE
        TYPE(pde), INTENT(INOUT) :: eqn
    END SUBROUTINE build_pde_prec


    MODULE SUBROUTINE free_pde_prec(eqn)
        IMPLICIT NONE
        TYPE(pde), INTENT(INOUT) :: eqn
    END SUBROUTINE free_pde_prec


    MODULE SUBROUTINE solve_pde_sys(eqn,b,x,iter,err)
        IMPLICIT NONE
        TYPE(pde),      INTENT(INOUT) :: eqn
        REAL(psb_dpk_), INTENT(IN)  :: b(:)
        REAL(psb_dpk_), INTENT(OUT) :: x(:)
        INTEGER,        INTENT(OUT) :: iter
        REAL(psb_dpk_), INTENT(OUT) :: err
    END SUBROUTINE solve_pde_sys

  END INTERFACE

END MODULE class_pde
