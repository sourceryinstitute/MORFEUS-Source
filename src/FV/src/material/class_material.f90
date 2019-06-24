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
! $Id: class_material.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    To be added...
!
MODULE class_material

    USE class_psblas
    USE MatLib, ONLY : MatProp, Set_Material_ids
    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: matptr                       ! Array of class_material
    PUBLIC :: material                       ! Class
    PUBLIC :: matlaw, check_temp, &          ! Other
        &    check_material_consistency

    ! String length constant
    INTEGER, PARAMETER :: name_len = 15

    TYPE material
        PRIVATE
        !
        ! From input file ...
        CHARACTER(len=name_len) :: name
        CHARACTER(len=name_len) :: mat_type
        INTEGER :: ilaw(4)
        INTEGER :: mat_id
        !
        ! ... and from DB
        CHARACTER(len=1) :: state
        REAL(psb_dpk_) :: dtemp, tmin, tmax
        REAL(psb_dpk_), ALLOCATABLE :: rho(:)
        REAL(psb_dpk_), ALLOCATABLE :: mu(:)
        REAL(psb_dpk_), ALLOCATABLE :: lambda(:)
        REAL(psb_dpk_), ALLOCATABLE :: sh(:)
    CONTAINS
        PROCEDURE :: create_material, free_material ! Constructor/Destructor
        PROCEDURE, PRIVATE :: get_material_name, get_material_id  ! Getter
        GENERIC, PUBLIC :: name_ => get_material_name
        GENERIC, PUBLIC :: mat_id_ => get_material_id
        PROCEDURE, PRIVATE :: nemo_material_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_material_sizeof
    END TYPE material

    TYPE matptr
        TYPE(material) :: mat
    END TYPE matptr


  ! ----- Generic Interfaces -----

  ! ----- Getters -----
  INTERFACE

    MODULE FUNCTION get_material_name(mat)
      !! Getter
        CHARACTER(len=name_len) :: get_material_name
        CLASS(material), INTENT(IN) :: mat
    END FUNCTION get_material_name

    MODULE FUNCTION get_material_id(mat)
      !! Getter
        INTEGER :: get_material_id
        CLASS(material), INTENT(IN) :: mat
    END FUNCTION get_material_id
  END INTERFACE


  INTERFACE matlaw
    !! Physical Properties Laws

    MODULE SUBROUTINE matlaw_v(mats,im,t,dim,f)
        USE class_dimensions
        USE tools_material
        USE tools_math
        TYPE(dimensions), INTENT(IN) :: dim
        INTEGER, INTENT(IN) :: im(:)
        CLASS(matptr),   INTENT(IN), OPTIONAL, TARGET :: mats(:)
        REAL(psb_dpk_), INTENT(IN) :: t(:)
        REAL(psb_dpk_), INTENT(OUT) :: f(:)
    END SUBROUTINE matlaw_v


    MODULE SUBROUTINE matlaw_s(mats,im,t,dim,f)
        USE class_dimensions
        USE tools_math
        TYPE(dimensions), INTENT(IN) :: dim
        CLASS(matptr),   INTENT(IN), OPTIONAL, TARGET :: mats(:)
        INTEGER, INTENT(IN) :: im
        REAL(psb_dpk_), INTENT(IN) :: t
        REAL(psb_dpk_), INTENT(OUT) :: f
    END SUBROUTINE matlaw_s

    MODULE SUBROUTINE matlaw_fast_s(mat, t, property, f)
        CLASS(material), INTENT(IN) :: mat
        REAL(psb_dpk_), INTENT(IN) :: t
        REAL(psb_dpk_), INTENT(OUT) :: f
        CHARACTER(LEN=*), INTENT(IN) :: property
    END SUBROUTINE matlaw_fast_s

  END INTERFACE matlaw

  ! Check Procedures
  INTERFACE check_temp

    MODULE SUBROUTINE check_temp_s(t,mat)
        REAL(psb_dpk_), INTENT(IN) :: t
        TYPE(material), INTENT(IN) :: mat
    END SUBROUTINE check_temp_s

    MODULE SUBROUTINE check_temp_v(t,mat)
        REAL(psb_dpk_), INTENT(IN) :: t(:)
        TYPE(material), INTENT(IN) :: mat
    END SUBROUTINE check_temp_v

  END INTERFACE check_temp

    LOGICAL, PARAMETER :: debug = .TRUE.
    INTEGER, PARAMETER :: nlen = 80

  INTERFACE
    MODULE FUNCTION nemo_material_sizeof(mat)
        USE psb_base_mod
        CLASS(material), INTENT(IN) :: mat
        INTEGER(kind=nemo_int_long_)   :: nemo_material_sizeof
    END FUNCTION nemo_material_sizeof

    ! ----- ILAW's Reference table -----
    ! ILAW(*)    = 1 -> constant
    ! ILAW(*)    = 2 -> piecewise linear
    ! ILAW(irho) = 3 -> perfect gas rho = f(T)
    ! ILAW(irho) = 4 -> perfect gas rho = f(p,T)

    ! ----- Constructors -----

    MODULE SUBROUTINE create_material(mat,input_file,sec)
      !! Global Constructor
        USE tools_material
        CLASS(material), INTENT(OUT) :: mat
        CHARACTER(len=*), INTENT(IN) :: input_file
        CHARACTER(len=*), INTENT(IN) :: sec
    END SUBROUTINE create_material


    ! ----- Destructor -----

    MODULE SUBROUTINE free_material(mat)
        CLASS(material), INTENT(INOUT) :: mat
    END SUBROUTINE free_material

    MODULE SUBROUTINE check_material_consistency(mat1,mat2,WHERE)
        TYPE(matptr), POINTER :: mat1(:), mat2(:)
        CHARACTER(len=*), INTENT(IN) :: WHERE
    END SUBROUTINE check_material_consistency

    ! ----- Debug -----

    MODULE SUBROUTINE debug_material(mat)
        USE tools_material
        TYPE(material), INTENT(IN) :: mat
    END SUBROUTINE debug_material

  END INTERFACE

END MODULE class_material
