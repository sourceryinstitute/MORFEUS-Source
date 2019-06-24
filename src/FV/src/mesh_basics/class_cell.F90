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
! $Id: class_cell.F90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    Provides cell class functionality.
!
! Includes:
!    CELL class              (all private)
!    CELL_                   basic elemental constructor
!    ALLOC_CELL              construct array pointer for numerous cells
!    FREE_CELL               free array pointers
!    BCAST_CELL              broadcast cells from node 0 to other processes
!    G2L_CELL                global to local reallocation of cells
!    GET_CELL_NV             get number of vertices of cell (elemental)
!    GET_CELL_GEO            get character code indicating type of cell
!    GET_CELLS_TYPE          concatenate lists of cells by type
!

MODULE class_cell

    USE class_psblas

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: cell                           !! Class
    PUBLIC :: cell_, alloc_cell, free_cell   !! Constructor/destructor
    PUBLIC :: bcast_cell, g2l_cell, l2g_cell !! Parallel ops.
    PUBLIC :: get_cells_type      ! Getters
    PUBLIC :: itri_, iqua_, &                !! Named constants
        &    itet_, ipyr_, ipri_, ihex_

    INTEGER, PARAMETER :: nlen = 3


    TYPE cell
        PRIVATE
        INTEGER :: nv  !! number of vertices
        INTEGER :: nf  !! number of faces
        INTEGER :: group !! group ID of cell
        CHARACTER(len=nlen) :: geo  !! abbreviation for kind of cell
    CONTAINS
        PROCEDURE, PRIVATE :: get_cell_nv, get_cell_geo, get_cell_group  !Getters
        GENERIC, PUBLIC :: nv_ => get_cell_nv
        GENERIC, PUBLIC :: geo_ => get_cell_geo
        GENERIC, PUBLIC :: group_ => get_cell_group
        PROCEDURE, PRIVATE :: nemo_cell_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_cell_sizeof
    END TYPE cell

    ! ----- Named Constants -----

    INTEGER, PARAMETER :: itri_ = 1, iqua_ = 2  !! triangle, quadrilateral
    INTEGER, PARAMETER :: itet_ = 3, ipyr_ = 4  !! tetrahedral, pyramid
    INTEGER, PARAMETER :: ipri_ = 5, ihex_ = 6  !! prism, hexahedron

  INTERFACE

    ELEMENTAL MODULE FUNCTION nemo_cell_sizeof(cll)
        USE psb_base_mod
        IMPLICIT NONE
        CLASS(cell), INTENT(IN) :: cll
        INTEGER(kind=nemo_int_long_)   :: nemo_cell_sizeof
    END FUNCTION nemo_cell_sizeof
      ! ----- Constructors -----
  
    ELEMENTAL MODULE FUNCTION cell_(nv,nf,group,geo)
        IMPLICIT NONE
        TYPE(cell) :: cell_
        INTEGER, INTENT(IN) :: nv, nf, group
        CHARACTER(len=nlen), INTENT(IN) :: geo
    END FUNCTION cell_


    MODULE SUBROUTINE alloc_cell(cells,n)
      !! Array constructor
        IMPLICIT NONE
        TYPE(cell), ALLOCATABLE :: cells(:)
        INTEGER, INTENT(IN)     :: n
    END SUBROUTINE alloc_cell


    ! ----- Destructor -----

    MODULE SUBROUTINE free_cell(cells)
        IMPLICIT NONE
        TYPE(cell), ALLOCATABLE  :: cells(:)
    END SUBROUTINE free_cell


    ! ----- Parallel Ops. -----

    MODULE SUBROUTINE bcast_cell(cells)
        IMPLICIT NONE
        TYPE(cell), ALLOCATABLE :: cells(:)
    END SUBROUTINE bcast_cell

    MODULE SUBROUTINE g2l_cell(cells,desc_c)
        USE psb_base_mod
        IMPLICIT NONE
        TYPE(cell), ALLOCATABLE :: cells(:)
        TYPE(psb_desc_type), INTENT(IN) :: desc_c
    END SUBROUTINE g2l_cell

    MODULE SUBROUTINE l2g_cell(cells_loc,cells_glob,desc_c)
        !! WARNING! The global results is allocated only on P0. After its usage
        !! it must be deallocated in the calling unit by means of the statement:
        !! "if(associated(glob_res)) deallocate(glob_res)"
        USE psb_base_mod
        IMPLICIT NONE
        TYPE(cell), ALLOCATABLE  :: cells_loc(:)
        TYPE(cell), ALLOCATABLE  :: cells_glob(:)
        TYPE(psb_desc_type), INTENT(IN) :: desc_c
    END SUBROUTINE l2g_cell

  ! ----- Generic Interfaces -----

  ! ----- Getters -----
    ELEMENTAL MODULE FUNCTION get_cell_nv(c)
        IMPLICIT NONE
        INTEGER :: get_cell_nv
        CLASS(cell), INTENT(IN) :: c
    END FUNCTION get_cell_nv

    MODULE FUNCTION get_cell_geo(c)
        IMPLICIT NONE
        CHARACTER(len=nlen) :: get_cell_geo
        CLASS(cell), INTENT(IN) :: c
    END FUNCTION get_cell_geo

    MODULE FUNCTION get_cell_group(c)
        IMPLICIT NONE
        INTEGER :: get_cell_group
        CLASS(cell), INTENT(IN) :: c
    END FUNCTION get_cell_group

    MODULE SUBROUTINE get_cells_type(cells,nctype,ictype,desc)
        USE psb_base_mod
        IMPLICIT NONE
        TYPE(cell), INTENT(IN) :: cells(:)     !!array of cells structs
        INTEGER, ALLOCATABLE, INTENT(OUT) :: nctype(:)   !!count of each type
        INTEGER, ALLOCATABLE, INTENT(OUT)  :: ictype(:)  !!array of cell id's sorted by type
        TYPE(psb_desc_type), INTENT(IN), OPTIONAL :: desc
    END SUBROUTINE get_cells_type

  END INTERFACE

END MODULE class_cell
