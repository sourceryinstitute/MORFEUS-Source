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
! $Id: class_mesh.F90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Mesh class
!
! Provides:
!    MESH                   class with public access.
!    CREATE_MESH            constructor, using file I/O. Also performs the domain
!                           partitioning and cell renumbering to reduce bandwidth.
!    IMPORT_MESH            opens file and reads appropriate format & broadcasts.
!    FREE_MESH              destructor.
!    G2L_MESH               global to local reallocation of mesh structures.
!    CHECK_MESH_UNUSED_EL   checks for unused vertices, faces, & cells in C2C, V2C
!                           F2C, V2F connectivity members.
!    CHECK_MESH_CONSISTENCY compares two meshes and checks their consistency.
!
MODULE class_mesh
    !! Define and manipulate data describing the discretization of space into connected finite-volume cells and surfaces
    USE psb_base_mod
    USE class_psblas, ONLY : nemo_int_long_
    USE class_cell, ONLY : cell
    USE class_connectivity, ONLY : connectivity
    USE class_face, ONLY : face
    USE class_least_squares, ONLY : least_squares
    USE class_vector, ONLY : vector
    USE class_keytable, ONLY : keytable
    USE class_surface, ONLY : surface
    USE class_vertex, ONLY : vertex
    USE grid_interface, ONLY : grid

    IMPLICIT NONE

    PRIVATE ! Default
    PUBLIC :: mesh                      ! Class
    PUBLIC :: check_mesh_consistency

    INTEGER, PARAMETER :: nlen = 80

    TYPE, EXTENDS(grid) :: mesh
      !! Encapsulate mesh connectivity, component parts, metrics, surfaces, and linear algebraic descriptors
        LOGICAL :: set = .FALSE.   ! Indicates if the mesh been created yet
        CHARACTER(len=nlen) :: id  ! Mesh ID
        INTEGER :: nbc             ! Number of BCs
        INTEGER :: ngp             ! Number of element groups
        INTEGER :: ncd             ! Dimensionality of the mesh (2d or 3d)

        ! Connectivity data
        TYPE(connectivity) :: v2c, v2f, f2c, c2g         ! From mesh importing
        TYPE(connectivity) :: v2v, f2f, c2c              ! Adjacency graphs
        TYPE(keytable)     :: ov2c_sup,c2ov_sup          ! Supplemental v2c/c2v connectivity
        TYPE(keytable)     :: ov2f_sup,f2ov_sup          ! Supplemental v2f/f2v connectivity
        TYPE(connectivity) :: v2b, f2b                   ! Boundary related

        ! Mesh sub-elements
        TYPE(vertex), ALLOCATABLE :: verts(:)            ! Vertex coordinates
        TYPE(face), ALLOCATABLE :: faces(:)              ! Faces description
        TYPE(cell), ALLOCATABLE :: cells(:)              ! Cells description

        ! PSBLAS descriptors
        TYPE(psb_desc_type) :: desc_v     ! Vertices
        TYPE(psb_desc_type) :: desc_f     ! Faces
        TYPE(psb_desc_type) :: desc_c     ! Cells

        ! Face-related metrics (see below for details)
        REAL(psb_dpk_), ALLOCATABLE :: area(:)
        REAL(psb_dpk_), ALLOCATABLE :: dist(:)
        REAL(psb_dpk_), ALLOCATABLE :: interp(:)
        TYPE(vector), ALLOCATABLE :: face_cntr(:)
        TYPE(vector), ALLOCATABLE :: af(:)
        TYPE(vector), ALLOCATABLE :: df(:)

        ! Cell-related metrics (see below for details)
        REAL(psb_dpk_), ALLOCATABLE :: vol(:)
        TYPE(vector), ALLOCATABLE :: cell_cntr(:)

        ! Metrics for cell-centered Least Squares Regression
        TYPE(least_squares), ALLOCATABLE :: lsr(:)

        ! Surface geometry and location
        TYPE(surface), ALLOCATABLE :: surf(:)
    CONTAINS
        PROCEDURE :: create_mesh          !! Constructor
        PROCEDURE :: free_mesh            !! Destructor
        PROCEDURE :: check_mesh_unused_el !! Check routines
        PROCEDURE, PRIVATE :: nemo_mesh_sizeof
        GENERIC, PUBLIC :: nemo_sizeof => nemo_mesh_sizeof
        !        PROCEDURE :: check_mesh_consistency
    END TYPE mesh

    ! Metrics description:
    ! AREA:      face area
    ! DIST:      center-to-center distance across the face
    ! INTERP:    interpolation weighting from master to face
    ! FACE_CNTR: face centroid
    ! AF:        normal area vector
    ! DF:        distance vector from master to slave
    ! VOL:       cell volume
    ! CELL_CNTR: cell centroid

    ! IMPORTANT! Class attributes are left PUBLIC because of:
    ! - faster access
    ! - no need to encapsulate the implementation

    INTERFACE

        !! ----- Constructors -----

        MODULE SUBROUTINE create_mesh(msh,input_file,sec)
            !! Global constructor
            IMPLICIT NONE
            CLASS(mesh),      INTENT(OUT) :: msh
            CHARACTER(len=*), INTENT(IN)  :: input_file
            CHARACTER(len=*), INTENT(IN)  :: sec
        END SUBROUTINE

        MODULE SUBROUTINE free_mesh(msh)
            !! ----- Destructor -----
            IMPLICIT NONE
            CLASS(mesh), INTENT(INOUT) :: msh
        END SUBROUTINE

        !! ----- Check Operations -----

        MODULE SUBROUTINE check_mesh_unused_el(msh)
            !! Scans through numerous connectivities, ensuring that in each one,
            !! all faces & cells are referenced at least once.
            IMPLICIT NONE
            CLASS(mesh), INTENT(IN) :: msh
        END SUBROUTINE check_mesh_unused_el

        MODULE FUNCTION nemo_mesh_sizeof(msh)
            IMPLICIT NONE
            CLASS(mesh), INTENT(IN) :: msh
            INTEGER(kind=nemo_int_long_)   :: nemo_mesh_sizeof
        END FUNCTION nemo_mesh_sizeof

        MODULE SUBROUTINE check_mesh_consistency(msh1,msh2,WHERE)
            !! Checks the consistency of two meshes: MSH1 and MSH2
            IMPLICIT NONE
            TYPE(mesh), POINTER :: msh1, msh2
            CHARACTER(len=*), INTENT(IN) :: WHERE
        END SUBROUTINE check_mesh_consistency

    END INTERFACE

END MODULE class_mesh
