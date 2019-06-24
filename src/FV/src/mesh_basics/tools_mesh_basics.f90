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
! $Id: tools_mesh_basics.f90 2469 2007-10-08 10:34:43Z sfilippo $
!
! Description:
!    Provides interfaces to geometric routines
!
! Provides:
!    geom_face   calculates face normal, area, and centroid
!    geom_cell   calculates volume and center of mass of cells by dividing into tets
!    geom_diff   find cell to cell distances, face interpolation factors
!    geom_tet    calculates volume and center of mass of tetrahedra

MODULE tools_mesh_basics

    IMPLICIT NONE

    ! ----- Geometry -----

    INTERFACE
        MODULE SUBROUTINE geom_face(verts,v2f,ncd, &
            & face_cntr,af,area)
            USE class_connectivity
            USE class_vector
            USE class_vertex
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            TYPE(vertex),       INTENT(IN), ALLOCATABLE  :: verts(:)
            TYPE(connectivity), INTENT(IN) :: v2f
            INTEGER,            INTENT(IN) :: ncd
            TYPE(vector),       INTENT(OUT), ALLOCATABLE :: face_cntr(:)
            TYPE(vector),       INTENT(OUT), ALLOCATABLE :: af(:)
            REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: area(:)
        END SUBROUTINE geom_face
    END INTERFACE

    INTERFACE
        MODULE SUBROUTINE geom_cell(verts,faces,cells,v2f,v2c,f2c,ncd, &
            & cell_cntr,vol,quiet)
            USE class_cell
            USE class_connectivity
            USE class_face
            USE class_vector
            USE class_vertex
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            TYPE(vertex),       INTENT(IN), ALLOCATABLE :: verts(:)
            TYPE(face),         INTENT(IN), ALLOCATABLE :: faces(:)
            TYPE(cell),         INTENT(IN), ALLOCATABLE :: cells(:)
            TYPE(connectivity), INTENT(IN) :: v2f, v2c, f2c
            INTEGER,            INTENT(IN) :: ncd
            TYPE(vector),       INTENT(OUT), ALLOCATABLE :: cell_cntr(:)
            REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: vol(:)
            LOGICAL,            INTENT(IN),  OPTIONAL    :: quiet

        END SUBROUTINE geom_cell
    END INTERFACE

    INTERFACE
        MODULE SUBROUTINE geom_diff(faces,f2b,face_cntr,af,cell_cntr, &
            & df,dist,int_fact)
            USE class_connectivity
            USE class_face
            USE class_vector
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            TYPE(face),         INTENT(IN), ALLOCATABLE  :: faces(:)
            TYPE(connectivity), INTENT(IN) :: f2b
            TYPE(vector),       INTENT(IN), ALLOCATABLE :: face_cntr(:)
            TYPE(vector),       INTENT(IN), ALLOCATABLE :: af(:)
            TYPE(vector),       INTENT(IN), ALLOCATABLE :: cell_cntr(:)
            TYPE(vector),       INTENT(OUT), ALLOCATABLE :: df(:)
            REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: dist(:)
            REAL(psb_dpk_),   INTENT(OUT), ALLOCATABLE :: int_fact(:)
        END SUBROUTINE geom_diff
    END INTERFACE

    INTERFACE
        MODULE FUNCTION geom_tet_center(v1,v2,v3,v4)
            USE class_vector
            USE class_vertex
            IMPLICIT NONE
            TYPE(vector) :: geom_tet_center
            TYPE(vertex), INTENT(IN) :: v1, v2, v3, v4
        END FUNCTION geom_tet_center
    END INTERFACE

    INTERFACE
        MODULE SUBROUTINE geom_tet_dihedral_angle(af,largest,smallest)
            USE class_vector
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            TYPE(vector), INTENT(IN) :: af(4)
            REAL(psb_dpk_), INTENT(OUT) :: largest, smallest
        END SUBROUTINE geom_tet_dihedral_angle
    END INTERFACE

    INTERFACE
        MODULE SUBROUTINE geom_hex_dihedral_angle(af,adjacent,largest,smallest)
            USE class_vector
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            TYPE(vector), INTENT(IN) :: af(6)
            INTEGER, INTENT(IN)      :: adjacent(12,2)
            REAL(psb_dpk_), INTENT(OUT) :: largest, smallest
        END SUBROUTINE  geom_hex_dihedral_angle
    END INTERFACE

    INTERFACE
        MODULE FUNCTION geom_tet_quality(v1,v2,v3,v4,vol)
            USE class_vertex
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            REAL(psb_dpk_) :: geom_tet_quality
            TYPE(vertex), INTENT(IN) :: v1, v2, v3, v4
            REAL(psb_dpk_), INTENT(IN) :: vol
        END FUNCTION geom_tet_quality
    END INTERFACE

    INTERFACE
        MODULE FUNCTION geom_hex_quality(verts, vol)
            USE class_vertex
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            REAL(psb_dpk_) :: geom_hex_quality
            TYPE(vertex), INTENT(IN) :: verts(8)
            REAL(psb_dpk_), INTENT(IN) :: vol
        END FUNCTION geom_hex_quality
    END INTERFACE

    INTERFACE
        MODULE FUNCTION geom_tet_volume(v1,v2,v3,v4)
            USE class_vertex
            USE class_psblas, ONLY : psb_dpk_
            IMPLICIT NONE
            REAL(psb_dpk_) :: geom_tet_volume
            TYPE(vertex), INTENT(IN) :: v1, v2, v3, v4
        END FUNCTION geom_tet_volume
    END INTERFACE


    ! ----- Named Constants -----

    INTEGER, PARAMETER :: iunknown_  = 0 ! unknown
    INTEGER, PARAMETER :: iplane_    = 1 ! plane
    INTEGER, PARAMETER :: icylinder_ = 2 ! cylinder
    INTEGER, PARAMETER :: isphere_   = 3 ! sphere

END MODULE tools_mesh_basics
