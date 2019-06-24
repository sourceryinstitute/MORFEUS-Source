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
! $Id: tools_mesh.f90 8157 2014-10-09 13:02:44Z sfilippo $
!
! Description:
!    Interfaces for import mesh data subroutines
!
! Provides:
!     interface for rd_inp_mesh             reading mesh-related inputs
!     interface for rd_gambit_mesh          reads gambit mesh file
!     interface for rd_cgns_mesh            reads cgns mesh file
!
MODULE tools_mesh
    USE class_psblas, ONLY : psb_dpk_, psb_desc_type
!    USE tools_mesh_basics
    USE class_connectivity, ONLY : connectivity
    IMPLICIT NONE

    PUBLIC

    INTERFACE

        ! ----- Reading: Input Parameters -----
        MODULE SUBROUTINE rd_inp_mesh(input_file,sec,&
            &mesh_file,scale,irenum,ipart,nswpref,mtx_pat)
            IMPLICIT NONE
            INTEGER, PARAMETER :: nlen = 80
            CHARACTER(len=*), INTENT(IN) :: input_file
            CHARACTER(len=*), INTENT(IN) :: sec
            CHARACTER(len=nlen), INTENT(INOUT) :: mesh_file
            REAL(psb_dpk_), INTENT(INOUT) :: scale
            INTEGER, INTENT(INOUT) :: irenum
            INTEGER, INTENT(INOUT) :: ipart
            INTEGER, INTENT(INOUT) :: nswpref
            LOGICAL, INTENT(INOUT) :: mtx_pat
        END SUBROUTINE rd_inp_mesh

    ! ----- Reading: Import mesh Gambit -----

        MODULE SUBROUTINE rd_gambit_mesh(mesh_file, mesh_id, nbc, ncd, &
            & verts, faces, cells, v2f, v2c, f2c, c2g)
            USE class_cell, ONLY : cell
            USE class_face, ONLY : face
            USE class_vertex, ONLY : vertex
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN) :: mesh_file
            CHARACTER(len=*), INTENT(OUT) :: mesh_id
            INTEGER, INTENT(OUT) :: nbc
            INTEGER, INTENT(OUT) :: ncd
            TYPE(vertex), ALLOCATABLE :: verts(:)
            TYPE(face), ALLOCATABLE :: faces(:)
            TYPE(cell), ALLOCATABLE :: cells(:)
            TYPE(connectivity), INTENT(OUT) :: v2f, v2c, f2c, c2g
        END SUBROUTINE rd_gambit_mesh

    ! ----- Reading: Import mesh CGNS -----

        MODULE SUBROUTINE rd_cgns_mesh(mesh_file, mesh_id, nbc, ncd, &
            & verts, faces, cells, v2f, v2c, f2c, c2g)
            USE class_cell, ONLY : cell
            USE class_face, ONLY : face
            USE class_vertex, ONLY : vertex
            IMPLICIT NONE
            CHARACTER(len=*), INTENT(IN) :: mesh_file
            CHARACTER(len=*), INTENT(OUT) :: mesh_id
            INTEGER, INTENT(OUT) :: nbc
            INTEGER, INTENT(OUT) :: ncd
            TYPE(vertex), INTENT(OUT), ALLOCATABLE :: verts(:)
            TYPE(face), INTENT(OUT), ALLOCATABLE :: faces(:)
            TYPE(cell), INTENT(OUT), ALLOCATABLE :: cells(:)
            TYPE(connectivity), INTENT(OUT) :: v2f, v2c, f2c, c2g
        END SUBROUTINE rd_cgns_mesh

    ! ----- Computational Routines -----

        MODULE SUBROUTINE supplement_v2c(v2c,desc_v,ov2c_suppl,c2ov_suppl)
            USE class_keytable, ONLY : keytable
            USE psb_base_mod, ONLY : psb_desc_type
            IMPLICIT NONE
            TYPE(connectivity), INTENT(IN) :: v2c
            TYPE(psb_desc_type), INTENT(IN) :: desc_v
            TYPE(keytable), INTENT(OUT) ::  ov2c_suppl, c2ov_suppl
        END SUBROUTINE supplement_v2c

        MODULE SUBROUTINE supplement_v2f(v2f, faces, desc_v,ov2f_suppl,f2ov_suppl)
            USE class_face, ONLY : face
            USE class_keytable, ONLY : keytable
            USE psb_base_mod, ONLY : psb_desc_type
            IMPLICIT NONE
            TYPE(connectivity), INTENT(INOUT)  :: v2f
            TYPE(psb_desc_type), INTENT(IN) :: desc_v
            TYPE(face)                      :: faces(:)
            TYPE(keytable), INTENT(OUT)     ::  ov2f_suppl, f2ov_suppl
        END SUBROUTINE supplement_v2f

        MODULE SUBROUTINE cmp_mesh_v2v(nverts,v2c,v2v)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nverts
            TYPE(connectivity), INTENT(IN) :: v2c
            TYPE(connectivity), INTENT(OUT) :: v2v
        END SUBROUTINE cmp_mesh_v2v

        MODULE SUBROUTINE cmp_mesh_v2ve(ncd,v2f,v2v)
            IMPLICIT NONE
            INTEGER,            INTENT(IN) :: ncd
            TYPE(connectivity), INTENT(IN) :: v2f
            TYPE(connectivity), INTENT(OUT) :: v2v
        END SUBROUTINE cmp_mesh_v2ve

        MODULE SUBROUTINE cmp_mesh_f2f(nfaces,f2c,f2f)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: nfaces
            TYPE(connectivity), INTENT(IN) :: f2c
            TYPE(connectivity), INTENT(OUT) :: f2f
        END SUBROUTINE cmp_mesh_f2f

        MODULE SUBROUTINE cmp_mesh_c2c(faces,f2c,c2c)
            USE class_face, ONLY : face
            IMPLICIT NONE
            TYPE(face), INTENT(IN) :: faces(:)
            TYPE(connectivity), INTENT(IN) :: f2c
            TYPE(connectivity), INTENT(OUT) :: c2c
        END SUBROUTINE cmp_mesh_c2c

        MODULE SUBROUTINE cmp_mesh_f2b(faces,nbc,f2b)
            USE class_face, ONLY : face
            IMPLICIT NONE
            TYPE(face), INTENT(IN) :: faces(:)
            INTEGER, INTENT(IN) :: nbc
            TYPE(connectivity), INTENT(OUT) :: f2b
        END SUBROUTINE cmp_mesh_f2b

        MODULE SUBROUTINE cmp_mesh_v2b(v2f,faces,nbc,v2b)
            USE class_face
            IMPLICIT NONE
            TYPE(connectivity), INTENT(IN) :: v2f
            TYPE(face), INTENT(IN) :: faces(:)
            INTEGER, INTENT(IN) :: nbc
            TYPE(connectivity), INTENT(OUT) :: v2b
        END SUBROUTINE cmp_mesh_v2b

        MODULE SUBROUTINE cmp_mesh_v2e(ncd,v2f,v2e)
            !USE class_connectivity, ONLY : connectivity
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: ncd
            TYPE(connectivity), INTENT(IN) :: v2f
            TYPE(connectivity), INTENT(OUT) :: v2e
        END SUBROUTINE cmp_mesh_v2e

        MODULE SUBROUTINE cmp_mesh_renum(irenum,cells,faces,c2c,f2c,v2c,c2g)
            USE class_cell, ONLY : cell
            !USE class_connectivity, ONLY : connectivity
            USE class_face, ONLY : face
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: irenum
            TYPE(cell), INTENT(INOUT) :: cells(:)
            TYPE(face), INTENT(INOUT) :: faces(:)
            TYPE(connectivity), INTENT(INOUT) :: c2c, f2c, v2c, c2g
        END SUBROUTINE cmp_mesh_renum

        MODULE SUBROUTINE cmp_mesh_part(ipart,nswpref,c2c)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: ipart, nswpref
            TYPE(connectivity), INTENT(IN) :: c2c
        END SUBROUTINE cmp_mesh_part

        MODULE SUBROUTINE cmp_mesh_desc(v2v,v2c,f2f,f2c,c2c,desc_v,desc_f,desc_c)
            IMPLICIT NONE
            !TYPE(connectivity), INTENT(INOUT) :: v2c!, v2v, f2f, f2c, c2c
            !! Note: IP 5/28/2019 - Had to change the variables below to INTENT(INOUT) rather than (IN) due to the
            !!                      fact that they get sent into procedures with the INTENT(INOUT) or INTENT(OUT).
            TYPE(connectivity), INTENT(INOUT) :: v2c, v2v, f2f, f2c, c2c
            TYPE(psb_desc_type), INTENT(OUT) :: desc_v, desc_f, desc_c
        END SUBROUTINE cmp_mesh_desc

        MODULE SUBROUTINE cmp_moving_surf(nbc,v2b,verts,surf)
            USE class_surface, ONLY : surface
            USE class_vertex, ONLY : vertex
            IMPLICIT NONE
            ! Parameters
            INTEGER, INTENT(IN)            :: nbc
            TYPE(connectivity), INTENT(IN) :: v2b
            TYPE(vertex),INTENT(IN)        :: verts(:)
            TYPE(surface),ALLOCATABLE, INTENT(OUT) :: surf(:)
        END SUBROUTINE cmp_moving_surf

    END INTERFACE

END MODULE tools_mesh
