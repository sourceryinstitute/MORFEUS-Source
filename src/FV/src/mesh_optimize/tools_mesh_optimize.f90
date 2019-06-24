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
! $Id: tools_mesh_optimize.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    to be added...
!
MODULE tools_mesh_optimize
    USE class_psblas
    IMPLICIT NONE

    INTERFACE

        MODULE SUBROUTINE smooth_mesh(msh,bc,surface_iter,interior_iter)
            USE class_iterating, ONLY : iterating
            USE class_mesh, ONLY : mesh
            USE class_bc, ONLY : bc_poly
            IMPLICIT NONE
            TYPE(mesh), INTENT(INOUT) :: msh
            TYPE(bc_poly), INTENT(IN) :: bc(:)
            TYPE(iterating), INTENT(INOUT):: surface_iter
            TYPE(iterating), INTENT(INOUT):: interior_iter

        END SUBROUTINE smooth_mesh

        MODULE SUBROUTINE mobile_verts(msh,bc,c2v,shared_flag,unconstrained, &
            & n_unconstrained,constrained, n_constrained, all_tets, mixed)
            USE class_bc, ONLY : bc_poly
            USE class_connectivity, ONLY : connectivity
            USE class_mesh, ONLY : mesh
            IMPLICIT NONE
            TYPE(mesh), INTENT(IN)    :: msh
            TYPE(bc_poly), INTENT(IN) :: bc(:)
            TYPE(connectivity), INTENT(IN) ::   c2v
            INTEGER, INTENT(OUT) :: unconstrained(:)
            INTEGER, INTENT(OUT) :: constrained(:)
            LOGICAL, INTENT(IN)  :: shared_flag(:)
            INTEGER, INTENT(OUT) :: n_unconstrained
            INTEGER, INTENT(OUT) :: n_constrained
            LOGICAL, INTENT(INOUT) :: all_tets(:)
            LOGICAL, INTENT(INOUT) :: mixed(:)

        END SUBROUTINE mobile_verts

        MODULE SUBROUTINE smooth_interior_vtx(iv,msh,c2v,shared_flag,all_tets)
            USE class_connectivity, ONLY : connectivity
            USE class_mesh, ONLY : mesh
            IMPLICIT NONE

            INTEGER,INTENT(IN)             :: iv             ! id # of this vertex
            TYPE(mesh), INTENT(INOUT)      :: msh            ! the mesh
            TYPE(connectivity), INTENT(IN) :: c2v            ! given vertices, find the cells
            LOGICAL, INTENT(IN)            :: shared_flag(:) ! flags if a vertex is shared
            LOGICAL, INTENT(IN)            :: all_tets(:)    ! is this vertex used only by tets?

        END SUBROUTINE smooth_interior_vtx

        MODULE SUBROUTINE check_right_handed(msh,shared,shared_flag,tangled,all_tets)
            USE class_mesh, ONLY : mesh
            IMPLICIT NONE
            TYPE(mesh), INTENT(IN) :: msh
            INTEGER, ALLOCATABLE, INTENT(IN) :: shared(:)
            LOGICAL, ALLOCATABLE, INTENT(IN) :: shared_flag(:)
            INTEGER, INTENT(OUT) :: tangled
            LOGICAL, INTENT(IN)    :: all_tets(:)     !true if only tets use this vertex
        END SUBROUTINE check_right_handed

        MODULE SUBROUTINE smooth_surf_vtx(iv,ib,msh,f2v,shared_flag,tangled,all_tets)
            USE class_connectivity, ONLY : connectivity
            USE class_mesh,         ONLY : mesh
            IMPLICIT NONE
            INTEGER,INTENT(IN)             :: iv             ! id # of this vertex
            INTEGER,INTENT(IN)             :: ib             ! id # of the boundary
            TYPE(mesh), INTENT(INOUT)      :: msh            ! the mesh
            TYPE(connectivity), INTENT(IN) :: f2v            ! given vertices, find the cells
            LOGICAL, INTENT(IN)            :: shared_flag(:) ! flags if a vertex is shared
            INTEGER, INTENT(IN)            :: tangled  ! is there a cell in the mesh that has
            LOGICAL, INTENT(IN)            :: all_tets(:)
        END SUBROUTINE smooth_surf_vtx

        MODULE SUBROUTINE laplacian_smooth(desc_v, v2v, n_unconstrained, unconstrained, verts, mixed)
            USE psb_base_mod!,       ONLY : psb_desc_type
            USE psb_desc_mod
            USE class_connectivity, ONLY : connectivity
            USE class_vertex,       ONLY : vertex
            IMPLICIT NONE
            TYPE(psb_desc_type),INTENT(INOUT) :: desc_v      ! Vertices
              !! Note: Had to change to INTENT(INOUT) rather than INTENT(IN) due to the procedures
              !!       being called in the subroutine having the INTENT(OUT) or INTENT(INOUT). IP - 6/6/2019
            TYPE(connectivity), INTENT(INOUT) :: v2v         ! given vertices, find the neighbors
            INTEGER,           INTENT(IN)  :: n_unconstrained  !# of interior, mobile vertices
            INTEGER,           INTENT(IN)  :: unconstrained(:) !interior, mobile vertices
            TYPE(vertex), ALLOCATABLE, INTENT(INOUT) :: verts(:)       ! Vertex coordinates
            LOGICAL, INTENT(IN)            :: mixed(:)
        END SUBROUTINE laplacian_smooth

        SUBROUTINE optimize_vertex_rand(msh,c2v,iv)

            ! right now, only for serial use

            USE class_connectivity, ONLY : connectivity
            USE class_mesh, ONLY : mesh

            IMPLICIT NONE

            ! Variable parameters

            TYPE(mesh), INTENT(INOUT) :: msh      ! mesh structure
            TYPE(connectivity), INTENT(IN) :: c2v ! given a vertex, finds connected cells
            INTEGER, INTENT(IN) :: iv             ! index of the vertex to be moved

        END SUBROUTINE optimize_vertex_rand

    ! ----- Routines for interfacing with OptMS -----

        FUNCTION initoptms(dims,technique,functionID)
            IMPLICIT NONE
            INTEGER :: initoptms
            INTEGER, INTENT(IN) :: dims
            INTEGER, INTENT(IN) :: technique
            INTEGER, INTENT(IN) :: functionID
        END FUNCTION initoptms

        FUNCTION initoptms2d(dims,technique,functionID)
            IMPLICIT NONE
            INTEGER :: initoptms2d
            INTEGER, INTENT(IN) :: dims
            INTEGER, INTENT(IN) :: technique
            INTEGER, INTENT(IN) :: functionID
        END FUNCTION initoptms2d

        FUNCTION freeoptms()
            IMPLICIT NONE
            INTEGER :: freeoptms
        END FUNCTION freeoptms

        FUNCTION freeoptms2d()
            IMPLICIT NONE
            INTEGER :: freeoptms2d
        END FUNCTION freeoptms2d

        FUNCTION right_handed(p1,p2,p3,p4)
            USE class_psblas, ONLY : psb_dpk_
            INTEGER :: right_handed
            REAL(psb_dpk_) :: p1(3)
            REAL(psb_dpk_) :: p2(3)
            REAL(psb_dpk_) :: p3(3)
            REAL(psb_dpk_) :: p4(3)
        END FUNCTION right_handed

        FUNCTION right_handed2d(p1,p2,p3)
            USE class_psblas, ONLY : psb_dpk_
            INTEGER :: right_handed2d
            REAL(psb_dpk_) :: p1(2)
            REAL(psb_dpk_) :: p2(2)
            REAL(psb_dpk_) :: p3(2)
        END FUNCTION right_handed2d

        FUNCTION call_smooth (num_incident_vtx,num_incident_tet,free_pos, &
            & tet_pos,tet_verts,tangled)
            USE class_psblas, ONLY : psb_dpk_
            INTEGER :: call_smooth
            INTEGER :: num_incident_vtx
            INTEGER :: num_incident_tet
            REAL(psb_dpk_) :: free_pos(3)
            REAL(psb_dpk_) :: tet_pos(3,num_incident_vtx)
            INTEGER :: tet_verts(3,num_incident_tet)
            INTEGER :: tangled
        END FUNCTION call_smooth

        ! C function prototype:
        !
        ! int call_smooth2d(int *num_incident_vtx, int *num_incident_tri,
        !                   double *free_pos,double incident_vtx[][2],
        !                    int vtx_connectivity[][2],int *tangled);

        FUNCTION call_smooth2d (num_incident_vtx,num_incident_tri,free_pos, &
            & tri_pos,tri_verts,tangled) BIND(C)
            USE class_psblas, ONLY : psb_dpk_
            USE iso_c_binding, ONLY : c_int, c_ptr, c_double
            INTEGER(c_int) :: call_smooth2d
            TYPE(c_ptr), VALUE :: num_incident_vtx
            TYPE(c_ptr), VALUE :: num_incident_tri
            REAL(c_double) :: free_pos(2)
            REAL(c_double) :: tri_pos(2,*)   ! tri_pos(2,num_incident_vtx)
            INTEGER(c_int) :: tri_verts(2,*) ! tri_verts(2,num_incident_tri)
            TYPE(c_ptr), VALUE :: tangled
        END FUNCTION call_smooth2d

    END INTERFACE

END MODULE tools_mesh_optimize
