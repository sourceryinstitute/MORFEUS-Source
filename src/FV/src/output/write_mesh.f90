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
! $Id: write_mesh.f90 9099 2015-04-24 15:34:00Z sfilippo $
!
SUBMODULE (tools_output) write_mesh_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE write_mesh
            USE class_psblas
            USE class_cell
            USE class_connectivity
            USE class_face
            USE class_iterating
            USE class_mesh
            USE class_output
            USE class_vertex
            USE tools_output_basics

            IMPLICIT NONE
            !
            INTEGER :: info, err_act
            INTEGER :: icontxt, mypnum
            INTEGER :: i, ic, ig, ncells, ngc, ngroups
            INTEGER, POINTER :: ic2g(:) => NULL()
            INTEGER, ALLOCATABLE  :: igroup(:), iproc(:)
            INTEGER, ALLOCATABLE :: i_loc(:)
            CHARACTER(len=32) :: path
            TYPE(cell), ALLOCATABLE :: cells(:)
            TYPE(face), ALLOCATABLE :: faces(:)
            TYPE(vertex), ALLOCATABLE :: verts(:)
            TYPE(connectivity) :: v2f, v2c, f2c

            ! Sets error handling for PSBLAS-2 routines
            CALL psb_erractionsave(err_act)

            mypnum  = mypnum_()
            icontxt = icontxt_()

            ! Sets output path
            IF(PRESENT(iter)) CALL out%set_output_path(iter)
            path = out%path_()

            ! Global number of cells
            ncells = psb_cd_get_global_cols(msh%desc_c)

            CALL psb_geall(i_loc,msh%desc_c,info)
            CALL psb_check_error(info,'write_mesh','psb_geall',icontxt)

            ALLOCATE(igroup(ncells),iproc(ncells),stat=info)
            IF(info /= 0) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            ! Gathers MESH components on P0
            CALL l2g_vertex(msh%verts,verts,msh%desc_v)
            CALL l2g_face(msh%faces,faces,msh%desc_f,msh%desc_c)
            CALL l2g_cell(msh%cells,cells,msh%desc_c)
            CALL l2g_conn(msh%v2f,v2f,msh%desc_v,msh%desc_f)
            CALL l2g_conn(msh%v2c,v2c,msh%desc_v,msh%desc_c)
            CALL l2g_conn(msh%f2c,f2c,msh%desc_f,msh%desc_c)

            ! Gathers processor IDs
            i_loc(:) = mypnum
            CALL psb_gather(iproc,i_loc,msh%desc_c,info,root=0)
            CALL psb_check_error(info,'write_mesh','psb_gather',icontxt)

            ! Gathers group IDs
            i_loc(:) = 0
            ngroups = msh%c2g%nel_()

            DO ig = 1, ngroups
                CALL msh%c2g%get_ith_conn(ic2g,ig)
                ngc = SIZE(ic2g) ! number of group cells
                DO i = 1, ngc
                    ic = ic2g(i)
                    i_loc(ic) = ig
                END DO
            END DO

            CALL psb_gather(igroup,i_loc,msh%desc_c,info,root=0)
            CALL psb_check_error(info,'write_mesh','psb_gather',icontxt)

            IF(mypnum == 0) THEN
                SELECT CASE(out%fmt_())
                CASE(vtk_)
                    !        call wr_vtk_mesh(msh%ncd,verts,cells,v2c,iproc,path)
                CASE default
                    WRITE(*,200)
                END SELECT
            END IF

            ! Frees Memory
            NULLIFY(ic2g)
            DEALLOCATE(igroup,iproc)

            CALL psb_gefree(i_loc,msh%desc_c,info)
            CALL psb_check_error(info,'write_mesh','psb_gefree',icontxt)

            CALL free_conn(f2c)
            CALL free_conn(v2c)
            CALL free_conn(v2f)
            IF(ALLOCATED(cells)) CALL free_cell(cells)
            IF(ALLOCATED(faces)) CALL free_face(faces)
            IF(ALLOCATED(verts)) CALL free_vertex(verts)


            ! ----- Normal termination -----
            CALL psb_erractionrestore(err_act)

100         FORMAT(' ERROR! Memory allocation failure in WRITE_MESH')
200         FORMAT(' ERROR! Unsupported output format in WRITE_MESH')

        END PROCEDURE write_mesh

END SUBMODULE write_mesh_implementation
