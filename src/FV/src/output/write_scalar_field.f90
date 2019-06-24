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
! $Id$
!
! Description:
!    To be added...
!
SUBMODULE (tools_output) write_scalar_field_implementation
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE  write_scalar_field
            USE class_psblas, ONLY : psb_dpk_, psb_erractionsave, mypnum_, icontxt_, psb_cd_get_global_cols, psb_gather
            IMPLICIT NONE
            !
            INTEGER :: err_act, icontxt, info, mypnum, ncells_glob
            REAL(psb_dpk_), ALLOCATABLE :: x_loc(:)
            REAL(psb_dpk_), ALLOCATABLE :: x_glob(:)
            CHARACTER(len=32) :: path
            TYPE(cell), ALLOCATABLE :: cells_glob(:)
            TYPE(mesh), POINTER :: msh => NULL()

            ! Sets error handling for PSBLAS-2 routines
            CALL psb_erractionsave(err_act)

            mypnum  = mypnum_()
            icontxt = icontxt_()

            ! Sets output path
            IF(PRESENT(iter)) CALL out%set_output_path(iter)
            path = out%path_()

        !!$  msh   => msh_(fld)
            CALL get_mesh(fld,msh)
            CALL get_x(fld,x_loc)

            ! Is FLD cell-centered?
            IF(fld%on_faces_()) THEN
                WRITE(*,100)
                CALL abort_psblas
            END IF

            ! Global number of cells
            ncells_glob = psb_cd_get_global_cols(msh%desc_c)

            ALLOCATE(x_glob(ncells_glob),stat=info)
            IF(info /= 0) THEN
                WRITE(*,200)
                CALL abort_psblas
            END IF

            ! Gathers cell-centered values
            CALL psb_gather(x_glob,x_loc,msh%desc_c,info,root=0)
            CALL psb_check_error(info,'write_scalar_field','psb_gather',icontxt)

            ! Local-to-global reallocation of CELL object
            CALL l2g_cell(msh%cells,cells_glob,msh%desc_c)

            IF(mypnum == 0) THEN
                SELECT CASE(out%fmt_())
                    !     case(vtk_)
                    !        call wr_vtk_field(field,x_glob,msh%ncd,path)
                CASE DEFAULT
                    WRITE(*,300)
                    CALL abort_psblas
                END SELECT
                WRITE(*,'()')
            END IF


            IF(ALLOCATED(cells_glob)) CALL free_cell(cells_glob)
            DEALLOCATE(x_glob,x_loc)
            NULLIFY(msh)


            ! ----- Normal Termination -----
            CALL psb_erractionrestore(err_act)

100         FORMAT(' ERROR! Face-centered field in WRITE_SCALAR_FIELD')
200         FORMAT(' ERROR! Memory allocation failure in WRITE_SCALAR_FIELD')
300         FORMAT(' ERROR! Unsupported output format in WRITE_SCALAR_FIELD')

        END PROCEDURE write_scalar_field

END SUBMODULE write_scalar_field_implementation
