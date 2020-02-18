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

SUBMODULE(class_output) class_output_procedures
    USE class_iterating
    IMPLICIT NONE

CONTAINS

    MODULE PROCEDURE nemo_output_sizeof
        USE class_psblas, ONLY : nemo_sizeof_int
        IMPLICIT NONE

        nemo_output_sizeof = nemo_sizeof_int + LEN(obj%basepath)  + LEN(obj%path)

    END PROCEDURE nemo_output_sizeof

    ! ----- Constructor -----

    MODULE PROCEDURE create_output
        USE tools_input
        USE tools_output_basics, ONLY : csv_, vtk_, cgns_, exodus_
        USE class_psblas, ONLY : abort_psblas, mypnum_
        USE json_module
        USE class_vtk_output, ONLY : vtk_output_
        IMPLICIT NONE

        CHARACTER(LEN=10) :: proc_id
        CHARACTER(LEN=80) :: output_sec
        CHARACTER(KIND=json_CK,LEN=:),ALLOCATABLE :: cval
        LOGICAL :: found
        TYPE(json_file) :: nemo_json

!        IF (mypnum_() /= 0) RETURN

        CALL open_file(input_file,nemo_json)
        output_sec = 'MORFEUS_FV.'//TRIM(sec)

        ! Gets format (i.e., 'vtk, 'csv', 'exodus', 'hdf5', etc)
        CALL nemo_json%get(TRIM(output_sec)//'.format', cval, found)
        IF (.NOT.found) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF
        SELECT CASE (cval)
        CASE ('csv')
            ALLOCATE(output::out)
            out%fmt = csv_
        CASE ('vtk')
            ALLOCATE(vtk_output_::out)
            out%fmt = vtk_
        CASE ('cgns')
            ALLOCATE(output::out)
            out%fmt = cgns_
        CASE ('exodus')
            ALLOCATE(output::out)
            out%fmt = exodus_
        CASE DEFAULT
            WRITE(*,105) cval
            CALL abort_psblas
        END SELECT

        ! Gets path basename
        CALL nemo_json%get(TRIM(output_sec)//'.base-path', cval, found)
        IF (.NOT.found) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        WRITE(proc_id,'(i10)') mypnum_()
        out%basepath  = cval // '_image_' // TRIM(ADJUSTL(proc_id))! // '_'
        ! Sets initial path
        out%path = out%basepath

100     FORMAT('Missing OUTPUT parameters')
105     FORMAT(/,'Unsupported output type: ',a,/)

    END PROCEDURE create_output

    ! ----- Getters -----

    MODULE PROCEDURE fmt_
        IMPLICIT NONE

        fmt_ = out%fmt

    END PROCEDURE fmt_

    MODULE PROCEDURE path_
        USE tools_output_basics, ONLY : csv_, vtk_, exodus_, cgns_
        IMPLICIT NONE

        SELECT CASE(out%fmt)
        CASE (csv_)
            path = TRIM(out%path) // '.csv'
        CASE (vtk_)
            path = TRIM(out%path) !! Do no extension. vtkmofo handles this based on cell type
        CASE (exodus_)
            path = TRIM(out%path) // '.e'
        CASE (cgns_)
            path = TRIM(out%path) // '.csv'
        END SELECT

    END PROCEDURE path_

    ! ----- Setters -----

    MODULE PROCEDURE set_output_path_h
        IMPLICIT NONE

        out%path = TRIM(path)

    END PROCEDURE set_output_path_h

    MODULE PROCEDURE set_output_path_iter
        USE tools_output_basics, ONLY : itoh
        IMPLICIT NONE
        INTEGER :: it, ndigits

        ndigits = INT(LOG10(REAL(iter%nmax_()))) + 1
        it = iter%current_iteration()

        out%path = TRIM(out%basepath)//itoh(it,ndigits)

    END PROCEDURE set_output_path_iter

    MODULE PROCEDURE write_output
        IMPLICIT NONE

    END PROCEDURE write_output

    MODULE PROCEDURE get_scalar_field
        USE class_psblas
        USE class_cell
        USE class_iterating
        USE class_mesh
        USE class_output
        USE class_scalar_field
        USE tools_output_basics
        IMPLICIT NONE
        !
        INTEGER :: err_act, icontxt, info, ncells_glob
        REAL(psb_dpk_), ALLOCATABLE :: x_loc(:)
        TYPE(cell), ALLOCATABLE :: cells_glob(:)
        TYPE(mesh), POINTER :: msh => NULL()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

    !!$  msh   => msh_(fld)
!        CALL fld%get_mesh(msh)
!        CALL fld%get_x(x_loc)

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

        CALL psb_check_error(info,'set_scalar_field','psb_gather',icontxt)

        ! Local-to-global reallocation of CELL object
        CALL l2g_cell(msh%cells,cells_glob,msh%desc_c)

        IF (ALLOCATED(cells_glob)) CALL free_cell(cells_glob)
        IF (ALLOCATED(x_loc)) DEALLOCATE(x_loc)
        IF (ASSOCIATED(msh)) NULLIFY(msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Face-centered field in set_scalar_field')
200     FORMAT(' ERROR! Memory allocation failure in set_scalar_field')
!300     FORMAT(' ERROR! Unsupported output format in set_scalar_field')

    END PROCEDURE get_scalar_field

    MODULE PROCEDURE get_vector_field
        USE class_psblas
        USE class_cell
        USE class_vector!, ONLY : vector, ASSIGNMENT(=)
        USE class_iterating
        USE class_mesh
        USE class_output
        USE class_vector_field
        USE tools_output_basics
        USE class_scalar_field
        IMPLICIT NONE
        !
        INTEGER :: err_act, icontxt, info, mypnum, ncells_glob, j
        TYPE(vector), ALLOCATABLE :: x_loc_v(:)
        TYPE(vector), ALLOCATABLE :: x_glob_v(:)
        TYPE(cell), ALLOCATABLE :: cells_glob(:)
        TYPE(mesh), POINTER :: msh => NULL()

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        mypnum  = mypnum_()
        icontxt = icontxt_()

    !!$  msh   => msh_(fld)
!        CALL fld%get_mesh(msh)
!        CALL fld%get_x(x_loc_v)

        ! Is FLD cell-centered?
        IF(fld%on_faces_()) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Global number of cells
        ncells_glob = psb_cd_get_global_cols(msh%desc_c)

        ALLOCATE(x_glob_v(ncells_glob),stat=info)
        IF(info /= 0) THEN
            WRITE(*,200)
            CALL abort_psblas
        END IF

        ! Gathers cell-centered values
        CALL l2g_vector(x_loc_v,x_glob_v,msh%desc_c)
    !!$  call psb_check_error(info,'set_vector_field','psb_gather',icontxt)

        ! Local-to-global reallocation of CELL object
        CALL l2g_cell(msh%cells,cells_glob,msh%desc_c)

        !if(mypnum == 0) then
        !   call wr_vtk_field(field,x_glob_v,msh%ncd,path)
        !end if
        DO j = 1, SIZE(x_glob_v)
            x_glob(1:3,j) = x_glob_v(j)
        END DO

        IF (ALLOCATED(cells_glob)) CALL free_cell(cells_glob)
        IF (ALLOCATED(x_loc_v)) DEALLOCATE(x_loc_v)
        IF (ALLOCATED(x_glob_v)) DEALLOCATE(x_glob_v)
        IF (ASSOCIATED(msh)) NULLIFY(msh)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Face-centered field in set_vector_field')
200     FORMAT(' ERROR! Memory allocation failure in set_vector_field')

    END PROCEDURE get_vector_field

END SUBMODULE class_output_procedures
