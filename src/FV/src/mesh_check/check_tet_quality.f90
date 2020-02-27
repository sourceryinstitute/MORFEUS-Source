!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
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
! $Id: check_tet_quality.f90 3093 2008-04-22 14:51:09Z sfilippo $
!
! Description:
!    Dumps information helpful in debugging a bad cell, including writing a
!    file in DX or VTK format for visualization.
!    WARNING! Only tet cells are currently supported.
!
SUBMODULE (tools_mesh_check) tools_mesh_check_tet
    IMPLICIT NONE

    CONTAINS

        MODULE PROCEDURE check_tet_quality
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_mesh
        USE class_vector
        USE class_vertex
        USE tools_math
        USE tools_mesh_basics, ONLY : geom_tet_quality, geom_tet_dihedral_angle
        USE tools_output_basics

        IMPLICIT NONE
        !
        INTEGER :: fmt_, ic_glob
        INTEGER :: i, iv, iv1, iv2, iv3, iv4
        INTEGER, POINTER :: iv2c(:) => NULL()
        INTEGER, POINTER :: if2c(:) => NULL()
        CHARACTER(len=20) :: file_name
        REAL(psb_dpk_) :: max_angle, min_angle
        TYPE(vector)     :: tet_af(4)

        ! Gets cell global index
        ic_glob = loc_to_glob_(msh%desc_c,ic)

        ! Get vertex indices
        CALL msh%v2c%get_ith_conn(iv2c,ic)
        iv1 = iv2c(1)
        iv2 = iv2c(2)
        iv3 = iv2c(3)
        iv4 = iv2c(4)

        ! Get face indices
        CALL msh%f2c%get_ith_conn(if2c,ic)

        DO i = 1,4
            tet_af(i) = msh%af(if2c(i))
        ENDDO

        ! Evaluates max and min dihedral angle
        CALL geom_tet_dihedral_angle( tet_af, max_angle, min_angle)

        ! Dumps log message
        WRITE(*,100) 'Debugging cell ',ic,' on proc. ', mypnum_()
        WRITE(*,200) '- Global ID: ', ic_glob
        WRITE(*,300) '- Cell type: ', msh%cells(ic)%geo_()
        WRITE(*,400) '- Expecting ', SIZE(iv2c), ' vertices'
        DO i = 1, 4
            iv = iv2c(i)
            WRITE(*,500) '- Coordinates of vertex ',i,': ', &
                &       msh%verts(iv)%x_(), &
                &       msh%verts(iv)%y_(), &
                &       msh%verts(iv)%z_()
        END DO
        WRITE(*,600) '- Quality: ', geom_tet_quality(&
            & msh%verts(iv1),msh%verts(iv2),msh%verts(iv3),msh%verts(iv4), &
            & msh%vol(ic))
        WRITE(*,600) '- Volume: ', msh%vol(ic)
        WRITE(*,700) '- Maximum dihedral angle: ', &
            & 180.0 * max_angle / pi
        WRITE(*,700) '- Mininum dihedral angle: ', &
            & 180.0 * min_angle / pi
        WRITE(*,'()')


        ! Dumps on file
        IF(PRESENT(fmt)) THEN
            fmt_ = fmt
        ELSE
            fmt_ = vtk_
        END IF

        file_name='cell_'//itoh(ic,INT(LOG10(REAL(ic))) + 1)//'p'//itoh(mypnum_(),2)

        SELECT CASE(fmt_)
        CASE(vtk_)
            !! TODO: Put in functionality for tets w/ VTKMOFO
            !!
            !! Old call:
            !!
            !! call wr_vtk_tet_cell(trim(file_name)//'.vtk',&
            !!   & (/msh%verts(iv1),msh%verts(iv2),msh%verts(iv3),msh%verts(iv4)/))
            !! interface
            !!    subroutine wr_vtk_tet_cell(file_name,verts)
            !!      use class_vertex
            !!      character(len=*), intent(in) :: file_name
            !!      type(vertex),     intent(in) :: verts(4)
            !!    end subroutine wr_vtk_tet_cell
            !! end interface

        CASE DEFAULT
            WRITE(*,900)
            CALL abort_psblas
        END SELECT

        NULLIFY(if2c,iv2c)

100     FORMAT(1x,a,i5,a,i2)
200     FORMAT(1x,a,i5)
300     FORMAT(1x,a,a3)
400     FORMAT(1x,a,i1,a)
500     FORMAT(1x,a,i1,a,3(e13.6,2x))
600     FORMAT(1x,a,e13.6)
700     FORMAT(1x,a,f8.4)
900     FORMAT(' ERROR! Unsupported output format in CHECK_TET_QUALITY')

        END PROCEDURE check_tet_quality

END SUBMODULE tools_mesh_check_tet
