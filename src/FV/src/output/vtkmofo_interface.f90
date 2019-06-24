!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
MODULE vtkmofo_io
    USE class_psblas, ONLY : psb_dpk_
    USE class_mesh,   ONLY : mesh
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_scalar_field, ONLY : scalar_field
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    USE class_vector_field, ONLY : vector_field
      !! An Intel 18.0.5 bug precludes putting this in the interface bodies
    IMPLICIT NONE
    !! author: Ian Porter, NRC
    !! date: 01/23/2019
    !!
    !! This module contains the routines necessary to interface morfeus w/ vtkmofo
    !!
    PRIVATE
    PUBLIC :: write_vtk_morfeus

    INTERFACE

        MODULE SUBROUTINE write_vtk_morfeus (msh, sfield, sname, vfield, vname, out, iter)
            USE class_iterating, ONLY : iterating
            USE class_output,    ONLY : output
            IMPLICIT NONE
            !! author: Ian Porter, NRC
            !! date: 01/23/2019
            !!
            !! This subroutine translates morfeus data structure
            !! into a vtkmofo data structure and writes the vtk file
            !!
            TYPE(mesh),                       INTENT(IN)           :: msh          !! DT of mesh info
            TYPE(scalar_field), DIMENSION(:), INTENT(IN), OPTIONAL :: sfield       !! DT of scalar info
            CHARACTER(LEN=*),   DIMENSION(:), INTENT(IN), OPTIONAL :: sname        !! Scalar names
            TYPE(vector_field), DIMENSION(:), INTENT(IN), OPTIONAL :: vfield       !! DT of vector info
            CHARACTER(LEN=*),   DIMENSION(:), INTENT(IN), OPTIONAL :: vname        !! Vector names
            TYPE(output),                     INTENT(INOUT)        :: out          !! DT of output file info
            TYPE(iterating),                  INTENT(IN), OPTIONAL :: iter         !! DT of iteration info

        END SUBROUTINE write_vtk_morfeus

        MODULE SUBROUTINE write_vtkmofo (ncd, points, cell_ids, v2cconn, icverts, &
            &                              iproc, sfield, sname, file_name)
            IMPLICIT NONE
            !! author: Ian Porter, NRC
            !! date: 01/23/2019
            !!
            !! This subroutine translates morfeus data structure
            !! into a vtkmofo data structure and writes the vtk file
            !!
            INTEGER,                          INTENT(IN) :: ncd          !!
            CHARACTER(LEN=*),                 INTENT(IN) :: file_name    !! Output file name
            INTEGER,          DIMENSION(:),   INTENT(IN) :: cell_ids     !! ID of each cell
            INTEGER,          DIMENSION(:),   INTENT(IN) :: icverts      !! # of vertices for each cell
            INTEGER,          DIMENSION(:),   INTENT(IN) :: v2cconn      !! For each cell, lists vertex ids
            INTEGER,          DIMENSION(:),   INTENT(IN) :: iproc        !! Processor ID
            REAL(psb_dpk_),   DIMENSION(:,:), INTENT(IN) :: points       !! Vertex x, y and z positions
            REAL(psb_dpk_),   DIMENSION(:),   INTENT(IN) :: sfield       !! Scalar fields
            CHARACTER(LEN=*), DIMENSION(:),   INTENT(IN) :: sname        !! Scalar name

        END SUBROUTINE write_vtkmofo

    END INTERFACE

END MODULE vtkmofo_io
