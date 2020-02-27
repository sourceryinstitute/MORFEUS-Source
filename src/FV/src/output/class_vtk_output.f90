!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
MODULE class_vtk_output
    USE class_output, ONLY : output
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
    PUBLIC :: vtk_output_

    TYPE, EXTENDS(output) :: vtk_output_
        !! DT for writing vtk files
    CONTAINS
        PROCEDURE :: write_output => write_vtk_morfeus
    END TYPE vtk_output_

    INTERFACE

        MODULE SUBROUTINE write_vtk_morfeus (out, msh, sfield, vfield, iter)
            USE class_iterating, ONLY : iterating
            IMPLICIT NONE
            !! author: Ian Porter, GSE
            !! date: 01/04/2020
            !!
            !! This subroutine translates morfeus data structure
            !! into a vtkmofo data structure and writes the vtk file
            !!
            CLASS(vtk_output_),               INTENT(INOUT)        :: out          !! DT of output file info
            TYPE(mesh),                       INTENT(IN)           :: msh          !! DT of mesh info
            TYPE(scalar_field), DIMENSION(:), INTENT(IN), OPTIONAL :: sfield       !! DT of scalar info
            TYPE(vector_field), DIMENSION(:), INTENT(IN), OPTIONAL :: vfield       !! DT of vector info
            TYPE(iterating),                  INTENT(IN), OPTIONAL :: iter         !! DT of iteration info

        END SUBROUTINE write_vtk_morfeus

    END INTERFACE

END MODULE class_vtk_output
