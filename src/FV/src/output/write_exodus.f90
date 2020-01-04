!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!

MODULE write_exodus

  USE iso_fortran_env, ONLY : r8k => real64
  USE class_mesh, ONLY : mesh

  IMPLICIT NONE
  INCLUDE 'exodusII.inc'

  !! author: Hari Radhakrishnan, GSE
  !! date 12/12/2019
  !!
  !! This module implements the routines for writing an exodus file with mesh and results

  CHARACTER(LEN=*), PARAMETER :: filename = "output.e"

CONTAINS

  SUBROUTINE write_exo_morfeus()

    INTEGER(kind = int64) :: ierr, titlelen
    INTEGER(kind = int64) :: exodus_file_id, num_dims, num_nodes, num_elems,  &
      num_elem_blks, num_node_sets, num_side_sets

    INTEGER(kind = int64) :: cpu_ws, io_ws, num_props, prop_value

    ! Open exodusII file for writing

    exo_file_id = exopen("output.exo", EXWRITE, cpu_ws, io_ws, vers, ierr)



  END SUBROUTINE write_exo_morfeus

END MODULE
