!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
PROGRAM T_shape_test
    !! author: Damian Rouson
    !! date: 05/03/2018
    !!
    !! Write a T-shaped, unstructured-grid geometry defined in VTK voxels
    USE kind_parameters, ONLY : i4k, r8k
    USE vtk_datasets,    ONLY : unstruct_grid
    IMPLICIT NONE
    TYPE (unstruct_grid)  :: t_shape
    INTEGER(i4k), PARAMETER :: n_points = 24, n_cells = 5

    BLOCK
        !! Initialize grid
        USE vtk_cells, ONLY : voxel, vtkcell_list
        TYPE(voxel), DIMENSION(n_cells) :: voxel_cells     !! Voxel cell type
        TYPE(vtkcell_list), DIMENSION(n_cells) :: cell_list       !! Full list of all cells
        INTEGER(i4k) :: i
        INTEGER(i4k), PARAMETER :: n_spacedim=3
        REAL(r8k), DIMENSION(*,*), PARAMETER :: points = RESHAPE ( &
            & [ 0.5_r8k, 0.0_r8k, 0.0_r8k, &
            &   1.0_r8k, 0.0_r8k, 0.0_r8k, &
            &   0.5_r8k, 0.5_r8k, 0.0_r8k, &
            &   1.0_r8k, 0.5_r8k, 0.0_r8k, &
            &   0.5_r8k, 0.0_r8k, 0.5_r8k, &
            &   1.0_r8k, 0.0_r8k, 0.5_r8k, &
            &   0.5_r8k, 0.5_r8k, 0.5_r8k, &
            &   1.0_r8k, 0.5_r8k, 0.5_r8k, &
            &   0.0_r8k, 0.0_r8k, 1.0_r8k, &
            &   0.5_r8k, 0.0_r8k, 1.0_r8k, &
            &   1.0_r8k, 0.0_r8k, 1.0_r8k, &
            &   1.5_r8k, 0.0_r8k, 1.0_r8k, &
            &   0.0_r8k, 0.5_r8k, 1.0_r8k, &
            &   0.5_r8k, 0.5_r8k, 1.0_r8k, &
            &   1.0_r8k, 0.5_r8k, 1.0_r8k, &
            &   1.5_r8k, 0.5_r8k, 1.0_r8k, &
            &   0.0_r8k, 0.0_r8k, 1.5_r8k, &
            &   0.5_r8k, 0.0_r8k, 1.5_r8k, &
            &   1.0_r8k, 0.0_r8k, 1.5_r8k, &
            &   1.5_r8k, 0.0_r8k, 1.5_r8k, &
            &   0.0_r8k, 0.5_r8k, 1.5_r8k, &
            &   0.5_r8k, 0.5_r8k, 1.5_r8k, &
            &   1.0_r8k, 0.5_r8k, 1.5_r8k, &
            &   1.5_r8k, 0.5_r8k, 1.5_r8k ], [n_spacedim, n_points] )

        CALL voxel_cells(1)%setup ( [ 0, 1, 2, 3, 4, 5, 6, 7 ] )
        CALL voxel_cells(2)%setup ( [ 4, 5, 6, 7, 9, 10, 13, 14 ] )
        CALL voxel_cells(3)%setup ( [ 8, 9, 12, 13, 16, 17, 20, 21 ] )
        CALL voxel_cells(4)%setup ( [ 9, 10, 13, 14, 17, 18, 21, 22 ] )
        CALL voxel_cells(5)%setup ( [ 10, 11, 14, 15, 18, 19, 22, 23 ] )

        DO i=1,SIZE(voxel_cells)                                 !! workaround gfortran 8.3.0 bug
            ALLOCATE(cell_list(i)%cell, source=voxel_cells(i))   !! Alternative: cell_list(i)%cell = voxel_cells(i)
        END DO

        CALL t_shape%init (points=points, cell_list=cell_list)
    END BLOCK

    BLOCK
        !! Defne scalar quantities and write grid
        USE vtk, ONLY : vtk_serial_write
        USE vtk_attributes, ONLY : attributes
        TYPE (attributes) :: cell_vals_to_write, point_vals_to_write
        INTEGER(i4k) :: j
        INTEGER(i4k), DIMENSION(*), PARAMETER ::  cell_ID = [ (j,j=1,n_cells ) ]
        INTEGER(i4k), DIMENSION(*), PARAMETER :: point_ID = [ (j,j=1,n_points) ]

        CALL define_scalar(  cell_vals_to_write, REAL( cell_ID, r8k),  'Cell_ID' )
        CALL define_scalar( point_vals_to_write, REAL(point_ID, r8k), 'Point_ID' )

        CALL vtk_serial_write (                  &
            t_shape,                             &
            celldatasets=[cell_vals_to_write],   &
            pointdatasets=[point_vals_to_write], &
            filename='test-write-voxels',        &
            multiple_io=.FALSE.)
    END BLOCK

    WRITE(*,*) 'Test passed.'

CONTAINS

    SUBROUTINE define_scalar( s, vals, dataname )
        USE vtk_attributes, ONLY : scalar, attributes
        CLASS(attributes), INTENT(INOUT) :: s
        REAL(r8k),         INTENT(IN)    :: vals(:)
        CHARACTER(LEN=*),  INTENT(IN)    :: dataname

        IF (.NOT. ALLOCATED(s%attribute)) THEN
            ALLOCATE(scalar::s%attribute)
        END IF

        CALL s%attribute%init (dataname, numcomp=1, real1d=vals)
    END SUBROUTINE

END PROGRAM T_shape_test
