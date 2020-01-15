!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
PROGRAM Slab_VTK_output
    USE kind_parameters, ONLY : i4k, r8k
    USE vtk, ONLY : vtk_serial_write
    USE vtk_attributes, ONLY : scalar, attributes
    USE vtk_datasets, ONLY : struct_grid
    IMPLICIT NONE
    !! author: Damian Rouson and Ian Porter
    !! date: 03/22/2019
    !!       11/25/2019 Modified by IP to remove DTIO due to vtkmofo doing the full file handling
    !!
    !! This tests output of a slab geometry as a VTK structured grid
    !!
    INTEGER(i4k), PARAMETER :: n_params_to_write = 1
    TYPE vtk_obj
      CHARACTER(LEN=LEN('test-write-slab')) :: filename = 'test-write-slab'
      TYPE(struct_grid) :: grid
      TYPE(attributes), DIMENSION(n_params_to_write) :: vals_to_write
    END TYPE vtk_obj
    TYPE (vtk_obj) :: slab
    INTEGER(i4k) :: i, j, k
    REAL(r8k), DIMENSION(*), PARAMETER :: x_vals = &
      & [ 0.00E+00_r8k, 8.03E-04_r8k, 1.51E-03_r8k, 2.12E-03_r8k, 2.64E-03_r8k, &
      &   3.08E-03_r8k, 3.45E-03_r8k, 3.75E-03_r8k, 3.99E-03_r8k, 4.18E-03_r8k, &
      &   4.32E-03_r8k, 4.42E-03_r8k, 4.49E-03_r8k, 4.53E-03_r8k, 4.56E-03_r8k, &
      &   4.56E-03_r8k, 4.56E-03_r8k, 4.65E-03_r8k, 5.38E-03_r8k ]
    REAL(r8k), DIMENSION(*), PARAMETER :: y_vals = &
      & [ 0.00E+00_r8k , 2.*maxval(x_vals) ]
    REAL(r8k), DIMENSION(*), PARAMETER :: z_vals = &
      & [ 2.50E-03_r8k, 5.00E-03_r8k, 7.50E-03_r8k, 1.00E-03_r8k ]
    INTEGER, PARAMETER :: n_x=SIZE(x_vals), n_y=size(y_vals), n_z=size(z_vals)
    REAL(r8k), DIMENSION(n_x*n_y*n_z) :: temperature
    REAL(r8k), DIMENSION(1:3,n_x*n_y*n_z) :: points

    DO CONCURRENT( i = 1:n_x, j = 1:n_y, k = 1:n_z )
      ASSOCIATE( cnt => (k-1)*n_x*n_y + (j-1)*n_x + i )
        points(:,cnt) = [x_vals(i), y_vals(j), z_vals(k)]
        temperature(cnt) = x_vals(i)
      END ASSOCIATE
    END DO

    CALL slab%grid%init (dims=[ n_x, n_y, n_z ], points=points)

    DO i = 1, n_params_to_write
        IF (.NOT. ALLOCATED(slab%vals_to_write(i)%attribute))THEN
            ALLOCATE(scalar::slab%vals_to_write(i)%attribute)
        END IF
        CALL slab%vals_to_write(i)%attribute%init ('Temperature_(K)     ' , numcomp=1, real1d=temperature)
    END DO

    !! Invoke vtk_serial_write with this object's components as arguments
    CALL vtk_serial_write (filename=slab%filename, geometry=slab%grid, pointdatasets=slab%vals_to_write)

    WRITE(*,*) 'Test passed.'

END PROGRAM Slab_VTK_output
