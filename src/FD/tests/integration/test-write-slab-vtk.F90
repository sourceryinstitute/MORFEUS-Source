!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
MODULE vtk_dtio_interface
  !! author: Damian Rouson
  !! date: 04/01/2019
  !!
  !! Encapsulate data for writing a VTK structured grid using derived-type input/output
  !!

  USE kind_parameters, ONLY : i4k
  USE vtk_datasets,   ONLY : struct_grid
  USE vtk_attributes, ONLY : attributes
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: vtk_dtio, n_params_to_write

  INTEGER(i4k), PARAMETER     :: n_params_to_write = 1

  TYPE vtk_dtio
    !! vtk_legacy_write required arguments
    CHARACTER(LEN=LEN('slab.vtk')) :: filename = 'slab.vtk'
    TYPE (struct_grid) grid
    TYPE (attributes), DIMENSION(n_params_to_write) :: vals_to_write
  CONTAINS
    PROCEDURE :: write_formatted
#ifdef HAVE_UDDTIO
    GENERIC :: WRITE(FORMATTED)=>write_formatted
#endif
  END TYPE

  INTERFACE

    MODULE SUBROUTINE write_formatted (this,unit,iotype, v_list, iostat, iomsg)
      !! Write a vtk_dtio object via user-defined derived type output wrapping vtk_legacy_write
      IMPLICIT NONE
      CLASS(vtk_dtio), INTENT(IN) ::this
      INTEGER, INTENT(IN) :: unit, v_list(:)
      CHARACTER (LEN=*), INTENT(IN) :: iotype
      INTEGER, INTENT(OUT) :: iostat
      CHARACTER(LEN=*), INTENT(INOUT) :: iomsg
    END SUBROUTINE

  END INTERFACE

END MODULE

SUBMODULE(vtk_dtio_interface) vtk_dtio_implementation
  IMPLICIT NONE
CONTAINS
  MODULE PROCEDURE write_formatted
    USE vtk, ONLY : vtk_legacy_write
    !! Invoke vtk_legacy_write with this object's components as arguments
    CALL vtk_legacy_write (unit=unit, geometry=this%grid, filename=this%filename, pointdatasets=this%vals_to_write)
  END PROCEDURE write_formatted
END SUBMODULE vtk_dtio_implementation

PROGRAM Slab_VTK_output
    USE kind_parameters, only : i4k, r8k
    USE vtk_attributes, ONLY : scalar
    USE vtk_dtio_interface, only : vtk_dtio, n_params_to_write
    IMPLICIT NONE
    !! author: Damian Rouson and Ian Porter
    !! date: 03/22/2019
    !!
    !! This tests output of a slab geometry as a VTK structured grid
    !!
    TYPE(vtk_dtio) slab
    INTEGER(i4k)                :: i, j, k
    INTEGER(i4k),     PARAMETER :: unit = 20
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

    OPEN(unit,FILE=slab%filename)
#ifdef HAVE_UDDTIO
    WRITE(unit,'(DT)') slab
#else
    BLOCK
      INTEGER :: io_status
      CHARACTER(LEN=132) :: io_message
      CALL slab%write_formatted(unit, iotype='DT', v_list=[INTEGER::], iostat=io_status, iomsg=io_message)
    END BLOCK
#endif

    WRITE(*,*) 'Test passed.'

END PROGRAM Slab_VTK_output
