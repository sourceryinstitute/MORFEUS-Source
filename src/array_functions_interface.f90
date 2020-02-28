!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
MODULE array_functions_interface
  !! author: Damian Rouson
  !! date: 04/25/2019
  !!
  !! Functionally pure array utilities
  USE kind_parameters, ONLY : r8k
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: OPERATOR(.catColumns.)
  PUBLIC :: OPERATOR(.catRows.)
  PUBLIC :: OPERATOR(.columnVectors.)
    !! Because the Fortran standard requires that user-defined operator dummy arguments have the INTENT(IN)
    !! attribute, exposing only the operator and not the function names communicates more information in the
    !! public interface and in code using this interface.

    INTERFACE OPERATOR(.columnVectors.)
        MODULE PROCEDURE :: column_vectors
    END INTERFACE

    INTERFACE OPERATOR(.catColumns.)
        MODULE PROCEDURE :: concatenate_columns
    END INTERFACE

    INTERFACE OPERATOR(.catRows.)
        MODULE PROCEDURE :: concatenate_rows
    END INTERFACE

    INTERFACE

        PURE MODULE FUNCTION column_vectors(vector_field) RESULT(array_of_3D_column_vectors)
            !! Result is array of 3D column vectors of dimension (space_dim,nx*ny*nz) reshaped from vector-field argument
            !! of dimension (nx,ny,nz,space_dim)
            IMPLICIT NONE
            REAL(r8k), DIMENSION(:,:,:,:), INTENT(IN) :: vector_field
            REAL(r8k), DIMENSION(:,:), ALLOCATABLE ::  array_of_3D_column_vectors
        END FUNCTION

        PURE MODULE FUNCTION concatenate_columns(a, b) RESULT(concatenated)
            !! Result contains the concatenation of the columns of argument a with the columns of argument b
            IMPLICIT NONE
            REAL(r8k), DIMENSION(:,:), INTENT(IN) :: a, b
            REAL(r8k), DIMENSION(:,:), ALLOCATABLE :: concatenated
        END FUNCTION

        PURE MODULE FUNCTION concatenate_rows(a, b) RESULT(concatenated)
            !! Result contains the concatenation of the rows of argument a with the rows of argument b
            IMPLICIT NONE
            REAL(r8k), DIMENSION(:,:), INTENT(IN) :: a, b
            REAL(r8k), DIMENSION(:,:), ALLOCATABLE :: concatenated
        END FUNCTION

    END INTERFACE

END MODULE
