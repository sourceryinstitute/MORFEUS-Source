!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
SUBMODULE(array_functions_interface) array_functions_implementation
  USE assertions_interface, ONLY : assert, assertions
  IMPLICIT NONE
CONTAINS

  MODULE PROCEDURE column_vectors
      INTEGER i, j, k

      ASSOCIATE( n => shape(vector_field) )
        IF (assertions) CALL assert(SIZE(n)==4, "3D vector field input")
        ALLOCATE( array_of_3D_column_vectors( n(4), PRODUCT(n(1:3)) ) )
        DO CONCURRENT( i=1:n(1), j=1:n(2), k=1:n(3) )
          ASSOCIATE( id => (k-1)*PRODUCT(n(1:2)) + (j-1)*n(1) + i )
            array_of_3D_column_vectors(:,id) =  vector_field(i,j,k,:)
          END ASSOCIATE
        END DO
      END ASSOCIATE

  END PROCEDURE

  MODULE PROCEDURE concatenate_columns
    !! Using RESHAPE rather than manipulating array elements directly frees the compiler to decide the particular order of array
    !! element references that best exploits the given platform.  Alternatively, DO CONCURRENT could instead free the compiler
    !! to order element accesses however is best. Trade-off: RESHAPE requires the creation of temporary array results but RESHAPE
    !! is likely to have more mature compiler support than DO CONCURRENT.  If this code turns out to be a critical performance
    !! bottleneck, try replacing this implementation with element-by-element copying using DO CONCURRENT.
    ASSOCIATE(rows=>SIZE(a,1))
    ASSOCIATE(cols=>SIZE(a,2)+SIZE(b,2))
    ASSOCIATE(a_unrolled=>RESHAPE(a,[SIZE(a)]))
    ASSOCIATE(b_unrolled=>RESHAPE(b,[SIZE(b)]))
      IF (assertions) CALL assert( rows==SIZE(b,1), "array_functions: compatible shapes")
      concatenated = RESHAPE( [a_unrolled, b_unrolled ],[rows, cols] )
    END ASSOCIATE; END ASSOCIATE; END ASSOCIATE; END ASSOCIATE
  END PROCEDURE

  MODULE PROCEDURE concatenate_rows
    !! For simplicity, this implementation invokes concatenate_columns at the cost of TRANSPOSE creating additional temporaries.
    !! If this code turns out to be a critical performance bottleneck, try replacing this implementation with element-by-element
    !! copying using DO CONCURRENT.
    concatenated = transpose( concatenate_columns(transpose(a),transpose(b)) )
  END PROCEDURE

END SUBMODULE
