!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
PROGRAM main
  USE Kinds, ONLY : r8k
  USE assertions_interface, ONLY : assert
  USE array_functions_interface, ONLY : OPERATOR(.catRows.), OPERATOR(.catColumns.), OPERATOR(.columnVectors.)
  IMPLICIT NONE

  REAL(r8k), PARAMETER, DIMENSION(*) :: u = [1.,2.,3.], v = [4.,5.,6.]
    !! Define two 3D vectors as test data
  INTEGER, PARAMETER :: num_vectors = SIZE( ["u","v"], 1 ), space_dimension=3

  test_catcolumns: BLOCK
    !! Store two 3D vectors u and v in columns of lhs and rhs, respectively and test the concatenation of the columns
    REAL(r8k), PARAMETER, DIMENSION(space_dimension,0) :: empty=RESHAPE( [REAL(r8k)::], [space_dimension,0])
      !! Array of zero 3D vectors
    REAL(r8k), PARAMETER, DIMENSION(*,*) :: lhs = RESHAPE( [u,v] , SHAPE=[SIZE(u,1), num_vectors] )
      !! Store u and v in columns of lhs
    REAL(r8k), PARAMETER, DIMENSION(*,*) :: rhs = 2.*lhs
      !! Store 2u and 2v in columns of rhs
    REAL(r8k), PARAMETER, DIMENSION(*,*) :: expected_result = RESHAPE( [ [u,v], 2.*[u,v] ], SHAPE=[SIZE(u,1), 2*num_vectors] )
      !! Expected result of concatenating the columns of lhs and rhs
    CALL assert( (empty .catColumns. (lhs .catColumns. rhs)) == expected_result ,"columns concatenated")
  END BLOCK test_catcolumns

  test_catrows: BLOCK
    !! Store two 3D vectors u and v in rows of lhs and rhs, respectively and test the concatenation of the rows
    INTEGER, PARAMETER :: space_dimension = SIZE(u,1 )
    REAL(r8k), DIMENSION(num_vectors,space_dimension) :: lhs , rhs
    REAL(r8k), DIMENSION(2*num_vectors,space_dimension) :: expected_result

    lhs(1,1:SIZE(u,1))  = u  !! Store u in 1st row of lhs
    lhs(2,1:SIZE(u,1))  = v  !! Store v in 2nd row of lhs
    rhs = 2.*lhs !! Store 2u and 2v in 1st and 2nd rows, respectively, of rhs
    expected_result(1:2,1:SIZE(u,1)) = lhs
    expected_result(3:4,1:SIZE(v,1)) = rhs

    CALL assert( (lhs .catRows. rhs) == expected_result ,"rows concatenated")
  END BLOCK test_catrows

  test_column_vectors: BLOCK
    INTEGER, PARAMETER :: n(4)=[2,2,2,3] !! [nx,ny,nz,nspace_dim]
    REAL(r8k), DIMENSION(n(1),n(2),n(3),n(4)) :: vector_field
    LOGICAL vectors_match(n(4),PRODUCT(n(1:3)))
    INTEGER i, j, k

    vector_field(1,1,1,:) = [0.,0.,0.]
    vector_field(2,1,1,:) = [1.,0.,0.]
    vector_field(1,2,1,:) = [0.,1.,0.]
    vector_field(2,2,1,:) = [1.,1.,0.]

    vector_field(1,1,2,:) = [0.,0.,1.]
    vector_field(2,1,2,:) = [1.,0.,1.]
    vector_field(1,2,2,:) = [0.,1.,1.]
    vector_field(2,2,2,:) = [1.,1.,1.]

    ASSOCIATE( colvec => .columnVectors. vector_field )
      CALL assert(SIZE(colvec)==SIZE(vector_field),"colummn-vector array size preserved")
      DO CONCURRENT(i=1:n(1), j=1:n(2), k=1:n(3))
        ASSOCIATE( id => (k-1)*PRODUCT(n(1:2)) + (j-1)*n(1) + i )
          vectors_match(:,id) = ( colvec(:,id) == vector_field(i,j,k,:) )
        END ASSOCIATE
      END DO
    END ASSOCIATE
    CALL assert( ALL(vectors_match),"colVectors extracted vector field")

  END BLOCK test_column_vectors

  print *,"Test passed."
END PROGRAM
