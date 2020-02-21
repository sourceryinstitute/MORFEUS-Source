!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program test_keytable
  !! label: Morfeus-FV
  !!
  !! Test keytable class behavior

    use assertions_interface, only : assert
    use class_keytable, only : keytable

    implicit none
    integer, allocatable :: row1(:), row2(:), row3(:)
    integer, pointer :: row1T(:), row2T(:), row3T(:)
    type(keytable) :: testTable

    row1 = [1,2,3]
    row2 = [4,5,6,7]
    row3 = [8,9,10,11,12,13]

    call testTable%alloc_keytable(1,3)
    call testTable%set_kt_row(1,row1)
    call testTable%set_kt_row(2,row2)
    call testTable%set_kt_row(3,row3)

    call testTable%get_kt_row(1,row1T)
    call testTable%get_kt_row(2,row2T)
    call testTable%get_kt_row(3,row3T)

    call assert( row1T == row1, "test_keytable: testTable%row(1) == row1" )
    call assert( row2T == row2, "test_keytable: testTable%row(2) == row2" )
    call assert( row3T == row3, "test_keytable: testTable%row(3) == row3" )
    call assert( testTable%get_row_size(1) == size(row1), "test_keytable: testTable%get_row_size(1) == size(row1)" )
    call assert( testTable%get_row_size(2) == size(row2), "test_keytable: testTable%get_row_size(2) == size(row2)" )
    call assert( testTable%get_row_size(3) == size(row3), "test_keytable: testTable%get_row_size(3) == size(row3)" )

    print *, "Test passed."

end program test_keytable
