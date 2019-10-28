!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module adi_mod
  implicit none

  integer ,parameter :: digits=8  ! num. digits of kind
  integer ,parameter :: decades=9 ! num. representable decades
  integer ,parameter :: rkind = selected_real_kind(digits)
  integer ,parameter :: ikind = selected_int_kind(decades)

  type grid_block
    real(rkind), dimension(8,3) :: v
    real(rkind), dimension(8) :: T, ddx2, ddy2, ddz2
  end type grid_block

  type(grid_block)  ,dimension(:,:,:) ,allocatable  :: global_grid_block
  integer(ikind)  :: i,j,k
  integer(ikind) ,parameter :: nx=21, ny=21, nz=21

contains

  pure function position_vectors(nx,ny,nz) result(vector_field)
    integer(ikind), intent(in) :: nx, ny, nz
    integer(ikind) :: ii, jj, kk
    integer(ikind), parameter :: components=3
    real(rkind), parameter :: dx=0.75/(5.0), dy=0.75/(5.0)
    real(rkind), parameter :: dx_s=0.25/5.0, dy_s=0.25/5.0
    real(rkind), dimension(:,:,:,:) ,allocatable  :: vector_field

    allocate( vector_field(nx,ny,nz,components) )

    associate( dz=>10.0/(nz-1) )
      do concurrent( ii=1:nx, jj=1:ny, kk=1:nz )
        vector_field(i,j,k,1) = merge(-0.75+(i-6)*dx, merge( -1.0+(i-1)*dx_s, 0.75+(i-16)*dx_s, i<6), i>=6 .and. i<=16)
        vector_field(i,j,k,2) = merge(-0.75+(j-6)*dy, merge( -1.0+(j-1)*dy_s, 0.75+(j-16)*dy_s, j<6), j>=6 .and. j<=16)
        vector_field(i,j,k,3) = (k-1)*dz
      end do
    end associate
  end function

  subroutine write_component(grid_blocks, component, file_unit, string)
    type(grid_block), intent(in), dimension(:,:,:) :: grid_blocks
    character(len=*), intent(in) :: component
    character(len=*), intent(in), optional :: string
    integer, intent(in) :: file_unit
    integer(ikind) :: i, j, k, nv, ndimension
    associate( n=>shape(grid_blocks) )
    do k=1, n(3)
      do j=1, n(2)
        do i=1, n(1)
          do nv=1,merge(1,8,component=="string")
            select case(component)
              case("T")
                write (file_unit,*) grid_blocks(i,j,k)%T(nv)
              case("v")
                write(file_unit, *) (grid_blocks(i,j,k)%v(nv,ndimension), ndimension=1,3)
              case("ddx2")
                write (file_unit,*) grid_blocks(i,j,k)%ddx2(nv)
              case("string")
                write (file_unit,*) string
              case default
                error stop "write_component: unrecognized component"
            end select
          end do
        end do
      end do
    end do
    end associate
  end subroutine

  subroutine output_result(grid_blocks)
    type(grid_block), intent(in), dimension(:,:,:) :: grid_blocks
    integer file_unit

    open(newunit=file_unit, file="spatial-derivatives.vtk")
    write(file_unit, '(a)') '# vtk DataFile Version 3.0'
    write(file_unit, '(a)') '3d-htc voxels'
    write(file_unit, '(a)') 'ASCII'
    write(file_unit, '(a)') 'DATASET UNSTRUCTURED_GRID'
    write(file_unit, '(a,i5,a)') "POINTS ", 8*(nx-1)*(ny-1)*(nz-1), " double"
    call write_component( grid_blocks, "v", file_unit)

    block
      integer i, j, k, nv
      write(file_unit,'(a,i5,3x,i5)') "CELLS ", (nx-1)*(ny-1)*(nz-1), 9*(nx-1)*(ny-1)*(nz-1)
      associate( n=>shape(grid_blocks) )
        do k=1, n(3)
          do j=1, n(2)
            do i=1, n(1)
              associate( grid_block_id => 8*(i+(j-1)*n(1)+(k-1)*n(1)*n(2)-1) )
                write(file_unit,*) "8", (grid_block_id+nv-1, nv=1,8)
              end associate
            end do
          end do
        end do
      end associate
    end block

    write(file_unit,'(a,I5)') "CELL_TYPES ", (nx-1)*(ny-1)*(nz-1)
    call write_component( grid_blocks, "string", file_unit, string="11")

    write(file_unit,'(a, I5)') 'POINT_DATA ', 8*(nx-1)*(ny-1)*(nz-1)
    write(file_unit,'(a)') 'SCALARS Temperature double 1'
    write(file_unit,'(a)') 'LOOKUP_TABLE default'
    call write_component( grid_blocks, "T", file_unit)

    write(file_unit,'(a)') 'SCALARS ddx2_Temperature double 1'
    write(file_unit,'(a)') 'LOOKUP_TABLE default'
    call write_component( grid_blocks, "ddx2", file_unit)

    close(file_unit)

    block
      integer i, error_file_unit
      open(newunit=error_file_unit, file="error-order.dat")
      do i=2, nx-2
        write(error_file_unit, *) &
          (grid_blocks(i,2,2)%v(1,1)-grid_blocks(i-1,2,2)%v(1,1))**2, grid_blocks(i,2,2)%ddx2(1)-2.0*grid_blocks(i,2,2)%v(1,2)
      end do
      close(error_file_unit)
    end block
  end subroutine output_result

end module adi_mod

program main
  use adi_mod
  implicit none

  real(rkind), dimension(8) :: dx_m, dy_m, dz_m
  real(rkind), dimension(8) :: dx_f, dy_f, dz_f
  real(rkind), dimension(8) :: dx_b, dy_b, dz_b
  real(rkind)               :: avg_err, minimum_dx
  integer(ikind)            :: nv

  associate( p => position_vectors(nx,ny,nz) )

    allocate(global_grid_block(nx-1,ny-1,nz-1))

    associate( x=>(p(:,:,:,1)), y=>(p(:,:,:,2)), z=>(p(:,:,:,3)) )
      do concurrent(k=1:nz-1, j=1:ny-1, i=1:nx-1)
        associate( v=>global_grid_block(i,j,k)%v )
          v(1,:) = [x(i,j,k)      , y(i,j,k)      , z(i,j,k)      ]
          v(2,:) = [x(i+1,j,k)    , y(i+1,j,k)    , z(i+1,j,k)    ]
          v(3,:) = [x(i,j+1,k)    , y(i,j+1,k)    , z(i,j+1,k)    ]
          v(4,:) = [x(i+1,j+1,k)  , y(i+1,j+1,k)  , z(i+1,j+1,k)  ]
          v(5,:) = [x(i,j,k+1)    , y(i,j,k+1)    , z(i,j,k+1)    ]
          v(6,:) = [x(i+1,j,k+1)  , y(i+1,j,k+1)  , z(i+1,j,k+1)  ]
          v(7,:) = [x(i,j+1,k+1)  , y(i,j+1,k+1)  , z(i,j+1,k+1)  ]
          v(8,:) = [x(i+1,j+1,k+1), y(i+1,j+1,k+1), z(i+1,j+1,k+1)]
          associate(T=>global_grid_block(i,j,k)%T, ddx2=>global_grid_block(i,j,k)%ddx2 )
            T(:) = v(:,1)**2 * v(:,2) ! T=x^2*y => d^2T/dx^2=y
            ddx2(:) = 0.0
          end associate
        end associate
      end do
    end associate
  end associate

!  compute_derivatives: block
!    real(rkind), dimension(8) :: dx_m, dy_m, dz_m
!    real(rkind), dimension(8) :: dx_f, dy_f, dz_f
!    real(rkind), dimension(8) :: dx_b, dy_b, dz_b

    do concurrent( k=2:nz-2, j=2:ny-2, i=2:nx-2)

      !calculate ddx2
      associate(v=>global_grid_block(i,j,k)%v , ddx2=>global_grid_block(i,j,k)%ddx2, &
        T=>global_grid_block(i,j,k)%T,  ddy2=>global_grid_block(i,j,k)%ddy2, ddz2=>global_grid_block(i,j,k)%ddz2 )

        associate(nv=>[1,3,5,7])
          dx_m(nv) = 0.5*(v(nv+1,1) - global_grid_block(i-1,j,k)%v(nv,1))
          dx_f(nv) =      v(nv+1,1) - v(nv,1)
          dx_b(nv) =      v(nv,1)   - global_grid_block(i-1,j,k)%v(nv,1)
          ddx2(nv) = (T(nv+1) - T(nv)  )/(dx_f(nv)*dx_m(nv)) - &
          &          (T(nv)   - global_grid_block(i-1,j,k)%T(nv))/(dx_b(nv)*dx_m(nv))
        end associate
        associate(nv=>[2,4,6,8])
          dx_m(nv) = 0.5*(global_grid_block(i+1,j,k)%v(nv,1) - v(nv-1,1))
          dx_f(nv) =      global_grid_block(i+1,j,k)%v(nv,1) - v(nv,1)
          dx_b(nv) =      global_grid_block(i,j,k)%v(nv,1)   - v(nv-1,1)
          ddx2(nv) = ( global_grid_block(i+1,j,k)%T(nv) - T(nv) )/(dx_f(nv)*dx_m(nv)) - (T(nv) - T(nv-1))/(dx_b(nv)*dx_m(nv))
        end associate

        !calculate ddy2
        associate(nv=>[1,2,5,6])
          dy_m(nv)=0.5*(v(nv+2,2) - global_grid_block(i,j-1,k)%v(nv,2))
          dy_f(nv)=v(nv+2,2) - v(nv,2)
          dy_b(nv)=v(nv,2) - global_grid_block(i,j-1,k)%v(nv,2)
          ddy2(nv)=(T(nv+2)-T(nv))/(dy_f(nv)*dy_m(nv)) - (T(nv)-global_grid_block(i,j-1,k)%T(nv))/(dy_b(nv)*dy_m(nv))
        end associate

        associate(nv=>[3,4,7,8])
          dy_m(nv) =0.5*(global_grid_block(i,j+1,k)%v(nv,2) - v(nv-2,2))
          dy_f(nv) =global_grid_block(i,j+1,k)%v(nv,2) - v(nv,2)
          dy_b(nv) = v(nv,2) - v(nv-2,2)
          ddy2(nv) = (global_grid_block(i,j+1,k)%T(nv) - T(nv))/(dy_f(nv)*dy_m(nv))- (T(nv) - T(nv-2))/(dy_b(nv)*dy_m(nv))
        end associate

        !calculate ddz2
        associate(nv=>[1,2,3,4])
          dz_m(nv)=0.5*(v(nv+4,3) - global_grid_block(i,j,k-1)%v(nv,3))
          dz_f(nv)=v(nv+4,3) - v(nv,3)
          dz_b(nv)=v(nv,3) - global_grid_block(i,j,k-1)%v(nv,3)
          ddz2(nv)=(T(nv+4)- T(nv))/(dz_f(nv)*dz_m(nv)) - (T(nv) - global_grid_block(i,j,k-1)%T(nv))/(dz_b(nv)*dz_m(nv))
        end associate

        associate(nv=>[5,6,7,8])
          dz_m(nv)=0.5*(global_grid_block(i,j,k+1)%v(nv,3)-v(nv-4,3))
          dz_f(nv)=global_grid_block(i,j,k+1)%v(nv,3)-v(nv,3)
          dz_b(nv)=v(nv,3)-v(nv-4,3)
          ddz2(nv)=(global_grid_block(i,j,k+1)%T(nv)-T(nv))/(dz_f(nv)*dz_m(nv))- (T(nv)-T(nv-4))/(dz_b(nv)*dz_m(nv))
        end associate
      end associate
    end do
!  end block compute_derivatives

!  error_calculation: block
!    real(rkind)            :: avg_err, minimum_dx
!    integer(ikind)         :: nv
    minimum_dx=0.25/5.0
    avg_err=0.0
    do concurrent(k=2:nz-2,j=2:ny-2,i=2:nx-2,nv=1:8)
      avg_err=avg_err+abs(global_grid_block(i,j,k)%ddx2(nv)-2.0*global_grid_block(i,j,k)%v(nv,2))
    end do
    avg_err=avg_err/((nz-3)*(ny-3)*(nx-3)*8)
    if (avg_err <= minimum_dx**2) print *, "Test passed."
!  end block error_calculation

  call output_result(global_grid_block)

end program main
