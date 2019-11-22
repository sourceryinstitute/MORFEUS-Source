!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  implicit none
  integer ,parameter :: digits=8  ! num. digits of kind
  integer ,parameter :: decades=9 ! num. representable decades
  integer ,parameter :: rkind = selected_real_kind(digits)
  integer ,parameter :: ikind = selected_int_kind(decades)

  type grid_block
    real(rkind), allocatable :: v(:,:,:,:)
    real(rkind), allocatable :: ddx2(:,:,:), ddy2(:,:,:), ddz2(:,:,:)
  end type grid_block

  integer(ikind)  :: i,j,k
  integer(ikind) ,parameter :: nx=21, ny=21, nz=21

  type(grid_block) :: global_grid_block
  allocate(global_grid_block%v(nx,ny,nz,4), global_grid_block%ddx2(nx,ny,nz), &
           global_grid_block%ddy2(nx,ny,nz), global_grid_block%ddz2(nx,ny,nz))
  global_grid_block%v(:,:,:,:)=0.0
  global_grid_block%ddx2(:,:,:)=0.0
  global_grid_block%ddy2(:,:,:)=0.0
  global_grid_block%ddz2(:,:,:)=0.0

  associate( p => position_vectors(nx,ny,nz) )
    global_grid_block%v(:,:,:,1) = p(:,:,:,1)
    global_grid_block%v(:,:,:,2) = p(:,:,:,2)
    global_grid_block%v(:,:,:,3) = p(:,:,:,3)
    global_grid_block%v(:,:,:,4) = global_grid_block%v(:,:,:,1)**2 * global_grid_block%v(:,:,:,2) ! T=x^2*y => d^2T/dx^2=y
  end associate

  call get_ddx2(global_grid_block)
  call get_ddy2(global_grid_block)
  call get_ddz2(global_grid_block)

  error_calculation: block
    real(rkind)            :: avg_err, minimum_dx
    minimum_dx=0.25/5.0
    avg_err=0.0
    do concurrent(k=1:nz,j=1:ny,i=2:nx-1)
      avg_err=avg_err+abs(global_grid_block%ddx2(i,j,k)-2.0*global_grid_block%v(i,j,k,2))
    end do
    avg_err=avg_err/((nz)*(ny)*(nx-2))
    if (avg_err <= minimum_dx**2) print *, "Test passed."
  end block error_calculation

  call output_result(global_grid_block)
contains

  pure function position_vectors(nx,ny,nz) result(vector_field)
    integer(ikind), intent(in) :: nx, ny, nz
    integer(ikind), parameter :: components=3
    real(rkind), parameter :: dx=0.75/(5.0), dy=0.75/(5.0)
    real(rkind), parameter :: dx_s=0.25/5.0, dy_s=0.25/5.0
    real(rkind), dimension(:,:,:,:) ,allocatable  :: vector_field

    allocate( vector_field(nx,ny,nz,components) )

    associate( dz=>10.0/(nz-1) )
      do concurrent( i=1:nx, j=1:ny, k=1:nz )
        vector_field(i,j,k,1) = merge(-0.75+(i-6)*dx, merge( -1.0+(i-1)*dx_s, 0.75+(i-16)*dx_s, i<6), i>=6 .and. i<=16)
        vector_field(i,j,k,2) = merge(-0.75+(j-6)*dy, merge( -1.0+(j-1)*dy_s, 0.75+(j-16)*dy_s, j<6), j>=6 .and. j<=16)
        vector_field(i,j,k,3) = (k-1)*dz
      end do
    end associate
  end function

  ! Below is a dummy conductivity function which assumed to a constant value, 1.0.
  ! User will need change this function if the conductivity is vary with Temetrature or space.
  ! kc in the deriveative calculation is assumed to be varing wih x for ddx2, varing with y for ddy2 and varing with z
  ! for ddz2;
  ! User need change kc term based on kc function in the derivative calculation.
  pure real(rkind) function kc(x)
    real(rkind), intent(in)   :: x
    kc=1.0
  end function

  subroutine get_ddx2(this)
    type(grid_block), intent(inout)     :: this
    integer(ikind)                      :: i, j, k
    real(rkind)                         :: dx_f, dx_b, dx_m
    associate( n=>shape(this%v) )
      do concurrent(k=1:n(3), j=1:n(2), i=2:n(1)-1)
        dx_m=0.5*(this%v(i+1,j,k,1)-this%v(i-1,j,k,1))
        dx_f=this%v(i+1,j,k,1)-this%v(i,j,k,1)
        dx_b=this%v(i,j,k,1)-this%v(i-1,j,k,1)
        this%ddx2(i,j,k)=kc(0.5*(this%v(i+1,j,k,1)+this%v(i,j,k,1)))*(this%v(i+1,j,k,4) - this%v(i,j,k,4))/(dx_f*dx_m) - &
                        kc(0.5*(this%v(i,j,k,1)+this%v(i-1,j,k,1)))*(this%v(i,j,k,4) - this%v(i-1,j,k,4))/(dx_b*dx_m)
      end do
    end associate
  end subroutine

  subroutine get_ddy2(this)
    type(grid_block), intent(inout)     :: this
    integer(ikind)                      :: i, j, k
    real(rkind)                         :: dy_f, dy_b, dy_m
    associate( n=>shape(this%v) )
      do concurrent(k=1:n(3), j=2:n(2)-1, i=1:n(1))
        dy_m=0.5*(this%v(i,j+1,k,2)-this%v(i,j-1,k,2))
        dy_f=this%v(i,j+1,k,2)-this%v(i,j,k,2)
        dy_b=this%v(i,j,k,2)-this%v(i,j-1,k,2)
        this%ddy2(i,j,k)=kc(0.5*(this%v(i,j+1,k,2)+this%v(i,j,k,2)))*(this%v(i,j+1,k,4) - this%v(i,j,k,4))/(dy_f*dy_m) - &
                        kc(0.5*(this%v(i,j,k,2)+this%v(i,j-1,k,2)))*(this%v(i,j,k,4) - this%v(i,j-1,k,4))/(dy_b*dy_m)
      end do
    end associate
  end subroutine

  subroutine get_ddz2(this)
    type(grid_block), intent(inout)     :: this
    integer(ikind)                      :: i, j, k
    real(rkind)                         :: dz_f, dz_b, dz_m
    associate( n=>shape(this%v) )
      do concurrent(k=2:n(3)-1, j=1:n(2), i=1:n(1))
        dz_m=0.5*(this%v(i,j,k+1,3)-this%v(i,j,k-1,3))
        dz_f=this%v(i,j,k+1,3)-this%v(i,j,k,3)
        dz_b=this%v(i,j,k,3)-this%v(i,j,k-1,3)
        this%ddz2(i,j,k)=kc(0.5*(this%v(i,j,k+1,3)+this%v(i,j,k,3)))*(this%v(i,j,k+1,4) - this%v(i,j,k,4))/(dz_f*dz_m) - &
                        kc(0.5*(this%v(i,j,k,3)+this%v(i,j,k-1,3)))*(this%v(i,j,k,4) - this%v(i,j,k-1,4))/(dz_b*dz_m)
      end do
    end associate
  end subroutine

  subroutine write_component(grid_blocks, component, file_unit, string)
    type(grid_block), intent(in) :: grid_blocks
    character(len=*), intent(in) :: component
    character(len=*), intent(in), optional :: string
    integer, intent(in) :: file_unit
    integer(ikind) :: i, j, k, ndimension
    associate( n=>shape(grid_blocks%v) )
    do k=1, n(3)-1
      do j=1, n(2)-1
        do i=1, n(1)-1
          select case(component)
            case("T")
              write (file_unit,*) grid_blocks%v(i,j,k,4)
              write (file_unit,*) grid_blocks%v(i+1,j,k,4)
              write (file_unit,*) grid_blocks%v(i,j+1,k,4)
              write (file_unit,*) grid_blocks%v(i+1,j+1,k,4)
              write (file_unit,*) grid_blocks%v(i,j,k+1,4)
              write (file_unit,*) grid_blocks%v(i+1,j,k+1,4)
              write (file_unit,*) grid_blocks%v(i,j+1,k+1,4)
              write (file_unit,*) grid_blocks%v(i+1,j+1,k+1,4)
            case("v")
              write(file_unit, *) (grid_blocks%v(i,j,k,ndimension), ndimension=1,3)
              write(file_unit, *) (grid_blocks%v(i+1,j,k,ndimension), ndimension=1,3)
              write(file_unit, *) (grid_blocks%v(i,j+1,k,ndimension), ndimension=1,3)
              write(file_unit, *) (grid_blocks%v(i+1,j+1,k,ndimension), ndimension=1,3)
              write(file_unit, *) (grid_blocks%v(i,j,k+1,ndimension), ndimension=1,3)
              write(file_unit, *) (grid_blocks%v(i+1,j,k+1,ndimension), ndimension=1,3)
              write(file_unit, *) (grid_blocks%v(i,j+1,k+1,ndimension), ndimension=1,3)
              write(file_unit, *) (grid_blocks%v(i+1,j+1,k+1,ndimension), ndimension=1,3)
            case("ddx2")
              write (file_unit,*) grid_blocks%ddx2(i,j,k)
              write (file_unit,*) grid_blocks%ddx2(i+1,j,k)
              write (file_unit,*) grid_blocks%ddx2(i,j+1,k)
              write (file_unit,*) grid_blocks%ddx2(i+1,j+1,k)
              write (file_unit,*) grid_blocks%ddx2(i,j,k+1)
              write (file_unit,*) grid_blocks%ddx2(i+1,j,k+1)
              write (file_unit,*) grid_blocks%ddx2(i,j+1,k+1)
              write (file_unit,*) grid_blocks%ddx2(i+1,j+1,k+1)
            case("string")
              write (file_unit,*) string
            case default
              error stop "write_component: unrecognized component"
          end select
        end do
      end do
    end do
    end associate
  end subroutine

  subroutine output_result(grid_blocks)
    type(grid_block), intent(in) :: grid_blocks
    integer file_unit

    open(newunit=file_unit, file="spatial-derivatives-new.vtk")
    write(file_unit, '(a)') '# vtk DataFile Version 3.0'
    write(file_unit, '(a)') '3d-htc voxels'
    write(file_unit, '(a)') 'ASCII'
    write(file_unit, '(a)') 'DATASET UNSTRUCTURED_GRID'
    write(file_unit, '(a,i5,a)') "POINTS ", 8*(nx-1)*(ny-1)*(nz-1), " double"
    call write_component( grid_blocks, "v", file_unit)

    block
      integer i, j, k, nv
      write(file_unit,'(a,i5,3x,i5)') "CELLS ", (nx-1)*(ny-1)*(nz-1), 9*(nx-1)*(ny-1)*(nz-1)
      associate( n=>shape(grid_blocks%v) )
        do k=1, n(3)-1
          do j=1, n(2)-1
            do i=1, n(1)-1
              associate( grid_block_id => 8*(i+(j-1)*(n(1)-1)+(k-1)*(n(1)-1)*(n(2)-1)-1) )
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
  end subroutine

end program
