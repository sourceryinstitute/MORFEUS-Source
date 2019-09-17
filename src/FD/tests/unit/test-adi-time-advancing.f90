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
    real(rkind), dimension(8,4) :: v
    real(rkind), dimension(8)   :: ddx2, ddy2, ddz2
    real(rkind), dimension(8)   :: T_analytical
  end type grid_block

  type(grid_block)  ,dimension(:,:,:) ,allocatable  :: global_grid_block
  integer(ikind)  :: i,j,k
  integer(ikind) ,parameter :: nx=41, ny=41, nz=41

  associate( p => position_vectors(nx,ny,nz) )

    allocate(global_grid_block(nx-1,ny-1,nz-1))

    associate( x=>(p(:,:,:,1)), y=>(p(:,:,:,2)), z=>(p(:,:,:,3)), T=>(p(:,:,:,4)) )
      do concurrent(k=1:nz-1, j=1:ny-1, i=1:nx-1)
        associate( v=>global_grid_block(i,j,k)%v )
          v(1,:) = [x(i,j,k)      , y(i,j,k)      , z(i,j,k)       , T(i,j,k)       ]
          v(2,:) = [x(i+1,j,k)    , y(i+1,j,k)    , z(i+1,j,k)     , T(i+1,j,k)     ]
          v(3,:) = [x(i,j+1,k)    , y(i,j+1,k)    , z(i,j+1,k)     , T(i,j+1,k)     ]
          v(4,:) = [x(i+1,j+1,k)  , y(i+1,j+1,k)  , z(i+1,j+1,k)   , T(i+1,j+1,k)   ]
          v(5,:) = [x(i,j,k+1)    , y(i,j,k+1)    , z(i,j,k+1)     , T(i,j,k+1)     ]
          v(6,:) = [x(i+1,j,k+1)  , y(i+1,j,k+1)  , z(i+1,j,k+1)   , T(i+1,j,k+1)   ]
          v(7,:) = [x(i,j+1,k+1)  , y(i,j+1,k+1)  , z(i,j+1,k+1)   , T(i,j+1,k+1)   ]
          v(8,:) = [x(i+1,j+1,k+1), y(i+1,j+1,k+1), z(i+1,j+1,k+1) , T(i+1,j+1,k+1) ]
          associate(ddx2=>global_grid_block(i,j,k)%ddx2, ddy2=>global_grid_block(i,j,k)%ddy2, &
                            ddz2=>global_grid_block(i,j,k)%ddz2 )
            ddx2(:) = 0.0
            ddy2(:) = 0.0
            ddz2(:) = 0.0
          end associate
        end associate
      end do
    end associate
  end associate

  time_advancing: block
    real(rkind)  :: dx_m, dy_m, dz_m
    real(rkind)  :: dx_f, dy_f, dz_f
    real(rkind)  :: dx_b, dy_b, dz_b
    real(rkind), dimension(:), allocatable :: a,b,c,d
    real(rkind)                            :: f,dt,t

    f=0.1
    dt=0.002
    t=0.0
    do while(t<=0.1)
    associate(n=>shape(global_grid_block))
      !x direction
      allocate(a(n(1)-1),b(n(1)-1),c(n(1)-1),d(n(1)-1))
      call get_ddy2(global_grid_block)
      call get_ddz2(global_grid_block)
      do k=2,n(3)
        do j=2,n(2)
          do i=2,n(1)
            associate(v=>global_grid_block(i,j,k)%v, ddy2=>global_grid_block(i,j,k)%ddy2, ddz2=>global_grid_block(i,j,k)%ddz2)
              dx_m = 0.5*(v(2,1) - global_grid_block(i-1,j,k)%v(1,1))
              dx_f =      v(2,1) - v(1,1)
              dx_b =      v(1,1)   - global_grid_block(i-1,j,k)%v(1,1)
              a(i-1)=-(3.0-2.0*f)/(dx_b*dx_m)
              b(i-1)=(3.0-2.0*f)/(dx_f*dx_m)+(3.0-2.0*f)/(dx_b*dx_m)+3.0/dt
              c(i-1)=-(3.0-2.0*f)/(dx_f*dx_m)
              d(i-1)=3.0*v(1,4)/dt+f*ddy2(1)+f*ddz2(1)
              if (i==2) d(i-1)=d(i-1)-a(i-1)*global_grid_block(i-1,j,k)%v(1,4)
              if (i==n(1)) d(i-1)=d(i-1)-c(i-1)*v(2,4)
            end associate
          end do
          associate(U=>tridiagonal_matrix_algorithm(a,b,c,d))
            do i=2,n(1)
              global_grid_block(i-1,j-1,k-1)%v(8,4)=U(i-1)
              global_grid_block(i-1,j-1,k)%v(4,4)=U(i-1)
              global_grid_block(i-1,j,k-1)%v(6,4)=U(i-1)
              global_grid_block(i-1,j,k)%v(2,4)=U(i-1)
              global_grid_block(i,j,k)%v(1,4)=U(i-1)
              global_grid_block(i,j-1,k-1)%v(7,4)=U(i-1)
              global_grid_block(i,j-1,k)%v(3,4)=U(i-1)
              global_grid_block(i,j,k-1)%v(5,4)=U(i-1)
            end do
          end associate
        end do
      end do
      deallocate(a,b,c,d)

      ! y direction
      allocate(a(n(2)-1),b(n(2)-1),c(n(2)-1),d(n(2)-1))
      call get_ddx2(global_grid_block)
      call get_ddz2(global_grid_block)
      do k=2,n(3)
        do i=2,n(1)
          do j=2,n(2)
            associate(v=>global_grid_block(i,j,k)%v, ddx2=>global_grid_block(i,j,k)%ddx2, ddz2=>global_grid_block(i,j,k)%ddz2)
              dy_m=0.5*(v(3,2) - global_grid_block(i,j-1,k)%v(1,2))
              dy_f=v(3,2) - v(1,2)
              dy_b=v(1,2) - global_grid_block(i,j-1,k)%v(1,2)
              a(j-1)=-(3.0-2.0*f)/(dy_b*dy_m)
              b(j-1)=(3.0-2.0*f)/(dy_f*dy_m)+(3.0-2.0*f)/(dy_b*dy_m)+3.0/dt
              c(j-1)=-(3.0-2.0*f)/(dy_f*dy_m)
              d(j-1)=3.0*v(1,4)/dt+f*ddx2(1)+f*ddz2(1)
              if (j==2) d(j-1)=d(j-1)-a(j-1)*global_grid_block(i,j-1,k)%v(1,4)
              if (j==n(2)) d(j-1)=d(j-1)-c(j-1)*v(3,4)
            end associate
          end do
          associate(U=>tridiagonal_matrix_algorithm(a,b,c,d))
          do j=2,n(2)
            global_grid_block(i-1,j-1,k-1)%v(8,4)=U(j-1)
            global_grid_block(i-1,j-1,k)%v(4,4)=U(j-1)
            global_grid_block(i-1,j,k-1)%v(6,4)=U(j-1)
            global_grid_block(i-1,j,k)%v(2,4)=U(j-1)
            global_grid_block(i,j,k)%v(1,4)=U(j-1)
            global_grid_block(i,j-1,k-1)%v(7,4)=U(j-1)
            global_grid_block(i,j-1,k)%v(3,4)=U(j-1)
            global_grid_block(i,j,k-1)%v(5,4)=U(j-1)
          end do
          end associate
        end do
      end do
      deallocate(a,b,c,d)

      ! z direction
      allocate(a(n(3)-1),b(n(3)-1),c(n(3)-1),d(n(3)-1))
      call get_ddx2(global_grid_block)
      call get_ddy2(global_grid_block)
      do j=2,n(2)
        do i=2,n(1)
          do k=2,n(3)
            associate(v=>global_grid_block(i,j,k)%v, ddx2=>global_grid_block(i,j,k)%ddx2, ddy2=>global_grid_block(i,j,k)%ddy2)
              dz_m=0.5*(v(5,3) - global_grid_block(i,j,k-1)%v(1,3))
              dz_f=v(5,3) - v(1,3)
              dz_b=v(1,3) - global_grid_block(i,j,k-1)%v(1,3)
              a(k-1)=-(3.0-2.0*f)/(dz_b*dz_m)
              b(k-1)=(3.0-2.0*f)/(dz_f*dz_m)+(3.0-2.0*f)/(dz_b*dz_m)+3.0/dt
              c(k-1)=-(3.0-2.0*f)/(dz_f*dz_m)
              d(k-1)=3.0*v(1,4)/dt+f*ddx2(1)+f*ddy2(1)
              if (k==2) d(k-1)=d(k-1)-a(k-1)*global_grid_block(i,j,k-1)%v(1,4)
              if (k==n(3)) d(k-1)=d(k-1)-c(k-1)*v(5,4)
            end associate
          end do
          associate(U=>tridiagonal_matrix_algorithm(a,b,c,d))
          do k=2,n(3)
            global_grid_block(i-1,j-1,k-1)%v(8,4)=U(k-1)
            global_grid_block(i-1,j-1,k)%v(4,4)=U(k-1)
            global_grid_block(i-1,j,k-1)%v(6,4)=U(k-1)
            global_grid_block(i-1,j,k)%v(2,4)=U(k-1)
            global_grid_block(i,j,k)%v(1,4)=U(k-1)
            global_grid_block(i,j-1,k-1)%v(7,4)=U(k-1)
            global_grid_block(i,j-1,k)%v(3,4)=U(k-1)
            global_grid_block(i,j,k-1)%v(5,4)=U(k-1)
          end do
          end associate
        end do
      end do
      deallocate(a,b,c,d)
    end associate
    t=t+dt
    end do
  end block time_advancing

  analytical_solution: block
    real(rkind), dimension(:,:,:), allocatable :: a_mnl, k_mnl
    integer(ikind)                             :: l,m,n
    real(rkind)                                :: pi, t

    pi=4.0*atan(1.0)
    t=0.1
    allocate(a_mnl(10,10,10), k_mnl(10,10,10))
    do concurrent(l=1:10, n=1:10, m=1:10)
      a_mnl(m,n,l)=-64.0/(pi**3*(2.0*m-1)*(2.0*n-1)*(2.0*l-1))*sin(0.5*(2.0*m-1)*pi)*sin(0.5*(2.0*n-1)*pi)*sin(0.5*(2.0*l-1)*pi)
      k_mnl(m,n,l)=(0.5*(2.0*m-1)*pi)**2+(0.5*(2.0*n-1)*pi)**2+(0.5*(2.0*l-1)*pi)**2
    end do
    do concurrent(k=1:nz-1,j=1:ny-1,i=1:nx-1)
      global_grid_block(i,j,k)%T_analytical(:)=2.0
    end do
    do concurrent(k=2:nz-1,j=2:ny-1,i=2:nx-1)
      do concurrent(l=1:10,n=1:10,m=1:10)
        global_grid_block(i,j,k)%T_analytical(1)=global_grid_block(i,j,k)%T_analytical(1)+a_mnl(m,n,l)*exp(-k_mnl(m,n,l)*t)* &
                                                  cos(0.5*(2.0*m-1)*pi*global_grid_block(i,j,k)%v(1,1))* &
                                                  cos(0.5*(2.0*n-1)*pi*global_grid_block(i,j,k)%v(1,2))* &
                                                  cos(0.5*(2.0*l-1)*pi*global_grid_block(i,j,k)%v(1,3))
      end do
    end do
    do concurrent(k=2:nz-1,j=2:ny-1, i=2:nx-1)
      global_grid_block(i-1,j-1,k-1)%T_analytical(8)=global_grid_block(i,j,k)%T_analytical(1)
      global_grid_block(i-1,j-1,k)%T_analytical(4)=global_grid_block(i,j,k)%T_analytical(1)
      global_grid_block(i-1,j,k-1)%T_analytical(6)=global_grid_block(i,j,k)%T_analytical(1)
      global_grid_block(i-1,j,k)%T_analytical(2)=global_grid_block(i,j,k)%T_analytical(1)
      global_grid_block(i,j,k)%T_analytical(1)=global_grid_block(i,j,k)%T_analytical(1)
      global_grid_block(i,j-1,k-1)%T_analytical(7)=global_grid_block(i,j,k)%T_analytical(1)
      global_grid_block(i,j-1,k)%T_analytical(3)=global_grid_block(i,j,k)%T_analytical(1)
      global_grid_block(i,j,k-1)%T_analytical(5)=global_grid_block(i,j,k)%T_analytical(1)
    end do
  end block analytical_solution

  error_calculation: block
    real(rkind)            :: avg_err_percentage
    real(rkind), parameter :: err_percentage=0.1
    integer(ikind)         :: nv
    avg_err_percentage=0.0
    do concurrent(k=1:nz-1,j=1:ny-1,i=1:nx-1,nv=1:8)
      avg_err_percentage=avg_err_percentage+abs((global_grid_block(i,j,k)%v(nv,4)-global_grid_block(i,j,k)%T_analytical(nv))/ &
              global_grid_block(i,j,k)%T_analytical(nv))
    end do
    avg_err_percentage=100*avg_err_percentage/((nz-1)*(ny-1)*(nx-1)*8)
    if (avg_err_percentage < 0.1) print *, "Test passed."
  end block error_calculation

  call output_result(global_grid_block)

contains

  pure function position_vectors(nx,ny,nz) result(vector_field)
    integer(ikind), intent(in) :: nx, ny, nz
    integer(ikind), parameter :: components=4
    real(rkind), parameter :: dx=0.6/(15.0), dy=0.6/(15.0)
    real(rkind), parameter :: dx_s=0.4/5.0, dy_s=0.4/5.0
    real(rkind), dimension(:,:,:,:) ,allocatable  :: vector_field

    allocate( vector_field(nx,ny,nz,components) )

    associate( dz=>2.0/(nz-1) )
      do concurrent( i=1:nx, j=1:ny, k=1:nz )
        vector_field(i,j,k,1) = merge(-0.6+(i-6)*dx, merge( -1.0+(i-1)*dx_s, 0.6+(i-36)*dx_s, i<6), i>=6 .and. i<=36)
        vector_field(i,j,k,2) = merge(-0.6+(j-6)*dy, merge( -1.0+(j-1)*dy_s, 0.6+(j-36)*dy_s, j<6), j>=6 .and. j<=36)
        vector_field(i,j,k,3) = -1.0+(k-1)*dz
        vector_field(i,j,k,4) = merge(2.0, 1.0, (i==1 .or. i==nx .or. j==1 .or. j==ny .or. k==1 .or. k==nz))
      end do
    end associate
  end function

  function tridiagonal_matrix_algorithm(a,b,c,d) result(x)
    real(rkind), dimension(:), intent(inout)  :: a,b,c,d
    real(rkind), dimension(:), allocatable :: x
    integer(ikind)                         :: i,n
    real(rkind)                            :: w
    n= size(b)
    allocate(x(n))
    do i=2,n
      w=a(i)/b(i-1)
      b(i)=b(i)-w*c(i-1)
      d(i)=d(i)-w*d(i-1)
    end do
    x(n)=d(n)/b(n)
    do i=n-1,1,-1
      x(i)=(d(i)-c(i)*x(i+1))/b(i)
    end do
  end function

  subroutine get_ddx2(this)
    type(grid_block), dimension(:,:,:), intent(inout)     :: this
    integer(ikind)                                        :: i, j, k, nv
    real(rkind), dimension(8)                             :: dx_f, dx_b, dx_m
    associate( n=>shape(this) )
    do concurrent(k=2:n(3)-1, j=2:n(2)-1, i=2:n(1)-1)
      associate(nv=>[1,3,5,7])
        dx_m(nv)=0.5*(this(i,j,k)%v(nv+1,1)-this(i-1,j,k)%v(nv,1))
        dx_f(nv)=this(i,j,k)%v(nv+1,1)-this(i,j,k)%v(nv,1)
        dx_b(nv)=this(i,j,k)%v(nv,1)-this(i-1,j,k)%v(nv,1)
        this(i,j,k)%ddx2(nv)=(this(i,j,k)%v(nv+1,4) - this(i,j,k)%v(nv,4))/(dx_f(nv)*dx_m(nv)) - &
                              (this(i,j,k)%v(nv,4) - this(i-1,j,k)%v(nv,4))/(dx_b(nv)*dx_m(nv))
      end associate
      associate(nv=>[2,4,6,8])
        dx_m(nv)=0.5*(this(i+1,j,k)%v(nv,1)-this(i,j,k)%v(nv-1,1))
        dx_f(nv)=this(i+1,j,k)%v(nv,1)-this(i,j,k)%v(nv,1)
        dx_b(nv)=this(i,j,k)%v(nv,1)-this(i,j,k)%v(nv-1,1)
        this(i,j,k)%ddx2(nv)=(this(i+1,j,k)%v(nv,4) - this(i,j,k)%v(nv,4))/(dx_f(nv)*dx_m(nv)) - &
                              (this(i,j,k)%v(nv,4) - this(i,j,k)%v(nv-1,4))/(dx_b(nv)*dx_m(nv))
      end associate
    end do
    end associate
  end subroutine

  subroutine get_ddy2(this)
    type(grid_block), dimension(:,:,:), intent(inout)     :: this
    integer(ikind)                                        :: i, j, k, nv
    real(rkind), dimension(8)                             :: dy_f, dy_b, dy_m
    associate( n=>shape(this) )
      do concurrent(k=2:n(3)-1, j=2:n(2)-1, i=2:n(1)-1)
        associate(nv=>[1,2,5,6])
          dy_m(nv)=0.5*(this(i,j,k)%v(nv+2,2) - this(i,j-1,k)%v(nv,2))
          dy_f(nv)=this(i,j,k)%v(nv+2,2) - this(i,j,k)%v(nv,2)
          dy_b(nv)=this(i,j,k)%v(nv,2) - this(i,j-1,k)%v(nv,2)
          this(i,j,k)%ddy2(nv)=(this(i,j,k)%v(nv+2,4) - this(i,j,k)%v(nv,4))/(dy_f(nv)*dy_m(nv)) - &
                                (this(i,j,k)%v(nv,4) - this(i,j-1,k)%v(nv,4))/(dy_b(nv)*dy_m(nv))
        end associate
        associate(nv=>[3,4,7,8])
          dy_m(nv) =0.5*(this(i,j+1,k)%v(nv,2) - this(i,j,k)%v(nv-2,2))
          dy_f(nv) =this(i,j+1,k)%v(nv,2) - this(i,j,k)%v(nv,2)
          dy_b(nv) =this(i,j,k)%v(nv,2) - this(i,j,k)%v(nv-2,2)
          this(i,j,k)%ddy2(nv)=(this(i,j+1,k)%v(nv,4) - this(i,j,k)%v(nv,4))/(dy_f(nv)*dy_m(nv)) - &
                                (this(i,j,k)%v(nv,4) - this(i,j,k)%v(nv-2,4))/(dy_b(nv)*dy_m(nv))
        end associate
      end do
    end associate
  end subroutine

  subroutine get_ddz2(this)
    type(grid_block), dimension(:,:,:), intent(inout)     :: this
    integer(ikind)                                        :: i, j, k, nv
    real(rkind), dimension(8)                             :: dz_f, dz_b, dz_m
    associate( n=>shape(this) )
      do concurrent(k=2:n(3)-1, j=2:n(2)-1, i=2:n(1)-1)
        associate(nv=>[1,2,3,4])
          dz_m(nv)=0.5*(this(i,j,k)%v(nv+4,3) - this(i,j,k-1)%v(nv,3))
          dz_f(nv)=this(i,j,k)%v(nv+4,3) - this(i,j,k)%v(nv,3)
          dz_b(nv)=this(i,j,k)%v(nv,3) - this(i,j,k-1)%v(nv,3)
          this(i,j,k)%ddz2(nv)=(this(i,j,k)%v(nv+4,4) - this(i,j,k)%v(nv,4))/(dz_f(nv)*dz_m(nv)) - &
                                (this(i,j,k)%v(nv,4) - this(i,j,k-1)%v(nv,4))/(dz_b(nv)*dz_m(nv))
        end associate
        associate(nv=>[5,6,7,8])
          dz_m(nv)=0.5*(this(i,j,k+1)%v(nv,3)-this(i,j,k)%v(nv-4,3))
          dz_f(nv)=this(i,j,k+1)%v(nv,3)-this(i,j,k)%v(nv,3)
          dz_b(nv)=this(i,j,k)%v(nv,3)-this(i,j,k)%v(nv-4,3)
          this(i,j,k)%ddz2(nv)=(this(i,j,k+1)%v(nv,4) - this(i,j,k)%v(nv,4))/(dz_f(nv)*dz_m(nv)) - &
                                (this(i,j,k)%v(nv,4) - this(i,j,k)%v(nv-4,4))/(dz_b(nv)*dz_m(nv))
        end associate
      end do
    end associate
  end subroutine

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
                write (file_unit,*) grid_blocks(i,j,k)%v(nv,4)
              case("v")
                write(file_unit, *) (grid_blocks(i,j,k)%v(nv,ndimension), ndimension=1,3)
              case("analytical_T")
                write (file_unit,*) grid_blocks(i,j,k)%T_analytical(nv)
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

    open(newunit=file_unit, file="3d-htc-grid.vtk")
    write(file_unit, '(a)') '# vtk DataFile Version 3.0'
    write(file_unit, '(a)') '3d-htc voxels'
    write(file_unit, '(a)') 'ASCII'
    write(file_unit, '(a)') 'DATASET UNSTRUCTURED_GRID'
    write(file_unit, '(a,i10,a)') "POINTS ", 8*(nx-1)*(ny-1)*(nz-1), " double"
    call write_component( grid_blocks, "v", file_unit)

    block
      integer i, j, k, nv
      write(file_unit,'(a,i10,3x,i10)') "CELLS ", (nx-1)*(ny-1)*(nz-1), 9*(nx-1)*(ny-1)*(nz-1)
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

    write(file_unit,'(a,I10)') "CELL_TYPES ", (nx-1)*(ny-1)*(nz-1)
    call write_component( grid_blocks, "string", file_unit, string="11")

    write(file_unit,'(a, I10)') 'POINT_DATA ', 8*(nx-1)*(ny-1)*(nz-1)
    write(file_unit,'(a)') 'SCALARS Temperature double 1'
    write(file_unit,'(a)') 'LOOKUP_TABLE default'
    call write_component( grid_blocks, "T", file_unit)

    write(file_unit,'(a)') 'SCALARS analytical_T double 1'
    write(file_unit,'(a)') 'LOOKUP_TABLE default'
    call write_component( grid_blocks, "analytical_T", file_unit)

    close(file_unit)

    block
      integer k, error_file_unit
      open(newunit=error_file_unit, file="error-order.dat")
      do k=1, nz-1
        write(error_file_unit, *) &
          grid_blocks(20,20,k)%v(1,4)-grid_blocks(20,20,k)%T_analytical(1)
      end do
      close(error_file_unit)
    end block
  end subroutine output_result

end program
