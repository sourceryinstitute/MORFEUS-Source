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
    real(rkind), allocatable :: v(:,:)
    real(rkind), allocatable :: ddr2(:)
    real(rkind), allocatable :: rho(:), cp(:)
    real(rkind), allocatable :: T_analytical(:)
  end type grid_block

  integer(ikind)  :: i
  integer(ikind) ,parameter :: nr=101

  real(rkind), parameter :: h=230, Tb=293.15

  type(grid_block) :: global_grid_block
  allocate(global_grid_block%v(nr,2), global_grid_block%ddr2(nr))
  allocate(global_grid_block%rho(nr), global_grid_block%cp(nr), global_grid_block%T_analytical(nr))
  global_grid_block%v(:,:)=0.0
  global_grid_block%ddr2(:)=0.0

  associate( p => position_vectors(nr) )
    global_grid_block%v(:,1) = p(:,1)
  end associate
  global_grid_block%v(:,2) = 1073.15

  time_advancing: block
    real(rkind)  :: dr_m, dz_m
    real(rkind)  :: dr_f, dz_f
    real(rkind)  :: dr_b, dz_b
    real(rkind), dimension(:), allocatable :: a,b,c,d
    real(rkind)                            :: dt,t,e, rf, rb

    dt=0.1
    t=0.0
    do while(t<1000)
      associate(n=>shape(global_grid_block%v))
        allocate(a(n(1)),b(n(1)),c(n(1)),d(n(1)))
        call get_rho(global_grid_block)
        call get_cp(global_grid_block)
        do i=1,n(1)
          e=1.0/(global_grid_block%rho(i)*global_grid_block%cp(i))
          if (i==1) then
            dr_f=global_grid_block%v(i+1,1) - global_grid_block%v(i,1)
            rf=0.5*(global_grid_block%v(i+1,1) + global_grid_block%v(i,1))
            b(i)=6.0*e*kc(rf)/(dr_f**2)+1.0/dt
            c(i)=-6.0*e*kc(rf)/(dr_f**2)
            d(i)=global_grid_block%v(i,2)/dt
          else if (i==n(1)) then
            dr_b=global_grid_block%v(i,1) - global_grid_block%v(i-1,1)
            dr_m = dr_b
            rb=0.5*(global_grid_block%v(i,1) + global_grid_block%v(i-1,1))
            rf=global_grid_block%v(i,1)+0.5*dr_b
            a(i)=-1.0*e*kc(rb)*rb**2/(dr_b*dr_m*global_grid_block%v(i,1)**2)
            b(i)=e*h*rf**2/(dr_m*global_grid_block%v(i,1)**2)+ &
                  e*kc(rb)*rb**2/(dr_b*dr_m*global_grid_block%v(i,1)**2)+1.0/dt
            d(i)=global_grid_block%v(i,2)/dt+e*h*Tb*rf**2/(dr_m*global_grid_block%v(i,1)**2)
          else
            dr_f=global_grid_block%v(i+1,1) - global_grid_block%v(i,1)
            rf=0.5*(global_grid_block%v(i+1,1) + global_grid_block%v(i,1))
            dr_b=global_grid_block%v(i,1) - global_grid_block%v(i-1,1)
            dr_m = 0.5*(global_grid_block%v(i+1,1) - global_grid_block%v(i-1,1))
            rb=0.5*(global_grid_block%v(i,1) + global_grid_block%v(i-1,1))
            a(i)=-1.0*e*kc(rb)*rb**2/(dr_b*dr_m*global_grid_block%v(i,1)**2)
            b(i)=e*kc(rf)*rf**2/(dr_f*dr_m*global_grid_block%v(i,1)**2)+ &
                  e*kc(rb)*rb**2/(dr_b*dr_m*global_grid_block%v(i,1)**2)+1.0/dt
            c(i)=-1.0*e*kc(rf)*rf**2/(dr_f*dr_m*global_grid_block%v(i,1)**2)
            d(i)=global_grid_block%v(i,2)/dt
          end if
        end do

        associate(U=>tridiagonal_matrix_algorithm(a,b,c,d))
          do i=1,n(1)
            global_grid_block%v(i,2)=U(i)
          end do
        end associate
        deallocate(a,b,c,d)
      end associate
      t=t+dt
    end do
  end block time_advancing

  analytical_solution: block
    real(rkind), dimension(20)             :: mu
    real(rkind)                            :: t,r_0, R, pi
    real(rkind)                            :: T0, T_inf
    integer(ikind)                         :: i,n

    T0=1073.15
    T_inf=293.15
    t=1000.0
    r_0=0.2
    pi=4.0*atan(1.0)
    do i=1, 20
      mu(i)=(2*i-1)*pi/2.0
    end do

    do i=1,nr
      R=global_grid_block%v(i,1)/r_0
      global_grid_block%T_analytical(i)=0.0
      if (i==1) then
        do n=1, 20
          global_grid_block%T_analytical(i)=global_grid_block%T_analytical(i)+ &
                2.0*(sin(mu(n))-mu(n)*cos(mu(n)))/(mu(n)-sin(mu(n))*cos(mu(n)))* &
                exp(-1.0e-5*(mu(n)/r_0)**2*t)
        end do
        global_grid_block%T_analytical(i)=(T0-T_inf)*global_grid_block%T_analytical(i)+T_inf
      else
        do n=1, 6
          global_grid_block%T_analytical(i)=global_grid_block%T_analytical(i)+ &
                2.0*(sin(mu(n))-mu(n)*cos(mu(n)))/(mu(n)-sin(mu(n))*cos(mu(n)))* &
                sin(mu(n)*R)/(mu(n)*R)*exp(-1.0e-5*(mu(n)/r_0)**2*t)
        end do
        global_grid_block%T_analytical(i)=(T0-T_inf)*global_grid_block%T_analytical(i)+T_inf
      end if
    end do
  end block analytical_solution


  error_calculation: block
    real(rkind)            :: avg_err_percentage
    real(rkind), parameter :: err_percentage=0.2
    avg_err_percentage=0.0
    do i=1,nr
      avg_err_percentage=avg_err_percentage+100.0*(abs(global_grid_block%v(i,2)- &
              global_grid_block%T_analytical(i))/global_grid_block%T_analytical(i))
    end do
    avg_err_percentage=avg_err_percentage/nr
    if (avg_err_percentage <= err_percentage) print *, "Test passed."
  end block error_calculation

contains

  pure function position_vectors(nr) result(vector_field)
    integer(ikind), intent(in) :: nr
    integer(ikind), parameter :: components=1
    real(rkind)               :: dx
    real(rkind), dimension(:,:) ,allocatable  :: vector_field
    integer(ikind)            :: i

    dx=1.0/(nr-1)
    allocate( vector_field(nr,components) )
    do i=1,nr
      vector_field(i,1) = sqrt((i-1)*dx)*0.2
      !vector_field(i,1) = (i-1)*dx
    end do
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

  pure real(rkind) function kc(x)
    real(rkind), intent(in)   :: x
    !kc=5.0*exp(3.0*x)
    kc=46.0
  end function

  subroutine get_cp(this)
    type(grid_block), intent(inout)     :: this
    integer(ikind)                      :: i
    associate( n=>shape(this%v) )
      do i=1, n(1)
        !this%cp(i)=exp(3.0*this%v(i,1))
        this%cp(i)=4600.0
      end do
    end associate
  end subroutine

  subroutine get_rho(this)
    type(grid_block), intent(inout)     :: this
    integer(ikind)                      :: i
    associate( n=>shape(this%v) )
      do i=1,n(1)
        this%rho(i)=1000.0
      end do
    end associate
  end subroutine

  subroutine get_ddr2(this)
    type(grid_block), intent(inout)     :: this
    integer(ikind)                      :: i
    real(rkind)                         :: dr_f, dr_b, dr_m
    real(rkind)                         :: rf, rb
    associate( n=>shape(this%v) )
      do i=2,n(1)-1
        dr_m=0.5*(this%v(i+1,1)-this%v(i-1,1))
        dr_f=this%v(i+1,1)-this%v(i,1)
        dr_b=this%v(i,1)-this%v(i-1,1)
        rf=0.5*(this%v(i+1,1)+this%v(i,1))
        rb=0.5*(this%v(i,1)+this%v(i-1,1))
        this%ddr2(i)=kc(rf)*rf**2*(this%v(i+1,2) - this%v(i,2))/(dr_f*dr_m) - &
                      kc(rb)*rb**2*(this%v(i,2) - this%v(i-1,2))/(dr_b*dr_m)
      end do
    end associate
  end subroutine

end program
