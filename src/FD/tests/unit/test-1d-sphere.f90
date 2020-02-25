!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module spherical_1D_solver_interface
  !! author: Xiaofeng Xu, Damian Rouson, and Karla Morris
  !! date: 2/24/2020
  !!
  !! Solve the 1D heat equation in spherically symmetric radial coordinates
  use kind_parameters, only : r8k, i4k
  implicit none

  private
  public :: grid_block

  type grid_block
    private
    real(r8k), allocatable :: v(:,:)
    real(r8k), allocatable :: ddr2(:)
    real(r8k), allocatable :: rho(:), cp(:)
    real(r8k), allocatable :: T_analytical(:)
  contains
    procedure :: set_v
    procedure :: set_ddr2_size
    procedure :: set_material_properties_size
    procedure :: set_expected_solution_size
    procedure :: set_rho
    procedure :: set_cp
    procedure :: run_test
  end type grid_block

  interface

    module subroutine set_v( this, nr, constants )
      implicit none
      class(grid_block), intent(inout) :: this
      integer, intent(in) :: nr
      real(r8k), intent(in) :: constants(:)
    end subroutine

    module subroutine set_ddr2_size(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine set_material_properties_size(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine set_expected_solution_size(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine set_rho(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine set_cp(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

    module subroutine run_test(this)
      implicit none
      class(grid_block), intent(inout) :: this
    end subroutine

  end interface

end module spherical_1D_solver_interface

submodule(spherical_1D_solver_interface) spherical_1D_solver_implementation
  use assertions_interface, only : assert
  implicit none

contains

  module procedure set_v

    integer i

    allocate(this%v(nr,2), source = constants(1))

    associate( dr => 1._r8k/(nr-1) )
      associate( radial_nodes => [( sqrt((i-1)*dr)*0.2, i = 1, nr )] )
        this%v(:,1) = radial_nodes
      end associate
    end associate
    this%v(:,2) = constants(2)

  end procedure

  module procedure set_ddr2_size
    call assert( allocated(this%v), "grid_block%set_ddr2_size: allocated(this%v)")
    associate( nr => size(this%v,1) )
      allocate(this%ddr2(nr))
    end associate
  end procedure

  module procedure set_material_properties_size
    call assert( allocated(this%v), "grid_block%set_material_properties_size: allocated(this%v)")
    associate( nr => size(this%v,1) )
      allocate(this%rho(nr), this%cp(nr))
    end associate
  end procedure

  module procedure set_expected_solution_size
    call assert( allocated(this%v), "grid_block%set_expected_solution_size: allocated(this%v)")
    associate( nr => size(this%v,1) )
      allocate( this%T_analytical(nr) )
    end associate
  end procedure

  module procedure run_test
    real(r8k), parameter :: h=230, Tb=293.15
    integer(i4k)  :: i

    call assert( allocated(this%v), "grid_block%run_test: allocated(this%v)")

    associate( nr => size(this%v,1) )

      time_advancing: block
        real(r8k)  :: dr_m, dz_m
        real(r8k)  :: dr_f, dz_f
        real(r8k)  :: dr_b, dz_b
        real(r8k), dimension(:), allocatable :: a,b,c,d
        real(r8k)                            :: dt,t,e, rf, rb

        dt=0.1
        t=0.0

        do while(t<1000)

          associate(nr=>size(this%v,1))

            allocate(a(nr),b(nr),c(nr),d(nr))

            call this%set_rho()
            call this%set_cp()

            i=1
            e=1.0/(this%rho(i)*this%cp(i))
            dr_f=this%v(i+1,1) - this%v(i,1)
            rf=0.5*(this%v(i+1,1) + this%v(i,1))
            b(i)=6.0*e*kc(rf)/(dr_f**2)+1.0/dt
            c(i)=-6.0*e*kc(rf)/(dr_f**2)
            d(i)=this%v(i,2)/dt

            do concurrent( i=2:nr-1 )
                e=1.0/(this%rho(i)*this%cp(i))
                dr_f=this%v(i+1,1) - this%v(i,1)
                rf=0.5*(this%v(i+1,1) + this%v(i,1))
                dr_b=this%v(i,1) - this%v(i-1,1)
                dr_m = 0.5*(this%v(i+1,1) - this%v(i-1,1))
                rb=0.5*(this%v(i,1) + this%v(i-1,1))
                a(i)=-1.0*e*kc(rb)*rb**2/(dr_b*dr_m*this%v(i,1)**2)
                b(i)=e*kc(rf)*rf**2/(dr_f*dr_m*this%v(i,1)**2)+ &
                      e*kc(rb)*rb**2/(dr_b*dr_m*this%v(i,1)**2)+1.0/dt
                c(i)=-1.0*e*kc(rf)*rf**2/(dr_f*dr_m*this%v(i,1)**2)
                d(i)=this%v(i,2)/dt
            end do

            i=nr
            e=1.0/(this%rho(i)*this%cp(i))
            dr_b=this%v(i,1) - this%v(i-1,1)
            dr_m = dr_b
            rb=0.5*(this%v(i,1) + this%v(i-1,1))
            rf=this%v(i,1)+0.5*dr_b
            a(i)=-1.0*e*kc(rb)*rb**2/(dr_b*dr_m*this%v(i,1)**2)
            b(i)=e*h*rf**2/(dr_m*this%v(i,1)**2)+ &
                  e*kc(rb)*rb**2/(dr_b*dr_m*this%v(i,1)**2)+1.0/dt
            d(i)=this%v(i,2)/dt+e*h*Tb*rf**2/(dr_m*this%v(i,1)**2)

            associate(U=>tridiagonal_matrix_algorithm(a,b,c,d))
              call assert(size(U)==nr, "size(U)")
              do i=1,nr
                this%v(i,2)=U(i)
              end do
            end associate
            deallocate(a,b,c,d)
          end associate
          t=t+dt
        end do
      end block time_advancing

      analytical_solution: block
        real(r8k), dimension(20)             :: mu
        real(r8k)                            :: t,r_0, R, pi
        real(r8k)                            :: T0, T_inf
        integer(i4k)                         :: i,n

        T0=1073.15
        T_inf=293.15
        t=1000.0
        r_0=0.2
        pi=4.0*atan(1.0)
        do i=1, 20
          mu(i)=(2*i-1)*pi/2.0
        end do

        do i=1,nr
          R=this%v(i,1)/r_0
          this%T_analytical(i)=0.0
          if (i==1) then
            do n=1, 20
              this%T_analytical(i)=this%T_analytical(i)+ &
                    2.0*(sin(mu(n))-mu(n)*cos(mu(n)))/(mu(n)-sin(mu(n))*cos(mu(n)))* &
                    exp(-1.0e-5*(mu(n)/r_0)**2*t)
            end do
            this%T_analytical(i)=(T0-T_inf)*this%T_analytical(i)+T_inf
          else
            do n=1, 6
              this%T_analytical(i)=this%T_analytical(i)+ &
                    2.0*(sin(mu(n))-mu(n)*cos(mu(n)))/(mu(n)-sin(mu(n))*cos(mu(n)))* &
                    sin(mu(n)*R)/(mu(n)*R)*exp(-1.0e-5*(mu(n)/r_0)**2*t)
            end do
            this%T_analytical(i)=(T0-T_inf)*this%T_analytical(i)+T_inf
          end if
        end do
      end block analytical_solution

    error_calculation: block
      real(r8k)            :: avg_err_percentage
      real(r8k), parameter :: err_percentage=0.2
      avg_err_percentage=0.0
      do i=1,nr
        avg_err_percentage=avg_err_percentage+100.0*(abs(this%v(i,2)- &
                this%T_analytical(i))/this%T_analytical(i))
      end do
      avg_err_percentage=avg_err_percentage/nr
      if (avg_err_percentage <= err_percentage) print *, "Test passed."
    end block error_calculation

    end associate

  end procedure run_test


  function tridiagonal_matrix_algorithm(a,b,c,d) result(x)
    real(r8k), dimension(:), intent(inout)  :: a,b,c,d
    real(r8k), dimension(:), allocatable :: x
    integer(i4k)                         :: i,n
    real(r8k)                            :: w
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

  pure real(r8k) function kc(x)
    real(r8k), intent(in)   :: x
    kc=46.0
  end function

  module procedure set_cp
    this%cp=4600.0
  end procedure

  module procedure set_rho
    this%rho=1000.0
  end procedure

  subroutine set_ddr2(this)
    class(grid_block), intent(inout)     :: this
    integer(i4k)                      :: i
    real(r8k)                         :: dr_f, dr_b, dr_m
    real(r8k)                         :: rf, rb
    associate( nr=>size(this%v,1) )
      do i=2,nr-1
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

end submodule spherical_1D_solver_implementation

program main
  use  spherical_1D_solver_interface, only : grid_block
  use kind_parameters, only : r8k
  implicit none

  type(grid_block) global_grid_block

  call global_grid_block%set_v( nr = 101, constants = [0._r8k, 1073.15_r8k] )
  call global_grid_block%set_ddr2_size()
  call global_grid_block%set_expected_solution_size()
  call global_grid_block%set_material_properties_size()
  call global_grid_block%run_test()

end program
