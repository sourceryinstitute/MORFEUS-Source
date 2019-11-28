!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
module spatial_derivative_procedures
  implicit none
  integer ,parameter :: digits=8  ! num. digits of kind
  integer ,parameter :: decades=9 ! num. representable decades
  integer ,parameter :: rkind = selected_real_kind(digits)
  integer ,parameter :: ikind = selected_int_kind(decades)

  type grid_block
    real(rkind), allocatable :: v(:,:,:,:)
    real(rkind), allocatable :: ddx2(:,:,:,:)
  end type grid_block

  integer(ikind)  :: i,j,k
  integer(ikind) ,parameter :: nx=21, ny=25, nz=30

  type(grid_block) :: global_grid_block

contains

  impure function position_vectors(nx,ny,nz) result(vector_field)
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

  subroutine set_scalar_fluxes(this)
    type(grid_block), intent(inout) :: this
    integer(ikind) i, j, k

    associate( n=>shape(this%v), s=>this%v(:,:,:,4) )
      associate( x=>this%v(:,:,:,1))
        do concurrent(k=1:n(3), j=1:n(2), i=2:n(1)-1)
          associate( &
            dx_m => 0.5_rkind*(x(i+1,j,k) - x(i-1,j,k)), &
            dx_f =>  x(i+1,j,k) - x(i,j,k), &
            dx_b =>  x(i,j,k)-x(i-1,j,k) )
            this%ddx2(i,j,k,1) = &
              kc(0.5_rkind*(x(i+1,j,k) + x(i,j,k))  )*(s(i+1,j,k) - s(i,j,k)  )/(dx_f*dx_m) - &
              kc(0.5_rkind*(x(i,j,k)   + x(i-1,j,k)))*(s(i,j,k)   - s(i-1,j,k))/(dx_b*dx_m)
          end associate
        end do
      end associate

      associate( y=>this%v(:,:,:,2) )
        do concurrent(k=1:n(3), j=2:n(2)-1, i=1:n(1))
          associate( &
            dy_m => 0.5_rkind*(y(i,j+1,k) - y(i,j-1,k)), &
            dy_f =>     y(i,j+1,k) - y(i,j,k ), &
            dy_b =>     y(i,j,k)   - y(i,j-1,k) )
            this%ddx2(i,j,k,2) = &
              kc(0.5_rkind*(y(i,j+1,k)+y(i,j,k)))*(s(i,j+1,k) - s(i,j,k))/(dy_f*dy_m) - &
              kc(0.5_rkind*(y(i,j,k)+y(i,j-1,k)))*(s(i,j,k) - s(i,j-1,k))/(dy_b*dy_m)
           end associate
        end do
      end associate

      associate( z=>this%v(:,:,:,3) )
        do concurrent(k=2:n(3)-1, j=1:n(2), i=1:n(1))
          associate( &
            dz_m => 0.5*(z(i,j,k+1)-z(i,j,k-1)), &
            dz_f => z(i,j,k+1)-z(i,j,k), &
            dz_b => z(i,j,k)-z(i,j,k-1) )
            this%ddx2(i,j,k,3) = &
              kc(0.5*(z(i,j,k+1)+z(i,j,k)))*(s(i,j,k+1) - s(i,j,k))/(dz_f*dz_m) - &
              kc(0.5*(z(i,j,k)+z(i,j,k-1)))*(s(i,j,k) - s(i,j,k-1))/(dz_b*dz_m)
          end associate
        end do
      end associate
    end associate
  end subroutine

end module spatial_derivative_procedures

program main
  use spatial_derivative_procedures
  implicit none
  real(rkind) ::  minimum_dx=0.25/5.0
  allocate(global_grid_block%v(nx,ny,nz,4), global_grid_block%ddx2(nx,ny,nz,3), source=0._rkind)

  associate( p => position_vectors(nx,ny,nz) )
    associate( x => p(:,:,:,1), y => p(:,:,:,2), z => p(:,:,:,3) )
      global_grid_block%v(:,:,:,1) = x
      global_grid_block%v(:,:,:,2) = y
      global_grid_block%v(:,:,:,3) = z
      global_grid_block%v(:,:,:,4) = x**2 * y ! T = x^2*y => (d^2/dx^2)x = 2x, (d^2/dy^2)y = 2y
    end associate
  end associate

  call set_scalar_fluxes(global_grid_block)

  associate( y => global_grid_block%v(2:nx-1,:,:,2) )
    associate( exact_answer => 2*y, approximation =>  global_grid_block%ddx2(2:nx-1,:,:,1) )
      associate( avg_err => sum ( abs( approximation - exact_answer ) )/ SIZE(exact_answer) )
        if (avg_err <= minimum_dx**2) print *, "Test passed."
      end associate
    end associate
  end associate

  end program main
