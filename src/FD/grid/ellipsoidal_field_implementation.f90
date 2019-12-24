!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(ellipsoidal_field_interface) ellipsoidal_field_implementation
  use kind_parameters, only : r8k
  use assertions_interface, only : assert, max_errmsg_len, assertions
  implicit none

contains

  module procedure evaluate
    use cartesian_grid_interface, only : cartesian_grid
    real(r8k), allocatable, dimension(:,:,:,:) :: position_vectors
    real(r8k), allocatable :: f_value(:,:,:)
    real(r8k), parameter ::  x_center = 3*(0.25E-01) + 1.E-01/2., x_max = 2*x_center
    real(r8k), parameter ::  y_center = 0.5E-01 + 2*(0.25E-01) + 3.E-01/2., y_max = 2*y_center
    real(r8k), parameter ::  z_center = 20.E-01/2., z_max=2*z_center
    integer, parameter :: success=0
    integer alloc_stat
    character(len=max_errmsg_len) error_message

    select type(grid_points)
      class is(cartesian_grid)
        position_vectors = grid_points%vectors()

        associate( x=>position_vectors(:,:,:,1), y=>position_vectors(:,:,:,2), z=>position_vectors(:,:,:,3) )
          associate( &
            r_sq=>((x-x_center)/(x_max-x_center))**2 + ((y-y_center)/(y_max-y_center))**2 + ((z-z_center)/(z_max-z_center))**2 )
            f_value = 400. - r_sq
          end associate
        end associate

        allocate( cartesian_grid :: f, stat=alloc_stat, errmsg=error_message)
        call assert(alloc_stat==success, "ellipsoidal_field%evaluate: f allocation (error message: "//error_message//")")

        select type(f)
          class is(cartesian_grid)
            call f%set_scalar( f_value )
          class default
            error stop "ellipsoidal_field%evaluate: unsupported class"
        end select

      class default
        error stop "ellipsoidal_field%evaluate: unsupported f class"
    end select

  end procedure

  module procedure laplacian
    use cartesian_grid_interface, only : cartesian_grid
    real(r8k), allocatable, dimension(:,:,:,:) :: position_vectors
    real(r8k), allocatable, dimension(:,:,:) :: laplacian_values
    real(r8k), parameter ::  x_center = 3*(0.25E-01) + 1.E-01/2., x_max = 2*x_center
    real(r8k), parameter ::  y_center = 0.5E-01 + 2*(0.25E-01) + 3.E-01/2., y_max = 2*y_center
    real(r8k), parameter ::  z_center = 20.E-01/2., z_max=2*z_center
    integer alloc_stat
    integer, parameter :: success=0
    character(len=max_errmsg_len) error_message

    select type(grid_points)
      class is(cartesian_grid)

        position_vectors = grid_points%vectors()

        associate( dim=>shape(position_vectors) )
          allocate( laplacian_values(dim(1),dim(2),dim(3)), stat=alloc_stat, errmsg=error_message, &
            source =  - ( (2./(x_max-x_center)) + (2./(y_max-y_center)) + (2./(z_max-z_center)) ) )
          call assert( alloc_stat==success, &
            "ellipsoidal_field%laplacian: allocate laplacian_values (error: "//adjustl(trim(error_message))//")" )
        end associate

        allocate( cartesian_grid :: laplacian_f, stat=alloc_stat, errmsg=error_message)
        call assert(alloc_stat==success, &
          "ellipsoidal_field%evaluate: allocate laplacian_f (error: "//adjustl(trim(error_message))//")")

        select type(laplacian_f)
          class is(cartesian_grid)
            call laplacian_f%set_scalar(laplacian_values)
          class default
            error stop "ellipsoidal_field%laplacian: unsupported laplacian_f type"
        end select

      class default
        error stop "ellipsoidal_field%laplacian: unsupported laplacian_f type"
    end select

  end procedure

end submodule ellipsoidal_field_implementation
