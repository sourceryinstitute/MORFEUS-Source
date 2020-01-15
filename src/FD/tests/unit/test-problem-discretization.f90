!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
program main
  !! Test the definition of a block-structured problem discretization
  use assertions_interface, only : assert, assertions
  use problem_discretization_interface, only :  problem_discretization
  use cartesian_grid_interface, only : cartesian_grid
  use kind_parameters, only : r8k
  implicit none

  type(problem_discretization) global_grid
    !! encapsulate the global grid structure
  type(cartesian_grid) prototype
    !! pass the cartesian_grid type
  integer, parameter :: space_dimensions=3
    !! maximum number of spatial dimensions
  integer, parameter :: lo_bound=1, up_bound=2
    !! array index values for defining lower & upper bounds for 1D spatial interval
  real(r8k), parameter :: global_domain(*,*) = reshape([0.,120., 0.,2.5, 0.,2.5],[up_bound,space_dimensions])
    !! overall rectangular domain boundaries
  integer, parameter :: num_structured_grids(*) = [240,5,5]
    !! number of blocks in each coordinate direction

  call global_grid%partition( num_structured_grids, prototype )
    !! partition the block-structured grid into blocks with connectivity implied by the supplied shape of the 3D array of blocks

  associate( my_blocks => global_grid%my_blocks() )

    block
      integer ijk(space_dimensions)
        !! indicial coordinates of this block
      real(r8k), allocatable, dimension(:,:,:) :: x, y, z
        !! coordinates of grid vertex locations
      real(r8k) block(lo_bound:up_bound,space_dimensions)
      integer m,n
      integer, parameter :: nx=11,ny=11,nz=11
        !! resolution within blocks

    !do concurrent( n = my_blocks(1) : my_blocks(2)) TODO: make concurrent after Intel supports co_sum
      do n = my_blocks(1) , my_blocks(2)
        ijk = global_grid%block_indicial_coordinates(n)

        do concurrent(m=1:space_dimensions)
          associate( block_width_m => (global_domain(up_bound,m) - global_domain(lo_bound,m))/num_structured_grids(m))
            block(:,m) = [ ijk(m)-1, ijk(m) ]*block_width_m
          end associate
        end do

        x = evenly_spaced_points(  block, [nx,ny,nz], direction=1 )
        y = evenly_spaced_points(  block, [nx,ny,nz], direction=2 )
        z = evenly_spaced_points(  block, [nx,ny,nz], direction=3 )

        call global_grid%set_vertices(x,y,z,block_identifier=n)
      end do
    end block
  end associate

  if (assertions) &
    call assert(load_balanced(global_grid%block_load(),product(num_structured_grids)),"test_problem_discretization: load balanced")

  print *,"Test passed"

contains

  pure function load_balanced(my_load, total_load) result(block_distribution)
    integer, intent(in) :: my_load, total_load
    logical :: block_distribution
    block_distribution = ( my_load == total_load/num_images() + merge(1,0,this_image()<=mod(total_load,num_images())) )
  end function

  pure function evenly_spaced_points( boundaries, resolution, direction ) result(grid_nodes)
    !! Define grid point coordinates with uniform spacing in the chosen block
    real(r8k), intent(in) :: boundaries(:,:)
      !! block boundaries of each coordinate direction
    integer, intent(in) :: resolution(:)
      !! number of grid points in each direction
    integer, intent(in) :: direction
      !! coordinate direction to define
    real(r8k), allocatable :: grid_nodes(:,:,:)
      !! grid node locations and spacing in each coordination direction
    real(r8k) dx(space_dimensions)

    integer ix,iy,iz

    allocate(grid_nodes(resolution(1),resolution(2),resolution(3)))

    dx = ( boundaries(up_bound,:) - boundaries(lo_bound,:) ) / resolution(:)

    associate( nx=>resolution(1), ny=>resolution(2), nz=>resolution(3) )

    select case(direction)
      case(1)
        do concurrent(iy=1:ny,iz=1:nz)
          grid_nodes(:,iy,iz) = [ boundaries(lo_bound,direction), (ix*dx(direction),ix=2,nx-1), boundaries(up_bound,direction) ]
        end do
      case(2)
        do concurrent(ix=1:nx,iz=1:nz)
          grid_nodes(ix,:,iz) = [ boundaries(lo_bound,direction), (iy*dx(direction),iy=2,ny-1), boundaries(up_bound,direction) ]
        end do
      case(3)
        do concurrent(ix=1:nx,iy=1:ny)
          grid_nodes(ix,iy,:) = [ boundaries(lo_bound,direction), (iz*dx(direction),iz=2,nz-1), boundaries(up_bound,direction) ]
        end do
      case default
        error stop "evenly_spaced_points: invalid direction"
    end select

    end associate

  end function

end program
