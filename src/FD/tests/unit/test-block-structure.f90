!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! Test the invertibility of mapping between 3D indicial coordinates and 1D identifiers
  !! for the regular grid blocks encapsulated inside a problem_discretization object.
  use assertions_interface, only : assert
  use problem_discretization_interface, only :  problem_discretization
  use cartesian_grid_interface, only : cartesian_grid
  implicit none

  type(problem_discretization) block_structured_grid
    !! encapsulate the global grid structure
  type(cartesian_grid) prototype
    !! pass the cartesian_grid type
  integer, parameter :: num_structured_grids(*) = [3,3,3]
    !! number of subdomains in each coordinate direction
  integer image
    !! implied do loop index

    associate( me => this_image() )
    associate( ni => num_images() )
    associate( num_blocks => product(num_structured_grids) )

    call assert( ni<=num_blocks, "test-problem-discretization-block-structure: enough blocks to distribute to images")

    call block_structured_grid%partition( num_structured_grids, prototype )
      !! partition the block-structured grid into subdomains with connectivity implied by the supplied shape array

      associate( remainder => mod(num_blocks,ni) )
      associate( quotient => num_blocks/ni )
      associate( my_first => 1 + sum([(quotient+overflow(image,remainder),image=1,me-1)]) )
      associate( my_last => my_first + quotient + overflow(me,remainder) - 1 )
      associate(nx=>num_structured_grids(1), ny=>num_structured_grids(2), nz=>num_structured_grids(3))
        block
          integer ix,iy,iz

          do concurrent( iz=1:nz, iy=1:ny, ix=1:nx )
            associate(i=>block_structured_grid%block_identifier([ix,iy,iz]))
              call assert( all( block_structured_grid%block_indicial_coordinates(i)==[ix,iy,iz]), "invertible 1D<->3D mapping" )
            end associate
          end do

        end block
    end associate; end associate; end associate; end associate; end associate; end associate; end associate

    sync all
    if (me==1) print *,"Test passed."
  end associate

contains

  pure function overflow(image,remainder) result(filler)
    integer, intent(in) :: image,remainder
    integer :: filler
    filler = merge(1,0,image<=remainder)
  end function

end program
