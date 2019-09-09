!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson
  !! date: 10/09/2018
  !!
  !! Test the problem_discretization creation and output of new problem_discretization to a json file

  use assertions_interface, only : assert, assertions
  use json_module, only : json_file, rk=>json_rk, ik=>json_ik
  use problem_discretization_interface, only :  problem_discretization
  implicit none

  type(problem_discretization) block_structured_grid
    !! encapsulate the global grid structure
  type(json_file) json_problem_description
  integer, parameter :: num_structured_grids(*) = [120,5,5]
    !! number of subdomains in each coordinate direction
  integer(ik), parameter :: space_dimensions=3, interval_extent=2
    !! number of spatial dimensions, array index extent for defining 1D spatial intervals
  real(rk), parameter:: global_domain(*,*)= reshape([real(rk):: 0.,120., 0.,2.5, 0.,2.5 ],[interval_extent,space_dimensions])
    !! overall rectangular domain boundaries
  integer, parameter :: max_block_identifier     =     999999999
  integer, parameter :: nx_local=11,ny_local=11,nz_local=11
  integer ix,iy,iz
  character(len=99) identifier_string

  call assert(num_images()==1,"single-image json file creation")
  call assert(product(num_structured_grids)<=max_block_identifier,"number of blocks is representable")

  call block_structured_grid%partition( num_structured_grids )
    !! partition the block-structured grid into subdomains with connectivity implied by the supplied shape array

  call json_problem_description%initialize()  ! specify whatever init options you want.
  call json_problem_description%add('geometry.time', 0.0_rk)
  call json_problem_description%add('geometry.domain.shape', shape(global_domain) )
  call json_problem_description%add('geometry.domain.sides', reshape(global_domain,[size(global_domain)]) )
  call json_problem_description%add('grid.resolution.block_array_shape', num_structured_grids)

  do ix=1,num_structured_grids(1)
    do iy=1,num_structured_grids(2)
      do iz=1,num_structured_grids(3)
        write (identifier_string,'(I0)') block_structured_grid%block_identifier([ix,iy,iz])
        call json_problem_description%add(trim('grid.resolution.id'//identifier_string), [nx_local,ny_local,nz_local])
  end do; end do; end do

  call assert( .not. json_problem_description%failed(), ".not. json_problem_description%failed()")

  call json_problem_description%print_file('fast_structured_grid.json')

  call assert( .not. json_problem_description%failed(), ".not. json_problem_description%failed()")

  print *,"Test passed."

end program
