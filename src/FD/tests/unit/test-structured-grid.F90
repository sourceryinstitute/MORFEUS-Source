!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
program test_structured_grid
  use kind_parameters, only : i4k, r8k
  use structured_grid_interface, only : structured_grid
  use cartesian_grid_interface, only : cartesian_grid
  use assertions_interface, only : assert, max_errmsg_len
#ifndef HAVE_COLLECTIVE_SUBROUTINES
  use emulated_intrinsics_interface, only : co_sum
#endif
  implicit none

  character(len=max_errmsg_len) :: alloc_error

  integer(i4k), parameter :: success=0
  integer(i4k) alloc_status
  integer(i4k), parameter :: nx=101,ny=11,nz=11
    !! Grid resolution in 3 coordinate directions
  real(r8k) :: x(nx,ny,nz)=1.,y(nx,ny,nz)=0.,z(nx,ny,nz)=-1.
    !! Vector components
  class(structured_grid), allocatable :: coordinate_plane
    !!
  type(cartesian_grid) prototype
    !! pass the cartesian_grid type

  allocate( coordinate_plane, stat=alloc_status, errmsg=alloc_error, mold=prototype )
  call assert( alloc_status==success, "test_structured_grid: allocation ("//alloc_error//")" )

  call coordinate_plane%set_vector_components(x,y,z)
  call assert(coordinate_plane%space_dimension()==3    ,"test_structured_grid: 3D")
  call assert(coordinate_plane%free_tensor_indices()==1,"test_structured_grid: 1 free tensor index (vector)")
  call assert(coordinate_plane%num_time_stamps()==1    ,"test_structured_grid: single time stamp")

  print *,"Test passed."

end program
