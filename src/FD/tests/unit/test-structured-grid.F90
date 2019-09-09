!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program test_structured_grid
  use structured_grid_interface, only : structured_grid
  use assertions_interface, only : assert
#ifndef HAVE_COLLECTIVE_SUBROUTINES
  use emulated_intrinsics_interface, only : co_sum
#endif
  implicit none

  type(structured_grid) :: coordinate_plane

  integer, parameter :: nx=101,ny=11,nz=11
    !! Grid resolution in 3 coordinate directions
  real :: x(nx,ny,nz)=1.,y(nx,ny,nz)=0.,z(nx,ny,nz)=-1.
    !! Vector components

  call coordinate_plane%set_vector_components(x,y,z)
  call assert(coordinate_plane%space_dimension()==3    ,"test_structured_grid: 3D")
  call assert(coordinate_plane%free_tensor_indices()==1,"test_structured_grid: 1 free tensor index (vector)")
  call assert(coordinate_plane%num_time_stamps()==1    ,"test_structured_grid: single time stamp")

  print *,"Test passed."

end program
