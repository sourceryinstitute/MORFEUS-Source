!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson
  !! date: 2019-06-21
  !! summary: test the initialization of a problem_discretization from json_file input
  use plate_3D_interface, only : plate_3D
  implicit none

  type(plate_3D) plate_3D_geometry

  call plate_3D_geometry%build( "3Dplate-low-resolution-layers.json" )

  print *,"Test passed."
end program
