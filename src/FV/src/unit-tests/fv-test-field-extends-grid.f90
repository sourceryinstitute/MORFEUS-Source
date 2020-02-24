!
!     (c) 2019-2020 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC) under contract
!     "Multi-Dimensional Physics Implementation into Fuel Analysis under Steady-state and Transients (FAST)",
!     contract # NRC-HQ-60-17-C-0007
!
program test_field_extends_grid
  !! author:  Damian Rouson
  !!
  !! Verify that the field child extends parent grid

  use grid_interface, only : grid
  use class_field, only : field
  use assertions_interface, only : assert
  implicit none

  type(grid) parent
  type(field) child

  call assert( extends_type_of(child, parent), "test_field_extends_grid: extends_type_of(child, parent)")

  print *,"Test passed."

end program test_field_extends_grid
