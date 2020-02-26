!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  !! author: Damian Rouson and Xiaofeng Xu
  !! date: 2/24/2020
  !!
  !! Test implicit time advancement of the unsteady, 1D spherical heat equation
  use spherical_1D_solver_interface, only : spherical_1D_solver
  use kind_parameters, only : r8k
  implicit none

  type(spherical_1D_solver) conducting_sphere

  call conducting_sphere%set_v( nr = 101, constants = [0._r8k, 1073.15_r8k] )
  call conducting_sphere%set_expected_solution_size()
  call conducting_sphere%set_material_properties_size()
  call conducting_sphere%time_advance_heat_equation( duration = 100._r8k )

end program
