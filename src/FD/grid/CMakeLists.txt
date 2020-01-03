#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
set( morfeus_fd_sources
  surfaces_interface.f90
  surfaces_implementation.f90
  curvilinear_grid_interface.f90
  units_implementation.f90
  units_interface.F90
  grid_implementation.f90
  grid_interface.F90
  problem_discretization_implementation.F90
  problem_discretization_interface.F90
  cartesian_grid_interface.f90
  cartesian_grid_implementation.f90
  structured_grid_implementation.f90
  structured_grid_interface.F90
  differentiable_field_interface.f90
  ellipsoidal_field_interface.f90
  ellipsoidal_field_implementation.f90
  ../geometry/geometry_interface.f90
  ../geometry/geometry_implementation.f90
  ../geometry/plate_3D_implementation.F90
  ../geometry/plate_3D_interface.F90
  ../geometry/block_metadata_interface.f90
  ../geometry/block_metadata_implementation.f90
)
add_library(morfeus_fd ${morfeus_fd_sources})
target_link_libraries(morfeus_fd PUBLIC
  morfeus_utils)
install(TARGETS morfeus_fd DESTINATION lib)

set_property(TARGET morfeus_fd
  PROPERTY
  FOLDER "MORFEUS-FD")
if(TARGET jsonfortran-static)
  target_link_libraries ( morfeus_fd
    PUBLIC jsonfortran-static )
else()
  message(FATAL_ERROR "Couldn't find JSON-Fortran! include_var=${jsonfortran_INCLUDE_DIRS}")
endif()
if(TARGET vtkmofo)
  target_link_libraries ( morfeus_fd
    PUBLIC vtkmofo)
else()
  message(FATAL_ERROR "Couldn't find VTKmofo! include_var=${vtkmofo_INCLUDE_DIRS}")
endif()
