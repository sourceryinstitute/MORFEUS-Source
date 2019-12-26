#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
set( base_names
  assertions
  units
  problem-discretization
  block-structure
  write-problem-discretization
  read-geometry
  array-functions
  spatial-derivatives
  adi-time-advancing
  spatial-derivatives-new
  halo-exchange
)
foreach(unit_test ${base_names})
  add_executable( test-${unit_test} test-${unit_test}.f90 )
  target_link_libraries( test-${unit_test} morfeus_fd  ) # also links in morfeus_utils through transitive dependencies
  set_property( TARGET test-${unit_test}
    PROPERTY FOLDER "MORFEUS-FD-Tests"
    )
endforeach()

set(requires_preprocessing
  structured-grid
  emulated-intrinsics
)
foreach(unit_test ${requires_preprocessing})
  add_executable( test-${unit_test} test-${unit_test}.F90 )
  target_link_libraries( test-${unit_test} morfeus_fd  ) # also links in morfeus_utils through transitive dependencies
  set_property( TARGET test-${unit_test}
    PROPERTY FOLDER "MORFEUS-FD-Tests"
  )
endforeach()
