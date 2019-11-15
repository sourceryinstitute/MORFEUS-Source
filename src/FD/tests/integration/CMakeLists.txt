#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
set( base_names
  plate-material-mapping
  discretize-3D-plate
  write-slab-vtk
  write-vtk-voxels
  differentiate-grid
)
foreach(integration_test ${base_names})
  add_executable( test-${integration_test} test-${integration_test}.F90 )
  target_link_libraries( test-${integration_test}
    morfeus_fd
    vtkmofo
  ) # also links in morfeus_utils through transitive dependencies
  set_property( TARGET test-${integration_test}
    PROPERTY FOLDER "MORFEUS-FD-Tests"
  )
endforeach()
