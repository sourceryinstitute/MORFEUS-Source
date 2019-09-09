#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
set( test_dir_list
  unit
  integration
)
foreach(test_dir ${test_dir_list})
  add_subdirectory( ${test_dir} )
endforeach()
