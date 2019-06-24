#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
foreach( prog test-rd_gambit_mesh.f90)
  add_morfeus_fv_exe( "${CMAKE_CURRENT_SOURCE_DIR}/${prog}" )
endforeach()
