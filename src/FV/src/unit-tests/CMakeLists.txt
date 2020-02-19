#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
foreach( prog test-rd_gambit_mesh.f90)
  add_morfeus_fv_exe( "${CMAKE_CURRENT_SOURCE_DIR}/${prog}" )
endforeach()

# Example of how to use ../../cmake/modules/AddParallelTest.cmake
# NOTE: The CMAKE_MODULE_PATH was modified in the project's top  level CMakeLists.txt
include(AddParallelTest)

add_parallel_test(NAME fv-unit-test # Name is optional
  PROCESSES 3                       # Processes is optional, will default to 1
  TEST_DIR ./paht/to/subdir         # Optional, path to directory where test should be run
  SOURCES fv-unit-test1.F90         # Specify one or more source files. If NAME is missing
  LIBS morfeus_fv
)
