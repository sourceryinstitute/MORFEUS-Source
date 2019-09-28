#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
set( morfeus_utils_sources
  assertions_interface.F90
  assertions_implementation.F90
  object_interface.f90
  object_implementation.f90
  co_object_interface.f90
  co_object_implementation.f90
  array_functions_interface.f90
  array_functions_implementation.f90
  string_functions_interface.f90
  string_functions_implementation.f90
  kind_parameters.f90
)
if(EMULATE_MISSING_INTRINSICS)
  set( morfeus_utils_sources
    ${morfeus_utils_sources}
    emulated_intrinsics_interface.F90
    emulated_intrinsics_implementation.F90
  )
endif()
add_library( morfeus_utils ${morfeus_utils_sources})
if(USE_ASSERTIONS)
  target_compile_definitions(morfeus_utils
    PUBLIC -DUSE_ASSERTIONS=.true.)
endif()
install(TARGETS morfeus_utils DESTINATION lib)
set_property(TARGET morfeus_utils
  PROPERTY
  FOLDER "MORFEUS-FD")
if(TARGET OpenCoarrays::caf_mpi_static)
  target_link_libraries(morfeus_utils
    PUBLIC OpenCoarrays::caf_mpi_static)
endif()
target_link_libraries(morfeus_utils
  PUBLIC utilities)
