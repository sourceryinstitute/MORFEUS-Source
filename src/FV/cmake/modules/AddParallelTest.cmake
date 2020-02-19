#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#

# Function to add MPI tests (or MPI-based coarray tests)
include_guard(DIRECTORY)
function(add_parallel_test)
  # Function to add MPI test
  #
  # Arguments:
  #   NAME: Name of the test and executable. Default value: basename without extension of the first source in SOURCES
  #   PROCESSES: Positive integer denoting the number of processes/images to spawn. Default value: 1
  #   TEST_DIR: The working directory for the test. Optional, default is whatever CMake defaults to
  #   SOURCES: A list of source files to compile, including the main program. The first entry is used for the name if NAME not given
  #   LIBS: One or more libraries to link the test executable against

  #  set(options "")
  set(oneValueArgs NAME PROCESSES TEST_DIR)
  set(multiValueArgs SOURCES LIBS)
  cmake_parse_arguments(APT "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT "${APT_SOURCES}")
    message( FATAL_ERROR "`add_parallel_test()` must be called with the SOURCES argument pointing to test sources!")
  endif()
  list(GET APT_SOURCES 0 _first_source)
  if(NOT "${APT_PROCESSES}")
    set(APT_PROCESSES 1)
  elseif("${APT_PROCESSES}" LESS 1)
    message( FATAL_ERROR
      "${APT_PROCESSES} processes requested in `add_parallel_test()`! An integer greater than or equal to 1 must be specified.")
  endif()
  if("${APT_NAME}")
    set(_test_name "${APT_NAME}")
  else()
    get_filename_component(_test_name "${_first_source}" NAME_WE)
    message( AUTHOR_WARNING "TEST name not explicitly passed to `add_parallel_test()`, using ${_test_name}")
  endif()

  add_executable("${_test_name}" ${APT_SOURCES})
  set_property(TARGET ${_test_name}
    PROPERTY FOLDER Morfeus_FV-Tests)

  foreach(lib IN LISTS APT_LIBS)
    target_link_libraries(${_test_name}
      PUBLIC ${lib})
  endforeach()

  if( MPI_FOUND )
    set(_test_launch
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${APT_PROCESSES} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${_test_name}> ${MPIEXEC_POSTFLAGS})
    target_link_libraries(${_test_name}
      PUBLIC MPI::MPI_Fortran)
  else()
    set(_test_launch $<TARGET_FILE:${_test_name}>)
  endif()

  add_test(NAME ${_test_name}
    COMMAND ${_test_launch})
  set_property(TEST ${name}
    APPEND
    PROPERTY LABELS "MORFEUS" "MORFEUS_FV" "integration-test")
  set_property(TEST ${name}
    APPEND
    PROPERTY PASS_REGULAR_EXPRESSION "[Tt]est [Pp]assed")
  set_property(TEST ${name}
    APPEND
    PROPERTY PROCESSORS ${APT_PROCESSES})
  if("${APT_TEST_DIR}")
    set_property(TEST ${name}
      APPEND
      PROPERTY WORKING_DIRECTORY ${APT_TEST_DIR})
  endif()
endfunction()
