#! author: Izaak Beekman, ParaTools, Inc, Damian Rouson, GSE LLC, Ian Porter, NRC
#! summary: Top level CMake file for building Morfeus-FD component
#! category: Morfeus-FD
#!
#! # Morfeus-FD Top Level CMake file
#!
#!
#!
#! @note
#! This file is currently being refactored
#!
#! ## Copyright Notice
#!
#!```
#!     (c) 2019 Guide Star Engineering, LLC
#!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#!```

cmake_minimum_required(VERSION 3.13.4)

# Name project and specify source languages
project(morfeus VERSION "0.1" LANGUAGES Fortran)

option(USE_ASSERTIONS "Turn on runtime assertion checking"
  FALSE)
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
  set (CMAKE_Fortran_FLAGS_DEBUG "-g -DUSE_ASSERTIONS=.true." CACHE
    STRING "Flags used by the compiler during debug builds."
    FORCE)
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -O2 -DNDEBUG -DUSE_ASSERTIONS=.true." CACHE
    STRING "Flags used by the compiler during release builds with debug info."
    FORCE)
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
  set (CMAKE_Fortran_FLAGS_DEBUG "-DUSE_ASSERTIONS=.true." CACHE
    STRING "Flags used by the compiler during debug builds."
    FORCE)
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -DNDEBUG -DUSE_ASSERTIONS=.true." CACHE
    STRING "Flags used by the compiler during release builds with debug info."
    FORCE)
endif()

  if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" )
    set(coarray_flags "-fcoarray=single")
  elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
    if ( "${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
      set(coarray_flags "${prefix}Qcoarray${eq}shared")
    else()
      set(coarray_flags "${prefix}coarray${eq}shared")
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG "${coarray_flags} ${CMAKE_Fortran_FLAGS_DEBUG}")
    set(CMAKE_Fortran_FLAGS_RELEASE "${coarray_flags} ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_MINSIZEREL "${coarray_flags} ${CMAKE_Fortran_FLAGS_MINSIZEREL}")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${coarray_flags} ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
  else()
    message(FATAL_ERROR "Only GFortran and Intel Fortran compilers supported.")
  endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
#install(TARGETS mod DESTINATION mod)


set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
set(CMAKE_REQUIRED_FLAGS "${coarray_flags}")
include(CheckFortranSourceCompiles)
check_fortran_source_compiles("
  program main
    implicit none
    integer :: i[*]
    i = this_image()
  end program
" Check_Simple_Coarray_Fortran_Source_Compiles
  SRC_EXT ".f90")
set (CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
unset(OLD_REQUIRED_FLAGS)

check_fortran_source_compiles("
  program main
  contains
  pure function foo() result(res)
    error stop 'Error stop is supported in pure functions (F2018)'
  end function
  end program
"
  HAVE_ERROR_STOP_IN_PURE
  SRC_EXT ".f90"
  )
if(HAVE_ERROR_STOP_IN_PURE)
  add_definitions(-DHAVE_ERROR_STOP_IN_PURE)
endif()

check_fortran_source_compiles("
  program main
  integer :: i
  i = 0
  error stop i
  end program
"
  HAVE_NON_CONSTANT_ERROR_STOP
  SRC_EXT ".f90")
if(HAVE_NON_CONSTANT_ERROR_STOP)
  add_definitions(-DHAVE_NON_CONSTANT_ERROR_STOP)
endif()

check_fortran_source_compiles("
  program main
  implicit none
  print *,findloc(array=[1,2], value=1, dim=1)
  end program
"
  HAVE_FINDLOC
  SRC_EXT ".f90")
if(HAVE_FINDLOC)
  add_definitions(-DHAVE_FINDLOC)
endif()

check_fortran_source_compiles("
  program main
  implicit none
  integer :: i=1
  call co_sum(i)
  call co_broadcast(i,source_image=1)
  end program
"
  HAVE_COLLECTIVE_SUBROUTINES
  SRC_EXT ".f90")
if(HAVE_COLLECTIVE_SUBROUTINES)
  add_definitions(-DHAVE_COLLECTIVE_SUBROUTINES)
endif()

if( (NOT HAVE_FINDLOC) OR (NOT HAVE_COLLECTIVE_SUBROUTINES))
  set(EMULATE_MISSING_INTRINSICS TRUE)
endif()

set( directories utilities grid tests ) # Order *may* matter here
foreach( directory ${directories} )
  add_subdirectory(${directory})
endforeach()

enable_testing()

# Add test fixtures for setup and teardown
set(MORFEUS_FD_INPUTS
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/3Dplate-low-resolution-layers.json"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/3Dplate-low-resolution-layers-material-map.json"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/3Dplate-high-resolution-layers.json"
  "${CMAKE_CURRENT_SOURCE_DIR}/tests/unit/3Dplate-low-resolution-halo.json"
)

#  "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/3Dplate-high-resolution.json"
add_test(NAME morfeus-fd-setup-fixture
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  COMMAND ${CMAKE_COMMAND} -E copy ${MORFEUS_FD_INPUTS} "${CMAKE_CURRENT_BINARY_DIR}")
set_tests_properties(morfeus-fd-setup-fixture
    PROPERTIES FIXTURES_SETUP MFD)

#------------------------------------------------------------------------------------------
# Define function for adding CAF tests, and ensuring proper flags are passed to MPI runtime
#------------------------------------------------------------------------------------------
function(add_caf_test name num_caf_img test_target)
  # Function to add MPI tests.
  if(TARGET ${test_target})
    get_target_property(min_test_imgs ${test_target} MIN_IMAGES)
  elseif(TARGET build_${test_target})
    get_target_property(min_test_imgs build_${test_target} MIN_IMAGES)
  endif()
  if(min_test_imgs)
    if(num_caf_img LESS min_test_imgs)
      message( FATAL_ERROR "Test ${name} requires ${min_test_imgs} but was only given ${num_caf_img}" )
    endif()
  endif()
  # Add a host file for OMPI
  if ( openmpi )
    set(test_parameters --hostfile ${CMAKE_BINARY_DIR}/hostfile)
  endif()
  if ( ((N_CPU LESS num_caf_img) OR (N_CPU EQUAL 0)) )
    message(STATUS "Test ${name} is oversubscribed: ${num_caf_img} CAF images requested with ${N_CPU} system processor available.")
    if ( openmpi )
      if (min_test_imgs)
        set( num_caf_img ${min_test_imgs} )
      elseif ( N_CPU LESS 2 )
        set( num_caf_img 2 )
      endif()
      set (test_parameters ${test_parameters} --oversubscribe)
    endif()
  endif()
  if( openmpi)
    set (test_parameters ${test_parameters} --allow-run-as-root)
  endif()
  set(test_parameters ${MPIEXEC_NUMPROC_FLAG} ${num_caf_img} ${test_parameters})
  if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
    add_test(NAME ${name} COMMAND "${test_target}" ${ARGN})
    set_property(TEST ${name}
      APPEND
      PROPERTY ENVIRONMENT FOR_COARRAY_NUM_IMAGES=${num_caf_img})
    message ( STATUS "Setting FOR_COARRAY_NUM_IMAGES=${num_caf_img} for test ${test}")
  else()
    add_test(NAME ${name} COMMAND ${MPIEXEC_EXECUTABLE} ${test_parameters} "${test_target}" ${ARGN})
  endif()
  set_tests_properties( ${name}
    PROPERTIES
    FIXTURES_REQUIRED MFD
    PASS_REGULAR_EXPRESSION "Test passed.")
  set_property(TEST ${name}
    APPEND
    PROPERTY LABELS "MORFEUS_FD" "MORFEUS")
  set_property(TEST ${name}
    APPEND
    PROPERTY PROCESSORS ${num_caf_img})
  set_property(TEST ${name}
    APPEND
    PROPERTY PROCESSOR_AFFINITY ON)

endfunction(add_caf_test)

# Single-image tests
set( base_names
  structured-grid
  write-problem-discretization
  read-geometry
  array-functions
  spatial-derivatives
  spatial-derivatives-new
  adi-time-advancing
  units
)
foreach(unit_test ${base_names})
  add_caf_test( fd-${unit_test} 1 ${CMAKE_CURRENT_BINARY_DIR}/tests/unit/test-${unit_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# Dual-image unit tests
set( base_names
  assertions
)
foreach(unit_test ${base_names})
  add_caf_test( fd-${unit_test} 2 ${CMAKE_CURRENT_BINARY_DIR}/tests/unit/test-${unit_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# Tri-image unit tests
if(EMULATE_MISSING_INTRINSICS)
set( base_names
  emulated-intrinsics
)
foreach(unit_test ${base_names})
  add_caf_test( fd-${unit_test} 3 ${CMAKE_CURRENT_BINARY_DIR}/tests/unit/test-${unit_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
endif()

# Quad-image unit tests
set( base_names
  problem-discretization
  block-structure
)
foreach(unit_test ${base_names})
  add_caf_test( fd-${unit_test} 4 ${CMAKE_CURRENT_BINARY_DIR}/tests/unit/test-${unit_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# 7-image unit tests
set( base_names
  halo-exchange
)
foreach(unit_test ${base_names})
  add_caf_test( fd-${unit_test} 7 ${CMAKE_CURRENT_BINARY_DIR}/tests/unit/test-${unit_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# Single-image integration tests
set( base_names
  plate-material-mapping
  discretize-3D-plate
  write-slab-vtk
  write-voxels-vtk
)
foreach(integration_test ${base_names})
  add_caf_test( fd-${integration_test} 1 ${CMAKE_CURRENT_BINARY_DIR}/tests/integration/test-${integration_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()

# Tri-image integration tests
set( base_names
  differentiate-grid
)
foreach(integration_test ${base_names})
  add_caf_test( fd-${integration_test} 3 ${CMAKE_CURRENT_BINARY_DIR}/tests/integration/test-${integration_test}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
