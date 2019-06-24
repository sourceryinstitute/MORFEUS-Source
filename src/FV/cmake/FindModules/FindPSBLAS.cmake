#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
# N.B.!: The order here is important: this dictates the order
#        passed to the linker on the command line; if it is wrong
#        you will likely get unresolved symbols!

if(WIN32)
  set(PSBLAS_SUFFIX_NAMES cbind util util_C krylov prec base)
  list(APPEND PSBLAS_SUFFIX_NAMES base_C)
else()
  set(PSBLAS_SUFFIX_NAMES cbind util krylov prec base)
endif()


foreach(suffix IN LISTS PSBLAS_SUFFIX_NAMES)
  LIST(APPEND PSBLAS_LIB_NAMES_WO_EXTENSION psb_${suffix})
endforeach()

# Prefer static libraries
if(NOT ${BUILD_SHARED_LIBS})
  if(WIN32)
    set(lib_suffix .lib)
  else()
    set(lib_suffix .a)
  endif()
endif()

# Look in certain places based on environment variables
if(DEFINED ENV{PSBLASDIR})
  if(NOT DEFINED PSBLAS_ROOT)
    set(PSBLAS_ROOT "$ENV{PSBLASDIR}")
  endif()
endif()

if( (DEFINED ENV{PSBLAS_ROOT}) OR (DEFINED PSBLAS_ROOT) )
  if( NOT DEFINED PSBLAS_ROOT)
    set(PSBLAS_ROOT "$ENV{PSBLAS_ROOT}")
  endif()
  set(PSBLAS_HINTS "${PSBLAS_ROOT}")
endif()

set(have_needed_libs TRUE)

foreach(lib IN LISTS PSBLAS_LIB_NAMES_WO_EXTENSION)
  message(STATUS "Looking for ${lib}${suffix}")
  find_library(PSBLAS_${lib} NAMES ${lib} # ${lib_suffix}
    HINTS
    ${PSBLAS_HINTS}
    PATHS
    /usr
    /usr/local
    /usr/local/opt/psblas
    PATH_SUFFIXES
    lib
    lib64
    DOC "Path to PSBLAS library ${lib}")
  if(NOT (PSBLAS_${lib} STREQUAL PSBLAS_psb_cbind))
    if( (NOT PSBLAS_${lib}) AND
	(NOT (PSBLAS_${lib} STREQUAL PSBLAS_psb_util_C))
	)
      set(have_needed_libs FALSE)
      message(WARNING "Failed to find ${lib}")
    else()
      message(STATUS "Adding ${PSBLAS_${lib}}")
      list(APPEND PSBLAS_Fortran_LIBRARIES ${PSBLAS_${lib}})
    endif()
  else()
    if(PSBLAS_${lib})
      set(PSBLAS_C_LIBRARY ${PSBLAS_${lib}})
    else()
      set(PSBLAS_C_LIBRARY NOTFOUND)
    endif()
  endif()
  mark_as_advanced(PSBLAS_${lib})
endforeach()

find_path(PSBLAS_Fortran_INCLUDE_DIR
  NAMES psb_psblas_mod.mod
  HINTS ${PSBLAS_HINTS}
  PATHS /usr/local /usr/local/opt/psblas
  PATH_SUFFIXES include lib modules include/debug include/release
  DOC "Location of PSBLAS Fortran module files."
  )

set(PSBLAS_FOUND "NO")
if(PSBLAS_Fortran_INCLUDE_DIR)
  if(have_needed_libs)
    set(PSBLAS_LIBRARIES ${PSBLAS_Fortran_LIBRARIES})
    set(PSBLAS_FOUND "YES")
  endif()
endif()

if(PSBLAS_FOUND)
  if(${PSBLAS_C_LIBRARY})
    list(INSERT PSBLAS_LIBRARIES 0 ${PSBLAS_C_LIBRARY})
  endif()
endif()

if(PSBLAS_FIND_REQUIRED AND NOT PSBLAS_FOUND)
  message(SEND_ERROR
    "Unable to find the requested PSBLAS libraries. Please point CMake to your
installation of PSBLAS using ccmake, cmake-gui or by passing the
-DPSBLAS_ROOT=/path/to/PSBLAS/installation flag on the command line")
endif()

# handle the QUIETLY and REQUIRED arguments and set PSBLAS_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PSBLAS DEFAULT_MSG PSBLAS_LIBRARIES PSBLAS_Fortran_INCLUDE_DIR)


mark_as_advanced(
  PSBLAS_Fortran_INCLUDE_DIR
  PSBLAS_LIBRARIES
  )
