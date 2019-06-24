#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under 
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
set(CGNS_LIB_NAMES cgns)

# Prefer static libraries
if(NOT BUILD_SHARED_LIBS)
  set(CGNS_LIB_NAMES libcgns.a ${CGNS_LIB_NAMES})
endif()

if( DEFINED ENV{CGNSDIR} )
  if( NOT DEFINED CGNS_ROOT )
    set(CGNS_ROOT "$ENV{CGNSDIR}")
  endif()
endif()

if( (DEFINED ENV{CGNS_ROOT}) OR (DEFINED CGNS_ROOT) )
  if( NOT DEFINED CGNS_ROOT)
    set(CGNS_ROOT "$ENV{CGNS_ROOT}")
  endif()
  set(CGNS_HINTS "${CGNS_ROOT}")
endif()

find_path(CGNS_INCLUDES cgnslib.h
  HINTS
  ${CGNS_ROOT}
  PATHS
  /usr/local/opt
  /usr/local
  /usr
  PATH_SUFFIXES
  include
  DOC
  "Path to the CGNS include files."
  )

if(CGNS_INCLUDES)
  foreach(include IN_LISTS CGNS_INCLUDES)
    get_filename_component(cgns_include_dir "${include}" DIRECTORY)
    get_filename_component(cgns_abs_include_dir "${cgns_include_dir}" ABSOLUTE)
    get_filename_component(new_cgns_hint "${include_dir}/.." ABSOLUTE )
    list(APPEND CGNS_HINTS "${new_cgns_hint}")
    break()
  endforeach()
endif()

if(CGNS_HINTS)
  list(REMOVE_DUPLICATES CGNS_HINTS)
endif()

find_path(CGNS_MODULE cgns.mod
  HINTS
  ${CGNS_HINTS}
  PATHS
  /usr/local/opt
  /usr/local
  /usr
  /usr/lib64
  PATH_SUFFIXES
  include
  gfortran/modules
  DOC
  "Location of the CGNS Fortran module file."
  )

list(APPEND CGNS_INCLUDES ${CGNS_MODULE})
list(REMOVE_DUPLICATES CGNS_INCLUDES)

find_library(CGNS_LIBRARY NAMES ${CGNS_LIB_NAMES}
  HINTS
  ${CGNS_HINTS}
  PATHS
  /usr/local/opt
  /usr/local
  /usr
  /usr
  PATH_SUFFIXES
  lib
  lib64
  DOC  "Path to the CGNS library"
)

set(CGNS_FOUND "NO")
if(CGNS_MODULE)
  if(CGNS_LIBRARY)
    set( CGNS_LIBRARIES ${CGNS_LIBRARY} )
    set( CGNS_FOUND "YES" )
  endif()
endif()

if(CGNS_FIND_REQUIRED AND NOT CGNS_FOUND)
  message(SEND_ERROR "Unable to find the requested CGNS libraries.")
endif()

# handle the QUIETLY and REQUIRED arguments and set CGNS_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CGNS DEFAULT_MSG CGNS_LIBRARY CGNS_INCLUDES)


mark_as_advanced(
  CGNS_INCLUDES
  CGNS_MODULE
  CGNS_LIBRARY
  )
