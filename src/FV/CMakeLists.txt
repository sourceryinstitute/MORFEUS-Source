#
#     (c) 2019 Guide Star Engineering, LLC
#     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
#     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
#     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
#
cmake_minimum_required(VERSION 3.11)
cmake_policy(VERSION 3.11.1...3.14.1)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS_ON)

project(MORFEUS_FV C Fortran)

#----------------------------------------------
# Determine system endian-ness and pointer size
#----------------------------------------------
include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)
if(${IS_BIG_ENDIAN})
  message( STATUS "System appears to be big endian.")
else()
  message( STATUS "System appears to be little endian.")
  add_definitions(-DNemoLittleEndian)
endif()
include(CheckIncludeFile)
CHECK_INCLUDE_FILE("endian.h" HAVE_ENDIAN_H)
if(${HAVE_ENDIAN_H})
  add_definitions(-DHAVE_ENDIAN)
endif()
include(CheckTypeSize)
CHECK_TYPE_SIZE("void *" VOID_P_SIZE LANGUAGE C)
if(${VOID_P_SIZE} EQUAL 8)
  add_definitions(-DPtr64Bits)
endif()
message(STATUS "Have 64bit pointers")

#-----------------------------------------
# Check for some Fortran compiler features
#-----------------------------------------
include(CheckFortranSourceCompiles)
CHECK_Fortran_SOURCE_COMPILES(
  "
integer, allocatable :: a(:), b(:)
allocate(a(5))
a = [1,2,3,4,5]
call move_alloc(from=a, to=b)
end
"
  HAVE_MOVE_ALLOC
  SRC_EXT f90
  )
if(HAVE_MOVE_ALLOC)
  add_definitions(-DHAVE_MOVE_ALLOC)
endif()
CHECK_Fortran_SOURCE_COMPILES(
  "integer, volatile :: i ; end"
  HAVE_VOLATILE
  SRC_EXT f90
  )
if(HAVE_VOLATILE)
  add_definitions(-DHAVE_VOLATILE)
endif()


#--------------------------------------------------------
# Tell build system where to find our custom find modules
#--------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindModules)

# Needed for mesh_optimize sources.
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_optimize/OptMS")


# Add option for enabling exodus
option(FAST_USE_EXODUS "Build with Exodus support (broken)" OFF)

# Make CGNS optional
option(FAST_USE_CGNS "Build with CGNS support (broken)" OFF)


set(SLATEC_srcs
  src/math/SLATEC/d1mach.f90
  src/math/SLATEC/dckder.f
  src/math/SLATEC/denorm.f
  src/math/SLATEC/dfdjc3.f
  src/math/SLATEC/dmpar.f
  src/math/SLATEC/dnls1e.f
  src/math/SLATEC/dnls1.f
  src/math/SLATEC/dqrfac.f
  src/math/SLATEC/dqrslv.f
  src/math/SLATEC/dwupdt.f
  src/math/SLATEC/fdump.f
  src/math/SLATEC/i1mach.f90
  src/math/SLATEC/j4save.f
  src/math/SLATEC/xercnt.f
  src/math/SLATEC/xerhlt.f
  src/math/SLATEC/xermsg.f
  src/math/SLATEC/xerprn.f
  src/math/SLATEC/xersve.f
  src/math/SLATEC/xgetua.f
  src/math/SLATEC/xsetf.f
  )
set(morfeus_fv_srcs
  src/basics/class_connectivity.f90
  src/basics/class_connectivity_procedures.f90
  src/basics/class_dimensions.f90
  src/basics/class_dimensions_procedures.f90
  src/basics/class_keytable.f90
  src/basics/class_keytable_procedures.f90
  src/basics/class_psblas.f90
  src/basics/class_psblas_procedures.f90
  src/basics/class_stopwatch.f90
  src/basics/class_stopwatch_procedures.f90
  src/basics/class_vector.f90
  src/basics/class_vector_procedures.f90
  src/basics/tools_psblas.f90
  src/basics/tools_psblas_procedures.f90
  src/basics/type_table.f90
  src/basics/type_table_procedures.f90
  src/bc/class_bc.f90
  src/bc/class_bc_procedures.f90
  src/bc/class_bc_math.f90
  src/bc/class_bc_math_procedures.f90
  src/bc/class_bc_wall.f90
  src/bc/class_bc_wall_procedures.f90
  src/bc/rd_inp_bc.f90
  src/bc/rd_inp_bc_math.f90
  src/bc/tools_bc.f90
  src/field/class_discretization.f90
  src/field/class_discretization_procedures.f90
  src/field/class_field.f90
  src/field/class_field_procedures.f90
  src/field/class_scalar_field.f90
  src/field/class_scalar_field_procedures.f90
  src/field/class_vector_field.f90
  src/field/class_vector_field_procedures.f90
  src/input/find_section.f90
  src/input/open_file.F90
  src/input/tools_input.f90
  src/material/class_material.f90
  src/material/class_material_procedures.f90
  src/material/load_material.f90
  src/material/mat_copper.f90
  src/material/mat_copper_procedures.f90
  src/material/mat_water.f90
  src/material/mat_water_procedures.f90
  src/material/rd_inp_material.f90
  src/material/tools_material.f90
  src/math/cart_to_polar_implementation.f90
  src/math/cholesky_implementation.f90
  src/math/class_iterating.f90
  src/math/class_iterating_procedures.f90
  src/math/isort_implementation.f90
  src/math/lin_interp_implementation.f90
  src/math/lu_implementation.f90
  src/math/psb_build_prec.f90
  src/math/psb_solve_sys_implementation.f90
  src/math/pwl_implementation.f90
  src/math/tools_math.f90
  src/mesh_basics/class_cell.F90
  src/mesh_basics/class_cell_procedures.F90
  src/mesh_basics/class_cylinder.f90
  src/mesh_basics/class_cylinder_procedures.f90
  src/mesh_basics/class_face.F90
  src/mesh_basics/class_face_procedures.F90
  src/mesh_basics/class_least_squares.f90
  src/mesh_basics/class_least_squares_procedures.f90
  src/mesh_basics/class_plane.f90
  src/mesh_basics/class_plane_procedures.f90
  src/mesh_basics/class_surface.f90
  src/mesh_basics/class_surface_procedures.f90
  src/mesh_basics/class_vertex.f90
  src/mesh_basics/geom_cell.f90
  src/mesh_basics/geom_diff.f90
  src/mesh_basics/geom_face.f90
  src/mesh_basics/geom_hex_dihedral_angle.f90
  src/mesh_basics/geom_hex_quality.f90
  src/mesh_basics/geom_tet_center.f90
  src/mesh_basics/geom_tet_dihedral_angle.f90
  src/mesh_basics/geom_tet_quality.f90
  src/mesh_basics/geom_tet_volume.f90
  src/mesh_basics/class_vertex_procedures.f90
  src/mesh_basics/tools_mesh_basics.f90
  src/mesh_check/check_mesh_quality.f90
  src/mesh_check/check_tet_quality.f90
  src/mesh_check/check_vertex_quality.f90
  src/mesh_check/tools_mesh_check.f90
  src/mesh_move/class_motion.f90
  src/mesh_move/class_motion_procedures.f90
  src/mesh_move/rd_inp_motion_law.f90
  src/mesh_move/tools_mesh_move.f90
  src/mesh_optimize/check_right_handed.f90
  src/mesh_optimize/laplacian_smooth.f90
  src/mesh_optimize/mobile_verts.f90
  src/mesh_optimize/smooth_interior_vtx.f90
  src/mesh_optimize/smooth_mesh.f90
  src/mesh_optimize/smooth_surf_vtx.f90
  src/mesh_optimize/tools_mesh_optimize.f90
  src/mesh/class_mesh.F90
  src/mesh/class_mesh_procedures.F90
  src/mesh/cmp_mesh_desc.f90
  src/mesh/cmp_mesh_implementation.f90
  src/mesh/part_block.f90
  src/mesh/part_block_procedures.f90
  src/mesh/part_graph.f90
  src/mesh/part_graph_procedures.f90
  src/mesh/part_random.f90
  src/mesh/part_random_procedures.f90
  src/mesh/rd_cgns_mesh.F90
#  src/mesh/rd_exodus_mesh.f90 # Commented out b/c functionality does not yet work
  src/mesh/rd_gambit_implementation.f90
  src/mesh/rd_inp_mesh_implementation.f90
  src/mesh/renum.F90
  src/mesh/renum_procedures.F90
  src/mesh/supplement_implementation.f90
  src/mesh/tools_mesh.f90
  src/mesh/tools_part.f90
  src/mesh/tools_part_procedures.f90
  src/modules/tools_nemo.f90
  src/operators/flux_pde_div.f90
  src/operators/op_d2dt2.f90
  src/operators/op_d2dt2_procedures.f90
  src/operators/op_ddt.f90
  src/operators/op_ddt_procedures.f90
  src/operators/op_div.f90
  src/operators/op_field.f90
  src/operators/op_grad.f90
  src/operators/op_laplacian.f90
  src/operators/op_laplacian_procedures.f90
  src/operators/op_source.f90
  src/operators/rhie_chow.f90
  src/operators/scalar_field_grad.f90
  src/operators/scalar_fld_vector_mul.f90
  src/operators/scalar_pde_d2dt2.f90
  src/operators/scalar_pde_ddt.f90
  src/operators/scalar_pde_div.f90
  src/operators/scalar_pde_laplacian.f90
  src/operators/scalar_pde_source.f90
  src/operators/scalar_vector_fld_mul.f90
  src/operators/tools_operators.f90
  src/operators/tools_operators_procedures.f90
  src/operators/vector_field_flux.f90
  src/operators/vector_field_grad.f90
  src/operators/vector_pde_d2dt2.f90
  src/operators/vector_pde_ddt.f90
  src/operators/vector_pde_div.f90
  src/operators/vector_pde_grad.f90
  src/operators/vector_pde_laplacian.f90
#  src/operators/vector_pde_source.f90  Not currently used
  src/output/class_output.f90
  src/output/class_output_procedures.f90
  src/output/htoi.f90
  src/output/itoh.f90
  src/output/tools_output_basics.f90
  src/output/tools_output.f90
  src/output/vtkmofo_interface.f90
  src/output/vtkmofo_implementation.f90
  src/output/wr_mtx_matrix.f90
  src/output/wr_mtx_pattern.f90
  src/output/wr_mtx_vector.f90
  src/output/write_mesh.f90
  src/output/write_scalar_field.f90
  src/output/write_vector_field.f90
  src/pde/class_pde.f90
  src/pde/class_pde_procedures.f90
  src/pde/class_scalar_pde.f90
  src/pde/class_scalar_pde_procedures.f90
  src/pde/class_scalar_source.f90
  src/pde/class_scalar_source_procedures.f90
  src/pde/class_vector_pde.f90
  src/pde/class_vector_pde_procedures.f90
  )
set(morfeus_fv_C_srcs
  src/mesh_optimize/call_smooth2d.h
  src/mesh_optimize/call_smooth2d.c
  src/mesh_optimize/call_smooth.c
  src/mesh_optimize/freeoptms2d.c
  src/mesh_optimize/freeoptms.c
  src/mesh_optimize/initoptms2d.c
  src/mesh_optimize/initoptms.c
  src/mesh_optimize/OptMS/SMactive.c
  src/mesh_optimize/OptMS/SMequil.c
  src/mesh_optimize/OptMS/SMerror.c
  src/mesh_optimize/OptMS/SMerror.h
  src/mesh_optimize/OptMS/SMfree.c
  src/mesh_optimize/OptMS/SMfunc2d.c
  src/mesh_optimize/OptMS/SMfunc3d.c
  src/mesh_optimize/OptMS/SMfunction.c
  src/mesh_optimize/OptMS/SMgrad3d.c
  src/mesh_optimize/OptMS/SMinit.c
  src/mesh_optimize/OptMS/SMintrinsic.c
  src/mesh_optimize/OptMS/SMintrinsic.h
  src/mesh_optimize/OptMS/SMlaplace.c
  src/mesh_optimize/OptMS/SMlog_reg.c
  src/mesh_optimize/OptMS/SMmatrix.c
  src/mesh_optimize/OptMS/SMoptimize.c
  src/mesh_optimize/OptMS/SMorient.c
  src/mesh_optimize/OptMS/SMsearch.c
  src/mesh_optimize/OptMS/SMsmooth.c
  src/mesh_optimize/OptMS/SMsmooth.h
  src/mesh_optimize/OptMS/SMstats.c
  src/mesh_optimize/OptMS/SMuntangle.c
  src/mesh_optimize/OptMS/SMuserFunc.c
  src/mesh_optimize/OptMS/SMvalid.c
  src/mesh_optimize/OptMS/SMwritesmooth.c
  src/mesh_optimize/right_handed2d.c
  src/mesh_optimize/right_handed.c
  )

if(FAST_USE_EXODUS)
  list(APPEND morfeus_fv_srcs
    src/mesh/rd_exodus_mesh.f90)
endif()

#----------------------
# Find pre-req packages
#----------------------

# Tell system where to find SI VM installed files
string(TOLOWER "${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION}" SI_postfix)

#----------------------------------------------------------
# Define a helper macro to automate the repetitive tasks
#----------------------------------------------------------

macro( find_morfeus_fv_prereq pkg_name )
  # This macro will find the morfeus_fv pre-requisites if:
  #
  #  1. They are installed into the "normal location", i.e., some system path, as package
  #     managers usually do.
  #  2. They are installed in the Sourcery Institute Linux VM, i.e., in:
  #         `/opt/<pkg-name-lowercase>/<pkg-ver>/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION}`
  #     This ${CMAKE_Fortran_COMPILER_ID} should always be gnu in this case, and currently the version is 8.1.0
  #  3. The USER_LOC_VAR and/or USER_LOC_ENV_VAR keyword arguments are passed to this macro defining the
  #     CMake variable the user should set or the environment variable they should set to find the package in
  #     question, and that CMake is called by the user providing that variable.
  #  4. If no USER_LOC_VAR is passed, a value of ${pkg_name}_ROOT will be assumed, which the user can use to
  #     specify the path to the package.
  #
  # Keyword options (don't take arguments):
  #
  #  * REQUIRED: If passed configuration will fail if the package is not found
  #  * CONFIG: All calls to find_package(...) using/assuming the use of Find modules by default. Set this to
  #    useConfig packages instead.
  #
  # One argument optional Keyword arguments:
  #
  #  * SI_VM_LOC_GLOB: A glob string to find where the package is installed on the SI VM. E.g.,
  #    `/opt/psblas/3.5.0/gnu/8.1.0`. If this is not specified and environment variables don't say where
  #    to look, then the SI VM path will be assumed
  #  * HELP_MSG: A `STATUS` or `WARNING` help message to display to the user if the package is not found
  #    and there are extra tricks to finding it
  #  * USER_LOC_VAR: Name of CMake variable that may be passed by the user to specify the package location
  #    This takes precedence over other locations when used by the user. If no value is specified it is
  #    assumed to be `${pkg_name}_ROOT`
  #  * USER_LOC_ENV_VAR: Name of environment variable that may be set by the user to specify the package
  #    location.
  #
  # Multi argument optional keyword arguments:
  #
  #  * COMPONENTS: Components of package that we need
  #

  # parse arguments
  set(options REQUIRED CONFIG)
  set(oneValueArgs SI_VM_LOC_GLOB HELP_MSG USER_LOC_VAR USER_LOC_ENV_VAR)
  set(multiValueArgs COMPONENTS)
  cmake_parse_arguments( FIND_MFV_REQS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  unset(options)
  unset(oneValueArgs)
  unset(multiValueArgs)

  string(TOLOWER "${pkg_name}" lower_pkg)

  if(FIND_MFV_REQS_REQUIRED)
    set( _required REQUIRED )
    set( _msg_type WARNING )
  else()
    unset( _required )
    set( _msg_type STATUS )
  endif()
  if( "${FIND_MFV_REQS_USER_LOC_VAR}" STREQUAL "" )
    # Call to macro didn't specify a CMake variable name to use for specifying package location
    # Assume a form of ${pkg_name}_ROOT for CMake variable specifying Find Module package locations
    set( FIND_MFV_REQS_USER_LOC_VAR ${pkg_name}_ROOT )
  endif()
  if( NOT (IS_DIRECTORY "${${FIND_MFV_REQS_USER_LOC_VAR}}" OR IS_DIRECTORY "${${FIND_MFV_REQS_USER_LOC_ENV_VAR}}" ) )
    # No user specified package root via CMake or Environment
    if( NOT FIND_MFV_REQS_SI_VM_LOC_GLOB )
      # Macro call did not specify a SI VM Glob location for finding the package
      # Call to macro did not specify a SI VM GLOB
      set( FIND_MFV_REQS_SI_VM_LOC_GLOB "/opt/${lower_pkg}/*/${SI_postfix}")
    endif()
  endif()
  if(FIND_MFV_REQS_SI_VM_LOC_GLOB)
    file( GLOB _${lower_pkg}_root
      LIST_DIRECTORIES true
      "${FIND_MFV_REQS_SI_VM_LOC_GLOB}")
    if( IS_DIRECTORY ${_${lower_pkg}_root} )
      message(STATUS "preferred ${pkg_name} root: ${_${lower_pkg}_root}")
      if(FIND_MFV_REQS_USER_LOC_VAR)
  if( NOT IS_DIRECTORY "${${FIND_MFV_REQS_USER_LOC_VAR}}" )
    set( ${FIND_MFV_REQS_USER_LOC_VAR} "${_${lower_pkg}_root}" )
  endif()
      endif()
      if(FIND_MFV_REQS_USER_LOC_ENV_VAR)
  if( NOT IS_DIRECTORY "$ENV{${FIND_MFV_REQS_USER_LOC_ENV_VAR}}" )
    set( ENV{${FIND_MFV_REQS_USER_LOC_ENV_VAR}} )
  endif()
      endif()
    endif()
  endif()
  if( NOT "${${FIND_MFV_REQS_USER_LOC_VAR}}" STREQUAL "" )
    if( CMAKE_LIBRARY_PATH )
      set( OLD_CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH})
      set( CMAKE_LIBRARY_PATH ${${FIND_MFV_REQS_USER_LOC_VAR}}/lib ${CMAKE_LIBRARY_PATH})
    else()
      set( CMAKE_LIBRARY_PATH ${${FIND_MFV_REQS_USER_LOC_VAR}}/lib )
    endif()
  endif()
  if(FIND_MFV_REQS_COMPONENTS)
    set( _components ${FIND_MFV_REQS_COMPONENTS})
  endif()
  if(NOT FIND_MFV_REQS_CONFIG )
    # Find module based package
    find_package(${pkg_name}
      ${_required} COMPONENTS "${_components}")
  else()
    message ( STATUS "Looking in: ${${FIND_MFV_REQS_USER_LOC_VAR}}" )
    if( NOT "${${FIND_MFV_REQS_USER_LOC_VAR}}" STREQUAL "" )
      set( _paths ${${FIND_MFV_REQS_USER_LOC_VAR}} )
    endif()
    # Config based package (Package provides its own info on how to find it)
    find_package(${pkg_name}
      ${_required}
      CONFIG
      PATHS ${_paths})
  endif()
  if( OLD_CMAKE_LIBRARY_PATH )
    set( CMAKE_LIBRARY_PATH ${OLD_CMAKE_LIBRARY_PATH})
    unset(OLD_CMAKE_LIBRARY_PATH)
  else()
    unset(CMAKE_LIBRARY_PATH)
  endif()
  unset(_paths)
  unset(_components)
  unset(_${lower_pkg}_root)
  unset(_required)
  unset(_msg_type)
  unset(lower_pkg)
endmacro()

#-------------------------------
# Actually find the dependencies
#-------------------------------

# FindMPI has quirks & unique requirements, so we can't use defaults
if(NOT WIN32)
find_morfeus_fv_prereq( MPI
  SI_VM_LOC_GLOB "/opt/mpich/*/${SI_postfix}"
  USER_LOC_VAR MPI_DIR
  USER_LOC_ENV_VAR MPI_DIR)
if(MPI_FOUND)
  list(REMOVE_DUPLICATES MPI_Fortran_INCLUDE_PATH)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
endif()
else()
  set(MPI_FOUND NO)
endif()

# METIS will be made optional
if(NOT APPLE)
  set(BLA_STATIC ON)
endif()

foreach ( _PKG PSBLAS METIS LAPACK BLAS )
  find_morfeus_fv_prereq( ${_PKG} REQUIRED )
endforeach()
include_directories(${PSBLAS_Fortran_INCLUDE_DIR})
include_directories(${METIS_INCLUDES})

if(${FAST_USE_CGNS})
  find_morfeus_fv_prereq( METIS REQUIRED )
  include_directories(${CGNS_INCLUDES})
  add_definitions(-DHAVE_CGNS)

  # HDF5 Find module has some quirky environment vs CMake variable logic
  if(DEFINED ENV{HDF5_ROOT})
    set(HDF5_ROOT $ENV{HDF5_ROOT})
  endif()
  if(DEFINED HDF5_ROOT)
    set(ENV{HDF5_ROOT} ${HDF5_ROOT})
  endif()

  # HDF5 is OPTIONAL, but may be required depending on how other prerequisites were built
  if(WIN32)
    set(HDF5_USE_STATIC_LIBRARIES ON)
  endif()
  find_morfeus_fv_prereq( HDF5
    HELP_MSG "If HDF5 is not found try setting HDF5_ROOT to point to installation"
    COMPONENTS
    Fortran
    Fortran_HL )
  if(HDF5_FOUND)
    add_definitions(${HDF5_Fortran_DEFINITIONS})
    include_directories(${HDF5_Fortran_INCLUDE_DIRS})
  endif()
endif()

if(FAST_USE_EXODUS)
  # NetCDF will be pulled in through SEACASExodus transitive target properties
  find_morfeus_fv_prereq(SEACASExodus_for
    REQUIRED
  CONFIG SI_VM_LOC_GLOB "/opt/seacas-exodus/${SI_postfix}"
  USER_LOC_VAR SEACASExodus_ROOT)
  include_directories(${SEACASExodus_for_INCLUDE_DIRS})
  message(${SEACASExodus_for_INCLUDE_DIRS})
  add_definitions(-DHAVE_EXODUS)
endif()

#----------------------------------------
# Set build dir and install dir locations
#----------------------------------------
set( MORFEUSFV_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}" )
set( MORFEUSFV_MODULE_DIRECTORY "${MORFEUSFV_BUILD_ROOT}" )

#-----------------------------------
# Add the actual MORFEUS_FV library target
#-----------------------------------
add_library(morfeus_fv_C STATIC ${morfeus_fv_C_srcs})
add_library(morfeus_fv STATIC ${morfeus_fv_srcs})
target_compile_options(morfeus_fv
  PUBLIC ${FAST_Fortran_FLAGS})
target_link_libraries(morfeus_fv
  PUBLIC morfeus_fv_C)
if(TARGET morfeus_fd)
  # Support FV dependencies on FD
  target_link_libraries(morfeus_fv
    PUBLIC morfeus_fd)
endif()
target_compile_definitions(morfeus_fv_C
  PRIVATE ${Fortran_INTROSPECTION_MANGLING})
if(TARGET MatLib)
  target_link_libraries(morfeus_fv
    PUBLIC MatLib)
else()
  message( FATAL_ERROR
    "Please run CMake from the top level FAST directory; the MatLib library is now required by MORFEUS_FV.")
endif()

if(TARGET vtkmofo)
  target_link_libraries(morfeus_fv
    PUBLIC vtkmofo)
endif()

if(TARGET jsonfortran-static)
  target_link_libraries(morfeus_fv
    PUBLIC jsonfortran-static)
endif()

if(${FAST_USE_CGNS})
  target_link_libraries(morfeus_fv
    PUBLIC ${CGNS_LIBRARIES})
endif()

target_link_libraries(morfeus_fv
  PUBLIC ${PSBLAS_LIBRARIES}
  PUBLIC ${METIS_LIBRARIES}
  PUBLIC ${LAPACK_LIBRARIES}
  PUBLIC ${BLAS_LIBRARIES}
  )
if(FAST_USE_EXODUS)
  target_link_libraries(morfeus_fv
    PUBLIC ${SEACASExodus_for_LIBRARIES})
endif()
if(HDF5_FOUND)
  target_link_libraries(morfeus_fv
    PUBLIC ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
endif()
if(MPI_Fortran_FOUND)
  target_link_libraries(morfeus_fv
    PUBLIC ${MPI_Fortran_LIBRARIES})
endif()
if(NETCDF_FOUND)
  target_link_libraries(morfeus_fv
    PUBLIC ${NETCDF_F90_LIBRARIES})
endif()
if(MPI_FOUND)
  target_compile_options(morfeus_fv PUBLIC
    ${MPI_Fortran_COMPILE_FLAGS})
  target_link_libraries(morfeus_fv
    PUBLIC ${MPI_Fortran_LINK_FLAGS})
endif()

foreach(lib morfeus_fv morfeus_fv_C)
  target_include_directories(${lib}
    PUBLIC "${CMAKE_Fortran_MODULE_DIRECTORY}"
    PUBLIC "${MORFEUSFV_MODULE_DIRECTORY}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_optimize/OptMS")
  set_property( TARGET ${lib}
    PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
  set_property( TARGET ${lib}
    PROPERTY LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
  set_property( TARGET ${lib}
    PROPERTY FOLDER Libraries/Morfeus_FV)
endforeach()

install(TARGETS morfeus_fv morfeus_fv_C
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include/morfeus_fv)

install(DIRECTORY "${MORFEUSFV_MODULE_DIRECTORY}" DESTINATION include
  FILES_MATCHING PATTERN "*.mod")


function( add_morfeus_fv_exe main_src )
  # Add an executable, linking it against libmorfeus_fv and setting the
  # appropriate properties
  #
  # inputs:
  #    main_src: The Fortran source file containing the program statement/declaration
  #              This should include the file extension (.f90, .f, .F90 etc.)

  get_filename_component( _prog_name "${main_src}" NAME_WE )
  add_executable(${_prog_name} ${main_src})
  target_link_libraries(${_prog_name}
    PRIVATE morfeus_fv)
  set_property(TARGET ${_prog_name}
    PROPERTY RUNTIME_OUTPUT_DIRECTORY "${MORFEUSFV_BUILD_ROOT}/bin")
  set_property(TARGET ${_prog_name}
    PROPERTY FOLDER Morfeus_FV-Tests)
endfunction()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/applications")
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/applications")
endif()
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src/unit-tests")

#------------------------------------------
# Add portable unistall command to makefile
#------------------------------------------
# Adapted from the CMake Wiki FAQ
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake"
    @ONLY)

add_custom_target ( uninstall
  COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake" )

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)
# See JSON-Fortran's CMakeLists.txt file to find out how to get the check target to depend
# on the test executables


enable_testing()

# Test directories
set( unit_test_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/unit")
set( integ_test_dir "${CMAKE_CURRENT_BINARY_DIR}/tests/integration")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/applications")
  set( example_mesh_dir "${CMAKE_CURRENT_SOURCE_DIR}/applications/examples/mesh")
  set( example_inputs_dir "${CMAKE_CURRENT_SOURCE_DIR}/applications/examples/input")
endif()

function(add_morfeus_fv_test name partitions test_target mesh_file input_file workdir)
  # This function will create a test directory for each test, where the specified input file
  # will be copied as `nemo.inp` and `@var_name@` variables will be substituted by values
  # defined by CMake. This allows the input file to be configured to use the mesh in the local
  # directory that we just copied in, and will let multiple different meshes be used with the
  # same input file.
  #
  # inputs:
  #
  #  * name: Test name for CTest
  #  * partitions: Number of MPI_RANKS/CoArray images---the domain will be partitioned into
  #                this many subdomains
  #  * test_target: Target name of the test executable. Usually this is the same as the
  #                 executable name.
  #  * mesh_file: Mesh file for the test. This will be copied to the local test directory in
  #                the build tree.
  #  * input_file: Input file copied and renamed into the test directory. Any `@var@` strings
  #                will be configured, substituting the value of the CMake `${var}` variable,
  #                if define. Use `@mesh_path@` and `@mesh_basename@` in the input file to
  #                allow different meshes.
  #  * workdir: Working directory to run the test in

  # Get mesh file basename
  get_filename_component( mesh_basename ${mesh_file} NAME )
  set( mesh_path "./")
  # First copy and configure the input file and the mesh file to the test working directory
  configure_file( "${input_file}" "${workdir}/nemo.inp"
    @ONLY)
  configure_file( "${mesh_file}" "${workdir}"
    COPYONLY)

  if(TARGET ${test_target})
    get_target_property(min_test_partitions ${test_target} MIN_PARTITIONS)
    get_target_property(build_bin_dir ${test_target} RUNTIME_OUTPUT_DIRECTORY)
  elseif(TARGET build_${test_target})
    get_target_property(min_test_partitions build_${test_target} MIN_PARTITIONS)
    get_target_property(build_bin_dir build_${test_target} RUNTIME_OUTPUT_DIRECTORY)
  endif()
  set( _loc_partitions ${partitions})
  if(min_test_partitions)
    if(partitions LESS min_test_partitions)
      message( FATAL_ERROR "Test ${name} requires ${min_test_partitions} partitions but was only given ${partitions} partitions" )
    endif()
  endif()
  if ( ((N_CPU LESS partitions) OR (N_CPU EQUAL 0)) )
    message(STATUS "Test ${name} is oversubscribed: ${partitions} partitions requested with ${N_CPU} system processor available.")
    if ( OPENMPI )
      if (min_test_partitions)
        set( _loc_partitions ${min_test_partitions} )
      elseif ( N_CPU LESS 2 )
        set( _loc_partitions 2 )
      endif()
      set (pre_flags --oversubscribe)
      message( STATUS
        "Open-MPI detected, passing --oversubscribe for oversubscribed test, ${name}, with ${_loc_partitions} ranks/partitions." )
    endif()
  endif()
  if( MPI_FOUND )
    set(pre_flags ${MPIEXEC_NUMPROC_FLAG} ${_loc_partitions} ${pre_flags} ${MPIEXEC_PREFLAGS})
    set(test_command ${MPIEXEC_EXECUTABLE} ${pre_flags} ${build_bin_dir}/${test_target} ${MPIEXEC_POSTFLAGS} ${ARGN})
    math(EXPR high_rank "${_loc_partitions} - 1")
    set(PASS_REGEX "Process *${high_rank} of *${_loc_partitions}: Normal Termination")
  else()
    set(test_command ${build_bin_dir}/${test_target} ${ARGN})
    set(PASS_REGEX "Process  0 of *1: Normal Termination")
  endif()

  add_test( NAME ${name} COMMAND ${test_command}
    WORKING_DIRECTORY ${workdir})
  set_property(TEST ${name} PROPERTY LABELS "MORFEUS_FV" "integration-test")
  set_property(TEST ${name} PROPERTY PASS_REGULAR_EXPRESSION "${PASS_REGEX}")
endfunction()

# FAST heat conduction multiple zones
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/applications")
  add_morfeus_fv_test( unsteady_heat_pipe 6 unsteady_fast
    applications/fast_heat/cases/test_unsteady/heat_pipe.neu
    applications/fast_heat/cases/test_unsteady/heat_pipe.inp
    ${integ_test_dir}/unsteady_heat_pipe)

  if(${FAST_USE_CGNS})
    # Steady piston conduction test
    # Should produce:
    #   a_sparse.mtx b_rhs.mtx out_nemo.vtk time.dat
    add_morfeus_fv_test( steady_piston 4 steady_conduction
      ${example_mesh_dir}/piston.cgns
      ${example_inputs_dir}/piston.inp
      ${integ_test_dir}/steady_piston)
  endif()
endif()

# # Steady conduction cylinder hybrid gambit neutral file
# add_morfeus_fv_test( steady_cylinder_hyb 4 steady_conduction
#   ${example_mesh_dir}/cylinder_hyb.neu
#   ${example_inputs_dir}/steady_conduction_3d.inp
#   ${integ_test_dir}/steady_cylinder_hyb)

# # Steady conduction cylinder tet only gambit neutral file
# add_morfeus_fv_test( steady_cylinder_tet 4 steady_conduction
#   ${example_mesh_dir}/cylinder_tet.neu
#   ${example_inputs_dir}/steady_conduction_3d.inp
#   ${integ_test_dir}/steady_cylinder_tet)

# # Unsteady slab example
# # Should produce:
# #   out_nemo0{0..9}0.vtk out_nemo100.vtk
# add_morfeus_fv_test( unsteady_slab 4 unsteady_conduction
#   ${example_mesh_dir}/slab.neu
#   ${example_inputs_dir}/unsteady_conduction_3d.inp
#   ${integ_test_dir}/unsteady_slab)

# # Add check mesh test for sphere gambit file
# add_morfeus_fv_test( check_sphere_tet_mesh 4 check_mesh
#   ${example_mesh_dir}/sphere_tet.neu
#   ${example_inputs_dir}/check_mesh_3D.inp
#   ${integ_test_dir}/check_sphere_tet_mesh)

# # Add check mesh test for concentric cylinder gambit file
# add_morfeus_fv_test( check_concentric_cylinder_mesh 4 check_mesh
#   applications/fast_heat/cases/test_unsteady/fast2.neu
#   ${example_inputs_dir}/check_mesh_3D.inp
#   ${integ_test_dir}/check_concentric_cylinder_mesh)

# # Add improve mesh test for concentric cylinder gambit file
# # This test does not pass reliably... may be due to geometry or initialization issue
# add_morfeus_fv_test( improve_concentric_cylinder_mesh 1 improve_mesh
#   applications/fast_heat/cases/test_unsteady/fast2.neu
#   applications/fast_heat/cases/improve_mesh/nemo.inp
#   ${integ_test_dir}/improve_concentric_cylinder_mesh)

# # Add rd_gambit_mesh unit test for sphere gambit file
# add_morfeus_fv_test( test_rd_gambit_mesh 4 test-rd_gambit_mesh
#   applications/fast_heat/cases/test_unsteady/fast2.neu
#   ${example_inputs_dir}/check_mesh_3D.inp
#   ${integ_test_dir}/test-rd_gambit_mesh)
