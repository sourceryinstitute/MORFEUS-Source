# This file controls the superbuild of all dependencies
#
####################### !!!!!!!!!WARNING!!!!!!!! ########################
#                                                                       #
# If you update downloaded software, make sure to update (and include)  #
# a cryptographicaly secure hash of the file to download! Also, if      #
# possible, an ssl encrypted URL should be used (https://...) and if    #
# the developers provide GPG/cryptographic signatures or code signing,  #
# then it would be a good idea to verify those __*BEFORE*__ you add the #
# new archive __*AND*__ its hash to this file.                          #
#                                                                       #
# The cryptographic hash ensures everyone gets the same file you got,   #
# and that no corruption or man-in-the-middle (MITM) attack occurred    #
# during the download.                                                  #
#                                                                       #
# The ssl encrypted URL also helps to prevent MITM attacks.             #
#                                                                       #
# Verifying the GPG signature or other cryptographic signature or code  #
# signing helps to ensure that the software you are updating was        #
# actually provided by the developers of that software, and that an     #
# attacker has not compromised the server hosting the software and      #
# replaced the official version with a malicious one.                   #
#                                                                       #
# At a very minimum, any software that is not "vendored" or included    #
# via a git submodule __*MUST*__ include a cryptographically secure     #
# hash. (SHA256 or better; md5 and SHA1 are not secure. SHA512 should   #
# be relatively future proof and SHA256 is the current gold standard.)  #
#                                                                       #
#########################################################################

message(STATUS
  "Configuring and building third party libraries (TPLs) required by MORFEUS.")

include (ExternalProject)

set_property (DIRECTORY PROPERTY EP_PREFIX ${TPL_DIR})
# See https://cmake.org/cmake/help/v3.11/module/ExternalProject.html
# This property causes ExternalProject_Add to create the following directory structure:
# ```
# TMP_DIR      = <prefix>/tmp
# STAMP_DIR    = <prefix>/src/<name>-stamp
# DOWNLOAD_DIR = <prefix>/src
# SOURCE_DIR   = <prefix>/src/<name>
# BINARY_DIR   = <prefix>/src/<name>-build
# INSTALL_DIR  = <prefix>
# ```

set (DEPENDENCIES)
set (EXTRA_CMAKE_ARGS)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

list(GET BLAS_LIBRARIES 0 FIRST_BLAS_LIB)
list(GET LAPACK_LIBRARIES 0 FIRST_LAPACK_LIB)
get_filename_component(BLAS_DIR_HOME ${FIRST_BLAS_LIB}
	DIRECTORY
	CACHE)
get_filename_component(LAPACK_DIR_HOME ${FIRST_LAPACK_LIB}
	DIRECTORY
	CACHE)
message(STATUS "BLAS found in: ${BLAS_DIR_HOME}")
message(STATUS "LAPACK found in: ${LAPACK_DIR_HOME}")

# Pass common options to all dependencies:
set (DEFAULT_CACHE_ARGS
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DCMAKE_INSTALL_PREFIX:PATH=${TPL_DIR}
  -DCMAKE_PREFIX_PATH:PATH=${TPL_DIR};${BLAS_DIR_HOME};${LAPACK_DIR_HOME})

if(DEFINED CMAKE_BUILD_TYPE AND (NOT ${CMAKE_BUILD_TYPE} STREQUAL ""))
  list(APPEND DEFAULT_CACHE_ARGS -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})
endif()

file(TO_NATIVE_PATH "${TPL_DIR}" NATIVE_TPL_DIR)

if( CMAKE_SYSTEM_NAME MATCHES "[Ll]inux|[Dd]arwin" )
  list(APPEND DEPENDENCIES hdf5)
  list(APPEND EXTRA_CMAKE_ARGS -DHDF5_ROOT:PATH=${TPL_DIR})
  set(HDF5_ROOT ${TPL_DIR})
  set(ENV{HDF5_ROOT} ${TPL_DIR})
  ExternalProject_Add( hdf5
    URL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.6/src/hdf5-1.10.6.tar.gz"
    URL_HASH SHA256=5f9a3ee85db4ea1d3b1fa9159352aebc2af72732fc2f58c96a3f0768dba0e9aa
    TLS_VERIFY ON
    SOURCE_DIR ${TPL_DIR}/src/hdf5
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND cd ${TPL_DIR}/src/hdf5 && autoreconf -fiv
    COMMAND ${TPL_DIR}/src/hdf5/configure --enable-build-mode=production --enable-optimization=debug --enable-dependency-tracking --enable-fortran --enable-hl --enable-static --enable-static-exec --disable-parallel --disable-preadwrite --disable-java --with-pic --with-default-api-version=v18 --prefix=${TPL_DIR}
    BUILD_IN_SOURCE ON
    TEST_BEFORE_INSTALL ON
    TEST_COMMAND make check
    INSTALL_COMMAND make install
    TEST_EXCLUDE_FROM_MAIN ON
    )
  ExternalProject_Add_StepTargets(hdf5 test)
endif()


list(APPEND DEPENDENCIES netcdf4)
ExternalProject_Add( netcdf4
  DEPENDS hdf5
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/netcdf4
  CMAKE_CACHE_ARGS ${DEFAULT_CACHE_ARGS} -DBUILD_TESTING:BOOL=OFF -DENABLE_TESTS:BOOL=OFF -DENABLE_NETCDF_4:BOOL=ON -DENABLE_DOXYGEN:BOOL=OFF -DENABLE_DAP:BOOL=OFF -DENABLE_EXAMPLES:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=OFF -DENABLE_V2_API:BOOL=OFF -DENABLE_CDF5:BOOL=ON -DHDF5_INCLUDE_DIR:PATH=${TPL_DIR}/include
  BUILD_IN_SOURCE OFF
  TEST_BEFORE_INSTALL ON)

if( CMAKE_SYSTEM_NAME MATCHES "[Ll]inux|[Dd]arwin" )
  list(APPEND DEPENDENCIES exodus)
  ExternalProject_Add( exodus
    DEPENDS netcdf4 hdf5
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/exodus
    BUILD_IN_SOURCE OFF
    CONFIGURE_COMMAND INSTALL_PATH=${TPL_DIR} ACCESS=${CMAKE_SOURCE_DIR}/exodus SHARED=OFF COMPILER=EXTERNAL FORTRAN=ON STATIC=ON ${CMAKE_SOURCE_DIR}/exodus/cmake-exodus
    TEST_BEFORE_INSTALL ON
    TEST_COMMAND ${CMAKE_COMMAND} --build . --target test
    INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install
    )
endif()

# ----
# hdf5
# ----
if(MORFEUS_USE_CGNS)
  option(UNIX_BUILD_HDF5_FROM_SOURCE "Always build hdf5 from source on unix-like platforms" OFF)
  set(HDF5_VER 1.10.1) # Update URLs/hashes if you change this
  list(APPEND DEPENDENCIES hdf5)
  list(APPEND EXTRA_CMAKE_ARGS -DHDF5_ROOT:PATH=${NATIVE_TPL_DIR})
  set(HDF5_ROOT ${NATIVE_TPL_DIR})
  set(env{HDF5_ROOT} ${NATIVE_TPL_DIR})
  # Fetch Linux binary, build from source on macOS
  if(UNIX_BUILD_HDF5_FROM_SOURCE)
    # Compile hdf5 from source
    ExternalProject_Add( hdf5
      URL https://gamma.hdfgroup.org/ftp/pub/outgoing/QATEST/hdf5110/source/CMake-hdf5-${HDF5_VER}.tar.gz
      URL_HASH SHA256=b8b8cb7b99e83cd9dfdf2cbebaf253349f2f4435868c9a65dfd93584d9c63c06
      TLS_VERIFY ON
      BUILD_IN_SOURCE ON
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ${CMAKE_CTEST_COMMAND} -S HDF5config.cmake,BUILD_GENERATOR=Unix,INSTALLDIR=${TPL_DIR},HDF5_BUILD_FORTRAN=ON,BUILD_SHARED_LIBS=OFF,LOCAL_SKIP_TEST=ON,BUILD_TESTING=OFF -V -O hdf5.log
      INSTALL_COMMAND make -C build install
      TEST_COMMAND ""
      )
  elseif(APPLE)
    ExternalProject_Add( hdf5
      URL https://gamma.hdfgroup.org/ftp/pub/outgoing/QATEST/hdf5110/binaries/unix/hdf5-${HDF5_VER}-Std-osx1011_64.tar.gz
      URL_HASH SHA256=d484a9db37f5929674e67456f55f5bd37a267776d63fc0592d756ce77f015670
      TLS_VERIFY ON
      BUILD_IN_SOURCE ON
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ${TPL_DIR}/src/hdf5/HDF5-${HDF5_VER}-Darwin.sh --prefix=${TPL_DIR} --skip-license --exclude-subdir
      INSTALL_COMMAND bash -c "cp -r ${TPL_DIR}/HDF_Group/HDF5/${HDF5_VER}/* ${TPL_DIR}"
      TEST_COMMAND ""
      )
  elseif(UNIX)
    message("Warning: Installing pre-compiled HDF5 binaries intended for Linux! (RHEL/Centos7)")
    ExternalProject_Add( hdf5
      URL https://gamma.hdfgroup.org/ftp/pub/outgoing/QATEST/hdf5110/binaries/unix/hdf5-${HDF5_VER}-Std-centos7_64.tar.gz
      URL_HASH SHA256=2041a34af5047f684eef15f3d23726aeb3273b0dc063995e52442428bfeae82a
      TLS_VERIFY ON
      BUILD_IN_SOURCE ON
      CONFIGURE_COMMAND ""
      BUILD_COMMAND cp -r . ${TPL_DIR}
      TEST_COMMAND ""
      INSTALL_COMMAND ""
      )
  elseif(WIN32) # Windows
    # https://stackoverflow.com/questions/8839978/install-msi-with-msiexec-in-a-specific-directory
    # ^ Discussion on how to specify installation path for MSI installers TL;DR property name depends
    # on package; for hdf5 this is INSTALL_ROOT.
    file(TO_NATIVE_PATH "${TPL_DIR}" WIN_TPL_DIR)
    ExternalProject_Add( hdf5
      URL https://gamma.hdfgroup.org/ftp/pub/outgoing/QATEST/hdf5110/binaries/windows/hdf5-${HDF5_VER}-Std-win10_64-vs14.zip
      URL_HASH SHA256=e9514d057303bd90109696039554db7cad745702c8e442c91c30654894d20b14
      TLS_VERIFY ON
      BUILD_IN_SOURCE ON
      CONFIGURE_COMMAND ""
      BUILD_COMMAND msiexec /i HDF5-${HDF5_VER}-win64.msi INSTALL_ROOT=${WIN_TPL_DIR} /quiet /qn /norestart /lv install.log
      TEST_COMMAND ""
      INSTALL_COMMAND ""
      )
  endif()
endif()

# ------
# NetCDF
# ------
if(MORFEUS_USE_EXODUS)
list(APPEND DEPENDENCIES netCDF-4)
ExternalProject_Add( netCDF-4
  URL https://github.com/Unidata/netcdf-c/archive/v4.6.1.tar.gz
  URL_HASH SHA256=a2fabf27c72a5ee746e3843e1debbaad37cd035767eaede2045371322211eebb
  TLS_VERIFY ON
  CMAKE_CACHE_ARGS ${DEFAULT_CACHE_ARGS} -DENABLE_NETCDF_4:BOOL=ON -DENABLE_DAP_REMOTE_TESTS:BOOL=OFF -DENABLE_DAP:BOOL=OFF -DENABLE_EXAMPLES:BOOL=OFF -DBUILD_UTILITIES:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=OFF -DENABLE_TESTS:BOOL=ON
  BUILD_IN_SOURCE OFF
  TEST_BEFORE_INSTALL ON)
endif()

# ----
# CGNS
# ----
if(${MORFEUS_USE_CGNS})
  list(APPEND DEPENDENCIES CGNS)
  # See http://cgns.github.io/FAQs.html (Installation Questions--Outdated as of 2018-09-14)
  # See also https://github.com/CGNS/CGNS/blob/4858d3852453770d2de63bcf8b0e5568daa362e4/appveyor.yml
  list(APPEND CGNS_CONFIG
    ${DEFAULT_CACHE_ARGS}
    -DCGNS_ENABLE_FORTRAN:BOOL=ON
    -DCGNS_BUILD_SHARED:BOOL=OFF
    -DCGNS_USE_SHARED:BOOL=OFF
    -DCGNS_ENABLE_64BIT:BOOL=ON
    -DCGNS_BUILD_CGNSTOOLS:BOOL=OFF
    -DCGNS_ENABLE_TESTS:BOOL=ON
    -DCGNS_ENABLE_HDF5:BOOL=ON
    -DCGNS_ENABLE_PARALLEL:BOOL=OFF
    -DCGNS_ENABLE_LFS:BOOL=OFF
    -DCGNS_ENABLE_SCOPING:BOOL=OFF
    -DHDF5_NEED_ZLIB:BOOL=ON
    )
  if(APPLE)
    list(APPEND CGNS_CONFIG
      -DCMAKE_STATIC_LINKER_FLAGS:FILEPATH=${TPL_DIR}/lib/libhdf5.a)
  endif()
  message(STATUS "Flags passed during CGNS configuration:
    ${CGNS_CONFIG}")
  ExternalProject_Add( CGNS
    DEPENDS hdf5
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/CGNS
    CMAKE_CACHE_ARGS ${CGNS_CONFIG}
    BUILD_IN_SOURCE OFF
    TEST_BEFORE_INSTALL ON)
  #  TEST_COMMAND "")
  list(APPEND EXTRA_CMAKE_ARGS -DCGNS_ROOT:PATH=${TPL_DIR})
endif()

if(MORFEUS_USE_EXODUS)
# -------------
# SEACAS-Exodus
# -------------
list(APPEND DEPENDENCIES SEACASExodus)
list(APPEND SEACAS_EXODUS_CONFIG
    ${DEFAULT_CACHE_ARGS}
    -DSEACASProj_ENABLE_SEACASExodus:BOOL=ON
    -DSEACASProj_ENABLE_SEACASExodus_for:BOOL=ON
    -DSEACASProj_ENABLE_SEACASExoIIv2for32:BOOL=ON
    -DSEACASProj_HIDE_DEPRECATED_CODE:BOOL=OFF
    -DSEACASProj_ENABLE_TESTS:BOOL=ON
    -DCMAKE_INSTALL_RPATH:PATH=${TPL_DIR}/lib
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DSEACASProj_ENABLE_Fortran:BOOL=ON
    -DTPL_ENABLE_Netcdf:BOOL=ON
    -DTPL_ENABLE_MPI:BOOL=OFF
    -DTPL_ENABLE_Pthread:BOOL=OFF
    -DSEACASExodus_ENABLE_THREADSAFE:BOOL=OFF
    -DNetCDF_ROOT:PATH=${TPL_DIR}
    -DHDF5_ROOT:PATH=${TPL_DIR}
    -DHDF5_NO_SYSTEM_PATHS:BOOL=ON
    -DSEACASProj_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON
    -DCMAKE_MODULE_PATH:PATH=${TPL_DIR}/src/SEACASExodus/cmake/tribits/common_tpls/find_modules
)

ExternalProject_Add( SEACASExodus
  DEPENDS netCDF-4 hdf5
  # GIT_REPOSITORY https://github.com/gsjaardema/seacas.git
  # GIT_TAG f32b74ac27df180823391f2670d624987247a9b7
  URL https://github.com/gsjaardema/seacas/archive/exodus.zip
  URL_HASH SHA256=99b981f745cce1f9d9ea5ef7ea33929044d71bd9774a8b981691b1d1ad66ac31
  TLS_VERIFY ON
  DOWNLOAD_NAME seacas-exodus.zip
  CMAKE_CACHE_ARGS ${SEACAS_EXODUS_CONFIG}
  BUILD_IN_SOURCE OFF
  TEST_BEFORE_INSTALL ON)
list(APPEND EXTRA_CMAKE_ARGS -DSEACASExodus_ROOT:PATH=${TPL_DIR})
endif()

# -----
# METIS
# -----
if(MORFEUS_USE_FV)
list(APPEND DEPENDENCIES METIS)
ExternalProject_Add( METIS
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/metis
  CMAKE_CACHE_ARGS -DBUILD_SHARED_LIBS:BOOL=OFF -DCMAKE_INSTALL_PREFIX:PATH=${TPL_DIR}
  BUILD_IN_SOURCE OFF
  TEST_BEFORE_INSTALL ON
  TEST_COMMAND "")
list(APPEND EXTRA_CMAKE_ARGS -DMETIS_ROOT:PATH=${TPL_DIR})
endif()

if(NOT WIN32)
  set(MPI_SKIP_GUESSING ON)
endif()
set(MPI_CXX_SKIP_MPICXX ON)
if(NOT WIN32)
  find_package(MPI)
endif()

if(NOT MPI_Fortran_FOUND)
  if(NOT WIN32)
    if(MORFEUS_USE_FV OR MORFEUS_USE_FD)
      ExternalProject_Add( MPICH
	URL "https://www.mpich.org/static/downloads/3.3/mpich-3.3.tar.gz"
	URL_HASH SHA256=329ee02fe6c3d101b6b30a7b6fb97ddf6e82b28844306771fa9dd8845108fa0b
	TLS_VERIFY ON
	SOURCE_DIR ${TPL_DIR}/src/MPICH
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ${TPL_DIR}/src/MPICH/configure --disable-cxx --enable-fortran=all --prefix=${TPL_DIR}
	BUILD_IN_SOURCE OFF
	TEST_BEFORE_INSTALL ON
	TEST_COMMAND make check
	)
      set(OC_DEPENDS MPICH)
      set(PSBLAS_DEPENDS MPICH)
      list(APPEND DEPENDENCIES MPICH)
    endif()
  endif()
elseif()
  set(OC_MPI_COMPILER="${MPI_Fortran_COMPILER}")
endif()

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL GNU)
  if(MORFEUS_USE_FD)
    ExternalProject_Add( OpenCoarrays
      DEPENDS "${OC_DEPENDS}"
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/OpenCoarrays
      CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${TPL_DIR} -DCAF_ENABLE_FAILED_IMAGES:BOOL=FALSE -DCMAKE_PREFIX_PATH:PATH=${TPL_DIR};${BLAS_DIR_HOME};${LAPACK_DIR_HOME}
      BUILD_IN_SOURCE OFF
      TEST_BEFORE_INSTALL ON
      )
    list(APPEND DEPENDENCIES OpenCoarrays)
  endif()
endif()

# PSBLAS
# Git as the source/upstream/download strategy fails on windows b/c
# Windows doesn't allow directories/files named aux
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL GNU AND MORFEUS_USE_FD)
  list(APPEND PSBLAS_DEPENDS OpenCoarrays)
  list(APPEND PSB_EXTRA_CFG_ARGS "-DPSBLAS_USE_OpenCoarrays:BOOL=ON")
endif()
if(MORFEUS_USE_FV)
  list(APPEND DEPENDENCIES PSBLAS)
  ExternalProject_Add( PSBLAS
    DEPENDS METIS ${PSBLAS_DEPENDS}
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/psblas
    CMAKE_CACHE_ARGS -DSHARED:BOOL=OFF ${DEFAULT_CACHE_ARGS} ${PSB_EXTRA_CFG_ARGS}
    BUILD_IN_SOURCE OFF
    TEST_BEFORE_INSTALL ON
    TEST_COMMAND "")
  list(APPEND EXTRA_CMAKE_ARGS -DPSBLAS_ROOT:PATH=${TPL_DIR})
endif()

# -------
# VTKmofo
# -------
list(APPEND DEPENDENCIES VTKmofo)
list(APPEND EXTRA_CMAKE_ARGS -DVTKmofo_ROOT:PATH=${TPL_DIR})
if((CMAKE_Fortran_COMPILER_ID MATCHES GNU) AND MORFEUS_USE_FD)
  set(VTKmofo_DEPENDS OpenCoarrays)
  set(VTKmofo_EXTRA_CFG_ARGS "-DVTKmofo_USE_OpenCoarrays:BOOL=ON")
endif()

ExternalProject_Add( VTKmofo
  DEPENDS "${VTKmofo_DEPENDS}"
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/vtkmofo
  CMAKE_CACHE_ARGS ${DEFAULT_CACHE_ARGS} ${VTKmofo_EXTRA_CFG_ARGS}
  BUILD_IN_SOURCE OFF
  TEST_BEFORE_INSTALL ON
  )

# ------------
# json-fortran
# ------------
if( ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL GNU) AND ("${CMAKE_Fortran_COMPILER_VERSION}" MATCHES 8.1) OR WIN32)
  message( STATUS "Skipping JSON-Fortran test jf_test_15 due to compiler bug in GFortran 8.1")
  set(JSON_SKIP_TEST_REGEX jf_test_15)
else()
  set(JSON_SKIP_TEST_REGEX non-existant-test)
endif()

list(APPEND DEPENDENCIES json-fortran)
list(APPEND EXTRA_CMAKE_ARGS -DJSON_Fortran_ROOT:PATH=${TPL_DIR})
if((CMAKE_Fortran_COMPILER_ID MATCHES GNU) AND MORFEUS_USE_FD)
  set(JSON_FORTRAN_DEPENDS OpenCoarrays)
  set(JSON_FORTRAN_EXTRA_CFG_ARGS "-DJSON_FORTRAN_USE_OpenCoarrays:BOOL=ON")
endif()
ExternalProject_Add( json-fortran
  DEPENDS "${JSON_FORTRAN_DEPENDS}"
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/json-fortran
  CMAKE_CACHE_ARGS ${DEFAULT_CACHE_ARGS} ${JSON_FORTRAN_EXTRA_CFG_ARGS} -DSKIP_DOC_GEN:BOOL=ON -DUSE_GNU_INSTALL_CONVENTION:BOOL=TRUE
  BUILD_IN_SOURCE OFF
  TEST_BEFORE_INSTALL ON
  # build_tests target depends on libjson_fortran and libjson_fortran-static
  BUILD_COMMAND ${CMAKE_COMMAND} --build . --target build_tests
  TEST_COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIG> -E ${JSON_SKIP_TEST_REGEX}
  )

if(DEFINED CACHE{MORFEUS_USE_FV})
  list(APPEND EXTRA_CMAKE_ARGS -DMORFEUS_USE_FV:BOOL=${MORFEUS_USE_FV})
endif()
if(DEFINED CACHE{MORFEUS_USE_FD})
  list(APPEND EXTRA_CMAKE_ARGS -DMORFEUS_USE_FD:BOOL=${MORFEUS_USE_FD})
endif()

ExternalProject_ADD( EP_MORFEUS
  DEPENDS ${DEPENDENCIES}
  SOURCE_DIR ${CMAKE_SOURCE_DIR}
  CMAKE_CACHE_ARGS ${DEFAULT_CACHE_ARGS} -DUSE_SUPERBUILD:BOOL=OFF ${EXTRA_CMAKE_ARGS}
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_BINARY_DIR}/MORFEUS
  )
ExternalProject_Add_StepTargets(EP_MORFEUS configure)
