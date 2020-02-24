cmake_minimum_required(VERSION 3.14.3)
cmake_policy(VERSION 3.14.3...3.16.2)

# Make project modules available
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

# High Sierra is the oldest supported macOS
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13") # sierra is oldest supported macOS

# Organize MSVS projects with folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Disallow in source builds
include(CheckOutOfSourceBuild)
check_out_of_source_build()

# Get information from git, if it repository
include(VersionHelper)
get_vcs_info()
message(STATUS "GIT_COMMIT: ${GIT_FULL_DESCRIBE}")

# Set default configuration types
include(SetConfigTypes)
set_config_types(DEFAULT_DEV_CONFIG Debug CONFIGS RuntimeDebug CodeCoverage)

# Use relative install rpaths and don't use build rpaths: Setup build directory to match install directory
include(RelocatableReproducibleHelper)

# Set some default compiler flags per compiler/language
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config_defaults.cmake)

# Set options for enabling Exodus; support currently broken
include(CMakeDependentOption)
option(MORFEUS_USE_FV "Build MORFEUS with finite volume PDE solver support" ON)
set(MORFEUS_FV_DIR "submodules/morfeus/src/FV")
option(MORFEUS_USE_FD "Build MORFEUS with finite difference PDE solver support" ON)
set(MORFEUS_FD_DIR "submodules/morfeus/src/FD")
CMAKE_DEPENDENT_OPTION(MORFEUS_USE_EXODUS "Build with Exodus support (broken)" OFF
  "MORFEUS_USE_FV" OFF) # Switch to on once exodus is sorted out

#----------------------------------------------------------------------------------------------------------
# SuperBuild: Configure, build (and sometimes test) third party libraries (TPLs) and then configure MORFEUS
#----------------------------------------------------------------------------------------------------------
# See https://github.com/Sarcasm/cmake-superbuild and
# https://blog.kitware.com/cmake-superbuilds-git-submodules/
option (USE_SUPERBUILD "Whether or not to perform the outer SuperBuild (config and build TPLS)" ON)
if(NOT APPLE)
  set(BLA_STATIC ON)
endif()
if (USE_SUPERBUILD)
  project (SUPERBUILD Fortran C CXX)
  # if(isMultiConfig)
  #   if(NOT "RuntimeDebug" IN_LIST CMAKE_CONFIGURATION_TYPES)
  #     list(APPEND CMAKE_CONFIGURATION_TYPES RuntimeDebug)
  #   endif()
  # endif()
  # add_missing_build_configs() # in PlatformAbstraction.cmake
  # set_property(GLOBAL APPEND PROPERTY DEBUG_CONFIGURATIONS RuntimeDebug)
  # execute the superbuild (this CMakeLists.txt will be called again without the
  # USE_SUPERBUILD option this time)
  set( TPL_DIR ${CMAKE_BINARY_DIR}/TPLs)
  include (cmake/SuperBuild.cmake)
  return() # Stop processing, recursively called again from SuperBuild.cmake
endif()

project(
  Morfeus
  VERSION 0.0.1
  DESCRIPTION "Multi-physics Object-oriented Reconfigurable Fluid Environment for Unified Simulations"
  HOMEPAGE_URL "https://github.com/sourceryinstitute/MORFEUS-Source"
  LANGUAGES Fortran C CXX
  )

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()


# add_subdirectory(src)

include(UninstallTarget)
