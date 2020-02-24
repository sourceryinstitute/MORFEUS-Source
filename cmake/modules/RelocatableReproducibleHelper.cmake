include_guard(GLOBAL)

# Use relative paths for the dynamic linker and setup build directory to mimic install directory
set(CMAKE_SKIP_RPATH
    ON
    CACHE BOOL "Don't add a build-dir rpath"
    )
set(CMAKE_BUILD_WITH_INSTALL_RPATH
    ON
    CACHE BOOL "Build using the install rpath"
    )
set(CMAKE_BUILD_RPATH_USE_ORIGIN
    ON
    CACHE BOOL "Use relative rpaths"
    )

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

if(NOT DEFINED ENV{SOURCE_DATE_EPOCH})
  if(${COMMIT_DATE_UNIX_TIME})
    message(WARNING "SOURCE_DATE_EPOCH was not defined in the build environment and the build will not be reproducible."
                    "\n`export SOURCE_DATE_EPOCH=${COMMIT_DATE_UNIX_TIME}`"
            )
  else()
    message(WARNING "SOURCE_DATE_EPOCH was not defined in the build environment and the build will not be reproducible."
                    "\n`export SOURCE_DATE_EPOCH=$(git show -s --format=%ct HEAD)` or\n`export SOURCE_DATE_EPOCH=$(date +%s)`\n"
            )
  endif()
endif()

if(DEFINED ENV{SOURCE_DATE_EPOCH})
  set(SOURCE_DATE_EPOCH "$ENV{SOURCE_DATE_EPOCH}")
elseif(DEFINED SOURCE_DATE_EPOCH OR DEFINED CACHE{SOURCE_DATE_EPOCH})
  set(ENV{SOURCE_DATE_EPOCH} ${SOURCE_DATE_EPOCH})
elseif(DEFINED COMMIT_DATE_UNIX_TIME)
  set(SOURCE_DATE_EPOCH ${COMMIT_DATE_UNIX_TIME})
  set(ENV{SOURCE_DATE_EPOCH} ${COMMIT_DATE_UNIX_TIME})
else()
  string(TIMESTAMP SOURCE_DATE_EPOCH %s)
  set(ENV{SOURCE_DATE_EPOCH} ${SOURCE_DATE_EPOCH})
endif()
