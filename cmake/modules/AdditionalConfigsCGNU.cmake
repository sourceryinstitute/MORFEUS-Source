if(CMAKE_C_COMPILER_VERSION VERSION_GREATER 8.0)
  set(CMAKE_C_FLAGS_RUNTIMEDEBUG
      "-g -fwrapv -fno-wrapv-pointer -fexceptions -fstack-check"
      CACHE STRING "Runtime debugging checks flags"
      )
else()
  set(CMAKE_C_FLAGS_RUNTIMEDEBUG
      "-g -fexceptions -fstack-check"
      CACHE STRING "Runtime debugging checks flags"
      )
endif()

set(CMAKE_C_FLAGS_CODECOVERAGE
    "-O0 -g --coverage"
    CACHE STRING "Code coverage flags"
    )
set(CMAKE_EXE_LINKER_FLAGS_CODECOVERAGE
    "--coverage"
    CACHE STRING "Exe linker flags"
    )
