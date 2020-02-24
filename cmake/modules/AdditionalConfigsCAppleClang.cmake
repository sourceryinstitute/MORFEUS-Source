set(CMAKE_C_FLAGS_RUNTIMEDEBUG
    "-g -fexceptions -fstack-check"
    CACHE STRING "Runtime debugging checks flags"
    )
# Clang sources don't want to link into Fortran executables The issue is that the right library and search path are not available to
# the fortran linker See https://stackoverflow.com/a/30313514/403516

set(CMAKE_C_FLAGS_CODECOVERAGE
    "-g"
    CACHE STRING "Code coverage flags"
    )
set(CMAKE_EXE_LINKER_FLAGS_CODECOVERAGE
    "--coverage"
    CACHE STRING "Exe linker flags"
    )
