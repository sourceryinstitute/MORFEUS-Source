set(CMAKE_Fortran_FLAGS_RUNTIMEDEBUG
    "-g -fcheck=all -ffpe-trap=invalid,zero,overflow"
    CACHE STRING "Runtime debugging checks flags"
    )

if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 8.0)
  set(CMAKE_EXE_LINKER_FLAGS_RUNTIMEDEBUG
      "-fcheck=all -ffpe-trap=invalid,zero,overflow -fwrapv -fno-wrapv-pointer -fexceptions -fstack-check"
      CACHE STRING "Exe linker flags"
      )
else()
  set(CMAKE_EXE_LINKER_FLAGS_RUNTIMEDEBUG
      "-fcheck=all -ffpe-trap=invalid,zero,overflow -fexceptions -fstack-check"
      CACHE STRING "Exe linker flags"
      )
endif()

set(CMAKE_Fortran_FLAGS_CODECOVERAGE
    "-O0 -g --coverage"
    CACHE STRING "Code coverage flags"
    )
set(CMAKE_EXE_LINKER_FLAGS_CODECOVERAGE
    "--coverage"
    CACHE STRING "Exe linker flags"
    )
