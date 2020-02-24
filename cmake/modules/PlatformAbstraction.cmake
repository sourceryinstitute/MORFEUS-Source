# GFortran build configs
set(GNU_Fortran_FLAGS_RUNTIMEDEBUG "-g -fcheck=all")
set(GNU_C_FLAGS_RUNTIMEDEBUG "-g -fstack-check")
set(GNU_EXE_LINKER_FLAGS "-fcheck=all")

# Intel compilers are different on Windows and *nix
include_guard(DIRECTORY)
if (WIN32)
  set(prefix "/")
  set(infix ":")
  set(Qf "Q")
  set(Q "Q")
  set(eq ":")
  set(colon ":")
  set(colon_ ":")
else()
  set(prefix "-")
  set( infix " ")
  set( Qf "f")
  set( Q "")
  set( eq "=")
  set( colon "")
  set( colon_ " ")
endif()

set(Intel_Fortran_FLAGS_RUNTIMEDEBUG "${prefix}check${colon_}bounds ${prefix}check${colon_}stack")
set(Intel_EXE_LINKER_FLAGS_RUNTIMEDEBUG "${prefix}check${colon_}bounds ${prefix}check${colon_}stack")

set(ADDITIONAL_FAST_BUILD_CONFIGS RUNTIMEDEBUG)
macro(add_missing_build_configs)
  foreach(_config ${ADDITIONAL_FAST_BUILD_CONFIGS})
    if("${CMAKE_Fortran_FLAGS_${_config}}" STREQUAL "")
      set(CMAKE_Fortran_FLAGS_${_config}
	"${${CMAKE_Fortran_COMPILER_ID}_Fortran_FLAGS_${_config}}"
	CACHE STRING "Additional RuntimeDebug flags" FORCE)
    endif()
    if("${CMAKE_C_FLAGS_${_config}}" STREQUAL "")
      set(CMAKE_C_FLAGS_${_config}
	"${${CMAKE_C_COMPILER_ID}_C_FLAGS_${_config}}"
	CACHE STRING "Additional RuntimeDebug flags" FORCE)
    endif()
    if("${CMAKE_EXE_LINKER_FLAGS_${_config}}" STREQUAL "")
      set(CMAKE_EXE_LINKER_FLAGS_${_config}
	"${${CMAKE_Fortran_COMPILER_ID}_EXE_LINKER_FLAGS_${_config}}"
	CACHE STRING "Additional RuntimeDebug flags" FORCE)
    endif()
  endforeach()
  if(MSVC)
    if("${CMAKE_C_FLAGS_RUNTIMEDEBUG}" STREQUAL "")
      set(CMAKE_C_FLAGS_RUNTIMEDEBUG "${CMAKE_C_FLAGS_DEBUG}")
    endif()
  endif()
endmacro()
