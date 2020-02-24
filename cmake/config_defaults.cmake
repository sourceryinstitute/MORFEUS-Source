if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # These are always the same whether on windows with *nix emulation or native *nix
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 8.1)
    set(Fortran_STD f2008ts)
    set(PREFIX_MAP_FLAG -fdebug-prefix-map)
  else()
    set(Fortran_STD f2018)
    set(PREFIX_MAP_FLAG -ffile-prefix-map)
  endif()

  set(GNU_Fortran_FLAGS
      -std=${Fortran_STD}
      -pedantic
      -fbacktrace
      -ffree-line-length-none
      -Wall
      -Wextra
      -fno-working-directory
      ${PREFIX_MAP_FLAG}=${CMAKE_SOURCE_DIR}=.
      )
  list(JOIN GNU_Fortran_FLAGS " " CMAKE_Fortran_FLAGS_INIT)

  include(AdditionalConfigsFortranGNU)
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU|AppleClang")
  if(CMAKE_C_COMPILER_ID MATCHES "AppleClang" OR CMAKE_C_COMPILER_VERSION VERSION_LESS 8.1)
    set(PREFIX_MAP_FLAG -fdebug-prefix-map)
  else()
    set(PREFIX_MAP_FLAG -ffile-prefix-map)
  endif()

  set(GNU_C_FLAGS -Wall -Wextra -Wpedantic ${PREFIX_MAP_FLAG}=${CMAKE_SOURCE_DIR}=.)
  if(NOT CMAKE_C_COMPILER_ID MATCHES "AppleClang")
    list(APPEND GNU_C_FLAGS -fno-working-directory)
  endif()
  list(JOIN GNU_C_FLAGS " " CMAKE_C_FLAGS_INIT)
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  include(AdditionalConfigsCGNU)
elseif(CMAKE_C_COMPILER_ID MATCHES "AppleClang")
  include(AdditionalConfigsCAppleClang)
endif()
