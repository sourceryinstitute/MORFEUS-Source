# Introspection CMake macros for Fortran language
include_guard(DIRECTORY)

#-------------------------------------
# Fortran name mangling intro-spection
#-------------------------------------
# returns Fortran_INTROSPECTION_MANGLING variable
macro(fortran_introspection_mangling)
  # Get capitalization and name mangling underscore conventions
  include(CapitalizeString)
  include(FortranCInterface)
  set(Fortran_INTROSPECTION_MANGLING NOTFOUND)
  capitalize_string(${FortranCInterface_GLOBAL__CASE} fc_case)
  message(STATUS "Name mangling capitalization: ${fc_case}")
  message(STATUS "Name mangling fortran global suffix underscore: ${FortranCInterface_GLOBAL__SUFFIX}")
  if(FortranCInterface_GLOBAL__SUFFIX STREQUAL "")
    set(Fortran_INTROSPECTION_MANGLING "${fc_case}Case")
  elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "_")
    set(Fortran_INTROSPECTION_MANGLING "${fc_case}Underscore")
  elseif(FortranCInterface_GLOBAL__SUFFIX STREQUAL "__")
    set(Fortran_INTROSPECTION_MANGLING "${fc_case}DoubleUnderscore")
  else()
    message( WARNING "Fortran name mangling suffix, \'${FortranCInterface_GLOBAL__SUFFIX}\', unknown to FortranIntrospection CMake module.")
  endif()
  message(STATUS "Fortran_INTROSPECTION_MANGLING = ${Fortran_INTROSPECTION_MANGLING}")
endmacro()

include(CheckFortranSourceCompiles)
#----------------------------------------------
# Test to ensure compiler works with submodules
#----------------------------------------------
# returns Fortran_INTROSPECTION_HAVE_SUBMODULES variable
macro(fortran_introspection_have_submodules)
  CHECK_Fortran_SOURCE_COMPILES("
    module foo_module
      implicit none
      interface
        module subroutine bar
          implicit none
        end subroutine
      end interface
    end module
    submodule(foo_module) foo_submodule
      implicit none
    contains
      module procedure bar
        print *,'Hello, world!'
      end procedure
    end submodule
    program main
      use foo_module, only : bar
      implicit none
      call bar
    end program
"   Fortran_INTROSPECTION_HAVE_SUBMODULES
    SRC_EXT ".f90")
  if(${Fortran_INTROSPECTION_HAVE_SUBMODULES})
    message(STATUS "Fortran compiler supports submodules")
  endif()
endmacro()

#-------------------------------------------------------------
# Test to ensure compiler allows error stop in pure procedures
#-------------------------------------------------------------
# returns Fortran_INTROSPECTION_HAVE_ERROR_STOP_IN_PURE variable
macro(fortran_introspection_have_error_stop_in_pure)
  check_fortran_source_compiles("
    program main
    contains
    pure function foo() result(res)
      error stop 'Error stop is supported in pure functions (F2018)'
    end function
    end program
"
    Fortran_INTROSPECTION_HAVE_ERROR_STOP_IN_PURE
    SRC_EXT ".f90"
    )
  if(Fortran_INTROSPECTION_HAVE_ERROR_STOP_IN_PURE)
    message(STATUS "Fortran compiler supports error stop in pure procedures")
  endif()
endmacro()

#-----------------------------------------------------
# Test to ensure error stop allows non-constant values
#-----------------------------------------------------
# returns Fortran_INTROSPECTION_HAVE_NONCONSTANT_ERROR_STOP
macro(fortran_introspection_have_nonconstant_error_stop)
  check_fortran_source_compiles("
    program main
    integer :: i
    i = 0
    error stop i
    end program
"
    Fortran_INTROSPECTION_HAVE_NONCONSTANT_ERROR_STOP
    SRC_EXT ".f90")
  if(Fortran_INTROSPECTION_HAVE_NON_CONSTANT_ERROR_STOP)
    message(STATUS "Fortran compiler supports non constant error stop")
  endif()
endmacro()

#--------------------------------------------
# Test that a simple coarray program compiles
#--------------------------------------------
# returns Fortran_INTROSPECTION_HAVE_COARRAYS
macro(fortran_introspection_have_coarrays)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(coarray_flags "-fcoarray=single")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(CMAKE_SYSTEM_NAME MATCHES "Windows")
      set(coarray_flags "/Qcoarray:shared /standard-semantics")
    else()
      set(coarray_flags "-coarray=shared")
    endif()
  else()
    message(WARNING "Coarray flags unknown for current Fortran compiler (${CMAKE_Fortran_COMPILER_ID})")
  endif()
  set(OLD_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
  set(CMAKE_REQUIRED_FLAGS "${coarray_flags}")
  check_fortran_source_compiles("
    program main
      implicit none
      integer :: i[*]
      i = this_image()
    end program
" Fortran_INTROSPECTION_HAVE_COARRAYS
    SRC_EXT ".f90")
  if(Fortran_INTROSPECTION_HAVE_COARRAYS)
    message(STATUS "Fortran compiler supports coarrays")
  endif()
  set (CMAKE_REQUIRED_FLAGS ${OLD_REQUIRED_FLAGS})
  unset(OLD_REQUIRED_FLAGS)
endmacro()
