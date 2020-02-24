include_guard(GLOBAL)

# Set a default build type if none was specified
function(set_config_types)
  set(options STRICT)
  set(oneValueArgs DEFAULT_DEV_CONFIG DEFAULT_RELEASE_CONFIG)
  set(multiValueArgs CONFIGS)
  cmake_parse_arguments(my_conf "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(NOT my_conf_STRICT) # add default build types
    set(build_types "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif()

  # Handle strict case, CONFIGS required!
  if(my_conf_STRICT AND NOT my_conf_CONFIGS)
    message(FATAL_ERROR "List of valid configurations must be passed to `set_config_type()` when `STRICT` is set.")
  elseif(my_conf_CONFIGS)
    # If CONFIGS is present add them to possible builds list
    set(build_types ${build_types} ${my_conf_CONFIGS})
    list(REMOVE_DUPLICATES my_conf_CONFIGS)
  endif()

  # Default to release builds if not specified and not developer build
  if(NOT my_conf_DEFAULT_RELEASE_CONFIG)
    set(my_conf_DEFAULT_RELEASE_CONFIG "Release")
  endif()
  # Default to RelWithDebInfo for developer builds
  if(NOT my_conf_DEFAULT_DEV_CONFIG)
    set(my_conf_DEFAULT_DEV_CONFIG "Debug")
  endif()

  # Determine if default should be release or developer build
  set(default_build_type ${my_conf_DEFAULT_RELEASE_CONFIG})
  if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(default_build_type ${my_conf_DEFAULT_DEV_CONFIG})
    if(DEFINED ENV{CMAKE_CONFIG_TYPE} AND "$ENV{CI}")
      set(default_build_type $ENV{CMAKE_CONFIG_TYPE})
    endif()
  endif()

  # Check if we're using an IDE
  get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  set(IS_MULTI_CONFIG
      ${isMultiConfig}
      PARENT_SCOPE
      )

  # If use did not specify build type and not an IDE:
  if(NOT CMAKE_BUILD_TYPE AND NOT IS_MULTI_CONFIG)
    message(
      STATUS
        "Setting build type to '${default_build_type}', use `-DCMAKE_BUILD_TYPE=<type>` or set with ccmake/cmake-gui to override."
      )
    set(CMAKE_BUILD_TYPE
        "${default_build_type}"
        CACHE STRING "Choose the type of build." FORCE
        )
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${build_types})
  endif()

  if(my_conf_STRICT)
    if(NOT CMAKE_BUILD_TYPE IN_LIST build_types)
      message(FATAL_ERROR "Invalid build type: ${CMAKE_BUILD_TYPE}, strict enforcement is ON!")
    endif()
  endif()
endfunction()
