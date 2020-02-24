# Add portable unistall command to makefile, adapted from the CMake Wiki FAQ
configure_file("${CMAKE_CURRENT_LIST_DIR}/../uninstall.cmake.in" "${CMAKE_BINARY_DIR}/uninstall.cmake" @ONLY)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/uninstall.cmake")
