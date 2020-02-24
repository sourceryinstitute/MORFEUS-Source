include_guard(DIRECTORY)

function(get_vcs_info)
  # Sets the following variables in the parent scope: GIT_PROJECT_ROOT GIT_SUBPROJECT_ROOT GIT_VERSION: output of git dscribe
  # --abbreve=0 validated to be a tag version string GIT_FULL_DESCRIBE: output of git describe --always GIT_LAST_TAG: Same as
  # GIT_VERSION but without validation COMMIT_DATE_UNIX_TIME: Commit date in unix time suitable for setting SOURCE_DATE_EPOCH
  # GIT_COMMIT_HASH: The SHA1 hash of the HEAD (last) commit

  if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(GIT_PROJECT_ROOT
        ${CMAKE_SOURCE_DIR}
        PARENT_SCOPE
        )
    set(WORKING_GIT_ROOT ${CMAKE_SOURCE_DIR})
  endif()
  if((NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR) AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    set(GIT_SUBPROJECT_ROOT
        ${CMAKE_CURRENT_SOURCE_DIR}
        PARENT_SCOPE
        )
    set(WORKING_GIT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
  endif()

  if(WORKING_GIT_ROOT)
    find_package(Git QUIET)

    if(GIT_FOUND)
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --abbrev=0
        WORKING_DIRECTORY "${WORKING_GIT_ROOT}"
        RESULT_VARIABLE git_status
        OUTPUT_VARIABLE git_output
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
        )
      if((git_status STREQUAL "0") AND (git_output MATCHES "[vV]?[0-9]+\\.[0-9]+\\.[0-9]+(-rc[0-9]+)?"))
        set(GIT_VERSION "${git_output}")
      else()
        set(GIT_VERSION NOTFOUND)
      endif()
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --always
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE git_status
        OUTPUT_VARIABLE GIT_FULL_DESCRIBE
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
        )
      if(NOT (git_status STREQUAL "0"))
        set(GIT_FULL_DESCRIBE NOTFOUND)
      endif()
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --abbrev=0
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE git_status
        OUTPUT_VARIABLE GIT_LAST_TAG
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
        )
      if(NOT (git_status STREQUAL "0"))
        set(GIT_LAST_TAG NOTFOUND)
      endif()
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" show -s --format=%ct HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE git_status
        OUTPUT_VARIABLE COMMIT_DATE_UNIX_TIME
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
        )
      if(NOT (git_status STREQUAL "0"))
        set(COMMIT_DATE_UNIX_TIME NOTFOUND)
      endif()
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse --verify HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE git_status
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
        )
      if(NOT (git_status STREQUAL "0"))
        set(GIT_COMMIT_HASH NOTFOUND)
      endif()
      set(GIT_VERSION
          ${GIT_VERSION}
          PARENT_SCOPE
          )
      set(GIT_FULL_DESCRIBE
          ${GIT_FULL_DESCRIBE}
          PARENT_SCOPE
          )
      set(GIT_LAST_TAG
          ${GIT_LAST_TAG}
          PARENT_SCOPE
          )
      set(COMMIT_DATE_UNIX_TIME
          ${COMMIT_DATE_UNIX_TIME}
          PARENT_SCOPE
          )
      set(GIT_COMMIT_HASH
          ${GIT_COMMIT_HASH}
          PARENT_SCOPE
          )
    else()
      message(WARNING "Could not find git executable!")
    endif()
  endif()
endfunction()
