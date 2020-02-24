# Determine properties of the host build machine and parallel runtime being used
include_guard(DIRECTORY)

#------------------------------------------------------------
# Determine if we're using Open MPI possibly via CAF wrappers
#------------------------------------------------------------
macro(detect_parallel_runtime)
  cmake_host_system_information(RESULT N_CPU QUERY NUMBER_OF_PHYSICAL_CORES)
  cmake_host_system_information(RESULT HOST_NAME QUERY HOSTNAME)
  if($ENV{FC} MATCHES caf)
    message(STATUS "User passed caf wrapper as $FC")
    get_filename_component(_CAF "$ENV{FC}"
      PROGRAM PROGRAM_ARGS _CAF_ARGS)
    message(STATUS "CAF = ${_CAF}")
    get_filename_component(_CAF_BINDIR "${_CAF}" DIRECTORY)
    message(STATUS "_CAF_BINDIR = ${_CAF_BINDIR}")
    find_program( CAF_RUN
      cafrun
      HINTS "${_CAF_BINDIR}")
    message(STATUS "CAF_RUN = ${CAF_RUN}")
    if(EXISTS "${CAF_RUN}")
      execute_process(COMMAND "${CAF_RUN}" --show
        OUTPUT_VARIABLE CAFRUN_COMMAND)
      message(STATUS "CAFRUN_COMMAND = ${CAFRUN_COMMAND}")
      string(REPLACE " " ";" CAFRUN_COMMAND_LIST ${CAFRUN_COMMAND})
      list(GET CAFRUN_COMMAND_LIST 0 CAF_MPIEXEC)
      message(STATUS "CAF_MPIEXEC = ${CAF_MPIEXEC}")
    endif()
  endif()
  if(EXISTS "${CAF_MPIEXEC}")
    set(MPIEXEC "${CAF_MPIEXEC}")
  else()
    find_package(MPI)
  endif()
  execute_process(COMMAND ${MPIEXEC} --version
    OUTPUT_VARIABLE mpi_version_out)
  #  message(STATUS "MPI Version info: ${mpi_version_out}")
  if (mpi_version_out MATCHES "[Oo]pen[ -][Mm][Pp][Ii]")
    message( STATUS "OpenMPI detected:")
    set ( OPENMPI true)
    # Write out a host file because OMPI's mpiexec is dumb
    set(MPI_HOSTFILE_PATH "${CMAKE_BINARY_DIR}/hostfile")
    file(WRITE "${MPI_HOSTFILE_PATH}" "${HOST_NAME} slots=${N_CPU}\n")
  endif ()
endmacro()
