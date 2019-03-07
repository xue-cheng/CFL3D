cmake_minimum_required(VERSION 3.13)
project(CFL3D
  VERSION 6.7
  LANGUAGES Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

###########################################################################
# Options
option(SINGLE_PRECISION  "create single precision executables." OFF)
option(NO_REDIRECT "disallow redirected input file" OFF)
option(WITH_CGNS "build with CGNS support" OFF)
option(WITH_FASTIO "build with fast I/O option for MPI versions" ON)

option(BUILD_REAL "build real versions of CFL3D " ON)
option(BUILD_CMPLX "build complex versions of CFL3D " OFF)

option(BUILD_SEQ "build sequential version of CFL3D" ON)
option(BUILD_MPI "build MPI versions of CFL3D " ON)
option(BUILD_SPLITTER "build splitter" ON)
option(BUILD_PRECFL3D "build precfl3d" ON)
option(BUILD_RONNIE "build ronnie" OFF)
option(BUILD_PRERONNIE "build preronnie" OFF)
option(BUILD_MAGGIE "build maggie" OFF)
option(BUILD_TOOLS "build cfl3d tools" ON)

###########################################################################
# setting up compiler flags

if(NO_REDIRECT)
  add_compile_definitions(NOREDIRECT) 
endif(NO_REDIRECT)

# CGNS
if(WITH_CGNS)
  message(WARNING "Only CGNS 2 is supported, download URL:")
  message(
    WARNING 
    "https://github.com/CGNS/CGNS_archives/raw/master/cgnslib_2.5-5.tar.gz"
  )

  find_package(CGNS 2.4 REQUIRED)
	add_library(cgns STATIC IMPORTED)
  set_property(TARGET cgns PROPERTY
    IMPORTED_LOCATION
    ${CGNS_LIBRARY}
  )
  set_property(TARGET cgns PROPERTY
    IMPORTED_LINK_INTERFACE_LANGUAGES
    Fortran
  )
  set_property(TARGET cgns PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES 
    $<BUILD_INTERFACE:${CGNS_INCLUDE_DIR}>
  )
  set_property(TARGET cgns PROPERTY
    INTERFACE_COMPILE_DEFINITIONS
    CGNS
  )
endif(WITH_CGNS)

# MPI
if(BUILD_MPI)
  if(WIN32)
    message(FATAL_ERROR "Windows is not supported, CFL3D can be compiled in WSL.")
  else(WIN32)
    find_package(MPI 1.0 REQUIRED)
    if(MPI_Fortran_FOUND)
      add_library(mpi STATIC IMPORTED)
      set_property(TARGET mpi PROPERTY 
        INTERFACE_INCLUDE_DIRECTORIES ${MPI_Fortran_INCLUDE_PATH})
      #message("Choose the first element of libs: ${MPI_Fortran_LIBRARIES}")
      set(other_elems ${MPI_Fortran_LIBRARIES})
      list(GET other_elems 0 first_elem)
      list(REMOVE_AT other_elems 0)
      set_property(TARGET mpi PROPERTY 
        IMPORTED_LOCATION ${first_elem})
      set_property(TARGET mpi PROPERTY 
        INTERFACE_LINK_LIBRARIES ${other_elems})
      set_property(TARGET mpi APPEND PROPERTY 
        INTERFACE_COMPILE_DEFINITIONS "DIST_MPI")
      set_property(TARGET mpi APPEND PROPERTY 
        LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
      set_property(TARGET mpi APPEND PROPERTY 
        INTERFACE_COMPILE_OPTIONS "${MPI_Fortran_COMPILE_FLAGS}")
      ## FAST-io
      if(WITH_FASTIO)
        set_property(TARGET mpi APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS "FASTIO")
      endif(WITH_FASTIO)
  else(MPI_Fortran_FOUND)
    ## FAST-io
    if(WITH_FASTIO)
      # TURN OFF FASTIO
      set(WITH_FASTIO OFF CACHE BOOL "build with fast I/O option for MPI versions (N/A)" FORCE) 
    endif(WITH_FASTIO)
    set(BUILD_MPI OFF CACHE BOOL "build MPI versions of CFL3D" FORCE) # set the MPI N/A
  endif(MPI_Fortran_FOUND)
  endif(WIN32)
  
endif(BUILD_MPI)

# in LINUX system
if(UNIX)
	message("in Linux or Linux-like system")
	add_compile_definitions(LINUX)
endif(UNIX)


# Single Precision
if(NOT SINGLE_PRECISION)
  add_compile_definitions(DBLE_PRECSN)
endif(NOT SINGLE_PRECISION)


if(MINGW)
  add_link_options("LINKER:-Wl,--allow-multiple-definition")
else(MINGW)
	add_link_options("LINKER:-z,muldefs")
endif(MINGW)

## project compiler dependent settings
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  add_compile_definitions(INTEL)
  add_compile_options(-xHost -traceback -fpe0 -ip -fno-alias)
  # Single Precision
  if(NOT SINGLE_PRECISION)
    add_compile_options(-r8)
  endif(NOT SINGLE_PRECISION)
	# no `-w` option, do not turn warnings off.
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  add_compile_options( -march=native -Wno-align-commons -fbacktrace)
  # Single Precision
  if(NOT SINGLE_PRECISION)
    add_compile_options(-fdefault-real-8)
  endif(NOT SINGLE_PRECISION)
else()
	message(FATAL_ERROR "Unknown Compiler ID: " ${CMAKE_Fortran_COMPILER_ID})
endif()

# it seems always use P3D_SINGLE in original Install script
add_compile_definitions(P3D_SINGLE)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include CACHE PATH 
  "Where to put all the Fortran modules when built.")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib CACHE PATH
  "Where to put all the LIBRARY target files when built.")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin CACHE PATH 
  "Where to put all the RUNTIME target files when built.")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib CACHE PATH
  "Where to put all the ARCHIVE target files when built.") 

add_subdirectory(source) 