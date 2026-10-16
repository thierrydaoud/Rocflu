cmake_minimum_required(VERSION 3.20)

project(Rocflu LANGUAGES C CXX Fortran)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(RocfluFilelist)

# --- Options -----------------------------------------------------------------
option(ROCFLU_PICL      "Couple the ppiclF particle solver"                ON)
option(ROCFLU_SPEC      "Species transport (also needed by the JWL gas model, GASMODEL 7)" OFF)
option(ROCFLU_PLAG      "Rocpart Lagrangian particles"                     OFF)
option(ROCFLU_FOLDER    "Write output into per-kind folders"               ON)
option(ROCFLU_UTILITIES "Build the pre-processing utilities (map, part, init)" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release or Debug" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug)

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(WARNING "Only gfortran has been tested; the flags below are gfortran flags")
endif()

# --- Flags -------------------------------------------------------------------
# The flags are set per part of the code, so CMake's own per-build-type
# defaults (-O3 -DNDEBUG, -g) are switched off.
foreach(lang C CXX Fortran)
  foreach(config RELEASE DEBUG RELWITHDEBINFO MINSIZEREL)
    set(CMAKE_${lang}_FLAGS_${config} "")
  endforeach()
endforeach()

set(DEBUG_FLAGS -g -O0 -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(ROCFLU_Fortran_FLAGS ${DEBUG_FLAGS} -fallow-argument-mismatch)
  set(PPICLF_Fortran_FLAGS -w -std=legacy -fmax-errors=10 -fallow-argument-mismatch ${DEBUG_FLAGS})
  set(ROCFLU_CXX_FLAGS -g)
  set(GSLIB_C_FLAGS -g -O0)
else()
  set(ROCFLU_Fortran_FLAGS -O -fallow-argument-mismatch)
  # -fcheck=all costs about 2% and has caught real bugs, so it stays on.
  set(PPICLF_Fortran_FLAGS -w -O3 -std=legacy -g -fmax-errors=10
                           -fallow-argument-mismatch -fcheck=all -fbacktrace)
  set(ROCFLU_CXX_FLAGS "")
  set(GSLIB_C_FLAGS -O2)
endif()

# Preprocessor flags (the Filelists test the same names with ifdef).
set(ROCFLU_DEF_RFLU TRUE)
set(ROCFLU_DEF_PICL ${ROCFLU_PICL})
set(ROCFLU_DEF_SPEC ${ROCFLU_SPEC})
set(ROCFLU_DEF_PLAG ${ROCFLU_PLAG})
set(ROCFLU_DEF_FOLDER ${ROCFLU_FOLDER})
if(ROCFLU_PICL OR ROCFLU_SPEC OR ROCFLU_PLAG)
  set(ROCFLU_DEF_INRT TRUE)
endif()
set(ROCFLU_DEFINES "")
foreach(flag RFLU PICL PLAG SPEC INRT FOLDER)
  if(ROCFLU_DEF_${flag})
    list(APPEND ROCFLU_DEFINES ${flag})
  endif()
endforeach()

# --- Dependencies ------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS C CXX Fortran)

# METIS 5 (rflupart only): built from 3rd_party/ unless METIS_ROOT points at
# an existing install with lib/libmetis.a, lib/libGKlib.a and include/metis.h.
set(METIS_ROOT "" CACHE PATH "Existing METIS/GKlib install (empty: build the bundled copy)")
add_library(metis INTERFACE)
if(ROCFLU_UTILITIES)
  if(METIS_ROOT)
    find_library(METIS_LIBRARY metis PATHS ${METIS_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
    find_library(GKLIB_LIBRARY GKlib PATHS ${METIS_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
    target_link_libraries(metis INTERFACE ${METIS_LIBRARY} ${GKLIB_LIBRARY} m)
    target_include_directories(metis INTERFACE ${METIS_ROOT}/include)
  else()
    add_subdirectory(3rd_party)
    target_link_libraries(metis INTERFACE metis_bundled)
  endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Installed executables keep the library paths they were linked with (MPI)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# Applies the Rocflu flags to a target with Fortran and C++ sources.
function(rocflu_flags target)
  target_compile_definitions(${target} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${ROCFLU_DEFINES}>
    $<$<COMPILE_LANGUAGE:CXX>:FORTRAN_UNDERSCORE FORTRAN_LOWERCASE>)
  target_compile_options(${target} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${ROCFLU_Fortran_FLAGS}>
    $<$<COMPILE_LANGUAGE:CXX>:${ROCFLU_CXX_FLAGS}>)
  target_include_directories(${target} PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${PROJECT_SOURCE_DIR}/include>)
  set_target_properties(${target} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/${target})
endfunction()

# --- ppiclF ------------------------------------------------------------------
if(ROCFLU_PICL)
  add_subdirectory(libpicl)
endif()

# --- libflu: everything shared by the solver and the utilities --------------
set(FLU_DIRS modfloflu modflu libfloflu libflu rocflu)
if(ROCFLU_DEF_INRT)
  list(APPEND FLU_DIRS rocinteract)
endif()
if(ROCFLU_PLAG)
  list(APPEND FLU_DIRS rocpart)
endif()
if(ROCFLU_PICL)
  list(APPEND FLU_DIRS rocpicl)
endif()
if(ROCFLU_SPEC)
  list(APPEND FLU_DIRS rocspecies)
endif()

set(FLU_SOURCES "")
foreach(dir IN LISTS FLU_DIRS)
  rocflu_read_filelist(${PROJECT_SOURCE_DIR}/${dir} dir_sources)
  list(APPEND FLU_SOURCES ${dir_sources})
endforeach()

add_library(flu STATIC ${FLU_SOURCES})
rocflu_flags(flu)
target_include_directories(flu INTERFACE ${CMAKE_BINARY_DIR}/modules/flu)
# The ppiclF layout macros (PPICLF_USER.h, PPICLF_STD.h); the user's copy first.
target_include_directories(flu PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:${PROJECT_SOURCE_DIR}/libpicl/user_files>
  $<$<COMPILE_LANGUAGE:Fortran>:${PROJECT_SOURCE_DIR}/libpicl/source>)
target_link_libraries(flu PUBLIC MPI::MPI_Fortran)

# Links like the Makefiles: C++ linker, own objects, then libflu.
function(rocflu_executable name dir)
  rocflu_read_filelist(${PROJECT_SOURCE_DIR}/${dir} sources)
  add_executable(${name} ${sources})
  rocflu_flags(${name})
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(${name} PRIVATE flu ${ARGN} MPI::MPI_Fortran MPI::MPI_CXX)
  # OPTIONAL: installing after building only some targets copies those
  install(TARGETS ${name} RUNTIME DESTINATION bin OPTIONAL)
endfunction()

# --- Solver ------------------------------------------------------------------
if(ROCFLU_PICL)
  rocflu_executable(rflump simulation ppiclF)
else()
  rocflu_executable(rflump simulation)
endif()

# --- Pre-processing utilities ---------------------------------------------------------------
if(ROCFLU_UTILITIES)
  rocflu_executable(rflumap utilities/map)
  rocflu_executable(rflupart utilities/part metis)
  rocflu_executable(rfluinit utilities/init)
endif()

message(STATUS "Rocflu ${CMAKE_BUILD_TYPE}: defines ${ROCFLU_DEFINES}")
