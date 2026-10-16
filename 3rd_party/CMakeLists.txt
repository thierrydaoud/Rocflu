# METIS 5 and GKlib, built from the sources in this directory for rflupart.
#
# The defines and type widths are those of METIS's own build (make config):
# 32-bit idx_t and real_t, -O3, asserts off. Two of its flags are left out:
# -march=native (a library built on one node type could fail on another) and
# -Werror (a newer compiler's warnings must not break the build).

set(GKLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/GKlib)
set(METIS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/METIS)

include(CheckIncludeFile)
include(CheckFunctionExists)
check_include_file(execinfo.h HAVE_EXECINFO_H)
check_function_exists(getline HAVE_GETLINE)

set(GK_DEFINES LINUX _FILE_OFFSET_BITS=64 NDEBUG NDEBUG2)
if(HAVE_EXECINFO_H)
  list(APPEND GK_DEFINES HAVE_EXECINFO_H)
endif()
if(HAVE_GETLINE)
  list(APPEND GK_DEFINES HAVE_GETLINE)
endif()
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(GK_FLAGS -O3 -std=c99 -fno-strict-aliasing)
else()
  set(GK_FLAGS -O3)
endif()

# --- GKlib ---------------------------------------------------------------------
file(GLOB GKLIB_SOURCES ${GKLIB_DIR}/*.c)
add_library(GKlib STATIC ${GKLIB_SOURCES})
target_include_directories(GKlib PRIVATE ${GKLIB_DIR})
target_compile_definitions(GKlib PRIVATE ${GK_DEFINES})
target_compile_options(GKlib PRIVATE ${GK_FLAGS})
target_link_libraries(GKlib PUBLIC m)

# --- METIS ---------------------------------------------------------------------
# metis.h with the type widths in front, as METIS's make config writes it.
set(METIS_GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)
file(READ ${METIS_DIR}/include/metis.h metis_header)
file(WRITE ${METIS_GEN_DIR}/metis.h.tmp
     "#define IDXTYPEWIDTH 32\n#define REALTYPEWIDTH 32\n${metis_header}")
file(COPY_FILE ${METIS_GEN_DIR}/metis.h.tmp ${METIS_GEN_DIR}/metis.h ONLY_IF_DIFFERENT)
file(REMOVE ${METIS_GEN_DIR}/metis.h.tmp)

file(GLOB METIS_SOURCES ${METIS_DIR}/libmetis/*.c)
add_library(metis_bundled STATIC ${METIS_SOURCES})
set_target_properties(metis_bundled PROPERTIES OUTPUT_NAME metis)
target_include_directories(metis_bundled PRIVATE ${METIS_DIR}/libmetis ${GKLIB_DIR})
target_compile_definitions(metis_bundled PRIVATE ${GK_DEFINES})
target_include_directories(metis_bundled PUBLIC ${METIS_GEN_DIR})
target_compile_options(metis_bundled PRIVATE ${GK_FLAGS})
target_link_libraries(metis_bundled PUBLIC GKlib)
