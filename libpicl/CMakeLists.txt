# ppiclF and gslib (replaces libpicl/Makefile and gslib's install script).

set(PPICLF_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(PPICLF_USER_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/user_files)
set(GSLIB_DIR         ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/gslib)

# --- gslib -------------------------------------------------------------------
# The objects gslib's own Makefile puts in libgs.a, with the options ppiclF's
# install script passes (MPI, trailing underscore, pgslib_/pfgslib_ prefixes,
# user exit, no BLAS choice).
set(GSLIB_SOURCES gs sort sarray_transfer sarray_sort gs_local fail crystal
                  comm tensor fcrystal findpts findpts_local obbox poly lob_bnd
                  findpts_el_3 findpts_el_2)
list(TRANSFORM GSLIB_SOURCES PREPEND ${GSLIB_DIR}/gslib/src/)
list(TRANSFORM GSLIB_SOURCES APPEND .c)

add_library(gs STATIC ${GSLIB_SOURCES})
target_include_directories(gs PRIVATE ${GSLIB_DIR}/gslib/src)
target_compile_definitions(gs PRIVATE MPI UNDERSCORE GLOBAL_LONG_LONG
                           PREFIX=pgslib_ FPREFIX=pfgslib_ USE_USR_EXIT)
target_compile_options(gs PRIVATE ${GSLIB_C_FLAGS})
target_link_libraries(gs PUBLIC MPI::MPI_C)

# --- ppiclF ------------------------------------------------------------------
# Compiled as one file, as the Makefile did: the core and user files are
# concatenated behind the two headers that define the layout macros.
set(PPICLF_SOURCES
    ${PPICLF_SOURCE_DIR}/ppiclf_unsteady.f
    ${PPICLF_SOURCE_DIR}/ppiclf_data.f
    ${PPICLF_USER_DIR}/ppiclf_user.f
    ${PPICLF_USER_DIR}/ppiclf_user_subbinMap.f
    ${PPICLF_USER_DIR}/ppiclf_user_QuasiSteady.f
    ${PPICLF_USER_DIR}/ppiclf_user_Fluctuations.f
    ${PPICLF_USER_DIR}/ppiclf_user_AddedMass.f
    ${PPICLF_USER_DIR}/Added_Mass/ppiclf_user_AddedMass_functions.f
    ${PPICLF_USER_DIR}/Added_Mass/ppiclf_user_AddedMass_subroutines.f
    ${PPICLF_USER_DIR}/ppiclf_user_HeatTransfer.f
    ${PPICLF_USER_DIR}/ppiclf_user_Torque.f
    ${PPICLF_USER_DIR}/ppiclf_user_Lift.f
    ${PPICLF_USER_DIR}/ppiclf_user_ViscousUnsteady.f
    ${PPICLF_USER_DIR}/ppiclf_user_EvalNearestNeighbor.f
    ${PPICLF_USER_DIR}/ppiclf_user_OutputFields.f
    ${PPICLF_USER_DIR}/ppiclf_user_Layout.f
    ${PPICLF_SOURCE_DIR}/ppiclf_comm.f
    ${PPICLF_SOURCE_DIR}/ppiclf_op.f
    ${PPICLF_SOURCE_DIR}/ppiclf_io.f
    ${PPICLF_SOURCE_DIR}/ppiclf_solve.f
    ${PPICLF_SOURCE_DIR}/ppiclf_wall.f
    ${PPICLF_SOURCE_DIR}/ppiclf_mxm.f)

# The user's headers come first so they win over any copy in source/.
set(PPICLF_HEADERS
    ${PPICLF_SOURCE_DIR}/PPICLF_SOLN.h
    ${PPICLF_SOURCE_DIR}/PPICLF_GRID.h
    ${PPICLF_SOURCE_DIR}/PPICLF_OPT.h
    ${PPICLF_SOURCE_DIR}/PPICLF_PARALLEL.h
    ${PPICLF_USER_DIR}/PPICLF_USER_COMMON.h
    ${PPICLF_SOURCE_DIR}/PPICLF_GEOM.h)

set(PPICLF_GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)

# The helper scripts take their file lists '|'-separated: a ';' would split
# the command line.
string(REPLACE ";" "|" headers_arg "${PPICLF_HEADERS}")
string(REPLACE ";" "|" sources_arg "${PPICLF_SOURCES}")
file(MAKE_DIRECTORY ${PPICLF_GEN_DIR})

# PPICLF: the Fortran include file with all common blocks, made by running
# the C preprocessor over the headers and dropping its # lines.
add_custom_command(
  OUTPUT ${PPICLF_GEN_DIR}/PPICLF
  COMMAND ${CMAKE_COMMAND} -DCPP=${CMAKE_C_COMPILER}
          -DOUTPUT=${PPICLF_GEN_DIR}/PPICLF
          "-DINCLUDES=${PPICLF_USER_DIR}|${PPICLF_SOURCE_DIR}|${PPICLF_LAYOUT_DIR}"
          "-DHEADERS=${headers_arg}"
          -P ${PROJECT_SOURCE_DIR}/cmake/PpiclfCommonBlocks.cmake
  DEPENDS ${PPICLF_HEADERS} ${PROJECT_SOURCE_DIR}/cmake/PpiclfCommonBlocks.cmake
          ${PPICLF_USER_DIR}/PPICLF_USER.h ${PPICLF_SOURCE_DIR}/PPICLF_STD.h
          ${PPICLF_LAYOUT_DIR}/PPICLF_LAYOUT.h
  COMMENT "Generating ppiclF common-block include file PPICLF"
  VERBATIM)

add_custom_command(
  OUTPUT ${PPICLF_GEN_DIR}/ppiclf.f
  COMMAND ${CMAKE_COMMAND} -DOUTPUT=${PPICLF_GEN_DIR}/ppiclf.f
          "-DSOURCES=${sources_arg}"
          -P ${PROJECT_SOURCE_DIR}/cmake/PpiclfConcatenate.cmake
  DEPENDS ${PPICLF_SOURCES} ${PROJECT_SOURCE_DIR}/cmake/PpiclfConcatenate.cmake
  COMMENT "Concatenating ppiclF sources into ppiclf.f"
  VERBATIM)

add_library(ppiclF STATIC ${PPICLF_GEN_DIR}/ppiclf.f ${PPICLF_GEN_DIR}/PPICLF)
set_source_files_properties(${PPICLF_GEN_DIR}/PPICLF PROPERTIES HEADER_FILE_ONLY TRUE)
set_target_properties(ppiclF PROPERTIES
  Fortran_FORMAT FIXED
  Fortran_PREPROCESS ON
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/ppiclF)
target_include_directories(ppiclF PRIVATE
  ${PPICLF_USER_DIR} ${PPICLF_SOURCE_DIR} ${PPICLF_GEN_DIR} ${PPICLF_LAYOUT_DIR}
  ${GSLIB_DIR}/include)
target_compile_options(ppiclF PRIVATE ${PPICLF_Fortran_FLAGS})
target_link_libraries(ppiclF PUBLIC gs MPI::MPI_Fortran)
