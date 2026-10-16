# Script mode: builds ppiclF's include file PPICLF (the common blocks).
# Runs the C preprocessor over HEADERS, in order, and keeps every line that
# does not start with '#' (empty lines are dropped too), like
#   mpicc -E <headers> | grep "^[^#]"
# Inputs: CPP (C compiler), INCLUDES (dirs), HEADERS, OUTPUT; lists are
# separated by "|".

string(REPLACE "|" ";" INCLUDES "${INCLUDES}")
string(REPLACE "|" ";" HEADERS "${HEADERS}")

set(include_flags "")
foreach(dir IN LISTS INCLUDES)
  list(APPEND include_flags "-I${dir}")
endforeach()

set(text "")
foreach(header IN LISTS HEADERS)
  execute_process(
    COMMAND ${CPP} -E ${include_flags} ${header}
    OUTPUT_VARIABLE expanded
    RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "Preprocessing ${header} failed")
  endif()
  string(APPEND text "${expanded}")
endforeach()

# Drop the lines that start with '#' and the empty lines.
string(PREPEND text "\n")
string(REGEX REPLACE "\n#[^\n]*" "" text "${text}")
string(REGEX REPLACE "\n\n+" "\n" text "${text}")
string(REGEX REPLACE "^\n" "" text "${text}")

file(WRITE ${OUTPUT} "${text}")
