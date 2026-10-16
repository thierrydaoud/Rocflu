# Script mode: writes ppiclf.f, the two layout headers followed by every
# ppiclF source in SOURCES (order matters: modules first).
# Inputs: SOURCES ("|"-separated), OUTPUT.

string(REPLACE "|" ";" SOURCES "${SOURCES}")

set(text "#include \"PPICLF_USER.h\"\n#include \"PPICLF_STD.h\"\n")
foreach(source IN LISTS SOURCES)
  file(READ ${source} contents)
  string(APPEND text "${contents}")
endforeach()

file(WRITE ${OUTPUT} "${text}")
