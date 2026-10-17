#!/bin/bash
# rocflu.sh: build and test Rocflu + ppiclF (a small front end to CMake).
#
#   ./rocflu.sh build [options]     configure and compile (every particle slot)
#   ./rocflu.sh run <case>.inp ...  check the case, build for it, pre-process it
#                                   here and submit the simulation to Slurm
#   ./rocflu.sh validate <case>.inp check a case's input and files only
#   ./rocflu.sh test  [options]     run the regression suite (tests/run_tests.py)
#   ./rocflu.sh clean               remove build/
#   source ./rocflu.sh load         load the HiPerGator modules into your shell
#
# Run ./rocflu.sh <command> --help for the options of a command.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_ROOT="$ROOT/build"
# Modules loaded on HiPerGator; override with ROCFLU_MODULES="..." (empty: none).
DEFAULT_MODULES="gcc/14.2.0 openmpi/5.0.7"
ROCFLU_MODULES="${ROCFLU_MODULES-$DEFAULT_MODULES}"

say()  { echo "${PROGRESS:-}rocflu:${PLAIN:-} $*"; }
fail() { echo "${ERROR_MARK:-}rocflu: ERROR:${ERROR_PLAIN:-} $*" >&2; exit 1; }

# Loads the modules if the module command exists (it does on HiPerGator).
load_modules() {
  if [ -n "$ROCFLU_MODULES" ] && command -v module > /dev/null 2>&1; then
    module purge > /dev/null 2>&1
    # shellcheck disable=SC2086
    module load $ROCFLU_MODULES || return 1
    say "loaded modules: $ROCFLU_MODULES"
  fi
}

# 'load' changes the caller's shell, so it only works when sourced.
if [ "$1" == "load" ]; then
  if [ "${BASH_SOURCE[0]}" == "$0" ]; then
    echo "rocflu: 'load' must be sourced: source ./rocflu.sh load" >&2
    exit 1
  fi
  load_modules
  return
fi

# Colours only on a terminal, and not when NO_COLOR is set (no-color.org):
# logs, pipes and srun output stay plain. Set after 'load', so that sourcing
# the script does not leave them in the shell.
colours_on() { [ -t "$1" ] && [ -z "${NO_COLOR:-}" ] && [ "${TERM:-}" != dumb ]; }
if colours_on 1; then
  HEADING=$'\033[1m'; OPTION=$'\033[36m'; DEFAULT=$'\033[33m'
  PROGRESS=$'\033[32m'; DIM=$'\033[2m'; COMPLETED=$'\033[1;36m'; PLAIN=$'\033[0m'
fi
if colours_on 2; then
  ERROR_MARK=$'\033[1;31m'; DEFAULT_ERR=$'\033[33m'; ERROR_PLAIN=$'\033[0m'
fi

usage() {
  cat << EOF
${HEADING}rocflu.sh${PLAIN}: build, check, run and test Rocflu + ppiclF (a small front end to CMake).

${HEADING}Commands:${PLAIN}
  ${OPTION}./rocflu.sh build${PLAIN} [options]       configure and compile (every particle slot)
  ${OPTION}./rocflu.sh run${PLAIN} <case>.inp ...    check the case, build for it, pre-process it
                                    here and submit the simulation to Slurm
  ${OPTION}./rocflu.sh post${PLAIN} <case>.inp ...   convert the case's dumps (restart_data/) into
                                    ParaView files in post/ (rflupost)
  ${OPTION}./rocflu.sh validate${PLAIN} <case>.inp   check a case's input and files only
  ${OPTION}./rocflu.sh test${PLAIN} [options]        run the regression suite (tests/run_tests.py)
  ${OPTION}./rocflu.sh clean${PLAIN}                 remove build/
  ${OPTION}source ./rocflu.sh load${PLAIN}           load the HiPerGator modules into your shell

${HEADING}Environment:${PLAIN}
  ${OPTION}ROCFLU_MODULES${PLAIN}    modules loaded by build, run, test and load
                    ${DEFAULT}(default: $DEFAULT_MODULES)${PLAIN}; "" loads none
  ${OPTION}NO_COLOR${PLAIN}          set to anything to print without colours

Run ${OPTION}./rocflu.sh <command> --help${PLAIN} for the options of a command.
EOF
}

VALIDATOR="$ROOT/utilities/scripts/case_validator.py"

usage_validate() {
  cat << EOF
${HEADING}Usage:${PLAIN} ${OPTION}./rocflu.sh validate${PLAIN} <dir>/<case>.inp

Checks the case before it runs and reports every problem at once:
unknown, removed or misplaced keys in # PICL, values out of range, particle
models that do not work together, the time-stepping and viscosity input
particles need, points.dat (format, particles inside the fluid domain), and
files in the case directory that would stop the run (filein.vtk,
par00002.vtu from an earlier run). No options.
./rocflu.sh run does this first. See docs/InputFile.md.
EOF
}

cmd_validate() {
  case "$1" in -h|--help|"") usage_validate; [ -n "$1" ]; exit $? ;; esac
  [ $# == 1 ] || fail "validate takes one case file"
  [ -f "$1" ] || fail "no file $1"
  exec python3 "$VALIDATOR" "$1"
}

usage_build() {
  cat << EOF
${HEADING}Usage:${PLAIN} ${OPTION}./rocflu.sh build${PLAIN} [options]

${HEADING}Targets${PLAIN} ${DEFAULT}(default: all)${PLAIN}:
  ${OPTION}-t, --targets${PLAIN} T ...   one or more of
                          simulation    the solver rflump
                          pre_process   rflumap, rflupart and rfluinit
                          post_process  the post program rflupost
                          rflump rflumap rflupart rfluinit rflupost
${HEADING}Build type and physics:${PLAIN}
  ${OPTION}--debug${PLAIN}               -g -O0 -fcheck=all ${DEFAULT}(default: optimised release build)${PLAIN}
  ${OPTION}--spec, --plag${PLAIN}        species transport (needed by the mixture and JWL
                        gas models) and Rocpart Lagrangian particles:
                        ${DEFAULT}disabled for now${PLAIN}, not tested for a while
  ${OPTION}--no-picl${PLAIN}             without the ppiclF particle solver ${DEFAULT}(default: with ppiclF)${PLAIN}
  ${OPTION}--no-folder${PLAIN}           write all output in the run directory
                        ${DEFAULT}(default: one folder per kind of output)${PLAIN}
${HEADING}Other:${PLAIN}
  ${OPTION}-j, --jobs${PLAIN} N          compile with N processes ${DEFAULT}(default: 16)${PLAIN}
  ${OPTION}--metis-root${PLAIN} DIR      use an installed METIS 5
                        ${DEFAULT}(default: build the copy in 3rd_party/)${PLAIN}
  ${OPTION}--fresh${PLAIN}               reconfigure from scratch ${DEFAULT}(default: reuse the configuration)${PLAIN}
  ${OPTION}-v, --verbose${PLAIN}         print every compile command ${DEFAULT}(default: off)${PLAIN}

Each combination of build type and physics options has its own name (for
example release, debug-nopicl): it is compiled in build/staging/<name> and its
executables are copied to build/install/<name>/bin. Switching between
combinations recompiles nothing. build/bin points at the last one built.
EOF
}

cmd_build() {
  local targets=() jobs=16 debug=0 spec=0 plag=0 picl=1 folder=1 fresh=0
  local verbose=0 metis_root=""
  while [ $# -gt 0 ]; do
    case "$1" in
      -t|--targets) shift
        while [ $# -gt 0 ] && [ "${1:0:1}" != "-" ]; do targets+=("$1"); shift; done
        [ ${#targets[@]} -gt 0 ] || fail "-t needs at least one target"
        continue ;;
      -j|--jobs)    jobs="$2"; shift ;;
      -j*)          jobs="${1#-j}" ;;
      --debug)      debug=1 ;;
      --spec|--plag)
        fail "$1 is disabled for now: species transport and Rocpart have not been tested for a while (see CMakeLists.txt)" ;;
      --no-picl)    picl=0 ;;
      --no-folder)  folder=0 ;;
      --metis-root) metis_root="$(cd "$2" && pwd)" || fail "no directory $2"; shift ;;
      --fresh)      fresh=1 ;;
      -v|--verbose) verbose=1 ;;
      -h|--help)    usage_build; exit 0 ;;
      *)            usage_build; fail "unknown option $1" ;;
    esac
    shift
  done

  # Target names -> CMake targets
  local cmake_targets=()
  for t in "${targets[@]}"; do
    case "$t" in
      simulation)  cmake_targets+=(rflump) ;;
      pre_process) cmake_targets+=(rflumap rflupart rfluinit) ;;
      post_process) cmake_targets+=(rflupost) ;;
      rflump|rflumap|rflupart|rfluinit|rflupost) cmake_targets+=("$t") ;;
      *) fail "unknown target $t (see ./rocflu.sh build --help)" ;;
    esac
  done
  [ ${#cmake_targets[@]} -gt 0 ] || cmake_targets=(all)

  local name type
  if [ $debug == 1 ]; then name=debug; type=Debug; else name=release; type=Release; fi
  [ $spec == 1 ]   && name+="-spec"
  [ $plag == 1 ]   && name+="-plag"
  [ $picl == 0 ]   && name+="-nopicl"
  [ $folder == 0 ] && name+="-nofolder"
  [ -n "$metis_root" ] && name+="-metis"

  on_off() { [ "$1" == 1 ] && echo ON || echo OFF; }
  load_modules || fail "could not load modules: $ROCFLU_MODULES"
  build_config "$name" "$jobs" "$fresh" "$verbose" "${cmake_targets[*]}" \
    -DCMAKE_BUILD_TYPE=$type -DROCFLU_CASE= \
    -DROCFLU_PICL=$(on_off $picl) -DROCFLU_SPEC=$(on_off $spec) \
    -DROCFLU_PLAG=$(on_off $plag) -DROCFLU_FOLDER=$(on_off $folder) \
    -DMETIS_ROOT="$metis_root"
  ln -sfn "install/$name/bin" "$BUILD_ROOT/bin"
  say "executables in build/install/$name/bin (build/bin points there)"
}

# Configures, builds and installs one configuration:
#   build_config <name> <jobs> <fresh 0|1> <verbose 0|1> "<targets>" <cmake -D options...>
# CMake's files go to build/staging/<name>, the executables to
# build/install/<name>/bin. The compiler output goes to
# build/staging/<name>/build.log; the terminal shows a progress bar and a
# line per finished part (plain lines when the output is not a terminal,
# every line with -v).
build_config() {
  local name="$1" jobs="$2" fresh="$3" verbose="$4" targets="$5"
  shift 5
  local dir="$BUILD_ROOT/staging/$name"      # CMake's files and the objects
  local inst="$BUILD_ROOT/install/$name"     # the executables
  local log="$dir/build.log" start=$SECONDS
  command -v cmake > /dev/null || fail "cmake not found"

  [ "$fresh" == 1 ] && rm -rf "$dir" "$inst"
  mkdir -p "$dir"
  say "building ${OPTION:-}$name${PLAIN:-} ($jobs jobs), log: ${log#$ROOT/}"

  # Configure
  if ! cmake -S "$ROOT" -B "$dir" -DCMAKE_INSTALL_PREFIX="$inst" "$@" \
         > "$dir/configure.log" 2>&1; then
    build_failed configure "" "$dir/configure.log"
  fi
  checkpoint "configure" $((SECONDS - start))

  # Compile, reading CMake's lines as they come
  local extra=() line percent=0 file="" status="" programs=() seen=" "
  [ "$verbose" == 1 ] && extra=(--verbose)
  : > "$log"
  # shellcheck disable=SC2086
  while IFS= read -r line; do
    case "$line" in
      "__rocflu_status__ "*) status="${line#__rocflu_status__ }"; continue ;;
    esac
    echo "$line" >> "$log"
    [ "$verbose" == 1 ] && echo "$line"
    if [[ "$line" =~ ^\[\ *([0-9]+)%\]\ (Building|Linking|Generating)\ .*[\ /]([^ /]+)$ ]]; then
      percent="${BASH_REMATCH[1]}"; file="${BASH_REMATCH[3]%.o}"
    elif [[ "$line" =~ Built\ target\ ([A-Za-z0-9_]+) ]] \
         && [[ "$seen" != *" ${BASH_REMATCH[1]} "* ]]; then
      # make reports a target again for each target that needs it: once is enough
      seen+="${BASH_REMATCH[1]} "
      case "${BASH_REMATCH[1]}" in
        gs)            checkpoint "gslib" $((SECONDS - start)) ;;
        metis_bundled) checkpoint "METIS + GKlib" $((SECONDS - start)) ;;
        ppiclF)        checkpoint "ppiclF" $((SECONDS - start)) ;;
        flu)           checkpoint "Rocflu library" $((SECONDS - start)) ;;
        rflump|rflumap|rflupart|rfluinit|rflupost) programs+=("${BASH_REMATCH[1]}") ;;
      esac
    fi
    [ "$verbose" == 1 ] || progress_bar "$percent" "$file" $((SECONDS - start))
  done < <(cmake --build "$dir" -j "$jobs" --target $targets "${extra[@]}" 2>&1
           echo "__rocflu_status__ $?")
  clear_bar
  [ "$status" == 0 ] || build_failed compile "$file" "$log"
  [ ${#programs[@]} -gt 0 ] && checkpoint "${programs[*]}" ""

  cmake --install "$dir" > /dev/null || fail "copying the executables to build/install/$name failed"
  local warnings
  warnings=$(grep -c "Warning:" "$log")
  printf '%s✔ Build completed:%s %s in %d s, ' "${COMPLETED:-}" "${PLAIN:-}" "$name" $((SECONDS - start))
  if [ "$warnings" == 0 ]; then
    printf '%s0 warnings%s\n' "${PROGRESS:-}" "${PLAIN:-}"
  else
    printf '%s%d warnings%s (grep -B4 Warning: %s)\n' "${DEFAULT:-}" "$warnings" "${PLAIN:-}" "${log#$ROOT/}"
  fi
}

# Progress display for build_config. On a terminal the bar is redrawn in
# place; elsewhere only the checkpoint lines are printed.
on_terminal() { [ -t 1 ]; }
progress_bar() {   # <percent> <file> <seconds>
  on_terminal || return 0
  local percent="$1" width=20 filled bar=""
  filled=$((percent * width / 100))
  printf -v bar '%*s' "$filled" ''; bar="${bar// /█}"
  local rest; printf -v rest '%*s' $((width - filled)) ''; rest="${rest// /░}"
  printf '\r\033[K  [%s%s%s%s] %3d%%  %s%-28.28s %d:%02d%s' \
    "${OPTION:-}" "$bar" "${PLAIN:-}" "$rest" "$percent" "${DIM:-}" "$2" \
    $(($3 / 60)) $(($3 % 60)) "${PLAIN:-}"
}
clear_bar() { on_terminal && printf '\r\033[K'; return 0; }
checkpoint() {     # <what> <seconds or "">
  clear_bar
  if [ -n "$2" ]; then
    printf '  %s✓%s %-32s %s%3d s%s\n' "${PROGRESS:-}" "${PLAIN:-}" "$1" "${DIM:-}" "$2" "${PLAIN:-}"
  else
    printf '  %s✓%s %s\n' "${PROGRESS:-}" "${PLAIN:-}" "$1"
  fi
}

# Reports a failed configure or compile step in red, with the compiler's
# messages, and stops. <step> configure|compile, <file> last file started.
build_failed() {   # <step> <file> <log>
  local step="$1" file="$2" log="$3" part="" object=""
  clear_bar
  if [ "$step" == configure ]; then
    part="configure (CMake)"
  else
    # make names the object it could not build: [...dir/<target>.dir/.../x.o] Error 1
    object=$(grep -m1 -E '\] Error [0-9]+' "$log" | sed -E 's/.*[[ ]([^][ ]+)\] Error [0-9]+.*/\1/')
    case "$object" in
      *CMakeFiles/gs.dir/*)                 part="gslib" ;;
      *GKlib.dir/*|*metis_bundled.dir/*)    part="METIS + GKlib" ;;
      *ppiclF.dir/*)                        part="ppiclF" ;;
      *CMakeFiles/flu.dir/*)                part="Rocflu library" ;;
      *rflump.dir/*|bin/rflump)             part="rflump" ;;
      *rflumap.dir/*|bin/rflumap)           part="rflumap" ;;
      *rflupart.dir/*|bin/rflupart)         part="rflupart" ;;
      *rfluinit.dir/*|bin/rfluinit)         part="rfluinit" ;;
      *)                                    part="build" ;;
    esac
    [ -n "$object" ] && file="$(basename "${object%.o}")"
    [ "$file" == "$part" ] && file=""
  fi
  {
    printf '  %s✗ %s%s%s\n' "${ERROR_MARK:-}" "$part" "${file:+: $file}" "${ERROR_PLAIN:-}"
    echo
    # The compiler's and linker's messages with the lines before them
    # (gfortran prints file:line:column, the source line and a marker first)
    awk -v red="${ERROR_MARK:-}" -v plain="${ERROR_PLAIN:-}" '
      /(^|[^A-Za-z])(Error|Fatal Error|error):|undefined reference|CMake Error/ {
        if (shown >= 8) { more++; next }
        for (i = n - 4; i < n; i++) if (i >= 0 && !(i in printed)) { print "    " buf[i % 5]; printed[i] = 1 }
        print "    " red $0 plain; printed[n] = 1; shown++; after = 3
      }
      { buf[n % 5] = $0 }
      after > 0 && !(n in printed) { print "    " $0; printed[n] = 1; after-- }
      { n++ }
      END { if (more) printf "    ... and %d more messages\n", more }
    ' "$log" | grep -v "gmake\[[0-9]*\]: \*\*\*\|^    gmake: \*\*\*\|^    \[ *[0-9]*%\]\|^    -- " \
      | sed "s|$ROOT/||g"
    echo
    # What these messages usually mean in this code
    local hint=""
    if grep -q "has no IMPLICIT type" "$log"; then
      hint="a variable is not declared, or its declaration is inside an #ifdef that is off in this build"
    elif grep -q "undefined reference" "$log"; then
      hint="a routine is called but not compiled into this build: is its file in the directory's Filelist.txt, or behind an #ifdef that is off?"
    elif grep -q "Cannot open module file" "$log"; then
      hint="a module is used that this build does not compile (an #ifdef, or a missing Filelist.txt entry)"
    elif grep -q "Rank mismatch\|Type mismatch\|Shape mismatch" "$log"; then
      hint="a call does not match the routine's arguments; check declarations that differ between #ifdef branches"
    elif [ "$step" == configure ]; then
      hint="the build options or CMakeLists.txt; the message above names the line"
    fi
    [ -n "$hint" ] && printf '  %shint:%s %s\n\n' "${DEFAULT_ERR:-}" "${ERROR_PLAIN:-}" "$hint"
  } >&2
  fail "the $part step failed; the full output is in ${log#$ROOT/}"
}

clock() {   # <seconds>: m:ss, or h:mm:ss from an hour on
  if [ "$1" -ge 3600 ]; then
    printf '%d:%02d:%02d' $(($1 / 3600)) $(($1 % 3600 / 60)) $(($1 % 60))
  else
    printf '%d:%02d' $(($1 / 60)) $(($1 % 60))
  fi
}

# Runs a pre-processing tool in <dir>, its output in <dir>/<tool>.log. On a
# terminal a status line shows the step the tool is in (the last step its log
# started, without the set-up and clean-up ones) and, once it works region by
# region, the region count and the time left at the pace so far. Rocflu's
# output is written as it comes (by default gfortran holds it back when it
# goes to a file). Ctrl-C stops the tool too. Returns the tool's status.
run_tool() {   # <dir> <executable> <regions> <arguments...>
  local dir="$1" exe="$2" nregions="$3"
  shift 3
  local tool log
  tool=$(basename "$exe"); log="$dir/$tool.log"
  (cd "$dir" && GFORTRAN_UNBUFFERED_PRECONNECTED=y exec "$exe" "$@") > "$log" 2>&1 &
  local pid=$! start=$SECONDS
  trap 'kill $pid 2> /dev/null; clear_bar; fail "$tool stopped (Ctrl-C)"' INT TERM

  # The line Rocflu writes once per region it has done
  local marker=""
  case $tool in
    rflupart) marker="Adding virtual cells done." ;;
    rfluinit) marker="Writing fluid dump" ;;
  esac
  local frames=(⠋ ⠙ ⠹ ⠸ ⠼ ⠴ ⠦ ⠧ ⠇ ⠏) frame=0 step regions_done first=0 left text width
  # The line fits the terminal: a wrapped line could not be redrawn in place
  width=$(( $(tput cols 2> /dev/null || echo 80) - 25 ))
  [ "$width" -lt 20 ] && width=20
  while kill -0 $pid 2> /dev/null; do
    if on_terminal; then
      step=$(grep -E '^ [A-Z][^.]*\.\.\.$' "$log" \
             | grep -vE '^ (Nullifying|Destroying|Deallocating|Creating|Setting|Allocating|Opening|Closing|Printing)' \
             | tail -1 | sed 's/^ //; s/\.\.\.$//; s/ for inviscid higher-order stencil//')
      step="${step,}"
      text="${step:-starting}"
      regions_done=0
      [ -n "$marker" ] && regions_done=$(grep -c "$marker" "$log")
      if [ "$regions_done" -ge "$nregions" ]; then
        text="all $nregions regions done, $step"
      elif [ "$regions_done" -gt 0 ]; then
        # Pace from the end of the first region: what came before is the serial part
        [ "$first" == 0 ] && first=$SECONDS
        text="region $regions_done/$nregions"
        if [ "$regions_done" -gt 1 ]; then
          left=$(( (SECONDS - first) * (nregions - regions_done) / (regions_done - 1) ))
          [ "$left" -gt 0 ] && text+=", about $(clock $left) left"
        fi
        text+=": $step"
      fi
      printf '\r\033[K  %s%s%s %-9s %s%-*.*s%s %s' "${PROGRESS:-}" "${frames[frame]}" \
        "${PLAIN:-}" "$tool" "${DIM:-}" "$width" "$width" "$text" "${PLAIN:-}" \
        "$(clock $((SECONDS - start)))"
      frame=$(( (frame + 1) % ${#frames[@]} ))
    fi
    sleep 1
  done
  trap - INT TERM
  wait $pid
}

usage_run() {
  cat << EOF
${HEADING}Usage:${PLAIN} ${OPTION}./rocflu.sh run${PLAIN} <dir>/<case>.inp [options]

Checks the case (./rocflu.sh validate), builds the code for it (only the
particle slots its input needs), pre-processes the case here (rflumap,
rflupart, rfluinit, one process each: an error in the .inp, .bc or grid shows
up now, before anything is queued), then writes <dir>/<case>.sbatch and
submits the simulation (rflump) to Slurm. The case directory must hold
<case>.inp, <case>.bc, the grid and, with particles, points.dat.

${HEADING}Job:${PLAIN}
  ${OPTION}-N, --nodes${PLAIN} N         nodes ${DEFAULT}(default: 1)${PLAIN}
  ${OPTION}-n, --tasks-per-node${PLAIN} N
                        MPI ranks per node ${DEFAULT}(default: 1)${PLAIN}; the case is split
                        into N x n regions, one per rank
  ${OPTION}-p, --partition${PLAIN} P     Slurm partition ${DEFAULT}(default: hpg-default)${PLAIN}
  ${OPTION}-a, --account${PLAIN} A       Slurm account ${DEFAULT}(default: bala1s)${PLAIN}
  ${OPTION}-q, --qos${PLAIN} Q           primary: the account's QOS (bala1s)
                        burst: its burst QOS (bala1s-b)
                        ${DEFAULT}(default: primary)${PLAIN}
  ${OPTION}-w, --walltime${PLAIN} T      wall time, HH:MM:SS ${DEFAULT}(default: 01:00:00)${PLAIN}
  ${OPTION}--mem${PLAIN} M               memory per node, e.g. 64gb
                        ${DEFAULT}(default: 2gb per rank, --mem-per-cpu=2gb)${PLAIN}
  ${OPTION}-@, --email${PLAIN} ADDRESS   send all Slurm notifications to ADDRESS
                        ${DEFAULT}(default: no e-mail)${PLAIN}
${HEADING}What to do:${PLAIN}
  ${OPTION}-t, --targets${PLAIN} T       pre_process: build and pre-process only, no job
                        simulation: build and submit only (the case was
                        pre-processed before with the same N x n)
                        ${DEFAULT}(default: both)${PLAIN}
  ${OPTION}-e, --engine${PLAIN} E        batch: submit the job script to Slurm (sbatch)
                        interactive: run it here, on a node you have already
                        allocated (srun --pty bash or salloc); the output is
                        shown and saved in <case>.interactive.out
                        ${DEFAULT}(default: batch)${PLAIN}
  ${OPTION}--dry-run${PLAIN}             build, pre-process and write the job script, but
                        do not submit or run it ${DEFAULT}(default: submit)${PLAIN}
  ${OPTION}--no-post${PLAIN}             do not convert the dumps after the simulation
                        ${DEFAULT}(default: when the simulation succeeds, rflupost
                        converts its dumps into ParaView files in post/, one
                        task; ./rocflu.sh post does it later)${PLAIN}
${HEADING}Build:${PLAIN}
  ${OPTION}--debug${PLAIN}               debug build (-g -O0 -fcheck=all)
                        ${DEFAULT}(default: optimised release build)${PLAIN}
  ${OPTION}-j, --jobs${PLAIN} N          compile with N processes ${DEFAULT}(default: 16)${PLAIN}

The build is named after the slot groups (e.g. release-case-qsfluct): cases
with the same physics share it. build/bin is not changed (it stays the
build used by ./rocflu.sh test).
EOF
}

# Slurm QOS of -q: primary is the account's own QOS, burst its burst QOS
qos_name() {   # <primary|burst> <account>
  case "$1" in
    primary) echo "$2" ;;
    burst)   echo "$2-b" ;;
    *) fail "unknown QOS $1 (primary or burst)" ;;
  esac
}

# Build for a case: a quick configure finds the particle slot groups its
# input needs; the build is named after them. Sets name and bin.
build_for_case() {   # <case.inp> <debug 0|1> <jobs>
  local inp="$1" debug="$2" jobs="$3" case_name
  case_name="$(basename "$inp" .inp)"
  local type=Release prefix=release
  if [ "$debug" == 1 ]; then type=Debug; prefix=debug; fi
  local probe="$BUILD_ROOT/staging/.probe-$prefix" groups
  mkdir -p "$BUILD_ROOT/staging"
  groups=$(cmake -S "$ROOT" -B "$probe" -DCMAKE_BUILD_TYPE=$type -DROCFLU_CASE="$inp" 2>&1 \
           | sed -n 's/^-- Particle slot groups: //p')
  [ -n "$groups" ] || fail "could not read the slot groups of $inp (try: cmake -S $ROOT -B $probe -DROCFLU_CASE=$inp)"
  local physics
  physics=$(echo "$groups" | tr ';' '\n' | grep -vxE 'core|props|gas' | paste -sd- -)
  [ -n "$physics" ] || physics=base
  name="$prefix-case-$physics"
  say "case $case_name needs slot groups: $groups"
  build_config "$name" "$jobs" 0 0 "rflump rflumap rflupart rfluinit rflupost" \
    -DCMAKE_BUILD_TYPE=$type -DROCFLU_CASE="$inp"
  bin="$BUILD_ROOT/install/$name/bin"
  say "executables: build/install/$name/bin"
}

cmd_run() {
  local inp="" nodes=1 tasks=1 partition=hpg-default account=bala1s qos=primary
  local walltime=01:00:00 mem="" email="" targets=() dry=0 debug=0 jobs=16
  local engine=batch post=1
  while [ $# -gt 0 ]; do
    case "$1" in
      -N|--nodes)          nodes="$2"; shift ;;
      -n|--tasks-per-node) tasks="$2"; shift ;;
      -p|--partition)      partition="$2"; shift ;;
      -a|--account)        account="$2"; shift ;;
      -q|--qos)            qos="$2"; shift ;;
      -w|--walltime)       walltime="$2"; shift ;;
      --mem)               mem="$2"; shift ;;
      -@|--email)          email="$2"; shift ;;
      -t|--targets) shift
        while [ $# -gt 0 ] && [ "${1:0:1}" != "-" ]; do targets+=("$1"); shift; done
        [ ${#targets[@]} -gt 0 ] || fail "-t needs pre_process and/or simulation"
        continue ;;
      -e|--engine)         engine="$2"; shift ;;
      --dry-run)           dry=1 ;;
      --post)              post=1 ;;   # the default; kept for old commands
      --no-post)           post=0 ;;
      --debug)             debug=1 ;;
      -j|--jobs)           jobs="$2"; shift ;;
      -h|--help)           usage_run; exit 0 ;;
      -*)                  usage_run; fail "unknown option $1" ;;
      *) [ -z "$inp" ] || fail "more than one case given ($inp, $1)"; inp="$1" ;;
    esac
    shift
  done
  [ -n "$inp" ] || { usage_run; fail "no case given"; }
  local slurm_qos
  slurm_qos=$(qos_name "$qos" "$account") || exit 1
  [ -f "$inp" ] || fail "no file $inp"
  case "$inp" in *.inp) ;; *) fail "$inp: the case file must end in .inp" ;; esac
  [[ "$nodes" =~ ^[1-9][0-9]*$ ]] || fail "-N must be a positive integer"
  [[ "$tasks" =~ ^[1-9][0-9]*$ ]] || fail "-n must be a positive integer"
  case "$engine" in
    batch) ;;
    interactive)
      # Only on a node already allocated: never run a simulation on a login node
      [ -n "$SLURM_JOB_ID" ] || [ $dry == 1 ] || fail "-e interactive runs on a node you have allocated; get one first, e.g. srun -A $account -p $partition -n $((nodes * tasks)) -t $walltime --pty bash (or salloc), then run this again there" ;;
    *) fail "unknown engine $engine (batch or interactive)" ;;
  esac

  local pre=1 sim=1 t
  if [ ${#targets[@]} -gt 0 ]; then
    pre=0; sim=0
    for t in "${targets[@]}"; do
      case "$t" in
        pre_process) pre=1 ;;
        simulation)  sim=1 ;;
        *) fail "unknown target $t (pre_process or simulation)" ;;
      esac
    done
  fi

  local case_dir case_name nranks
  case_dir="$(cd "$(dirname "$inp")" && pwd)"
  case_name="$(basename "$inp" .inp)"
  inp="$case_dir/$case_name.inp"
  nranks=$((nodes * tasks))
  [ -f "$case_dir/$case_name.bc" ] || fail "no $case_name.bc in $case_dir"

  # 0. Check the case before building anything
  local check
  if check=$(python3 "$VALIDATOR" "$inp" 2>&1); then
    checkpoint "case checked (./rocflu.sh validate)" ""
  else
    printf '  %s✗ case check%s\n\n' "${ERROR_MARK:-}" "${ERROR_PLAIN:-}" >&2
    echo "$check" | grep -v "^case_validator: [0-9]* errors\? in" \
      | sed "s/^case_validator: ERROR: /    /" >&2
    echo >&2
    fail "the case has errors (above); nothing was built or submitted"
  fi

  load_modules || fail "could not load modules: $ROCFLU_MODULES"

  # 1. Build for the case
  local name bin
  build_for_case "$inp" "$debug" "$jobs"

  # 2. Pre-process here, one process each
  if [ $pre == 1 ]; then
    say "pre-processing $case_name into $nranks regions in $case_dir"
    local tool args start=$SECONDS
    for tool in rflumap rflupart rfluinit; do
      args=(-c "$case_name" -v 1)
      [ $tool == rflumap ] && args=(-c "$case_name" -m 1 -p "$nranks" -r "$nranks" -v 1)
      if ! run_tool "$case_dir" "$bin/$tool" "$nranks" "${args[@]}" \
         || grep -q "ERROR" "$case_dir/$tool.log"; then
        {
          printf '  %s✗ %s%s\n\n' "${ERROR_MARK:-}" "$tool" "${ERROR_PLAIN:-}"
          # Rocflu's error block (ERROR, Function, Called from), else the last lines
          if grep -q "ERROR" "$case_dir/$tool.log"; then
            grep -A4 "ERROR" "$case_dir/$tool.log" | head -12 \
              | sed "s/^/    /; s/\(.*ERROR.*\)/${ERROR_MARK:-}\1${ERROR_PLAIN:-}/"
          else
            tail -15 "$case_dir/$tool.log" | sed "s/^/    /"
          fi
          echo
        } >&2
        fail "$tool failed (log: $case_dir/$tool.log); nothing was submitted"
      fi
      checkpoint "$tool ($tool.log)" $((SECONDS - start))
    done
  fi
  if [ $sim == 0 ]; then
    say "pre-processing done; no job submitted (-t pre_process)"
    return
  fi

  # 3. Job script
  local script="$case_dir/$case_name.sbatch" build title path
  build="$name, commit $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  git -C "$ROOT" diff --quiet HEAD 2>/dev/null || build="$build with uncommitted changes"
  title="Rocflu+ppiclF case $case_name @ "
  path="$inp"
  [ $((${#title} + ${#path})) -gt 105 ] && path="...${path: -$((102 - ${#title}))}"
  title="$title$path"
  {
    echo "#!/bin/bash"
    echo "# Written by rocflu.sh run on $(date '+%Y-%m-%d %H:%M'): $case_name, $nodes node(s) x $tasks rank(s)"
    echo "#SBATCH --job-name=$case_name"
    echo "#SBATCH --account=$account"
    echo "#SBATCH --qos=$slurm_qos"
    echo "#SBATCH --partition=$partition"
    echo "#SBATCH --nodes=$nodes"
    echo "#SBATCH --ntasks-per-node=$tasks"
    if [ -n "$mem" ]; then
      echo "#SBATCH --mem=$mem"
    else
      echo "#SBATCH --mem-per-cpu=2gb"
    fi
    echo "#SBATCH --time=$walltime"
    echo "#SBATCH --output=$case_name.%j.out"
    if [ -n "$email" ]; then
      echo "#SBATCH --mail-user=$email"
      echo "#SBATCH --mail-type=ALL"
    fi
    echo
    [ -n "$ROCFLU_MODULES" ] && echo "module purge; module load $ROCFLU_MODULES"
    echo "cd $case_dir"
    # Header box of the run log (rflump prints the Setup box and the footer)
    cat << EOF

rule() { printf '+%s+\n' "\$(printf -- '-%.0s' \$(seq 107))"; }
row()  { printf '| * %-15s%-35.35s* %-15s%-35.35s  |\n' "\$1" "\$2" "\$3" "\$4"; }
wide() { printf '| * %-15s%-88.88s |\n' "\$1" "\$2"; }
modules=\$(echo "\${LOADEDMODULES:-N/A}" | tr ':' ' ')
rule
printf '| %-105.105s |\n' "$title"
rule
row "Start-time" "\$(date +%T)"                    "Start-date" "\$(date +%F)"
row "Partition"  "\${SLURM_JOB_PARTITION:-N/A}"    "Walltime"   "$walltime"
row "Account"    "\${SLURM_JOB_ACCOUNT:-N/A}"      "QOS"        "\${SLURM_JOB_QOS:-$slurm_qos} ($qos)"
row "Nodes"      "\${SLURM_JOB_NUM_NODES:-1} (\${SLURM_NTASKS:-$nranks} ranks)" "Engine" "$engine"
row "Job Name"   "\${SLURM_JOB_NAME:-$case_name}"  "Job ID"     "\${SLURM_JOB_ID:-N/A}"
wide "Email"     "${email:-N/A}"
wide "Build"     "$build"
wide "Modules"   "\$modules"
rule
EOF
    echo
    # In an interactive session (srun --pty bash) the session's own step holds
    # the job's CPUs and memory: --overlap shares them, --mem=0 takes the
    # job's memory (without both this srun waits forever or is refused)
    local overlap=""
    [ $engine == interactive ] && overlap="--overlap --mem=0 "
    echo "srun --mpi=pmix ${overlap}-n $nranks $bin/rflump -c $case_name -v 1"
    if [ $post == 1 ]; then
      # Only after a simulation that succeeded
      echo "status=\$?"
      echo "[ \$status == 0 ] || exit \$status"
      echo "srun --mpi=pmix ${overlap}-n 1 $bin/rflupost -c $case_name -v 0"
    fi
  } > "$script"
  say "job script: $script"

  # 4. Submit, or run here
  if [ $dry == 1 ]; then
    say "dry run: not submitted (submit with: sbatch $script)"
    return
  fi
  if [ $engine == interactive ]; then
    say "running $case_name here (job $SLURM_JOB_ID); output also in $case_dir/$case_name.interactive.out"
    (cd "$case_dir" && bash "$script" 2>&1 | tee "$case_name.interactive.out"; exit "${PIPESTATUS[0]}") \
      || fail "the run failed (output: $case_dir/$case_name.interactive.out)"
    return
  fi
  command -v sbatch > /dev/null || fail "sbatch not found"
  local out
  out=$(cd "$case_dir" && sbatch "$script" 2>&1) || { echo "$out" >&2; fail "sbatch failed"; }
  checkpoint "$out (output: $case_dir/$case_name.<job id>.out)" ""
}

usage_post() {
  cat << EOF
${HEADING}Usage:${PLAIN} ${OPTION}./rocflu.sh post${PLAIN} <dir>/<case>.inp [options]

Converts the dumps of the case (<dir>/restart_data/) into ParaView files in
<dir>/post/ with rflupost: open <dir>/post/collection.pvd. The fluid fields
are chosen with the *_WRT keys of # POST, the particle fields with those of
# PICL. It uses the build of the case (as ./rocflu.sh run) and can run while
the simulation is still running.

${HEADING}Steps:${PLAIN}
  ${OPTION}-s, --steps${PLAIN} S         new: the steps not yet in post/ ${DEFAULT}(default)${PLAIN}
                        all: every step (e.g. after changing the *_WRT keys)
                        last: the last step
                        N,M,...: the steps given
${HEADING}Where:${PLAIN}
  ${OPTION}-e, --engine${PLAIN} E        here: run on this node, one process ${DEFAULT}(default)${PLAIN}
                        batch: submit a one-task Slurm job (<case>.post.sbatch)
  ${OPTION}-p, --partition${PLAIN} P     Slurm partition for batch ${DEFAULT}(default: hpg-default)${PLAIN}
  ${OPTION}-a, --account${PLAIN} A       Slurm account for batch ${DEFAULT}(default: bala1s)${PLAIN}
  ${OPTION}-q, --qos${PLAIN} Q           for batch: primary (the account's QOS, bala1s) or
                        burst (bala1s-b) ${DEFAULT}(default: primary)${PLAIN}
  ${OPTION}-w, --walltime${PLAIN} T      wall time for batch, HH:MM:SS ${DEFAULT}(default: 00:30:00)${PLAIN}
  ${OPTION}--mem${PLAIN} M               memory for batch, e.g. 16gb ${DEFAULT}(default: 8gb)${PLAIN}
${HEADING}Build:${PLAIN}
  ${OPTION}--debug${PLAIN}               debug build ${DEFAULT}(default: release)${PLAIN}
  ${OPTION}-j, --jobs${PLAIN} N          compile with N processes ${DEFAULT}(default: 16)${PLAIN}
EOF
}

cmd_post() {
  local inp="" steps=new engine=here partition=hpg-default account=bala1s qos=primary
  local walltime=00:30:00 mem=8gb debug=0 jobs=16
  while [ $# -gt 0 ]; do
    case "$1" in
      -s|--steps)     steps="$2"; shift ;;
      -e|--engine)    engine="$2"; shift ;;
      -p|--partition) partition="$2"; shift ;;
      -a|--account)   account="$2"; shift ;;
      -q|--qos)       qos="$2"; shift ;;
      -w|--walltime)  walltime="$2"; shift ;;
      --mem)          mem="$2"; shift ;;
      --debug)        debug=1 ;;
      -j|--jobs)      jobs="$2"; shift ;;
      -h|--help)      usage_post; exit 0 ;;
      -*)             usage_post; fail "unknown option $1" ;;
      *) [ -z "$inp" ] || fail "more than one case given ($inp, $1)"; inp="$1" ;;
    esac
    shift
  done
  [ -n "$inp" ] || { usage_post; fail "no case given"; }
  local slurm_qos
  slurm_qos=$(qos_name "$qos" "$account") || exit 1
  [ -f "$inp" ] || fail "no file $inp"
  case "$inp" in *.inp) ;; *) fail "$inp: the case file must end in .inp" ;; esac
  case "$engine" in here|batch) ;; *) fail "unknown engine $engine (here or batch)" ;; esac
  [[ "$steps" =~ ^(new|all|last|[0-9]+(,[0-9]+)*)$ ]] \
    || fail "-s must be new, all, last or step numbers separated by commas"

  local case_dir case_name
  case_dir="$(cd "$(dirname "$inp")" && pwd)"
  case_name="$(basename "$inp" .inp)"
  inp="$case_dir/$case_name.inp"
  [ -d "$case_dir/restart_data" ] || fail "no restart_data/ in $case_dir: run the case first"

  load_modules || fail "could not load modules: $ROCFLU_MODULES"
  local name bin
  build_for_case "$inp" "$debug" "$jobs"

  if [ $engine == here ]; then
    say "converting the dumps of $case_name (steps: $steps) into $case_dir/post"
    (cd "$case_dir" && "$bin/rflupost" -c "$case_name" -s "$steps" -v 0) \
      || fail "rflupost failed"
    say "open $case_dir/post/collection.pvd in ParaView"
    return
  fi

  local script="$case_dir/$case_name.post.sbatch"
  {
    echo "#!/bin/bash"
    echo "# Written by rocflu.sh post on $(date '+%Y-%m-%d %H:%M'): $case_name, steps $steps"
    echo "#SBATCH --job-name=$case_name-post"
    echo "#SBATCH --account=$account"
    echo "#SBATCH --qos=$slurm_qos"
    echo "#SBATCH --partition=$partition"
    echo "#SBATCH --ntasks=1"
    echo "#SBATCH --mem=$mem"
    echo "#SBATCH --time=$walltime"
    echo "#SBATCH --output=$case_name.post.%j.out"
    echo
    [ -n "$ROCFLU_MODULES" ] && echo "module purge; module load $ROCFLU_MODULES"
    echo "cd $case_dir"
    echo "srun --mpi=pmix -n 1 $bin/rflupost -c $case_name -s $steps -v 0"
  } > "$script"
  command -v sbatch > /dev/null || fail "sbatch not found"
  local out
  out=$(cd "$case_dir" && sbatch "$script" 2>&1) || { echo "$out" >&2; fail "sbatch failed"; }
  checkpoint "$out (output: $case_dir/$case_name.post.<job id>.out)" ""
}

usage_test() {
  cat << EOF
${HEADING}Usage:${PLAIN} ${OPTION}./rocflu.sh test${PLAIN} [options]

Runs tests/run_tests.py with the executables in build/bin (the last build).
Options are passed on to run_tests.py:
  ${OPTION}-o${PLAIN} PATTERN ...        only the cases whose name contains a PATTERN
                        ${DEFAULT}(default: all cases)${PLAIN}
  ${OPTION}-j, --jobs${PLAIN} N          run N cases at once ${DEFAULT}(default: 8)${PLAIN}
  ${OPTION}--account${PLAIN} A           Slurm account for the solver runs
                        ${DEFAULT}(default: \$ROCPICL_TEST_ACCOUNT, else bala1s)${PLAIN}
  ${OPTION}--partition${PLAIN} P         Slurm partition ${DEFAULT}(default: \$ROCPICL_TEST_PARTITION, else hpg-dev)${PLAIN}
  ${OPTION}--local-ranks${PLAIN} N       run cases with at most N ranks here with mpirun
                        ${DEFAULT}(default: 0, every case through srun)${PLAIN}
  ${OPTION}--generate${PLAIN}            rewrite the golden files ${DEFAULT}(default: compare to them)${PLAIN}
  ${OPTION}--keep${PLAIN}                keep the run directories of passing cases
                        ${DEFAULT}(default: delete them)${PLAIN}
  ${OPTION}--list${PLAIN}                list the cases
See tests/README.md.
EOF
}

cmd_test() {
  case "$1" in -h|--help) usage_test; exit 0 ;; esac
  local exe
  for exe in rflumap rflupart rfluinit rflump rflupost; do
    [ -x "$BUILD_ROOT/bin/$exe" ] || fail "no $exe in build/bin ($(readlink "$BUILD_ROOT/bin" 2>/dev/null)): the last build had -t; build again with the same options and no -t"
  done
  load_modules || fail "could not load modules: $ROCFLU_MODULES"
  local sub=run
  if [ "$1" == "--list" ]; then sub=list; shift; fi
  cd "$ROOT/tests" || exit 1
  if [ $sub == list ]; then
    exec python3 run_tests.py list "$@"
  fi
  exec python3 run_tests.py run --exe "$BUILD_ROOT/bin" "$@"
}

case "$1" in
  build) shift; cmd_build "$@" ;;
  run)   shift; cmd_run "$@" ;;
  post)  shift; cmd_post "$@" ;;
  validate) shift; cmd_validate "$@" ;;
  test)  shift; cmd_test "$@" ;;
  clean) rm -rf "$BUILD_ROOT"; say "removed build/" ;;
  ""|-h|--help|help) usage ;;
  *) usage; fail "unknown command $1" ;;
esac
