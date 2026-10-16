#!/bin/bash
# rocflu.sh: build and test Rocflu + ppiclF (a small front end to CMake).
#
#   ./rocflu.sh build [options]     configure and compile
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

say()  { echo "rocflu: $*"; }
fail() { echo "rocflu: ERROR: $*" >&2; exit 1; }

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

usage() {
  sed -n '2,10p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

usage_build() {
  cat << EOF
Usage: ./rocflu.sh build [options]

Targets (default: all):
  -t, --targets T ...   one or more of
                          simulation    the solver rflump
                          pre_process   rflumap, rflupart and rfluinit
                          rflump rflumap rflupart rfluinit
Build type and physics:
  --debug               -g -O0 -fcheck=all instead of the optimised build
  --spec                species transport (needed by the JWL gas model, GASMODEL 7)
  --plag                Rocpart Lagrangian particles
  --no-picl             without the ppiclF particle solver
  --no-folder           write all output in the run directory
Other:
  -j, --jobs N          compile with N processes (default: 16)
  --metis-root DIR      use an installed METIS 5 instead of the bundled copy
  --fresh               reconfigure from scratch
  -v, --verbose         print every compile command

Each combination of build type and physics options has its own name (for
example release, debug-spec): it is compiled in build/staging/<name> and its
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
      --spec)       spec=1 ;;
      --plag)       plag=1 ;;
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
      rflump|rflumap|rflupart|rfluinit) cmake_targets+=("$t") ;;
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
  local dir="$BUILD_ROOT/staging/$name"      # CMake's files and the objects
  local inst="$BUILD_ROOT/install/$name"     # the executables

  on_off() { [ "$1" == 1 ] && echo ON || echo OFF; }
  load_modules || fail "could not load modules: $ROCFLU_MODULES"
  command -v cmake > /dev/null || fail "cmake not found"

  [ $fresh == 1 ] && rm -rf "$dir" "$inst"
  mkdir -p "$BUILD_ROOT/staging"
  say "configuring $name in build/staging/$name"
  cmake -S "$ROOT" -B "$dir" -DCMAKE_BUILD_TYPE=$type -DCMAKE_INSTALL_PREFIX="$inst" \
        -DROCFLU_PICL=$(on_off $picl) -DROCFLU_SPEC=$(on_off $spec) \
        -DROCFLU_PLAG=$(on_off $plag) -DROCFLU_FOLDER=$(on_off $folder) \
        -DMETIS_ROOT="$metis_root" > "$dir.configure.log" 2>&1 \
    || { cat "$dir.configure.log"; fail "configure failed"; }
  rm -f "$dir.configure.log"

  say "building ${cmake_targets[*]} with $jobs processes"
  local extra=()
  [ $verbose == 1 ] && extra=(--verbose)
  cmake --build "$dir" -j "$jobs" --target "${cmake_targets[@]}" "${extra[@]}" \
    || fail "build failed"

  cmake --install "$dir" > /dev/null || fail "copying the executables to build/install/$name failed"
  ln -sfn "install/$name/bin" "$BUILD_ROOT/bin"
  say "done: executables in build/install/$name/bin (build/bin points there)"
}

usage_test() {
  cat << EOF
Usage: ./rocflu.sh test [options]

Runs tests/run_tests.py with the executables in build/bin (the last build).
Options are passed on to run_tests.py, for example:
  -o PATTERN       only the cases whose name contains PATTERN
  -j N             run N cases at once
  --generate       rewrite the golden files
  --list           list the cases
See tests/README.md.
EOF
}

cmd_test() {
  case "$1" in -h|--help) usage_test; exit 0 ;; esac
  local exe
  for exe in rflumap rflupart rfluinit rflump; do
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
  test)  shift; cmd_test "$@" ;;
  clean) rm -rf "$BUILD_ROOT"; say "removed build/" ;;
  ""|-h|--help|help) usage ;;
  *) usage; fail "unknown command $1" ;;
esac
