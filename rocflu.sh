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
  PROGRESS=$'\033[32m'; PLAIN=$'\033[0m'
fi
if colours_on 2; then ERROR_MARK=$'\033[1;31m'; ERROR_PLAIN=$'\033[0m'; fi

usage() {
  cat << EOF
${HEADING}rocflu.sh${PLAIN}: build, check, run and test Rocflu + ppiclF (a small front end to CMake).

${HEADING}Commands:${PLAIN}
  ${OPTION}./rocflu.sh build${PLAIN} [options]       configure and compile (every particle slot)
  ${OPTION}./rocflu.sh run${PLAIN} <case>.inp ...    check the case, build for it, pre-process it
                                    here and submit the simulation to Slurm
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
                          rflump rflumap rflupart rfluinit
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

  on_off() { [ "$1" == 1 ] && echo ON || echo OFF; }
  load_modules || fail "could not load modules: $ROCFLU_MODULES"
  build_config "$name" "$jobs" "$fresh" "$verbose" "${cmake_targets[*]}" \
    -DCMAKE_BUILD_TYPE=$type -DROCFLU_CASE= \
    -DROCFLU_PICL=$(on_off $picl) -DROCFLU_SPEC=$(on_off $spec) \
    -DROCFLU_PLAG=$(on_off $plag) -DROCFLU_FOLDER=$(on_off $folder) \
    -DMETIS_ROOT="$metis_root"
  ln -sfn "install/$name/bin" "$BUILD_ROOT/bin"
  say "done: executables in build/install/$name/bin (build/bin points there)"
}

# Configures, builds and installs one configuration:
#   build_config <name> <jobs> <fresh 0|1> <verbose 0|1> "<targets>" <cmake -D options...>
# CMake's files go to build/staging/<name>, the executables to
# build/install/<name>/bin.
build_config() {
  local name="$1" jobs="$2" fresh="$3" verbose="$4" targets="$5"
  shift 5
  local dir="$BUILD_ROOT/staging/$name"      # CMake's files and the objects
  local inst="$BUILD_ROOT/install/$name"     # the executables
  command -v cmake > /dev/null || fail "cmake not found"

  [ "$fresh" == 1 ] && rm -rf "$dir" "$inst"
  mkdir -p "$BUILD_ROOT/staging"
  say "configuring $name in build/staging/$name"
  cmake -S "$ROOT" -B "$dir" -DCMAKE_INSTALL_PREFIX="$inst" "$@" \
        > "$dir.configure.log" 2>&1 \
    || { cat "$dir.configure.log"; fail "configure failed"; }
  rm -f "$dir.configure.log"

  say "building $targets with $jobs processes"
  local extra=()
  [ "$verbose" == 1 ] && extra=(--verbose)
  # shellcheck disable=SC2086
  cmake --build "$dir" -j "$jobs" --target $targets "${extra[@]}" \
    || fail "build failed"
  cmake --install "$dir" > /dev/null || fail "copying the executables to build/install/$name failed"
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
  ${OPTION}--dry-run${PLAIN}             build, pre-process and write the job script, but
                        do not submit it ${DEFAULT}(default: submit)${PLAIN}
${HEADING}Build:${PLAIN}
  ${OPTION}--debug${PLAIN}               debug build (-g -O0 -fcheck=all)
                        ${DEFAULT}(default: optimised release build)${PLAIN}
  ${OPTION}-j, --jobs${PLAIN} N          compile with N processes ${DEFAULT}(default: 16)${PLAIN}

The build is named after the slot groups (e.g. release-case-qsfluct): cases
with the same physics share it. build/bin is not changed (it stays the
build used by ./rocflu.sh test).
EOF
}

cmd_run() {
  local inp="" nodes=1 tasks=1 partition=hpg-default account=bala1s
  local walltime=01:00:00 mem="" email="" targets=() dry=0 debug=0 jobs=16
  while [ $# -gt 0 ]; do
    case "$1" in
      -N|--nodes)          nodes="$2"; shift ;;
      -n|--tasks-per-node) tasks="$2"; shift ;;
      -p|--partition)      partition="$2"; shift ;;
      -a|--account)        account="$2"; shift ;;
      -w|--walltime)       walltime="$2"; shift ;;
      --mem)               mem="$2"; shift ;;
      -@|--email)          email="$2"; shift ;;
      -t|--targets) shift
        while [ $# -gt 0 ] && [ "${1:0:1}" != "-" ]; do targets+=("$1"); shift; done
        [ ${#targets[@]} -gt 0 ] || fail "-t needs pre_process and/or simulation"
        continue ;;
      --dry-run)           dry=1 ;;
      --debug)             debug=1 ;;
      -j|--jobs)           jobs="$2"; shift ;;
      -h|--help)           usage_run; exit 0 ;;
      -*)                  usage_run; fail "unknown option $1" ;;
      *) [ -z "$inp" ] || fail "more than one case given ($inp, $1)"; inp="$1" ;;
    esac
    shift
  done
  [ -n "$inp" ] || { usage_run; fail "no case given"; }
  [ -f "$inp" ] || fail "no file $inp"
  case "$inp" in *.inp) ;; *) fail "$inp: the case file must end in .inp" ;; esac
  [[ "$nodes" =~ ^[1-9][0-9]*$ ]] || fail "-N must be a positive integer"
  [[ "$tasks" =~ ^[1-9][0-9]*$ ]] || fail "-n must be a positive integer"

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
  python3 "$VALIDATOR" "$inp" || fail "the case has errors; nothing was built or submitted"

  load_modules || fail "could not load modules: $ROCFLU_MODULES"

  # 1. Build for the case. A quick configure finds the slot groups its
  #    input needs; the build is named after them.
  local type=Release prefix=release
  if [ $debug == 1 ]; then type=Debug; prefix=debug; fi
  local probe="$BUILD_ROOT/staging/.probe-$prefix" groups
  mkdir -p "$BUILD_ROOT/staging"
  groups=$(cmake -S "$ROOT" -B "$probe" -DCMAKE_BUILD_TYPE=$type -DROCFLU_CASE="$inp" 2>&1 \
           | sed -n 's/^-- Particle slot groups: //p')
  [ -n "$groups" ] || fail "could not read the slot groups of $inp (try: cmake -S $ROOT -B $probe -DROCFLU_CASE=$inp)"
  local physics
  physics=$(echo "$groups" | tr ';' '\n' | grep -vxE 'core|props|gas' | paste -sd- -)
  [ -n "$physics" ] || physics=base
  local name="$prefix-case-$physics"
  say "case $case_name needs slot groups: $groups"
  build_config "$name" "$jobs" 0 0 "rflump rflumap rflupart rfluinit" \
    -DCMAKE_BUILD_TYPE=$type -DROCFLU_CASE="$inp"
  local bin="$BUILD_ROOT/install/$name/bin"
  say "executables: build/install/$name/bin"

  # 2. Pre-process here, one process each
  if [ $pre == 1 ]; then
    say "pre-processing $case_name into $nranks regions in $case_dir"
    local tool args
    for tool in rflumap rflupart rfluinit; do
      args=(-c "$case_name" -v 1)
      [ $tool == rflumap ] && args=(-c "$case_name" -m 1 -p "$nranks" -r "$nranks" -v 1)
      if ! (cd "$case_dir" && "$bin/$tool" "${args[@]}" > "$tool.log" 2>&1) \
         || grep -q "ERROR" "$case_dir/$tool.log"; then
        echo "---- last lines of $case_dir/$tool.log:" >&2
        tail -15 "$case_dir/$tool.log" >&2
        fail "$tool failed; nothing was submitted"
      fi
      say "  $tool done ($tool.log)"
    done
  fi
  if [ $sim == 0 ]; then
    say "pre-processing done; no job submitted (-t pre_process)"
    return
  fi

  # 3. Job script
  local script="$case_dir/$case_name.sbatch"
  {
    echo "#!/bin/bash"
    echo "# Written by rocflu.sh run on $(date '+%Y-%m-%d %H:%M'): $case_name, $nodes node(s) x $tasks rank(s)"
    echo "#SBATCH --job-name=$case_name"
    echo "#SBATCH --account=$account"
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
    echo "srun --mpi=pmix $bin/rflump -c $case_name -v 1"
  } > "$script"
  say "job script: $script"

  # 4. Submit
  if [ $dry == 1 ]; then
    say "dry run: not submitted (submit with: sbatch $script)"
    return
  fi
  command -v sbatch > /dev/null || fail "sbatch not found"
  local out
  out=$(cd "$case_dir" && sbatch "$script" 2>&1) || { echo "$out" >&2; fail "sbatch failed"; }
  say "$out (output: $case_dir/$case_name.<job id>.out)"
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
  run)   shift; cmd_run "$@" ;;
  validate) shift; cmd_validate "$@" ;;
  test)  shift; cmd_test "$@" ;;
  clean) rm -rf "$BUILD_ROOT"; say "removed build/" ;;
  ""|-h|--help|help) usage ;;
  *) usage; fail "unknown command $1" ;;
esac
