# Building Rocflu + ppiclF

The code is built with CMake. The script `rocflu.sh` at the top of the repository
wraps CMake, in the style of MFC's `mfc.sh`. You pick what to build (targets) and
how to build it (flags) on one command line, and the script does the rest: it loads
the modules, configures, compiles and keeps each configuration in its own directory.

The input file and `points.dat` are described in [InputFile.md](InputFile.md).

## Contents

- [Quick start](#quick-start)
- [What gets built](#what-gets-built)
- [`rocflu.sh build`](#rocflush-build)
- [`rocflu.sh test`, `clean` and `load`](#rocflush-test-clean-and-load)
- [Running a case: `rocflu.sh run`](#running-a-case-rocflush-run)
  - [Checking a case: `rocflu.sh validate`](#checking-a-case-rocflush-validate)
  - [Running the steps by hand](#running-the-steps-by-hand)
- [Using CMake directly](#using-cmake-directly)
- [How the build is organised](#how-the-build-is-organised)
  - [Adding a source file](#adding-a-source-file)
  - [Compiler flags](#compiler-flags)
  - [ppiclF](#ppiclf)
  - [METIS](#metis)
  - [Particle slots: `ppiclf_slots.txt`](#particle-slots-ppiclf_slotstxt)
- [Troubleshooting](#troubleshooting)

## Quick start

On HiPerGator:

```bash
./rocflu.sh build        # optimised build of everything, about 2-3 minutes
./rocflu.sh test         # regression suite (34 cases, about a minute)
./rocflu.sh run /path/to/case/rectshktb.inp -n 4   # build for the case, pre-process, submit
```

The executables are in `build/bin`. There is nothing to set up first: the script loads
`gcc/14.2.0` and `openmpi/5.0.7` itself, and METIS is built from the copy in `3rd_party/`.

You need CMake 3.20 or newer, a Fortran, C and C++ compiler (only gfortran has been
tested), and MPI. HiPerGator's default `cmake` (3.31) is new enough.

## What gets built

| Executable | Target name | Purpose |
|---|---|---|
| `rflump` | `simulation` | the flow solver, coupled to the ppiclF particle solver |
| `rflumap` | `pre_process` | maps the grid onto regions (one region per MPI rank) |
| `rflupart` | `pre_process` | partitions the grid into those regions (uses METIS) |
| `rfluinit` | `pre_process` | writes the initial flow solution |

`rflump` is built from `simulation/` (its main program and driver) and the library
`libflu`, which holds the rest of Rocflu (`modfloflu`, `modflu`, `libfloflu`,
`libflu`, `rocflu`, `rocinteract`, `rocpicl`, and the optional physics modules). The
pre-processing tools are built from `utilities/map`, `utilities/part` and
`utilities/init`, plus `libflu`.

## `rocflu.sh build`

```bash
./rocflu.sh build [options]
```

| Option | Meaning |
|---|---|
| `-t, --targets T ...` | build only these targets: `simulation`, `pre_process`, or executable names (`rflump rflumap rflupart rfluinit`). Default: all |
| `--debug` | debug build: `-g -O0 -fbacktrace -fcheck=all` everywhere |
| `--spec` | species transport (also needed by the JWL gas model, `GASMODEL 7`) |
| `--plag` | Rocpart Lagrangian particles |
| `--no-picl` | without the ppiclF particle solver |
| `--no-folder` | write all output files in the run directory instead of sub-folders |
| `-j, --jobs N` | compile with N processes (default 16) |
| `--metis-root DIR` | use an installed METIS 5 instead of the bundled copy |
| `--fresh` | delete this configuration's build directory and start again |
| `-v, --verbose` | print every compile command |

Examples:

```bash
./rocflu.sh build -t simulation -j 8     # only rflump
./rocflu.sh build --debug                # debug build of everything
./rocflu.sh build --debug --spec         # debug build with species
```

**One directory per configuration.** Each combination of build type and physics
options has a name: `release`, `debug`, `debug-spec`, `release-nopicl`, and so on.
Everything lives under `build/`:

```
build/
  bin -> install/release/bin     the last configuration built (used by rocflu.sh test)
  install/
    release/bin/                 rflump rflumap rflupart rfluinit
    debug-spec/bin/              ...
  staging/                       CMake's files and the compiled objects; no need to look here
    release/
    debug-spec/
```

Switching between configurations recompiles nothing that is already built. To run a
specific configuration, use its `install/` directory, for example
`build/install/debug/bin/rflump`. If you build only some targets (`-t simulation`),
only those are in that configuration's `bin/`; build the others too before running
`rocflu.sh test` with it.

**Rebuilding after an edit.** Run the same `./rocflu.sh build` line again. Only the
changed files and the files that depend on them are recompiled. A change to one
Rocflu routine takes about 10 seconds; a change to a ppiclF header about 25 seconds.

## `rocflu.sh test`, `clean` and `load`

```bash
./rocflu.sh test                    # run the suite with build/bin
./rocflu.sh test -o shocktube       # only the cases whose name contains "shocktube"
./rocflu.sh test --list             # list the cases
./rocflu.sh test --generate         # rewrite the goldens (see tests/README.md first)
./rocflu.sh clean                   # remove build/ (all configurations)
source ./rocflu.sh load             # load the modules into your own shell
```

`test` passes its options on to `tests/run_tests.py`; see
[`tests/README.md`](../tests/README.md). `load` is only needed when you run the
executables yourself, outside the script, because they need the same compiler
libraries they were built with.

The modules the script loads can be changed with the environment variable
`ROCFLU_MODULES`, for example `ROCFLU_MODULES="gcc/12.2.0 openmpi/4.1.6" ./rocflu.sh build`.
`ROCFLU_MODULES=""` loads nothing (for machines without `module`, or when your shell
is already set up).

## Running a case: `rocflu.sh run`

```bash
./rocflu.sh run /path/to/case/rectshktb.inp -N 1 -n 4 -p hpg-default -w 02:00:00 -@ you@ufl.edu
```

The case directory must hold `<case>.inp`, `<case>.bc`, the grid (`<case>.hyb.asc`) and,
with particles, `points.dat`. `run` then:

0. **Checks the case** ([`rocflu.sh validate`](#checking-a-case-rocflush-validate)). An
   error stops here, before anything is built.
1. **Builds the code for this case.** It reads the PICL and FLOWMODEL sections of the
   input and builds only the particle slots the models that are on need (see
   [Particle slots](#particle-slots-ppiclf_slotstxt)). The build is named after those
   slot groups, e.g. `build/install/release-case-qsfluct/bin`; every case with the same
   physics shares it, so a second case usually compiles nothing.
2. **Pre-processes the case here**, on the login node: `rflumap` (into N x n regions),
   `rflupart`, `rfluinit`, one process each, with their output in `rflumap.log`,
   `rflupart.log` and `rfluinit.log` in the case directory. A mistake in the input, the
   boundary conditions or the grid stops here, before anything is queued.
3. **Writes `<case>.sbatch`** in the case directory and **submits it**: the job runs
   only `rflump`, with `srun` on N nodes x n ranks per node. The solver's output goes to
   `<case>.<job id>.out`.

| Option | Meaning | Default |
|---|---|---|
| `-N, --nodes N` | nodes | 1 |
| `-n, --tasks-per-node N` | MPI ranks **per node**; the case is split into N x n regions | 1 |
| `-p, --partition P` | Slurm partition | `hpg-default` |
| `-a, --account A` | Slurm account | `bala1s` |
| `-w, --walltime T` | wall time `HH:MM:SS` | `01:00:00` |
| `--mem M` | memory per node, e.g. `64gb` | 2 GB per rank (`--mem-per-cpu=2gb`) |
| `-@, --email A` | send all Slurm notifications (start, end, failure) to A | none |
| `-t pre_process` | build and pre-process only, no job | |
| `-t simulation` | build and submit only; the case was pre-processed before with the same N x n | |
| `--dry-run` | build, pre-process and write the job script, but do not submit | |
| `--debug` | debug build of the case | |
| `-j, --jobs N` | compile with N processes | 16 |

`run` does not change `build/bin`, which stays the full build `rocflu.sh test` uses.
Pre-processing runs on the login node, which is fine for test-size grids; very large
grids (a million cells and more) can take long there (to be improved).

### Checking a case: `rocflu.sh validate`

```bash
./rocflu.sh validate /path/to/case/rectshktb.inp
```

`utilities/scripts/case_validator.py` reads the input file the way Rocflu does and
reports every problem at once, with the line it is on:

- keys in `# PICL` that Rocflu would ignore or misread: unknown keys (with the closest
  known key), removed keys (with what replaces them), a key after a space, a key given
  twice, and `QSFLUCTFILT` before `QSFLUCT` (Rocflu reads it as `QSFLUCT`);
- flags that are not whole numbers or out of range;
- particle models that do not work together, e.g. `QSFLUCT` without `QSFORCE`;
- with `USED 1`: `RKSCHEME 2` in `# TIMESTEP`, `VISCOSITY` in `# VISCMODEL`, no
  `AXIFLAG 1`, `DENSITY` and `SPECIFICHEAT` in `# PICL`, and no `filein.vtk` or
  `par00002.vtu` (which would start a particle restart) in the case directory;
- `points.dat`: `!` comment lines, the number of particles, then exactly that many
  lines of `x y z vx vy vz dp` (position, initial velocity, diameter) with `dp > 0`
  (lines in an older form, with a material name first or with `x y z dp` only, are
  reported as such), and every particle inside the fluid domain. The domain
  is the box around the grid's vertices, read from `<case>.hyb.asc` (CENTAUR ASCII,
  `GRIDSRC 0`; about a second for a million vertices). In a domain that is not a box
  (a cylinder, a wedge), a particle can be inside the box but outside the grid.

Settings of a model that is off (e.g. `KSP` without collisions) are ignored. `rocflu.sh
run` and the regression suite run the same checks. The validator's own tests:
`python3 tests/test_case_validator.py`.

### Running the steps by hand

With the modules loaded (`source ./rocflu.sh load`) and N MPI ranks, using a full build
(`build/bin`) or a case build (`build/install/<name>/bin`):

```bash
B=/path/to/repo/build/bin
$B/rflumap  -c <case> -m 1 -p N -r N -v 1
$B/rflupart -c <case> -v 1
$B/rfluinit -c <case> -v 1
srun -n N $B/rflump -c <case> -v 1      # or mpirun -n N
```

A case build stops at start-up if the input asks for a model whose slots it does not
have, with a message such as *Input needs ppiclF slots 'qsfluct' (QSFLUCT >= 1 or
PSEUDOTURB >= 1), not in this build: use ./rocflu.sh run*.

## Using CMake directly

`rocflu.sh` only runs CMake, so the same build can be done by hand:

```bash
module purge; module load gcc/14.2.0 openmpi/5.0.7
cmake -S . -B build/staging/mine -DCMAKE_BUILD_TYPE=Debug -DROCFLU_SPEC=ON
cmake --build build/staging/mine -j 16 --target rflump   # executables in build/staging/mine/bin
```

| CMake option | Default | `rocflu.sh` flag |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `--debug` (= `Debug`) |
| `ROCFLU_PICL` | `ON` | `--no-picl` |
| `ROCFLU_SPEC` | `OFF` | `--spec` |
| `ROCFLU_PLAG` | `OFF` | `--plag` |
| `ROCFLU_FOLDER` | `ON` | `--no-folder` |
| `ROCFLU_UTILITIES` | `ON` | (off: only `rflump`) |
| `METIS_ROOT` | empty (bundled) | `--metis-root` |
| `ROCFLU_CASE` | empty (every particle slot) | set by `rocflu.sh run` |

## How the build is organised

| File | Role |
|---|---|
| `CMakeLists.txt` | options, flags, `libflu`, the four executables |
| `cmake/RocfluFilelist.cmake` | reads the `Filelist.txt` of each source directory |
| `libpicl/CMakeLists.txt` | gslib and ppiclF |
| `cmake/PpiclfCommonBlocks.cmake`, `cmake/PpiclfConcatenate.cmake` | generate ppiclF's `PPICLF` include file and `ppiclf.f` |
| `3rd_party/CMakeLists.txt` | GKlib and METIS |

### Adding a source file

Rocflu's source lists are the `Filelist.txt` files in each directory. To add a routine, put the file in the right directory and add
its name to that directory's `Filelist.txt`. The next build notices the changed
Filelist and reconfigures by itself. `ifdef SPEC`-style blocks in a Filelist are
honoured, so a file listed under `ifdef SPEC` is compiled only with `--spec`.

A new ppiclF user file goes into the list `PPICLF_SOURCES` in
`libpicl/CMakeLists.txt`.

### Compiler flags

| Part | Release | Debug |
|---|---|---|
| Rocflu (Fortran) | `-O -fallow-argument-mismatch` | `-g -O0 -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow` |
| ppiclF | `-O3 -g -fcheck=all -fbacktrace -std=legacy -w` | as Rocflu's debug flags, with `-std=legacy -w` |
| Rocflu C++ (`main.C`, METIS interface) | none | `-g` |
| gslib | `-O2` | `-g -O0` |
| METIS, GKlib | `-O3 -std=c99` | `-O3 -std=c99` |

These are the flags the code was always built with (the make build used until
October 2026), and a release build gives results identical to it, byte for byte. `-fcheck=all` stays on in ppiclF in release builds: it costs about
2% and has caught real bugs. `-ffpe-trap` has no effect because the main program is C++
(see `simulation/rflump.F90` to trap floating-point errors).

Fortran `.F90` files are preprocessed by gfortran itself (`-cpp`). The defines are
`RFLU`, `FOLDER`, `INRT` (with any particle or species module), and `PICL`, `SPEC`,
`PLAG` from the options. The build prints no warnings; please keep it that way.

### ppiclF

ppiclF is compiled as a single file: CMake concatenates the core and user files into
`build/staging/<config>/libpicl/generated/ppiclf.f`, behind the headers `PPICLF_USER.h` and
`PPICLF_STD.h`, as the old ppiclF makefile did. It also generates the include file
`PPICLF` (all ppiclF common blocks) by running the C preprocessor over the
`PPICLF_*.h` headers. The user's headers in `libpicl/user_files/` take precedence over
the copies in `libpicl/source/`.

### METIS

`rflupart` partitions the grid with METIS 5. CMake builds METIS and its helper library
GKlib from `3rd_party/METIS` and `3rd_party/GKlib`, with the settings of METIS's own
build (32-bit integers and reals). The partitions are identical, file for file, to
those of the prebuilt libraries used before. To use another installation instead, pass
`--metis-root DIR` (or `-DMETIS_ROOT=DIR`), where `DIR` holds `lib/libmetis.a`,
`lib/libGKlib.a` and `include/metis.h`.

### Particle slots: `ppiclf_slots.txt`

ppiclF stores each particle's state in `ppiclf_y` (position, velocity, temperature,
spin) and its properties and the gas values at the particle in `ppiclf_rprop`; the
fields it projects to the gas are in `ppiclf_feedbk`. Every entry ("slot") has a name,
`PPICLF_JX`, `PPICLF_R_JRHOF`, `PPICLF_P_JFX`, ..., used in ppiclF and Rocflu.

The slots are listed in `libpicl/user_files/ppiclf_slots.txt`, in groups by the
physics that needs them, with a rule for each optional group:

```
rule   spin     COLLISIONS >= 2
rule   qsfluct  QSFLUCT >= 1 or PSEUDOTURB >= 1
rule   ptdiag   PSEUDOTURB >= 1 and OUTMODEL >= 1
...
y      core     JX JY JZ JVX JVY JVZ JT
y      spin     JOX JOY JOZ
rprop  qsfluct  R_FLUCTFX R_FLUCTFY R_FLUCTFZ
...
```

CMake numbers the slots of the groups a build includes and writes `PPICLF_LAYOUT.h`
into the build directory (`PPICLF_USER.h` includes it). A build without a case
(`rocflu.sh build`) includes every group; a case build (`rocflu.sh run`, or
`-DROCFLU_CASE=<case>.inp`) only the groups whose rule holds for that input, so a
drag-only case carries 7 values in `y` and 11 in `rprop` per particle instead of 10 and 53.
Slots a build does not have are defined as 0: code that only runs with a model is never
reached, and the few places that always write a slot are compiled out
(`#if PPICLF_JOX > 0`). A generated check in `PICL_InitSolver` stops a run whose input
needs a group the build lacks.

To add a slot, add its name to a line of `ppiclf_slots.txt` (a new group needs a `rule`
line too; groups without a rule are always built) and rebuild: CMake reconfigures by
itself and recompiles the files that use the slots.

## Troubleshooting

- **``GLIBCXX_3.4.32' not found``** when running an executable: the modules are not
  loaded in that shell. Run `source ./rocflu.sh load` (or use `./rocflu.sh test`,
  which loads them).
- **A build directory is in a strange state** (for example after changing compilers):
  `./rocflu.sh build --fresh ...` reconfigures that configuration from scratch;
  `./rocflu.sh clean` removes all of them.
- **To see the exact compile commands:** `./rocflu.sh build -v`, or
  `build/staging/<config>/compile_commands.json`.
