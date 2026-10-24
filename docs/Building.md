# Building Rocflu

The code is built with CMake. The script `rocflu.sh` at the top of the repository
wraps CMake, in the style of MFC's `mfc.sh`. You pick what to build (targets) and
how to build it (flags) on one command line, and the script does the rest: it loads
the modules, configures, compiles and keeps each configuration in its own directory.

Running a case (`rocflu.sh run`, `validate` and `post`) is described in
[Running.md](Running.md), the input file and `points.dat` in [InputFile.md](InputFile.md).

## Contents

- [Quick start](#quick-start)
- [What gets built](#what-gets-built)
- [`rocflu.sh build`](#rocflush-build)
- [`rocflu.sh test`, `clean` and `load`](#rocflush-test-clean-and-load)
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
./rocflu.sh build        # optimised build of everything (about 30 s on 16 cores)
./rocflu.sh test         # regression suite (36 cases, about a minute)
./rocflu.sh run /path/to/case/rectshktb.inp -n 4   # build for the case, pre-process, submit
./rocflu.sh post /path/to/case/rectshktb.inp       # dumps -> ParaView files in post/
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
| `rflupost` | `post_process` | converts the dumps in `restart_data/` into ParaView files in `post/` |

`rflump` is built from `simulation/` (its main program and driver) and the library
`libflu`, which holds the rest of Rocflu (`modfloflu`, `modflu`, `libfloflu`,
`libflu`, `rocflu`, `rocinteract`, `rocpicl`, and the optional physics modules). The
pre-processing tools are built from `utilities/map`, `utilities/part` and
`utilities/init`, the post program from `utilities/post`, each plus `libflu`.

## `rocflu.sh build`

```bash
./rocflu.sh build [options]
```

| Option | Meaning |
|---|---|
| `-t, --targets T ...` | build only these targets: `simulation`, `pre_process`, or executable names (`rflump rflumap rflupart rfluinit`). Default: all |
| `--debug` | debug build: `-g -O0 -fbacktrace -fcheck=all` everywhere |
| `--spec`, `--plag` | species transport (needed by the mixture and JWL gas models, `GASMODEL` 3-7) and Rocpart Lagrangian particles: **disabled for now**, see below |
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
./rocflu.sh build --debug --no-picl      # debug build without particles
```

**Species transport and Rocpart are disabled for now.** They have not been tested or
maintained for a while, so `--spec` and `--plag` stop with a message, and CMake refuses
`ROCFLU_SPEC=ON` and `ROCFLU_PLAG=ON`. Their sources and `#ifdef SPEC`/`PLAG` branches
are still in the code; to re-enable one, remove it from the check near the top of
`CMakeLists.txt`. Inputs that need species (`# SPECIES USED 1`, `GASMODEL` 3-7) stop in
pre-processing and are reported by the case validator.

**One directory per configuration.** Each combination of build type and physics
options has a name: `release`, `debug`, `release-nopicl`, `debug-nopicl`, and so on.
Everything lives under `build/`:

```
build/
  bin -> install/release/bin     the last configuration built (used by rocflu.sh test)
  install/
    release/bin/                 rflump rflumap rflupart rfluinit
    debug/bin/                   ...
  staging/                       CMake's files and the compiled objects; no need to look here
    release/
    debug/
```

Switching between configurations recompiles nothing that is already built. To run a
specific configuration, use its `install/` directory, for example
`build/install/debug/bin/rflump`. If you build only some targets (`-t simulation`),
only those are in that configuration's `bin/`; build the others too before running
`rocflu.sh test` with it.

**What a build shows.** On a terminal, a progress bar with the file being compiled, and a
line for each part as it finishes, then the summary:

```
rocflu: building release (16 jobs), log: build/staging/release/build.log
  ✓ configure                          3 s
  ✓ gslib                              4 s
  ✓ METIS + GKlib                      9 s
  ✓ ppiclF                            21 s
  ✓ Rocflu library                    24 s
  ✓ rflump rflumap rflupart rfluinit
✔ Build completed: release in 31 s, 0 warnings
```

The compiler's full output is in `build/staging/<name>/build.log` (`configure.log` for
CMake). If a step fails, the script names it in red with the file (`✗ Rocflu library:
RFLU_ModMPI.F90`), prints the compiler's or linker's messages with the source lines they
point at, and, for the usual cases (undeclared variable, undefined reference, missing
module, argument mismatch), a hint at what they mean in this code. When the output is
not a terminal (`srun`, a log file) the same lines are printed without the bar; `-v`
prints every line of the build instead. `./rocflu.sh run` shows the same, plus a line for
the case check, each pre-processing tool and the submitted job.

**Build times.** A build from scratch with `-j 16` (the default), on 16 cores of a
HiPerGator compute node (`srun -A bala1s -p hpg-dev -n 1 -c 16 --mem=16gb`), measured
on 745fb70 in a clean checkout:

| Command | Time | Warnings |
|---|---|---|
| `./rocflu.sh build` | 31 s | 0 |
| `./rocflu.sh build --no-picl` | 27 s | 0 |

The parallel build needs nothing built beforehand: CMake orders the Fortran modules
itself. On a login node, with fewer free cores, it takes longer; use `-j` with the number
of cores you have.

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

`./rocflu.sh --help` and `./rocflu.sh <command> --help` list every option with its
default. On a terminal the help and the script's messages are coloured (options cyan,
defaults yellow, errors red); in a log file or a pipe they are plain. Set `NO_COLOR=1` to
turn the colours off.

## Using CMake directly

`rocflu.sh` only runs CMake, so the same build can be done by hand:

```bash
module purge; module load gcc/14.2.0 openmpi/5.0.7
cmake -S . -B build/staging/mine -DCMAKE_BUILD_TYPE=Debug -DROCFLU_PICL=OFF
cmake --build build/staging/mine -j 16 --target rflump   # executables in build/staging/mine/bin
```

| CMake option | Default | `rocflu.sh` flag |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `--debug` (= `Debug`) |
| `ROCFLU_PICL` | `ON` | `--no-picl` |
| `ROCFLU_SPEC` | `OFF` (`ON` refused for now) | `--spec` (disabled) |
| `ROCFLU_PLAG` | `OFF` (`ON` refused for now) | `--plag` (disabled) |
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
honoured, so a file listed under `ifdef SPEC` is compiled only with `ROCFLU_SPEC` (disabled for now).

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
rule   ptdiag   PSEUDOTURB >= 1 and MODEL_WRT >= 1
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
