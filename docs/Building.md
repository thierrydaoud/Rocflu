# Building Rocflu + ppiclF

The code is built with CMake. The script `rocflu.sh` at the top of the repository
wraps CMake, in the style of MFC's `mfc.sh`. You pick what to build (targets) and
how to build it (flags) on one command line, and the script does the rest: it loads
the modules, configures, compiles and keeps each configuration in its own directory.

## Contents

- [Quick start](#quick-start)
- [What gets built](#what-gets-built)
- [`rocflu.sh build`](#rocflush-build)
- [`rocflu.sh test`, `clean` and `load`](#rocflush-test-clean-and-load)
- [Running a case](#running-a-case)
- [Using CMake directly](#using-cmake-directly)
- [How the build is organised](#how-the-build-is-organised)
  - [Adding a source file](#adding-a-source-file)
  - [Compiler flags](#compiler-flags)
  - [ppiclF](#ppiclf)
  - [METIS](#metis)
- [Troubleshooting](#troubleshooting)

## Quick start

On HiPerGator:

```bash
./rocflu.sh build        # optimised build of everything, about 2-3 minutes
./rocflu.sh test         # regression suite (35 cases, about a minute)
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
options gets its own directory: `build/release`, `build/debug`, `build/debug-spec`,
`build/release-nopicl`, and so on. Switching between them recompiles nothing that is
already built. `build/bin` is a link to the `bin/` of the last configuration built,
which is what `rocflu.sh test` and `tests/run_tests.py` use. To run a specific
configuration, point at its directory, for example `build/debug/bin/rflump`.

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

## Running a case

A case directory holds `<case>.inp`, `<case>.bc`, the grid `<case>.hyb.asc` and, with
particles, `points.dat`. With the modules loaded (`source ./rocflu.sh load`) and N MPI
ranks:

```bash
B=/path/to/repo/build/bin
$B/rflumap  -c <case> -m 1 -p N -r N -v 1
$B/rflupart -c <case> -v 1
$B/rfluinit -c <case> -v 1
srun -n N $B/rflump -c <case> -v 1      # or mpirun -n N
```

## Using CMake directly

`rocflu.sh` only runs CMake, so the same build can be done by hand:

```bash
module purge; module load gcc/14.2.0 openmpi/5.0.7
cmake -S . -B build/mine -DCMAKE_BUILD_TYPE=Debug -DROCFLU_SPEC=ON
cmake --build build/mine -j 16 --target rflump
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
`build/<config>/libpicl/generated/ppiclf.f`, behind the headers `PPICLF_USER.h` and
`PPICLF_STD.h`, as the old ppiclF makefile did. It also generates the include file
`PPICLF` (all ppiclF common blocks) by running the C preprocessor over the
`PPICLF_*.h` headers. The user's headers in `libpicl/user_files/` take precedence over
the copies in `libpicl/ppiclF/source/`.

### METIS

`rflupart` partitions the grid with METIS 5. CMake builds METIS and its helper library
GKlib from `3rd_party/METIS` and `3rd_party/GKlib`, with the settings of METIS's own
build (32-bit integers and reals). The partitions are identical, file for file, to
those of the prebuilt libraries used before. To use another installation instead, pass
`--metis-root DIR` (or `-DMETIS_ROOT=DIR`), where `DIR` holds `lib/libmetis.a`,
`lib/libGKlib.a` and `include/metis.h`.

## Troubleshooting

- **``GLIBCXX_3.4.32' not found``** when running an executable: the modules are not
  loaded in that shell. Run `source ./rocflu.sh load` (or use `./rocflu.sh test`,
  which loads them).
- **A build directory is in a strange state** (for example after changing compilers):
  `./rocflu.sh build --fresh ...` reconfigures that configuration from scratch;
  `./rocflu.sh clean` removes all of them.
- **To see the exact compile commands:** `./rocflu.sh build -v`, or
  `build/<config>/compile_commands.json`.
