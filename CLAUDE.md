# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A collection of C++ (plus some Fortran) libraries and executables for nuclear and hadronic physics calculations by Wim Cosyn. It covers Glauber/RMSGA final-state interactions, quasi-elastic knockout, deuteron and He3 inclusive and tagged DIS, rho/pion transparencies, neutrino QE scattering, and (optionally) deuteron GPDs via PARTONS. The repo has no test suite. The programs in `progs/` serve as both test drivers and production calculations. More background is in `manual/manual.tex` (LaTeX) and in the Doxygen comments in the headers.

## Build

Dependencies: ROOT (with `mathmore`, `minuit2`), GSL, a Fortran compiler, and CMake ≥ 3.6. C++17 is required. Release builds use `-O3 -march=native -fopenmp`. MPI and PARTONS are optional. PARTONS v5 (`elementary-utils`, `numa`, `partons`, no longer uses Qt) needs SFML 2.6, LHAPDF, APFEL++, CLN, Eigen3 and libxml2. Ubuntu 26.04 only packages SFML 3, which removed the `sf::Thread`/`sf::Mutex`/`sf::Lock` classes PARTONS uses, so SFML 2.6 is built from source. PARTONS and these dependencies are installed in `~/.local`, which CMake searches automatically when `PARTONS_ENABLED` is set; the PARTONS sources are cloned next to this repo in `~/Git`.

```bash
cd build
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DBUILD_PROGS=yes ..
make install -j4          # the VS Code default build task (cwd = build/)
make doc                  # Doxygen API docs -> doc/
```

- `-DBUILD_PROGS=yes`: executables in `progs/` are `EXCLUDE_FROM_ALL` unless this is set. Without it, build one program with `make <target>` (e.g. `make poldeut`), but it is not installed into `bin/`.
- Other options: `-DPARTONS_ENABLED=yes` (builds `GPD_Partons/` and `progs/GPD`), `-DMPI_ENABLED=yes` (adds `TMPI`), `-DSTATIC_BUILD=yes` (adds `*ST` static libs), and `-DCMAKE_BUILD_TYPE=Debug`.
- `make install` copies shared libs into `<repo>/lib` and executables into `<repo>/bin`. Binaries use an RPATH that points at `<repo>/lib` plus the directories of linked external libraries (e.g. `/opt/root/lib`, `~/.local/lib`), so run them from `bin/` after installing.
- Cuba (`Integrator/Cuba` → symlink to `Cuba-4.2.1`) is built as a CMake `ExternalProject` with autotools. It runs serially (`make -j1`) with `-std=gnu17`, because it does not compile as C23 under gcc ≥ 15. The resulting `lib/libcuba.a` is linked by path in `DIS`.
- If CMake cannot find ROOT, check `ROOTSYS` and the search path in `cmake/Modules/FindROOT.cmake`.

## Library layering

Libraries build in this order, and each links against the ones before it:

1. **Integrator** (`libIntegrator`): third-party integration code. `numint/` provides templated Romberg/trapezoid/Clenshaw-Curtis/Fubini routines. `adaptive/cubature.c` provides adaptive cubature. `numint/numint2Cuba.hpp` wraps the `numint::mdfunction` functors for Cuba.
2. **Tools** (CMake project and library named **MePhys**, `libMePhys`): general physics and math. It includes four-vectors, gamma matrices and spinors, kinematics classes (`TKinematics2to2/2to3`, `TElectronKinematics`, `LightConeKin*`), deuteron wavefunctions (`TDeuteron` plus the `T*PWF` parametrizations), nucleon/nuclear structure functions and PDFs (MSTW, CTEQ, Alekhin, GRSV, mostly Fortran), and `He3wf`.
3. **Glauber** (`libGlauber`): mean-field nuclei (`MeanFieldNucleus[Thick]`), nucleon EM and weak current operators, `FastParticle` (the ejected particle and its scattering parameters), and the FSI grid hierarchy described below.
4. **Knockout** (`libKnockout`): quasi-elastic A(e,e'p), 2N knockout, rho/pion transparency, and neutrino QE cross sections. These are built on Glauber grids. `Model`/`Cross` form the main A(e,e'p) pair.
5. **DIS** (`libDIS`): inclusive and tagged DIS on the deuteron and He3, the polarized deuteron (`Poldeut`), resonances, and the EMT light-front convolution. `DIS/evolution/` is a separate Fortran library, `libEvolution`, for QCD evolution.
6. **GPD_Partons** (`libGPD`, optional): deuteron vector and transversity GPD convolutions that use PARTONS.

To add a source file, list it explicitly in that library's `CMakeLists.txt`; there is no globbing. To add a program, create `progs/<dir>/CMakeLists.txt` and follow the pattern in `progs/poldeut/CMakeLists.txt`: `add_executable`, then `target_link_libraries`, then the `EXCLUDE_FROM_ALL` setting when `BUILD_PROGS` is off, then `install` to `${CMAKE_SOURCE_DIR}/bin`. Finally, register the directory in `progs/CMakeLists.txt`.

## Key cross-cutting conventions

- **`Tools/constants.hpp` is generated** from `Tools/constants.hpp.in` at configure time and is gitignored, so edit the `.in` file. The generated file defines the physical constants (masses in MeV, `HBARC`, etc.) as macros. It also defines `HOMEDIR` as the absolute path to `<repo>/share`.
- **The `share/` data directory** holds input data such as nuclear mean-field wavefunctions (`share/input/*.radial*.bound`), PDF grids, and correlation grids. Most physics classes take a `dir`/`homedir` string argument, usually `HOMEDIR`, and read their input relative to it.
- **The FSI grid hierarchy (Glauber).** `AbstractFsiGrid` → `AbstractFsiCTGrid`/`AbstractFsiGridThick` → … → concrete `GlauberGrid`, `GlauberGridThick`, `GlauberDecayGridThick`, `OneGlauberGrid`, `ROMEAGrid`, and `ClassGridThick_SCX`. A user adds `FastParticle`s, then calls `fillGrids()`/`updateGrids()`. The grids are expensive to compute, so they are cached as files under `share/grids/` (gitignored) and reused on later runs. If you change the physics or the grid parameters, delete the stale grid files.
- **ROOT dictionaries.** Each `T*` class in Tools (and `TMFSpinor` in Glauber) has a `*LinkDef.h`, and CMake generates ROOT dictionaries from them via `ROOT_GENERATE_DICTIONARY`. A new class of this kind needs a LinkDef entry in that library's CMakeLists.
- Integrands for `numint` are usually written as a struct functor that holds a pointer to the physics object, plus a static `exec` function. Integration precision is passed in as a parameter, and the integrator choice is often an `int` flag on the class.
- The code mixes MeV and GeV units (for example `MASSP` versus `MASSP_G` in `constants.hpp`), so check which unit a class expects.
- `bin/`, `lib/`, `build/`, and `doc/` are build outputs and are gitignored.
