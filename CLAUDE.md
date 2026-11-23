# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**nucqcd** ("QCD in nuclei", formerly `physics-code`): a collection of C++ (plus some Fortran) libraries and executables by Wim Cosyn for high-energy reactions on nuclei (JLab to EIC energies), used to learn about the strong force in nuclei. It covers Glauber/RMSGA final-state interactions, quasi-elastic knockout, deuteron and He3 inclusive and tagged DIS, rho/pion transparencies, neutrino QE scattering, and (optionally) deuteron GPDs via PARTONS. A Catch2/CTest test suite lives in `tests/` (see below). The programs in `progs/` serve as both drivers and production calculations. More background is in `manual/manual.tex` (LaTeX) and in the Doxygen comments in the headers. Not everything has Doxygen comments, especially the older code (mostly in `Tools`, e.g. the `TKinematics*`, `TSpinor` and `TDeuteron` classes). After `make doc`, `build/doxygen_warnings.log` lists what is missing.

## Build

Dependencies: ROOT 6 (with `minuit2`; `mathmore` is no longer needed), GSL, a Fortran compiler, and CMake ≥ 3.16. `make doc` needs Doxygen, and Graphviz (`dot`) for the class and include graphs; both are optional, and CMake turns the graphs off when `dot` is missing. The code is compiled as C++17 (ROOT 6.40 itself is built with C++23; ROOT's libraries are linked by path so that its C++ standard is not imposed). Release builds use `-O3 -march=native` (no OpenMP: the C/C++ code has no OpenMP directives). MPI, LHAPDF and PARTONS are optional. LHAPDF (found through its pkg-config file, also in `~/.local`; `-DLHAPDF_ENABLED=no` to skip it) adds the `"LHAPDF:<set>[:<polarized set>]"` structure functions of `NucleonStructure` (`Tools/LHAPDFStructure`); its tests need the sets CT18LO and NNPDFpol11_100 and are skipped without them. PARTONS v5 (`elementary-utils`, `numa`, `partons`, no longer uses Qt) needs SFML 2.6, LHAPDF, APFEL++, CLN, Eigen3 and libxml2. Ubuntu 26.04 only packages SFML 3, which removed the `sf::Thread`/`sf::Mutex`/`sf::Lock` classes PARTONS uses, so SFML 2.6 is built from source. PARTONS and these dependencies are installed in `~/.local`, which CMake searches automatically when `PARTONS_ENABLED` is set; the PARTONS sources are cloned next to this repo in `~/Git`. LHAPDF is built with its Python bindings (`import lhapdf` from the system `python3`; the module goes to the user site-packages, so a venv only sees it with `--system-site-packages`); PDF sets are installed with `lhapdf install <set>` into `~/.local/share/LHAPDF`. To set this up on a new machine, run `scripts/install_partons.sh` (needs ROOT installed for its final test build). PARTONS programs only start with `partons.properties` (and the `logger.properties` it points to) in the executable's directory; the script creates both in `GPD_Partons/bin`, where the `progs/GPD` programs are installed.

```bash
cd build
cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DBUILD_PROGS=yes ..
make install -j4          # the VS Code default build task (cwd = build/)
make doc                  # Doxygen API docs -> doc/, warnings -> build/doxygen_warnings.log
```

- `-DBUILD_PROGS=yes`: executables in `progs/` are `EXCLUDE_FROM_ALL` unless this is set. Without it, build one program with `make <target>` (e.g. `make poldeut`), but it is not installed into `bin/`.
- Other options: `-DLHAPDF_ENABLED=no` (LHAPDF is used when found), `-DPARTONS_ENABLED=yes` (builds `GPD_Partons/` and `progs/GPD`), `-DMPI_ENABLED=yes` (adds `TMPI`), `-DSTATIC_BUILD=yes` (adds `*ST` static libs), and `-DCMAKE_BUILD_TYPE=Debug`.
- `make install` copies shared libs into `<repo>/lib` and executables into `<repo>/bin` (the default `CMAKE_INSTALL_PREFIX` is the repo itself). Binaries use an RPATH that points at `<repo>/lib` plus the directories of linked external libraries (e.g. `/opt/root/lib`, `~/.local/lib`), so run them from `bin/` after installing.
- Cuba (`Integrator/Cuba` → symlink to `Cuba-4.2.1`) is built as a CMake `ExternalProject` with autotools into the build tree. It runs serially (`make -j1`) with `-std=gnu17`, because it does not compile as C23 under gcc ≥ 15. Link the INTERFACE target `cuba` to use it (provides `<Cuba/cuba.h>`, `libcuba.a` and the build dependency); `libcuba.a` is also installed into `lib/`.
- ROOT is found through its own CMake config (`ROOTConfig.cmake`): source `thisroot.sh` (sets `ROOTSYS`), or it is taken from `/opt/root`. GSL and MPI use CMake's built-in find modules; only the PARTONS libraries (`ElementaryUtils`, `NumA++`, `PARTONS`) use find modules in `GPD_Partons/cmake/Modules`, and SFML its own config.

## Roadmap

`docs/ROADMAP.md` holds the ordered cleanup and modernization plan and its status. Read it before cleanup or refactoring work, and update its status in the same PR as the work. `docs/HISTORY.md` records where the code came from (SVN, UGent GitHub, github.com), its contributors, the archived branches, and which programs produced which results.

## Tests

`tests/` holds a Catch2 v3 suite, fetched by CMake on the first configure, and is built when `BUILD_TESTING` is ON (the default). Details are in `tests/README.md`.

```bash
cd build && make -j8 && ctest -LE slow -j8    # unit + regression, ~30 s
ctest -L slow                                # Glauber grids computed from scratch
NUCQCD_UPDATE_REFERENCE=1 ctest -L "regression|slow"   # re-baseline, then review git diff tests/reference/
```

- Labels: `unit` (analytic or exact checks), `regression` (golden values in `tests/reference/*.ref`, captured from master on 2026-11-20), `slow`.
- A regression failure is a lead to investigate, not a verdict: the old value may itself be wrong. Re-baseline only when the change is understood.
- Tests use a sandbox share dir (`build/tests/share`, with an empty `grids/`) via `nucqcd_test::share_dir()`, never `share/grids/`; the test `main()` unsets `NUCQCD_GRID_DIR`.
- Bugs found by the suite are documented as `[known-bug]` tests (`[!mayfail]`, or hidden `[.]` when they crash).
- Test target names are prefixed `unit_`/`reg_`. `progs/attic/camille/test` already defines `test_*` targets.

## Library layering

Libraries build in this order, and each links against the ones before it:

1. **Integrator** (`libIntegrator`): third-party integration code. `numint/` provides templated Romberg/trapezoid/Clenshaw-Curtis/Fubini routines. `adaptive/cubature.c` provides adaptive cubature. `numint/numint2Cuba.hpp` wraps the `numint::mdfunction` functors for Cuba.
2. **Tools** (CMake project and library named **MePhys**, `libMePhys`): general physics and math. It includes four-vectors, gamma matrices and spinors, kinematics classes (`TKinematics2to2/2to3`, `TElectronKinematics`, `LightConeKin*`), deuteron wavefunctions (`TDeuteron` plus the `T*PWF` parametrizations), nucleon/nuclear structure functions and PDFs (MSTW, CTEQ, Alekhin, GRSV, mostly Fortran), and `He3wf`.
3. **Glauber** (`libGlauber`): mean-field nuclei (`MeanFieldNucleus[Thick]`), nucleon EM and weak current operators, `FastParticle` (the ejected particle and its scattering parameters), and the FSI grid hierarchy described below.
4. **Knockout** (`libKnockout`): quasi-elastic A(e,e'p), rho/pion transparency, and neutrino QE cross sections. These are built on Glauber grids. `Model`/`Cross` form the main A(e,e'p) pair.
5. **DIS** (`libDIS`): inclusive and tagged DIS on the deuteron and He3, the polarized deuteron (`Poldeut`), resonances, and the EMT light-front convolution. `DIS/evolution/` is a separate Fortran library, `libEvolution`, for QCD evolution.
6. **GPD_Partons** (`libGPD`, optional): deuteron vector and transversity GPD convolutions that use PARTONS.

All CMakeLists use the helpers in `cmake/NucQCD.cmake`:

- `nucqcd_add_library(<name> SOURCES ... LINK ... [LINKDEFS ...])`: shared library (version 1.0.0) with its directory as public include dir, installed into `lib/`; with `STATIC_BUILD` also `<name>ST`. Dependencies in `LINK` are public, so programs only need to name the highest library they use.
- `nucqcd_add_program(<name> [SOURCES ...] LINK ... [DESTINATION dir])` and `nucqcd_add_programs(LINK ... PROGRAMS a b c ...)`: executables (default source `<name>.cpp`), built by `make` and installed only with `BUILD_PROGS`.
- Interface targets for external dependencies: `root_libs`, `GSL::gsl`, `cuba`, `MPI::MPI_CXX`, `partons_libs`.

To add a source file, list it explicitly in that library's `CMakeLists.txt`; there is no globbing. To add a program, add a `nucqcd_add_program(s)` line to the `CMakeLists.txt` of its `progs/` subdirectory (new directories are registered in `progs/CMakeLists.txt`).

## Key cross-cutting conventions

- **`Tools/constants.hpp` is generated** from `Tools/constants.hpp.in` at configure time and is gitignored, so edit the `.in` file. The generated file defines the physical constants as `inline constexpr` variables in namespace `nucqcd` (`nucqcd::mass_p`, `nucqcd::hbarc`, ...; masses in MeV unless the name ends in `_GeV`). The old macro aliases (`MASSP`, `HBARC`, `HOMEDIR`, ...) were removed; only `PI` is left, a literal macro because `Integrator/numint/macros.hpp` defines it too. `nucqcd::share_dir` is the absolute path to `<repo>/share`.
- **The `share/` data directory** holds input data such as nuclear mean-field wavefunctions (`share/input/*.radial*.bound`), PDF grids, and correlation grids. Experiment data used for fits and comparisons (DEEPS, BONuS) are in `share/expdata/`. Most physics classes take a `dir`/`homedir` string argument, usually `nucqcd::share_dir`, and read their input relative to it.
- **The FSI grid hierarchy (Glauber).** `AbstractFsiGrid` → `AbstractFsiCTGrid`/`AbstractFsiGridThick` → … → concrete `GlauberGrid`, `GlauberGridThick`, `GlauberDecayGridThick`, `OneGlauberGrid`, `ROMEAGrid`, and `ClassGridThick_SCX`. A user adds `FastParticle`s, then calls `fillGrids()`/`updateGrids()`. The grids are expensive to compute, so they are cached as files under `share/grids/` (gitignored), or in `$NUCQCD_GRID_DIR` if that is set, and reused on later runs (`Tools/GridCache.hpp`; files are written atomically via a temporary file and a rename). The file name is the cache key; it bins the particle kinematics (10 MeV in p, 0.1 rad in θ and φ), and `addParticle` moves each particle to the centre of its bin (`FastParticle::gridBinCentre`), so grids do not depend on which kinematics came first. The key does not record the physics, so if you change the physics or the grid parameters, delete the stale grid files.
- **No ROOT dictionaries.** The `T*` classes in Tools (and `TMFSpinor`) were ROOT classes (`TObject`/`TNamed` bases, `ClassDef`, LinkDef dictionaries); that was removed (roadmap step 9, phase 1), since the code never used ROOT I/O. The libraries still use ROOT types (`TVector3`, `TLorentzVector`, `TString`, ...); see step 9 for the plan to remove those.
- Integrands for `numint` are usually written as a struct functor that holds a pointer to the physics object, plus a static `exec` function. Integration precision is passed in as a parameter, and the integrator choice is often an `int` flag on the class.
- **Errors** are exceptions: `std::invalid_argument`, `std::out_of_range`, `std::logic_error` or `std::runtime_error`, with the message built by `nucqcd::ErrorMessage` (`Tools/ErrorMessage.hpp`), e.g. `throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range: " << r);`. Do not call `exit()` in the libraries.
- The code mixes MeV and GeV units (for example `nucqcd::mass_p` versus `nucqcd::mass_p_GeV` in `constants.hpp`), so check which unit a class expects.
- `bin/`, `lib/`, `build/`, and `doc/` are build outputs and are gitignored.
