# Roadmap: cleanup and modernization of physics-code

This roadmap follows from an assessment of the repository done in November 2026. The steps are
ordered so that every change is protected by the test suite (`tests/`, see `tests/README.md`). Work
on one step per branch, run `ctest -LE slow` after each change, and update the status here in the
same PR.

## Ground rules

- **Tests first.** The test suite is the safety net for everything below.
- **Work on a branch.** Each step, or each part of a larger step, gets its own branch and PR.
- **Regression failures are leads, not verdicts.** The reference values in `tests/reference/` were
  captured from master on 2026-11-20 and may contain long-standing bugs. When a change moves them,
  find out which value is physically right. If the new value is right, re-baseline deliberately
  (`PHYSCODE_UPDATE_REFERENCE=1 ctest -L "regression|slow"`) and review `git diff tests/reference/`.
- **Bug fixes remove their `[known-bug]` test tag** (see `tests/README.md`).

## Status

| Step | Content | Status |
|---|---|---|
| 1 | Test suite | ✅ done (PR #2) |
| 2 | Repo hygiene, including branches | 🔶 in progress: merged branches deleted |
| 3 | Rename the repository | ⬜ (name to be chosen) |
| 4 | Bug fixes | ⬜ |
| 5 | Low-risk modernization | ⬜ |
| 6 | Structural refactors | ⬜ |
| 7 | Unbundle Cuba | ⬜ (can be done any time after step 2) |
| 8 | Repository README | ⬜ (last) |

---

## 1. Test suite ✅

The suite is a Catch2 v3 + CTest test suite in `tests/`:

- **Unit tests** check exact or analytic results.
- **Regression tests** compare against reference values captured from master.
- **Slow tests** compute Glauber grids.

Usage is described in `tests/README.md`.

## 2. Repo hygiene

This step needs no physics changes, so the tests must stay green throughout.

- [ ] **Remove junk files and stale build files:**
  - `todo` and `list of changes` (a 2010s log, with a space in the filename)
  - `.vscode/.cmaketools.json` and `.vscode/c_cpp_properties.json` (per-machine IDE state)
  - `progs/GPD/diDVCS_gen_grid.cpp.backup`, `DIS/evolution/kinds.mod` (a build artifact),
    `share/input/.directory`, `Integrator/Cuba-4.2.1/demo/demo-c.out`
  - `progs/2bodymom/src/libMeanFieldNucleusWrapper.so`
  - the hand-written makefiles that CMake doesn't use: `progs/deepsdis/makefile` and
    `progs/2bodymom/src/Makefile`
  - add `.gitignore` entries so these don't come back
- [ ] **Move `layout.xml`** next to `Doxyfile.in`, for instance in `cmake/` or a `doc/` source folder.
- [ ] **`share/`:**
  - 719 of its 843 committed files are `share/recmosh/`, a generated Moshinsky-bracket cache
    written by `Tools/recmosh.cpp:63`. Decide whether to generate it on demand instead of committing it.
  - Move the Fortran sources out of `share/input/` (`bystricky.f`, `fpn_chex.f`).
  - Separate the experiment data (`share/deepsdata/`, used only by `progs/deepsdis`) from the library data.
  - Resolve the two different `exc.dat` files (`share/input/` and `share/nuclei/`).
- [ ] **Duplicate class:** `progs/2bodymom/src/ClassGridThick_SCX.*` is a diverged copy of
  `Glauber/ClassGridThick_SCX.*`, and both get linked into the same program (an ODR violation).
  Merge them into one. `SCX_test.cpp` and `2bodymom.cpp` are also duplicated between `progs/` subdirectories.
- [ ] **`progs/` organization:** move the personal or one-off directories (`camille/`, `LDRD_test`,
  `twonucl`, `DQEinc`, `gammatest`, `cuba_test`, `integrator`) into `progs/attic/` or a separate
  repo, and collect the Python scripts that are mixed in with the C++.
- [ ] **Branches:**
  - [x] Deleted the merged branches `fix-ww-integral`, `pionT`, `brandon` and `partons_integration` (2026-11-20).
  - [ ] Review the unmerged branches: merge them, cherry-pick what's useful, or archive them. To
    archive, tag the tip (e.g. `archive/LDRD-branch-0.01`) and delete the branch; the commits stay
    reachable through the tag.
    - `Alan`: 6 commits, 2024-09 ("fixed b2 changes")
    - `evolution_integration`: 2 commits, 2018
    - `LDRD-branch-0.01`: 5 commits, 2014

## 3. Rename the repository

`physics-code` is too generic for what the repository is: nuclear and hadronic reaction theory,
including Glauber/RMSGA final-state interactions, deuteron and He3 structure, tagged and inclusive
DIS, and deuteron GPDs. Do the rename after the hygiene step, so the new name is in place before the
bug-fix PRs and long before the README.

- [ ] **Choose a name.** Candidates:
  - `hadronuc`: hadronic and nuclear physics. Short, and covers the whole scope.
  - `nucreact`: nuclear reactions. Describes most of the code (knockout, FSI, DIS, neutrino QE).
  - `glauber-lf`: the two main tools, Glauber FSI and light-front deuteron physics. Narrower.
  - `cosyn-nuclear`: tied to the author. Unambiguous, but less descriptive.
- [ ] **Rename on GitHub:** `gh repo rename <new>`, or use the repository settings. GitHub permanently
  redirects the old web and git URLs, as long as no new repo named `physics-code` is created under
  the account. Commit hashes are unchanged, so references such as `e9d1e6a` in other repos stay valid.
- [ ] **Update local clones:** `git remote set-url origin https://github.com/wcosyn/<new>.git`.
- [ ] **Update the name inside the repo:**
  - the CMake `project()` name and the header comments of `CMakeLists.txt` and `cmake/PhysicsCode.cmake`
  - `scripts/install_partons.sh`
  - usage comments in `progs/` (`~/physics-code/...`)
  - CLAUDE.md and this roadmap
  - the text in `tests/support/test_support.*` and the `.ref` headers

  Consider renaming `cmake/PhysicsCode.cmake` and the `physcode_*` helpers too.
- [ ] **Optional: rename the local folder** `~/Git/physics-code`. This requires:
  - a fresh build directory, because `HOMEDIR` in `Tools/constants.hpp` and the RPATHs are absolute paths
  - `make install` again
  - updating `physics-code.code-workspace` and the VS Code settings
  - moving the Claude Code project memory folder (`~/.claude/projects/-home-wcosyn-Git-physics-code`)
- [ ] **Update the docs in `tagged-poldeut`** (github.com/wcosyn/tagged-poldeut), on both `master` and
  the `photoproduction` branch. That branch exists only on GitHub; check it out with
  `git switch photoproduction`. Refer to the new name, and say "formerly known as `physics-code`" once,
  e.g. in the README and `docs/DEVELOPMENT_NOTES.md`. Mentions as of 2026-11-20:
  - `master`, 23 lines: `README.md`, `docs/DEVELOPMENT_NOTES.md`, `notebooks/build_notebook.py`,
    `notebooks/paper_figures.ipynb`, and docstrings and comments in `taggedpd/{constants,nucleon,wavefunctions}.py`
    and `tests/test_nucleon.py`
  - `photoproduction`, 44 lines: the same files, plus `taggedpd/{christy_bosted,inclusive}.py`,
    `notebooks/{build_inclusive_notebook.py,inclusive.ipynb}`, `tests/test_christy_bosted.py` and
    `tests/reference/fortran/README.md`

  Regenerate the notebooks from their `build_*.py` scripts rather than editing the `.ipynb` by hand.
  The commit hashes cited there stay valid.

## 4. Bug fixes

Fix one bug per commit or PR. Each fix that changes outputs re-baselines its reference values
deliberately and removes its `[known-bug]` tag.

- [ ] **He3 is broken.** Since commit `a458b63`, `Tools/getwavemod_tm.f:25` declares the file name as
  a `bind(C)` character of length 1. `He3wf` and `He3Cross` therefore open `'/'` and stop with a
  Fortran runtime error. Fix it with a C_CHAR array `datname(256)`, converted to a string before `OPEN`.
- [ ] **`WeakQECross` neutral-current constructor:** it never initializes `lepton`, but
  `getDiffWeakQECross` dereferences it, so nuclear NC cross sections segfault or return garbage.
- [ ] **NC prefactor:** `Knockout/WeakQECross.cpp:638` looks like it has a misplaced parenthesis
  (`E - ω·G_F/…` instead of `(E-ω)·G_F/…`; compare the charged-current branch at line 157). As written,
  the values are about 1e23 too large.
- [ ] **Constants in `Tools/constants.hpp.in`:**
  - `M_Z` is 10× too large; it is used on the NC line above.
  - `MASSn_G` evaluates in MeV instead of GeV.
  - `INVHBARC` is 1/197.328 while `HBARC` is 197.327.
  - The proton and neutron doc comments on `MASSN`/`MASSP` (and on the `_G` versions) are swapped.
  - There are two inconsistent mass sets, `MASSP`/`MASSN` and `M_P`/`M_N`.
- [ ] **`Poldeut::operator=` leaks** the old wave function (`DIS/Poldeut.cpp:48`).
- [ ] **Cuba Suave crashes** with flatness < 1 ("buffer overflow detected" or a segfault). Fix the
  parameters in `progs/cuba_test/*.cpp` and `progs/2bodymom/src/SCX_test.cpp`; the demo value is 25.
- [ ] **`progs/momdistr/momdistr.cpp`** Fourier-transforms the lower component F with l instead of l̄.
- [ ] **ROOT `.pcm` dictionary files** are not copied to `build/lib`. This only shows up in crash output.
- [ ] **Check how the structure functions reach `DeuteronCross` and `InclusiveCross`.** A mutation test
  of SLAC F2 was not caught by their tests, so they appear to use a different code path.

## 5. Low-risk modernization

- [ ] **Constants:** `constexpr` constants in a namespace, keeping the macros as aliases during the transition.
- [ ] **Ownership:** use smart pointers for owned members and rule-of-zero classes. About 17 classes own
  raw pointers without copy or assignment operators (e.g. `Knockout/Model`, `Knockout/Cross`,
  `DIS/InclusiveCross`, `Glauber/FsiCorrelator`), and `MeanFieldNucleus[Thick]` has no `operator=`.
- [ ] **Errors:** replace `exit()` in the libraries (about 295 calls) with exceptions, so errors can be
  tested and handled. For example, Christy-Bosted calls `exit()` outside its kinematic range.
- [ ] **Grid cache:**
  - Make the grid cache directory configurable independently of the input directory.
  - Revisit the cache key: the momentum is rounded to 10 MeV (`Glauber/AbstractFsiGrid.cpp:262`).
  - Make cache writes atomic.
- [ ] **Fortran interop:** stop filling COMMON blocks with `strcpy` into fixed 128-byte buffers
  (`Tools/NuclStructure.cpp`). Use `ISO_C_BINDING` and bounded copies.
- [ ] **`-fopenmp`:** remove it, or actually use OpenMP. No C/C++ code contains OpenMP directives.
- [ ] **Doxygen:** turn off `EXTRACT_ALL` so that `WARN_IF_UNDOCUMENTED` reports the gaps.

## 6. Structural refactors

Do these one at a time, running the full test suite (including `slow`) after each.

- [ ] **Split `Tools` (library `MePhys`)** into coherent parts:
  - math and Lorentz algebra
  - kinematics
  - deuteron wave functions
  - structure functions and PDFs
  - nuclear structure
  - MPI
- [ ] **Include paths:** prefix them (e.g. `<mephys/TDeuteron.h>`) and settle on one header extension.
  Consider an `include/` and `src/` split.
- [ ] **Glauber grid hierarchy:** deduplicate it. `GlauberGridThick` and `GlauberDecayGridThick` are
  largely copy-pasted, as are the {plain, Thick} × {CT, Decay} abstract classes.
- [ ] **Model selection:** use enums or factory registries instead of model-name strings (about 700
  string comparisons). A typo should fail at compile time, not with `exit()` at runtime.
- [ ] **Very long files and functions:** break them up, for example `Tools/TKinematics2to3.cpp`
  (6.8k lines, one function of about 1.3k lines) and `Glauber/NucleonEMOperator.cpp` (4.2k lines).

## 7. Unbundle Cuba

This is a build-system change only, and the `unit_cuba` tests verify it.

- [ ] Remove `Integrator/Cuba-4.2.1` and the `Integrator/Cuba` symlink. They are 75 of the 97 files in `Integrator/`.
- [ ] CMake first looks for an installed Cuba (in `CUBA_ROOT`, `~/.local` or the system). If none is
  found, it downloads the release tarball via `ExternalProject` (URL and SHA256 pinned) and builds
  it with the current workarounds (`-std=gnu17`, serial `make`).
- [ ] Keep the `cuba` INTERFACE target name, so nothing that links against it changes. Pin Cuba 4.2.x,
  because the `numint` wrappers rely on the 4.2 API.
- [ ] Optionally, install it into `~/.local` from `scripts/`, the way `scripts/install_partons.sh` does for PARTONS.

## 8. Repository README

Write this last, after the refactors, so it describes the final layout and uses the new name.

- [ ] Cover:
  - what the code computes, and the libraries and how they depend on each other
  - dependencies and how to install them (ROOT, GSL, Cuba, and optionally PARTONS and MPI)
  - how to build, test and run the main programs
  - the data in `share/`
  - which papers to cite

  CLAUDE.md and `manual/manual.tex` already contain parts of this.
