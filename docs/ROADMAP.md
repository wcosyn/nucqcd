# Roadmap: cleanup and modernization of nucqcd

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
  (`NUCQCD_UPDATE_REFERENCE=1 ctest -L "regression|slow"`) and review `git diff tests/reference/`.
- **Bug fixes remove their `[known-bug]` test tag** (see `tests/README.md`).

## Status

| Step | Content | Status |
|---|---|---|
| 1 | Test suite | ✅ done (PR #2) |
| 2 | Repo hygiene, including branches | ✅ done (PR #4) |
| 3 | Rename the repository | ✅ done (PR #7, tagged-poldeut PRs #1 and #2) |
| 4 | Bug fixes | 🔶 listed bugs fixed (PRs #8, #9); open: He3 normalization, `TKinematics` format check |
| 5 | Low-risk modernization | 🔶 `-fopenmp`, Doxygen and Fortran interop done |
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

## 2. Repo hygiene ✅

This step needed no physics changes, and the tests stayed green throughout. Done on 2026-11-20:

- [x] **Removed junk and stale build files:** `todo`, `list of changes`, `.vscode/.cmaketools.json`,
  `.vscode/c_cpp_properties.json`, `progs/GPD/diDVCS_gen_grid.cpp.backup`, `DIS/evolution/kinds.mod`,
  `share/input/.directory`, `Integrator/Cuba-4.2.1/demo/demo-c.out`, and the unused makefiles
  `progs/deepsdis/makefile` and `progs/2bodymom/src/Makefile`. `.gitignore` entries keep them out.
  `MeanFieldNucleusWrapper` is now installed only into `lib/` (not into the source tree), and
  `MeanFieldNucleus.py` loads it from there instead of a hard-coded path.
- [x] **Moved** `Doxyfile.in` and `layout.xml` to `cmake/doxygen/`.
- [x] **`share/`:**
  - `share/recmosh/` (719 files) is gone. `RecMosh` memoizes its recursion and computes all brackets
    on demand in milliseconds. The cache was wrong in 29 files (shells 2n1+l1+2n2+l2 ≥ 10), because
    its integer key collided once l or Λ reached 10. This changes results of `DoubleNModel`
    (`progs/attic/twonucl`) for those shells. `tests/unit/unit_moshinsky.cpp` checks orthonormality.
  - `bystricky.f`, `fpn_chex.f` and their `B_ampl_*.dat` tables moved to `docs/reference-code/nn_amplitudes/`.
  - Experiment data moved to `share/expdata/` (`deeps/`, read by `DeuteronCross::readin_deeps`, and `bonus/`).
  - `share/input/exc.dat` was removed. `share/nuclei/exc.dat` is a superset of it, and no code reads either file.
- [x] **Duplicate class:** the 2014 `progs/2bodymom/src/ClassGridThick_SCX` was renamed to
  `ClassGridThickNucleus_SCX`, which removes the ODR violation. It has a different API from the
  generic `Glauber/ClassGridThick_SCX`, so the two were not merged. `SCX_test.cpp` and `2bodymom.cpp` turned
  out to be different programs that share a file name, not duplicates.
- [x] **`progs/`:** `camille`, `LDRD_test`, `twonucl`, `DQEinc`, `gammatest`, `cuba_test` and
  `integrator` moved to `progs/attic/`. They are still built with `BUILD_PROGS`. The Python scripts are
  helpers for the program next to them (plotting or fitting its output), so they stay where they are.
- [x] **Branches:**
  - Deleted the merged branches `fix-ww-integral`, `pionT`, `brandon` and `partons_integration`.
  - Archived the unmerged branches as the tags `archive/Alan` (b1–b4 extraction in
    `progs/deepsdis/b1_calc.cpp`, 2024), `archive/evolution_integration` (2018, its change is already on
    master) and `archive/LDRD-branch-0.01` (2014), and deleted the branches. To revive one, run
    `git switch -c <name> archive/<name>`.

## 3. Rename the repository ✅

`physics-code` is too generic for what the repository is: nuclear and hadronic reaction theory,
including Glauber/RMSGA final-state interactions, deuteron and He3 structure, tagged and inclusive
DIS, and deuteron GPDs. Do the rename after the hygiene step, so the new name is in place before the
bug-fix PRs and long before the README.

- [x] **Name:** `nucqcd`, "QCD in nuclei": high-energy (JLab to EIC) reactions on nuclei, used to learn
  about the strong force in nuclei, as opposed to low-energy nuclear structure. Other candidates were
  `qcd-nuclei`, `hinuc`, `hadronuc`, `nucreact`, `glauber-lf` and `cosyn-nuclear`.
- [x] **Renamed on GitHub** (2026-11-20). GitHub permanently redirects the old web and git URLs, as long
  as no new repo named `physics-code` is created under the account. Commit hashes are unchanged.
- [x] **Updated the remote of the local clone:** `git remote set-url origin https://github.com/wcosyn/nucqcd.git`.
- [x] **Updated the name inside the repo:** the CMake `project()`, `cmake/PhysicsCode.cmake` →
  `cmake/NucQCD.cmake`, the `physcode_*` helpers → `nucqcd_*`, `PHYSCODE_*` variables → `NUCQCD_*`
  (including `NUCQCD_UPDATE_REFERENCE`), `nucqcd_test::` in the tests, the `.ref` headers,
  `scripts/install_partons.sh`, the usage comments in `progs/`, README, CLAUDE.md and HISTORY.md.
- [x] **Renamed the local folder** `~/Git/physics-code` → `~/Git/nucqcd` with `mv`, not a fresh clone: a clone would lose the
  untracked grid cache in `share/grids/` (367 MB), `share/gpd_deutgrids/` and the PARTONS configuration in `GPD_Partons/bin`. This requires:
  - a fresh build directory, because `HOMEDIR` in `Tools/constants.hpp` and the RPATHs are absolute paths
  - `make install` again
  - updating `physics-code.code-workspace` and the VS Code settings
  - moving the Claude Code project memory folder (`~/.claude/projects/-home-wcosyn-Git-physics-code`)
- [x] **Updated the docs in `tagged-poldeut`** (its PRs #1 and #2) (github.com/wcosyn/tagged-poldeut), on both `master` and
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
deliberately and removes its `[known-bug]` tag. The fixes below were done on 2026-11-20 in one PR,
one commit per bug; no `[known-bug]` tests remain.

- [x] **He3 was broken.** Since commit `a458b63` (2020), `Tools/getwavemod_tm.f` declared the file name
  as a `bind(C)` character of length 1, so `He3wf` and `He3Cross` opened `'/'` and stopped with a
  Fortran runtime error. It is now a `C_CHAR` array `datname(256)`, converted to a string before
  `OPEN`. `readwave` also frees its module arrays first, so a wave function can be read twice in one
  process (`He3Cross` for proton and then neutron crashed on the second read). The `He3Cross` test is
  now a regular regression test.
- [ ] **`TKinematics::ReadFormatString` never rejects a bad number of variables:** the check
  `(fNVars<2) && (fNVars>4)` is always false (it should be `||`), so format strings with 1 or more
  than 4 variables are not caught. Found during step 5.
- [ ] **He3 normalization (lead).** A 6D Vegas integral of Σ|ψ|² d³p₁d³p₂ of `He3wf("AV18")` gave
  0.918 ± 0.007 (Vegas did not reach its target precision), not 1. Check this with a better integral
  (Jacobi coordinates p₁₂, p₃, as used inside `getwavemod_tm.f`) before trusting He3 results. The
  normalization of the momentum grid in the `.form` files, the `INVHBARC³` conversion in `He3wf::getWF`
  and the sum of the three cyclic Faddeev components are the places to look.
- [x] **`WeakQECross` neutral current:** the NC constructor never initialized `lepton`, but
  `getDiffWeakQECross` dereferenced it, so nuclear NC cross sections segfaulted. NC now takes its
  kinematics from the electron object.
- [x] **NC prefactor:** `E - ω·G_F/…` is now `(E-ω)·G_F/…` (the values were about 1e23 too large).
  The prefactor was also a factor 2 larger than the CC one, although both currents use the same
  G_F/√2 normalization (the NC form factors carry the ½ of the Z coupling); with that fixed, NC/CC ≈ 0.24
  at the test kinematics, as expected. The NC branch was marked "still to be implemented" and had
  never been validated: the NC cross section has not been used in any calculation yet, so it is
  untested. Compare with a published NC calculation before using it (a note says so in the code).
- [x] **Constants in `Tools/constants.hpp.in`:** `M_Z` was 10× too large, `MASSn_G` was in MeV,
  `INVHBARC` was 1/197.328 while `HBARC` is 197.327 (`He3wf.cpp` had its own copy), the proton and
  neutron doc comments were swapped, and `M_P`/`M_N`/`M_NUCL` are now aliases of `MASSP`/`MASSN`/`MASSn`
  instead of a second set of values. Results moved by at most 4e-5 relative (more only near zeros).
- [x] **`Poldeut::operator=` leaked** the old wave function.
- [x] **Cuba Suave crashed** with flatness < 1 ("buffer overflow detected"), and for flatness 1–5 it
  silently returned wrong values with `fail = 1`. The cause is a bug in Cuba 4.2.1 itself:
  `src/suave/Fluct.c` capped the fluctuation with `MaxL` instead of `MinL`, so it overflowed. Patched in
  the bundled copy (see step 7); the programs' flatness values can stay as they are.
- [x] **`progs/momdistr/momdistr.cpp`** Fourier-transformed the lower component F with l instead of l̄.
- [x] **ROOT `.pcm` dictionary files** are now copied to `build/lib`, next to the libraries.
- [x] **Structure functions in `DeuteronCross` and `InclusiveCross`:** they use `NuclStructure` (via
  `DeuteronStructure`), a second wrapper around the same parametrizations as `NucleonStructure`. The
  mutation test had changed `NucleonStructure`, which these classes do not use; a mutation of
  `NuclStructure` is caught. A new test checks the two classes against each other. It found that
  `NuclStructure` still passed Q² instead of Q to MSTW (fixed in `NucleonStructure` in 2016,
  `2e05c63`), and that its copy constructor and assignment left the MSTW pointer uninitialized. Both
  are fixed. Merging the two classes is part of step 6.
  **Impact:** only results computed with the `"MSTW"` option through `NuclStructure` were affected;
  SLAC, CB, Alekhin and the other options were not. No program has `"MSTW"` in its code. The b1
  programs (`b1_calc`, `b1_calc_comp`) use `NucleonStructure`, so they were correct. The programs that
  take the structure-function name from the command line and pass it to `NuclStructure` are
  `deepsdis/bonus_F2nfit_allangles` (through `DeuteronCross`), the `incldis` programs (through
  `InclusiveCross` or directly), `threeres` and `increstest`, and `poldeut` (through `Poldeut`). To
  the author's knowledge, none of them was run with MSTW for published results.
- [x] **FSI integrals in `DeuteronMomDistr` were not converged.** The four 2D `cube_adaptive` FSI
  integrals (used by `DeuteronCross`) asked for `PREC` = 1e-4 with at most 2000 evaluations, and every
  call hit that cap. The results were accurate only to ~1e-5–1e-4 (worst for the `sigmaparam` FSI
  value) and depended on the CPU, so `reg_deuteron_tagged` failed on a second machine. They now use
  relative precision 1e-6 with at most 2e4 evaluations (error <~1e-6, about 13× slower); the FSI
  references moved by up to 1e-4 and were re-baselined, with a test tolerance of 5e-6. Published
  DEEPS/BONuS results were affected only at the 1e-4 level.
- [x] **Tests used installed libraries.** `CMAKE_BUILD_WITH_INSTALL_RPATH` gave the libraries in
  `build/lib` the RUNPATH `<repo>/lib`, so a test loaded e.g. `libDIS` from the build tree but its
  dependencies `libKnockout` and `libGlauber` from the last `make install`. Libraries now use a
  build-tree RPATH; `make install` still sets `<repo>/lib`.
- [x] **`progs/deepsdis/b1_calc.cpp` Azz angular factor:** the longitudinal tensor term used
  `1/4 + 3/4 cos θq` instead of d²₀₀(θq) = `1/4 + 3/4 cos 2θq` (the integrand `res[10]` already used
  `cos 2θq`). Fixed in PR #5. Found while reviewing the archived `Alan` branch, which had the same fix.

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
- [x] **Fortran interop:** file names reached Fortran through unbounded `strcpy`/`strcat` into fixed
  buffers: the 128-character COMMON block strings of the Alekhin and GRSV code (in both
  `NuclStructure` and `NucleonStructure`) and the 256-character He3 file name. A path longer than the
  buffer overwrote memory silently. `Tools/FortranString.hpp` now provides bounded copies that throw
  `std::length_error`: `to_fortran_string` (blank-padded, no NUL, as Fortran expects) and
  `to_c_string` (for the `bind(C)` He3 buffer). `TKinematics::ReadFormatString` checks the length
  of the format string against its 100-character buffer. The COMMON blocks were kept: moving them to
  `ISO_C_BINDING` would mean rewriting the third-party Fortran, since interoperable strings must be
  arrays of single characters. `unit_fortran_string` tests the copies and the length checks.
- [x] **`-fopenmp`:** removed. No C/C++ code contains OpenMP directives, so the flag had no effect.
  The Fortran in `DIS/evolution` has `!$omp` directives but was never built with OpenMP.
- [x] **Doxygen:** `EXTRACT_ALL` is off, so `WARN_IF_UNDOCUMENTED` reports the gaps, written to
  `build/doxygen_warnings.log` by `make doc`. Undocumented entities are no longer in the output.
  The Doxyfile template was upgraded with `doxygen -u` (Doxygen 1.15; no more obsolete-tag warnings),
  and CMake sets `HAVE_DOT` only when Graphviz's `dot` is found.
- [ ] **Doxygen gaps:** with `EXTRACT_ALL` off there are about 1400 warnings, about 1300 of them
  undocumented members, mostly in `Tools` (`TKinematics2to3.h` alone has 300). About 80 are real doc
  errors that were already there: `\param` names that do not match the arguments, unclosed `\f$`
  formulas, and member documentation that matches no declaration. Fix the errors first, then
  document the classes that programs use. The gaps are also a good task for Claude, one class or library per PR:
  it can write the missing comments from the code and the existing docs, and the warning log shows
  when a class is done. Units, conventions and physics meaning still need the author's review,
  because the code often does not state them. Where known, cite in each class the paper or note it
  implements; candidates come from `manual/manual.tex`, references already in the code comments
  and the author's publication list, each marked for the author to confirm. Also replace the Doxygen main page, which is still
  titled "Glauber ISI/FSI RMSGA C++ Project" and shares its comment with the Glauber group
  (`Glauber/AbstractFsiGrid.hpp`), together with the README (step 8).

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
- [ ] **Merge `NuclStructure` into `NucleonStructure`.** They wrap the same parametrizations and
  have already drifted apart once (the MSTW scale, step 4). `DeuteronStructure` uses `NuclStructure`
  and builds F1 for SLAC itself, where `NuclStructure` returns inf.
- [ ] **Very long files and functions:** break them up, for example `Tools/TKinematics2to3.cpp`
  (6.8k lines, one function of about 1.3k lines) and `Glauber/NucleonEMOperator.cpp` (4.2k lines).

## 7. Unbundle Cuba

This is a build-system change only, and the `unit_cuba` tests verify it.

- [ ] Remove `Integrator/Cuba-4.2.1` and the `Integrator/Cuba` symlink. They are 75 of the 97 files in `Integrator/`.
- [ ] CMake first looks for an installed Cuba (in `CUBA_ROOT`, `~/.local` or the system). If none is
  found, it downloads the release tarball via `ExternalProject` (URL and SHA256 pinned) and builds
  it with the current workarounds (`-std=gnu17`, serial `make`).
- [ ] The bundled Cuba 4.2.1 carries a patch to `src/suave/Fluct.c` (`MaxL` → `MinL`, step 4). The
  latest release, 4.2.2, still has the bug (checked 2026-11-21), so keep the patch as a file in the repo
  and apply it in the `ExternalProject` (`PATCH_COMMAND`) and in any install script. A system Cuba
  cannot be patched; consider reporting the bug upstream.
  `unit_cuba` tests Suave with a small flatness, so it catches this.
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
  - a link to `docs/HISTORY.md`

  CLAUDE.md and `manual/manual.tex` already contain parts of this.
