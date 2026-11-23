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
| 4 | Bug fixes | ✅ done: listed bugs (PRs #8, #9, #12); found later: stale FSI+CT grid (#18), kinematics (#19, #20, #33), DoubleNModel/Pair (#21, then moved to the attic in #22), π⁻ Glauber parameters (#28), TDeuteron errors (#31), He3 spin sum (#47) and factor A in `He3Cross` (#48), `TMPI` not compiling with MPI-3 (fixed in the attic with step 6, "Retire MPI"). New bugs found later are added here |
| 5 | Low-risk modernization | ✅ done: `-fopenmp`, Fortran interop (#11), constants (#13, #14), grid cache (#15–#17), ownership (#23), exceptions (#24), Doxygen (#10, #25–#30, #32, #34–#39: no warnings left) |
| 6 | Structural refactors | 🔶 done: constants macros, `NuclStructure` merge, code layout items 1–3 (`Tools` split, prefixed includes, `.hpp`), `Model`/`Cross` renamed; MPI retired; to do: namespace, Glauber grids, model selection, integrators, long files |
| 7 | Unbundle Cuba and cubature | ⬜ (can be done any time after step 2) |
| 8 | Repository README | ⬜ (last, after 9, 10 and 11) |
| 9 | Decouple from ROOT where possible | 🔶 phases 1 (#44), 2 (#45), 3 and 4a (vector types) done; 4b (conversion), 5, 6 to do |
| 10 | LHAPDF interface for the nucleon structure functions | ✅ done (#46), including MSTW through LHAPDF |
| 11 | Run bookkeeping (provenance) for all programs | 🔶 part 1 (design, `RunLog`) done; `Bookkeep` replaced in its 8 programs; `RunLog` in the other programs to do |

## Where we are / next up

*Updated 2026-11-23.*

**Done recently** (PRs #44–#54 and the phase 4a PR): ROOT decoupling phases 1–3; step 10 (LHAPDF);
He3 fixes; step 11 part 1 (`RunLog`); step 6 code layout: `Tools` split into Core, Kinematics,
WaveFunctions and StructureFunctions (#52), then includes as paths from the repository root, `.hpp`
as the only extension, and `Model`/`Cross` renamed to `EMQEHadronCurrent`/`EMQECross` (#53); MPI
retired (`TMPI` and `rhod`/`rhottest` to `progs/attic/mpi`, #54); ROOT phase 4a: own vector and
Lorentz types, bitwise equal to ROOT.

**Recommended next step:** ROOT phase 4b: convert the code to the new `nucqcd::Vector3`,
`LorentzVector`, `Rotation` and `LorentzRotation` (phase 4a added them, bitwise equal to ROOT),
library by library, starting with the Tools kinematics and spinors; add `namespace nucqcd` to the
classes touched where that stays small.

Phase 4b touches ~700 uses in the libraries.

**After that:** step 11 part 2 (`RunLog` in the other ~69 programs, `// ` headers to `# `) together with ROOT phase 6 (programs), and ROOT phase 5
(ROOT optional in the build).

**Open choices for the author:** the license (step 8, see the item there); the PDG constants update (changes values); whether to unbundle Cuba
(step 7); the README (step 8) comes last.

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
- [x] **`TKinematics::ReadFormatString` never rejected a bad number of variables:** the check
  `(fNVars<2) && (fNVars>4)` was always false, and its error branch used `assert`, which does nothing
  in release builds. Found during step 5. It now uses `||` and throws `std::invalid_argument`
  (tested in `unit_kinematics`). `TKinematics2to3` has its own parser and is not affected.
- [x] **He3 normalization** (done 2026-11-22). The norm of `He3wf` was 0.967–0.975 for the five
  parametrizations (the 2026-11 Vegas estimate 0.918 was not converged; a Gauss-Legendre integral in
  Jacobi coordinates p₁₂ = (p₁−p₂)/2, p₃ converges to 1e-5). The cause was nucqcd's wrapper, not the
  external Fortran code or data: `He3wf::getWF` took no ms3 and set ms3 = mA − ms1 − ms2 (0 when that
  was not ±1), but in the D and P waves the nucleon spins need not add up to mA. With ms3 summed over
  independently the norm is 0.9998–0.9999 for all five, the spin-summed density is rotation invariant
  (1e-7; before up to 30% off), and the Faddeev identity ⟨Ψ|Ψ⟩ = 3⟨ψ₃|Ψ⟩ holds (1e-4). The
  charge-dependent (T = 3/2) parts, left out by `He3wf`, change the norm by 1e-5. `getWF` now takes
  ms3 (the old 7-argument calls no longer compile), `He3Cross::getDensity` sums over it, and the
  Fortran m_A state is shared by all `He3wf` objects (the TVector3 version reset it on every call).
  **Impact:** `He3Cross` density and cross section +5.1% (n) and +6.2% (p) at the regression point
  (spectators at 200 MeV), re-baselined; momentum-dependent elsewhere (larger at high momenta, where
  the D wave dominates). Not used for published results (the author). New test `unit_he3wf` (norm,
  proton/neutron split, rotation invariance, antisymmetry).
- [x] **He3Cross: factor A = 3 missing** (done 2026-11-22, found with the above). Checked against the
  LDRD note "Tagged Deep Inelastic Scattering for general nucleus" (he3.tex, by the author; not in the
  repository). The note normalizes the wave function to 1 for each struck nucleon (its optical-theorem
  check gives Σ_N σ_hN), `He3wf` to 1 over all isospin arrangements (antisymmetric, isospin included),
  so each arrangement carries 1/A: S = A Σ|Ψ|². With the spectators not identified (the author's
  choice; for a struck proton the pn pair is then treated like the pp pair of a struck neutron,
  (A−1)!Θ = 1), the spectator momentum distribution now integrates to Z = 2 (struck proton) and N = 1
  (struck neutron); it gave 2/3 and 1/3. `getCross` otherwise matches the note's final formula. New
  `He3Cross::getMomentumDistribution` (rest frame), tested in `unit_he3wf`. **Impact:** `He3Cross`
  density and cross section ×3, re-baselined; not used for published results.
- [x] **Stale FSI+CT grid after `updateGrids()` (found during step 5).** `AbstractFsiGrid::fillGrids`
  sets `filledallgrid` when it computes the grids, but does not clear it when it reads the FSI grid
  from the cache. If a grid object first computes a grid and then, in `updateGrids()` for new
  kinematics, finds the new FSI grid in the cache, `AbstractFsiCTGrid::fillGrids` takes the old
  FSI+CT grid still in memory as the new one, uses it, and overwrites the new FSI+CT cache file with
  it. This affects Glauber+CT results in runs that mix computed and cached grids (the reproducer
  gave a CT grid 11% off, from the previous particle). **Fixed:** reading the FSI grid now clears
  `filledallgrid`, so the FSI+CT grid is read or computed for the new particle. The getters that check
  `filledallgrid` are all in classes with a CT grid, whose `fillGrids` sets it again; `ROMEAGrid` checks
  `filledgrid`, and `DistMomDistrGrid` has its own logic without this bug. Tested in a `slow` case
  of `reg_glauber_grid`. **FSI+CT grid files cached before this fix may be wrong** (those written by
  a run that had computed a grid before reading another one from the cache); delete them if in doubt.
- [x] **`TKinematics2to3WithLabAngles`: overlapping `sprintf` (found during step 5).** The format
  string of each solution was built with `sprintf(formatString, "%s...", formatString, ...)`, where
  the destination is also a source: undefined behaviour. With gcc -O3 the first four tokens were
  lost, and constructing the solutions failed with "Format string should provide 6 independent
  variables" (and exited), so `DoubleNCross` (`progs/attic/twonucl`) could not run. Fixed by copying
  the first tokens first (and using `snprintf`); it was the only such call in the code. Tested in
  `unit_kinematics` (12C(e,e'pn) kinematics: conservation laws, and the solutions are built).
- [x] **`TKinematics2to3WithLabAngles::GetSolution` has no masses for isospin channel 7 (found
  during step 5).** `UpdateKinematics` constructs the solutions as `TKinematics2to3` without passing
  the masses, so for channel 7 (masses given by the caller, used for A(e,e'pN)) their masses are 0
  and their vectors NaN. The object itself is correct, and no code calls `GetSolution`. **Fixed:**
  the masses `fMd, fMn, fMk, fMy` are passed (channels 1–6 set their own masses anyway). Tested in
  `unit_kinematics` (the solutions conserve energy and momentum and reproduce the input).
- [x] **`DoubleNModel` and `Pair` (found during step 5).** A(e,e'NN) with correlated pairs, unfinished
  since the first commit (2012) and never used for results: the 2012 Vegas logs of
  `progs/twonucl` did not converge (correlated −1.65 ± 1.18, uncorrelated 0.183 ± 0.054). Fixed:
  (1) the z component of r₂ in the integrand was φ₂ instead of cos θ₂; (2) with `pw = 1` the FSI
  grids were not created but still read (segfault); now the Glauber factors are 1; (3) `Pair`
  returned NaN at r₁ = ±r₂ (0/0 in the relative and c.m. coordinates, and u(r)/r at r = 0); (4) the
  azimuth of the relative coordinate in `Pair::getrelcoord` was computed from r₁ + r₂ (the c.m.
  azimuth) instead of r₁ − r₂. The headers now say the code is unfinished. Tested in `unit_pair`.
  Not fixed: the integral does not converge at the twonucl SRC kinematics, where the uncorrelated
  matrix element nearly vanishes. Observed on the way: the Paris deuteron u(r)/r and w(r)/r
  (`TDeuteron`) are numerically noisy below about 0.01 fm (w(r)/r = −0.09 at 10⁻⁴ fm), probably
  cancellations in the analytic form; this carries almost no weight in integrals over d³r.
  **Moved to the attic** afterwards (chosen by the author): `DoubleNModel`, `DoubleNCross` (from
  Knockout) and `Pair`, `Newcoef` (from MePhys) now live in `progs/attic/twonucl` and are compiled
  into `twonucl`; `unit_pair` compiles `pair.cpp` and `newcoef.cpp` from there (`nucqcd_add_test`
  takes `SOURCES`). `RecMosh` (Moshinsky brackets, tested in `unit_moshinsky`) stays in MePhys.
- [x] **π⁻ Glauber parameters (found while documenting `FastParticle`).** For `PI_MIN` the constructor
  takes the π⁺p tables for scattering off protons and the π⁺n tables for neutrons, exactly as for
  `PI_PLUS`. By isospin π⁻p = π⁺n, and until commit `1666202` (2012-04-23), when the file reads were
  replaced by static tables, the π⁻ case did swap them (the old call is still there, commented out).
  No library code or program creates a `PI_MIN` explicitly (`PionTCross` uses π⁺), so no known results
  are affected. **Fixed** (with the author's OK): `PI_MIN` swaps the tables, in the constructor and in
  `interpPionGlauberData`. Tested in `unit_grid_bins`. Cached π⁻ grids (`Ptype3` in the file name)
  were computed with the old parameters; the cache key contains σ, so most of them are not reused,
  but delete them to be sure.
- [x] **`TDeuteron::Wavefunction` error handling (found while documenting `TDeuteron`).**
  `CreateWavefunction(name)` prints a warning and returns a null pointer for an unknown name (e.g.
  "Nijmegen" instead of "nijm3"), so the caller crashes later; it should throw
  `std::invalid_argument`. `getResidu()` on a non-Yukawa wave function asserts and, with assertions
  off (release builds), falls off the end without a return value (undefined behaviour); it should
  throw `std::logic_error`. No behaviour change for valid input. **Fixed:** both throw, with the list
  of known names in the message. Tested in `unit_deuteron_wf`.
- [x] **`TKinematics2to3::GetPhik(kS)` ignores the axis convention (found while documenting).**
  For the standard convention it calls `GetPhik(fAC)` without `return`, falls through and returns
  the kaon azimuth in the `kY` convention whatever `fAC` is (`GetPhiy` and `GetPhin` do return). Wrong
  whenever `fAC` is not `kY`, e.g. `GetPhik()` and the variable "phik" (`GetVarByGetter`) with `kN` or
  `kK`. No results are affected: the only call in the code is `fPhik = GetPhik(kN)` inside the class
  (explicit convention), and nothing calls `GetPhik()`, `GetVar("phik")`, `GetVarArray("phik")` or
  `GetVarByGetter`; "phik" in a format string (A(e,e'pN) with `kN`) points directly at `fPhik`.
  **Fixed:** the `return` is added. Tested in `unit_kinematics`.
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

### Open leads (to investigate)

- [ ] **Rho transparency: t- and z-integrated results disagree (added 2026-11-23, from the author's
  memory of ~2012; details may be off).** The rho programs offer two ways to do the integrals, which
  should agree but gave results differing by orders of magnitude, while the Glauber/plane-wave ratios
  (the transparencies) were similar; the cause was never found. What the code and history show:
  - The two options are the `torz` flag of `progs/RhoTransparency/rhoset.cpp` (argument 1; `rhoset`
    is the program used for the actual calculations), which integrates either
    `RhoTCross::getCrosst` over ν and t (t from -0.4 to -0.1 GeV²) or `RhoTCross::getCrossz` over ν
    and z (0.9 to 1). `rhoset` can also do the integral with `rombergerN` (nested) or the 2D
    adaptive `numint` integrator, and `RhoDeuteron` (`rhodint`) has the same t/z pair.
  - The history (Sept.–Nov. 2012) records the problem: `374b867` "wrong mass I used in the fixed z
    integrations"; `58fd5ea` "finally figured out jacobians in RhoDeuteron, but still difference
    between z and t integrations, but goes like A....very puzzling"; `5fe7fae` "deuteron ... z and t
    integrations giving similar results and similar to the t integrations on large nuclei. still
    problem in the RhoTCross z integration, absolutely no idea now"; `67ae550` "corrected small error
    in t integrated jacobian". No later commit says the RhoTCross z problem was solved.
  - Lead (hypothesis, not checked): a difference that scales with A and cancels in the
    Glauber/plane-wave ratio points at a factor common to both, such as the Jacobian dt/dz or a
    flux/phase-space factor using the mass of the nucleus instead of the nucleon (or the reverse),
    rather than at the FSI.
  - Plan: run `rhoset` (and `rhodint` for the deuteron) with both options for one nucleus at a few
    kinematics, compare plane-wave cross sections first, check the t↔z Jacobian in `RhoTCross`
    analytically (and against `RhoDeuteron`, which was fixed), and check against a plane-wave case
    that can be integrated independently. Any fix changes values and needs the author's OK.

## 5. Low-risk modernization

- [x] **Constants:** `Tools/constants.hpp.in` defines `inline constexpr` variables in namespace
  `nucqcd` (also reachable as `nucqcd::constants::`), and the macros are aliases of them. Values are
  unchanged, so no reference moved. `PI` stays a literal, because `Integrator/numint/macros.hpp`
  defines the same macro; `HOMEDIR` stays a string literal (also `nucqcd::share_dir`). `M_W` is now a
  double instead of an int. The `MUP` doc said "anomalous" magnetic moment, but 2.79 is the full
  moment; fixed in the doc.
- [x] **Rho mass and degree conversions:** there were two rho masses, `MASSRHO` 775.49 (the PDG 2012
  value) and `M_RHO` 775.5 (unused). Both are now `nucqcd::mass_rho` = 775.26 MeV, the PDG 2022–2026
  average, so rho results move by about −0.03% (rho transparency and production in `FastParticle`,
  `RhoTCross`, `RhoDeuteron`). Earlier rho results used 775.49. No test covers the rho mass. Delete
  cached rho Glauber grids in `share/grids/`, since the file names do not record the mass.
  `DEGRTORAD`/`RADTODEGR` are now exactly π/180 and 180/π (they were accurate to about 1e-9); no
  reference moved.
- [x] **Ownership:** classes that own memory now either copy it correctly or cannot be copied. A
  build with the copy operations of all owning classes deleted showed that only `Cross` is copied
  (`vector<Cross>` in `progs/onenucl/JlabPAC16_calc`); `Cross` creates its `Model` per call and owns
  nothing in between, so that copy was already safe.
  - `std::unique_ptr` for the deuteron wave function of `DQEinclusive`, `InclusiveCross`,
    `InclusiveCrossRes` and `Deut_Conv_EMTLF` (movable, not copyable; their destructors are gone).
  - `std::vector<std::unique_ptr<…>>` for the grids of `PionTCross` and `RhoTCross`. These grids point
    to the member `nucleusthick`, so the classes are neither copyable nor movable (deleted copy
    operations). The unused `pdistgrid` member was removed.
  - `std::vector` for the arrays of `MeanFieldNucleus` (quantum numbers, F/G grids, Y_κ grids),
    `MeanFieldNucleusThick` (densities), `FsiCorrelator` (correlation grids) and `ClassGridThick_SCX`:
    rule of zero, so the hand-written copy constructors and destructors are gone and assignment,
    which was missing for `MeanFieldNucleus[Thick]`, now works. The getters of the 1D arrays still
    return `const` pointers; `getF`/`getG`/`getYkappa`/`getYminkappa` (unused outside the class) return
    `const` references to the vectors.
  - The Glauber grids (`AbstractFsiGrid` and everything derived from it, so also `Model`, which holds
    grids by value) are not copyable: copy operations deleted on `AbstractFsiGrid`. Their grid
    arrays (up to 5D `complex<double>` raw arrays, used in the interpolation and the binary cache
    I/O) stay raw pointers for now; converting them belongs with the grid-hierarchy refactor in
    step 6.
  Tested in `unit_ownership` (compile-time checks of which classes are copyable or movable, and
  copies and assignments of `MeanFieldNucleusThick` that outlive the original). Function-local
  `new`/`delete` pairs (scratch arrays in integrands) were left alone. Not covered: the Tools
  classes with hand-written copy operations (`TKinematics2to3`, the `T*PWF` wave functions,
  `GammaStructure`, ...), which were not checked line by line.
- [x] **Errors:** the libraries throw standard exceptions instead of printing to `std::cerr` and calling
  `exit(1)` (about 270 calls in Tools, Glauber, Knockout, GPD_Partons; DIS had none):
  `std::invalid_argument` for invalid options, names, format strings and quantum numbers (most of
  them, including the ~110 format-string checks of `TKinematics2to3`), `std::out_of_range` for values
  outside a grid or parametrization (r beyond the mean-field grid, Christy-Bosted above W² = 25 GeV²,
  matrix indices), `std::logic_error` for calls in the wrong order (grids read before they are filled
  or initialized twice), and `std::runtime_error` for input files that cannot be opened and failed
  calculations. The messages are unchanged; `Tools/ErrorMessage.hpp` (`nucqcd::ErrorMessage`) builds
  them with stream syntax. An uncaught exception still ends the program with the message (now via
  `std::terminate`, so the exit status is that of an abort rather than 1). Exceptions cannot pass
  through the C and Fortran integrators (Cuba, cubature, GSL), so inside integrands called from there
  they still end the program. Left as they were: third-party code (Cuba, cubature, MSTW (removed in step 10)), `TMPI`, and
  the `STOP`s in the Fortran code. Tested in `unit_errors`. Converted with a script and reviewed per
  message (CRLF line endings of `WeakQEHadronCurrent.cpp` kept).
- [x] **Grid cache:**
  - [x] **Cache directory:** `Tools/GridCache.hpp` provides `nucqcd::grid_cache_dir(share)`, which is
    `$NUCQCD_GRID_DIR` if that is set and `share/grids/` otherwise. `AbstractFsiGrid` (and so all
    Glauber grids), `DistMomDistrGrid` and `ClassGridThick_SCX` use it, so the cache can live on
    scratch space or per job while the input is still read from `share/`. The test `main()` unsets
    the variable. Grids for particles with user-set scattering parameters used to go to a directory
    `fixgrids` next to the cache, found by cutting 5 characters (`grids`) from the path; that
    directory was never created, so these grids were recomputed on every run. They now go to the
    cache directory itself; their file names (with `.beta`/`.eps` instead of `.p`) already differ
    from the others.
  - [x] **Atomic writes:** `nucqcd::write_file_atomically` writes a temporary file in the same
    directory and renames it, creating missing directories. A crash, a full disk or a concurrent job
    can no longer leave a truncated grid file. Used by all grid writers above and for the
    correlation grids of `FsiCorrelator` (which stay in `share/input/`, where precomputed ones are
    tracked in git). Readers also check the stream after reading, and recompute a grid whose file is
    incomplete (old files can still be truncated). Tested in `unit_grid_cache` and in a `slow` case of
    `reg_glauber_grid`. Not done for the deuteron GPD grids in `share/gpd_deutgrids/`
    (`GPD_Partons`), which are text files written point by point during the computation.
  - [x] **Cache key: grids at the bin centre** (option (a) of three, chosen by the author on
    2026-11-22; the others were exact keys, which recompute every new kinematics, and finer bins). The
    file name is the cache key and bins the kinematics: 10 MeV in p, 0.1 in cos θ and φ, 0.1 GeV² in
    the hard scale, 0.01 in the CT parameters. Until now the grid was computed at whatever kinematics
    came first in its bin and then reused (also within a run, by `updateGrids()`), so results
    depended on the contents of the cache and on the order of the calculations.
    `AbstractFsiGrid::addParticle` now stores the particle moved to the centre of its bin
    (`FastParticle::gridBinCentre`), so the grid in a file is the same whoever computed it. The
    bin widths are unchanged, but the bins are centred on the lattice (rounding instead of
    truncating p and the CT parameters), so round inputs such as p = 1500 MeV, Q² = 8 GeV² and
    lc_mod = 1 are bin centres and stay as they are. φ is taken in [-π, π] and the edge bin
    |φ| > 3.05 goes to π, so particles near the φ = π plane lie in it. When the binned cos θ is ±1, φ is set to 0:
    it is irrelevant along the z axis, and φ = 0 gives the in-plane symmetry and a single cache file. cos θ = 1
    stays exactly 1, so single particles near the z axis use the faster code for `getCosTheta() == 1.`. User-set scattering
    parameters are binned too (σ to 1 mb, β² and ε to 0.01), and their key now also contains the
    momentum, which the CT and decay grids depend on. File names end in `.key2.dat`, so grids cached
    with the old key are not reused (delete them).
    **Impact:** results change by the distance to the bin centre. The A(e,e'p) programs
    (`Knockout/Model`) put the proton along z, so only the momentum moves (by at most 5 MeV): the RMSGA
    cross sections in `reg_glauber_grid` moved by 4–6e-5 relative and were re-baselined. Particles
    off the z axis move further: cos θ bins of 0.1 put every θ below about 18° at θ = 0, so
    `reg_glauber_grid` now uses cos θ = 0.9 (a bin centre) instead of θ = 0.3 rad, and its grid
    references were re-baselined. Earlier results with off-axis particles used the angle of the first
    calculation in each bin. (The cos θ bins were replaced by θ bins in the next item.)
    Also fixed: the `FastParticle` copy constructor and assignment did not copy the decay width
    `Gamma`, which `setScatter` and the binning use. Tested in `unit_grid_bins` and in a `slow` case
    of `reg_glauber_grid` (two particles in the same bin give bitwise identical grids, computed or
    read from the cache).
  - [x] **Cache key: θ bins of 0.1 rad instead of cos θ bins of 0.1** (chosen by the author on
    2026-11-22, after a measurement). Bins of 0.1 in cos θ are narrow in angle near θ = 90° but wide
    near the z axis: the cos θ = 1 bin holds every θ below 18.2°, all computed at θ = 0. Bins of
    0.1 rad are never more than 0.05 rad (2.9°) from their centre. The integrals still run in cos θ;
    only the angle the grid is computed at and the key change. Below 0.05 rad the particle still goes
    on the z axis with cos θ exactly 1 (fast path, φ = 0); the edge bin θ > 3.05 goes to π. Keys
    now contain `.tht<θ/0.1>` and end in `.key3.dat` (the key2 files are not reused; delete them).
    **Measurement** (local scratch driver, not committed): distorted momentum distribution
    ρ_D(p_m) = Σ_m |∫d³r e^{-ip_m·r} ψ_m(r) F(r)|² for carbon s1/2 and p3/2, proton 1 GeV at
    θ = 5°–90°, p_m = 150 and 250 MeV parallel, perpendicular and out of plane. With cos θ bins the
    deviation from the exact angle reached 8% at θ = 15° (155% in the deep FSI minimum of the
    s-shell at 250 MeV perpendicular, 17% of the plane-wave value there); with θ bins 4% (86%,
    9.5%), and below 1% at most angles. **Impact:** only grids with particles off the z axis:
    `DoubleNModel` and `progs/2bodymom`. The single-nucleon knockout codes (`Model`, `Cross`,
    `RhoTCross`, `PionTCross`, `WeakQEHadronCurrent`, `DQEinclusive`) put the particle along z and are
    unchanged (the RMSGA references in `reg_glauber_grid` did not move). `reg_glauber_grid` now
    uses θ = 0.5 rad, and its grid references were re-baselined. Rotating each grid into the frame
    of its particle (no direction binning) is not straightforward for `DoubleNModel`: the
    knocked-out states (shell, m) that the grid divides out are quantized along the grid's z axis.
  - Optional: the file names do not record the physics (model version, constants), so stale grids must still
    be deleted by hand when the physics changes; the key version tag could be bumped for that.
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
- [x] **Doxygen gaps:** with `EXTRACT_ALL` off there are about 1400 warnings, about 1300 of them
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
  - [x] **Doc errors fixed** (2026-11-22): all warnings other than "not documented" are gone (1623 →
    1231 warnings, all undocumented members). Unclosed `\f$` formulas (`ClassGridThick_SCX`,
    `WeakQECross`, which also produced a bogus `#if` imbalance), stray backticks (`GlauberGridThick`),
    LaTeX outside formulas, `\param` names that matched no argument (`Deut_Conv_EMTLF`, `PionTCross`,
    `Poldeut`, `WeakQECross`), `\return` on `void` functions (`Cross`), unknown commands. Definitions in
    `.cpp` files that Doxygen could not match to their declaration because they wrote `string` or
    `complex<double>` (under `using namespace std`) where the header has `std::`: the signatures now
    use `std::`. In the Doxyfile, the ROOT macros `ClassImp` and `ClassDef` are predefined (`ClassDef`
    declares `Streamer`, as ROOT's does) and `Streamer` is excluded from the docs.
  - [x] **Undocumented members** (1231, done 2026-11-22; `make doc` now gives no warnings and no errors): Tools 1049 (`TKinematics2to3.h` 303, `TKinematics.h` 107,
    `TSpinor.h` 94, `TDeuteron.h` 66, `TInterpolatingWavefunction.h` 52, ...), Knockout 91, Glauber 58,
    DIS 33. One class or library per PR.
    - [x] DIS (`Resonance`, `He3Cross`, `InclusiveCrossRes`, `Deut_Conv_EMTLF`): no warnings left.
    - [x] Glauber: no warnings left.
    - [x] Knockout: no warnings left. Remaining (1049): all in Tools.
    - [x] `TDeuteron.h`, including the D-wave conventions (also stated in `Poldeut` and `Deut_Conv_EMTLF`).
    - [x] `TKinematics2to3.h` (303 members).
    - [x] `TKinematics.h` (107 members).
    - [x] `TSpinor.h` (94 members), with the spinor conventions.
    - [x] Deuteron wave function implementations: `TWavefunctionImplementation`, `TYukawaPWF`, `TGrossPWF`, `TCstPWF`, `TInterpolatingWavefunction` (150 members).
    - [x] `constants.hpp` (macro aliases); the third-party MSTW code (`mstwpdf.h/.cpp`) is excluded from the docs.
    - [x] Lorentz and Dirac algebra: `FourVector`, `Matrix`, `GammaStructure`, `TLorentzQuaternion` (132 members).
    - [x] The rest of Tools: `TElectronKinematics`, `TLeptonKinematics`, `TMPI`, `RecMosh`,
      `OldDeuteron` (same amplitudes as `TDeuteron`, checked; its radial functions absorb 1/sqrt(4 pi)
      and the D-wave sign, and their documented unit fm^3/2 was wrong: MeV^-3/2),
      `TKinematics2to2`, `TKinematics2to3WithLabAngles`, structure-function getters, and small helpers.
      Formulas in the docs can use `amsmath` (Doxyfile `EXTRA_PACKAGES`).
    - Still open, with the README (step 8): the Doxygen main page.

## 6. Structural refactors

Do these one at a time, running the full test suite (including `slow`) after each.

- [x] **Constants macros:** all uses of the macro aliases (about 1700, in 143 files: libraries,
  programs, tests) are replaced by the `nucqcd::` names, and the aliases are removed from
  `constants.hpp.in`, including `HOMEDIR` (now `nucqcd::share_dir`). Only `PI` is left, a literal,
  until `numint` drops its own. Done with a script that skips string literals and keeps line
  endings; `progs/GPD/diDVCS_gen_gridonly.cpp` keeps its own `#define ALPHA` (same value).
  `ClassGridThick_SCX` takes its share directory as `const char*` (it took `char*`, which only
  accepted the string literal `HOMEDIR` by a GCC extension). Values are unchanged: all tests,
  including `slow`, pass with the same references. `unit_constants` now fails to compile if an old
  alias comes back.
- [ ] **Update the constants to current PDG values** (split off from the item above): several values
  are older PDG values (e.g. the proton mass 938.27231 versus 938.27209 in PDG 2024); updating them
  changes results, so it needs the author's OK.
- [ ] **Code layout** (sharpened 2026-11-22 with the author). Headers next to their sources is fine and
  stays (no general `include/`/`src/` split: that pays off for libraries installed for outside
  projects, which nucqcd's are not, and moving every file would cost the `git blame` history). What
  does not follow good practice, in this order:
  1. ✅ **Split `Tools` (library `MePhys`)** (done 2026-11-23): the 36 headers plus Fortran are now
     four libraries, one per subdirectory (moved with `git mv`), each linking the ones before it:
     - `Tools/Core` (`libToolsCore`): `constants`, `ErrorMessage`, `FortranString`, `StringUtil`,
       `MathUtil`, `Utilfunctions`, `AuxFunction`, `GridCache`, `RunLog`, the Dirac/Lorentz algebra
       (`FourVector`, `Matrix`, `GammaStructure`, `TLorentzQuaternion`, `TSpinor`), `recmosh`, and
       `TMPI` (with MPI);
     - `Tools/Kinematics` (`libToolsKinematics`): `TKinematics`, `TKinematics2to2`,
       `TKinematics2to3(WithLabAngles)`, `TElectronKinematics`, `TLeptonKinematics`,
       `LightConeKin2to2/2to3`;
     - `Tools/WaveFunctions` (`libToolsWaveFunctions`, links Kinematics because `TDeuteron.h` uses
       `TKinematics2to3`): `TDeuteron` and its implementations, `OldDeuteron`, `He3wf` with the
       Fortran `getwavemod_tm.f`, `clebsch.f`, `spline.f` (only used by `getwavemod_tm.f`);
     - `Tools/StructureFunctions` (`libToolsStructureFunctions`, links Kinematics for
       `DeuteronStructure`, but not WaveFunctions): `NucleonStructure`, `NuclStructure`,
       `LHAPDFStructure`, `DeuteronStructure`, and the Fortran parametrizations.

     The libraries have a `Tools` prefix because `libCore.so` would clash with ROOT's `libCore.so`
     (the RPATH of the programs contains both `lib/` directories). `MePhys` is an INTERFACE target
     that links all four (`MePhysST` with `STATIC_BUILD`), so programs are unchanged. Glauber and
     Knockout link `ToolsWaveFunctions` (Glauber only for `TDeuteron::Wavefunction::ClebschGordan` in
     `TMFSpinor`; moving that function to Core would let Glauber link Core only), DIS and GPD_Partons
     link `ToolsStructureFunctions`, and each test links the smallest library it needs.
     `nucqcd_add_library` now links with `--no-undefined` (Linux), so a missing dependency in a `LINK`
     list is a link error; the whole tree (programs and PARTONS included) builds with it. The include
     lines are unchanged (every library directory is still a public include directory); they change
     in item 2. Not tested: `MPI_ENABLED` (no MPI on the development machine).
  2. ✅ **Prefixed includes, with the repository root as the include root** (done 2026-11-23).
     Every library directory used to be a public include directory, so code wrote
     `#include <Matrix.h>`, `<constants.hpp>`: generic names that could collide with headers of ROOT,
     PARTONS or LHAPDF, with the layering invisible. Now `nucqcd_add_library` exports only the
     repository root, and includes are paths from there (`#include <Tools/WaveFunctions/TDeuteron.hpp>`,
     `<Integrator/numint/numint.hpp>`). Scripted rewrite of 1075 `#include` lines in 290 files
     (libraries, programs including the attic, tests; line endings kept); headers local to a program
     keep their quoted includes. GPD_Partons was flattened (its `include/` and `src/` moved next to
     its `CMakeLists.txt`, as in the other libraries), so its headers are `<GPD_Partons/...>`; its
     unused `include/examples.h` (declarations of PARTONS example functions, never implemented or
     included) was removed. Found along the way: `numint.hpp` includes `<Cuba/cuba.h>`, which only
     worked because the Integrator directory (with the `Cuba` symlink) was on every include path;
     `Integrator` now links the `cuba` target publicly. `cubature.c` includes its own header as
     `"cubature.h"`. Doxygen shows the include paths from the repository root
     (`STRIP_FROM_INC_PATH`).
  3. ✅ **One header extension: `.hpp`** (done 2026-11-23, same PR): the 19 `.h` headers in Tools
     were renamed with `git mv`, and their `\file` tags and mentions updated. Only the third-party
     `Integrator/adaptive/cubature.h` and Cuba keep `.h`.
  4. **Namespace `nucqcd`.** 64 of the library headers declare their classes in the global
     namespace (`Matrix`, `Model`, `Cross`, `Resonance`, `FastParticle`, ...), the same collision
     risk at the symbol level. Add the namespace class by class, best together with ROOT phase 4,
     which touches the same classes; `using` declarations in a compatibility header keep the
     programs compiling in between.
- [ ] **Glauber grid hierarchy:** deduplicate it. `GlauberGridThick` and `GlauberDecayGridThick` are
  largely copy-pasted, as are the {plain, Thick} × {CT, Decay} abstract classes.
- [ ] **Speed of the Glauber grids (to investigate, added 2026-11-23).** From the author's experience,
  grids for particles that do not all move along the z axis are slow to compute (cached afterwards,
  but still). The `numint` integrators used there (`cube_adaptive`, `cube_romb`) are serial; only the
  explicit `numint::cuhre/vegas/...` wrappers use Cuba, which forks workers on one node
  (`CUBACORES`). First measure where the time goes, then compare: (a) filling the grid points in
  parallel threads (each point is an independent integral, but the integrand reads the point from
  members such as `r_hit`, so it has to be passed in explicitly first); (b) the vectorized
  `hcubature_v` (see the cubature item under step 7), which hands over batches of points that can be
  evaluated in threads; (c) Cuhre with `CUBACORES`. Best done together with the deduplication above.
  (b) and (c) change values within the precision, so they need the author's OK.
- [ ] **Remove the integrators `rombergerN` and `cube_romb`?** (added 2026-11-23). The author
  remembers the adaptive integrators doing better, especially for multidimensional integrands with
  non-trivial boundaries or cuts inside the integrand. **Benchmark this first; remove only if it
  holds.** Compare `rombergerN`, `cube_romb`, `cube_adaptive` (and Cuhre) on accuracy versus number of
  evaluations and time: (1) test integrands with known values in 2 and 3 dimensions: smooth ones, ones
  with a step or cut inside the domain, a non-rectangular domain mapped onto a box (as with the
  impact-parameter and z integrations in the Glauber phases), and a peaked one; (2) real integrands:
  a Glauber grid point for a particle along z and one off-axis, an A(e,e'p) amplitude from `Model`,
  and an `InclusiveCrossRes` integral, each against a high-precision reference. Record the results
  here (or in `docs/`), then decide. If they are removed: most uses are selected by the `integrator` flag of the grid
  and cross-section classes (0: `rombergerN` from `Tools/Core/Utilfunctions.hpp`, 1:
  `numint::cube_romb`, 2: adaptive `cube_adaptive`; the tests use 2). Those branches can go once the
  programs are checked for passing 0 or 1. Hard-coded uses have to be replaced, which changes values
  (author's OK per class): `rombergerN` in `InclusiveCrossRes` (6), `RhoDeuteron` (8), `RhoTCross`
  (4), `OneGlauberGrid` (4), `DistMomDistrGrid` (1) and the programs `rhodint`, `rhoset`;
  `cube_romb` in `Model.cpp:744` and `progs/momdistr/momdistr.cpp`. Then remove
  `Integrator/numint/cube_romb.cpp` (and `cube_tr`/`fubini`/`vfubini` if unused) and its tests in
  `unit_numint`, and simplify the `integrator` flag (an enum, see "Model selection" below).
- [x] **Better names for `Model` and `Cross`** (done 2026-11-23, in the include PR): `Knockout/Model`
  and `Knockout/Cross` are now `EMQEHadronCurrent` and `EMQECross` (files and classes, `git mv`),
  following `WeakQEHadronCurrent`/`WeakQECross` and `NucleonEMOperator`/`NucleonWeakOperator`; the
  class names were replaced in the 19 files that include them (library, programs including the
  attic, tests). Also fixed: the `\file` tag of `DQEinclusive.hpp` said `Cross.hpp`.
- [ ] **Model selection:** use enums or factory registries instead of model-name strings (about 700
  string comparisons). A typo should fail at compile time, not with `exit()` at runtime.
- [x] **Merge `NuclStructure` into `NucleonStructure`.** The parametrizations are now implemented
  once, in `NucleonStructure`; `NuclStructure` keeps its interface (about 20 programs,
  `DeuteronStructure`, `DeuteronCross`, `Poldeut` use it unchanged) and is a thin wrapper that stores
  the kinematics. `NuclStructure` takes x and W^2 independently (off-shell nucleons in
  `DeuteronStructure`), which `NucleonStructure` did not support: it now has `getF_xQW` (Christy-Bosted
  and SLAC use W^2, the others x; `getF_xQ` calls it with the on-shell W^2). The CB → SLAC switch above
  W^2 = 9 GeV^2 (and its effect on `getName`) is kept. Checked against the former stand-alone
  `NuclStructure` on 216 points (six parametrizations, proton and neutron, mostly off-shell): F2, F1
  and the name agree to 1e-13. **One value change (author's decision):** F1 for SLAC and CTEQ through
  `NuclStructure` was infinite (`1./0.`) and is now the `NucleonStructure` value
  F2(1+gamma^2)/(2x(1+R1998)); this only affects paths that gave inf/NaN (CTEQ in
  `DeuteronStructure`, programs asking for F1 with SLAC or CTEQ); `DeuteronStructure` with SLAC builds
  its own F1 and is unchanged. Also fixed: `NucleonStructure::operator=` leaked its MSTW object. About
  500 duplicated lines removed. Tested in `reg_nucleon_structure` (F1 now compared for all
  parametrizations, and a new off-shell case).
- [x] **Retire MPI (`TMPI`)** (done 2026-11-23, author's decision). `TMPI` (from strangecalc) only
  distributed independent kinematic points over ranks (rank-strided loop plus `GatherResults`); it
  did not parallelize integrals. Its use in the libraries was removed in 2014 (commit `24fffe1`: "I
  don't use it anyway"); only `rhod` and `rhottest` (progs/RhoTransparency) still used it, and
  neither was used for actual calculations. `TMPI` and those two programs moved to
  `progs/attic/mpi`, which is built only when CMake finds MPI (no option); the `MPI_ENABLED` option,
  the MPI links of `ToolsCore` and `Integrator`, and the Intel `-parallel` flag (automatic
  parallelization, oddly tied to `MPI_ENABLED`) were removed, as were the commented-out `TMPI`
  includes in `TKinematics2to3*.cpp`. The single-process programs `rhoset` (used for the actual
  calculations), `rhosingle` and `rhodint` do not use MPI; they are not one-to-one replacements
  (`rhod` scans an 11×101 grid in ν and t or z including coherent pieces, `rhodint` integrates over
  them; `rhottest` scans 32 points with `RhoTCross`), but the code is kept in the attic. **Bug found
  on the way:** `TMPI` no longer compiled with a current MPI, so `MPI_ENABLED=yes` was broken: it
  used `MPI_Type_struct` (removed in MPI-3.0, 2012) and called `dup`/`dup2`/`close` without
  `<unistd.h>`. Fixed in the attic copy (`MPI_Type_create_struct`, same arguments); it now compiles
  and links with Open MPI 5.0 (not run). `RunLog` keeps its rank detection from the launcher's
  environment (also covers `srun` without MPI).
- [ ] **Very long files and functions:** break them up, for example `Tools/Kinematics/TKinematics2to3.cpp`
  (6.8k lines, one function of about 1.3k lines) and `Glauber/NucleonEMOperator.cpp` (4.2k lines).

## 7. Unbundle Cuba (and cubature)

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

- [ ] **cubature (`Integrator/adaptive`, added 2026-11-23).** The bundled `cubature.c` is Steven G.
  Johnson's pre-1.0 `adapt_integrate` from October 2010, with a `minEval` argument that upstream does
  not have (probably added locally before the import; no convergence before `minEval` evaluations; the code passes 1e2 to 1e4). Upstream
  (github.com/stevengj/cubature) is now at 1.0.4: `adapt_integrate` became `hcubature`, plus
  `pcubature` (p-adaptive Clenshaw-Curtis), a choice of error norm for vector integrands, `size_t`
  limits, and a fix for a hang in `hcubature` for some integrands (1.0.4). Ubuntu packages it
  (`libcubature-dev` 1.0.4), but without `minEval`. Plan: fetch the pinned 1.0.4 release as for Cuba
  (CMake `FetchContent`/`ExternalProject`, SHA256) and keep the `minEval` change as a small patch
  file; port `numint2adapt.hpp` (`cube_adaptive`) to `hcubature` with `ERROR_INDIVIDUAL`, which
  should match the old convergence test. Run the regression and slow tests; any change in values
  needs to be understood and OK'd by the author. Also consider `pcubature` for smooth integrands,
  and `hcubature_v` for parallel evaluation (see "Speed of the Glauber grids", step 6).

## 8. Repository README

Write this last, after the refactors, so it describes the final layout and uses the new name.

- [ ] Cover:
  - what the code computes, and the libraries and how they depend on each other
  - dependencies and how to install them (ROOT, GSL, Cuba, and optionally LHAPDF with its PDF sets,
    PARTONS and MPI; what is still needed of ROOT after step 9)
  - the physics inputs that can be chosen by name: deuteron wave functions, nuclear mean fields, and the
    nucleon structure functions (`NucleonStructure`: CB, SLAC, Alekhin, CTEQ, HMRS, MSTW, `LHAPDF:<set>`,
    and combinations such as `CB+LHAPDF:CT18LO@4`, step 10)
  - how to build, test and run the main programs
  - the data in `share/`
  - which papers to cite
  - a link to `docs/HISTORY.md`

  CLAUDE.md and `manual/manual.tex` already contain parts of this.

- [ ] **License** (added 2026-11-23; the author's decision, best before the README). The repository
  has no LICENSE file; by default copyright law then gives others no right to use, copy or modify
  the code (this is a summary, not legal advice; the university may have rules too). Points to settle:
  - **What the dependencies allow.** GSL, PARTONS and LHAPDF are GPL-3, and the bundled
    `Integrator/adaptive/cubature.c` is GPL-2-or-later, so binaries that link them are covered by
    the GPL: nucqcd's own code can carry any GPL-3-compatible license (GPL-3.0-or-later is the
    simplest; a permissive license such as MIT or BSD-3 for nucqcd's own files also works, but the
    combined programs are still distributed under the GPL). Cuba is LGPL-3.
  - **Code derived from ROOT** (LGPL-2.1-or-later): `Tools/Core/Vector3`, `LorentzVector`,
    `Rotation` (phase 4a) follow ROOT's TVector3 & co. closely, and `nucqcd::solve_3x3` in
    `Tools/Core/MathUtil.cpp` copies `TDecompLU`'s algorithm (phase 2). LGPL-2.1+ code may be used
    in a GPL-3 project; keep ROOT's copyright notice and license reference in those files (to add),
    or rewrite them independently if a permissive license is chosen.
  - **Bundled code without a stated license**: `Integrator/numint` (Klaas Vantournhout, Ghent
    University and GSI), the QCD evolution code in `DIS/evolution`, and the Fortran
    parametrizations in `Tools/StructureFunctions` and `Tools/WaveFunctions` (Alekhin, Christy-Bosted,
    CTEQ, GRSV, HMRS, R1990/R1998, the He3 wave function code by Misak Sargsian). Check with the
    authors or their published terms whether they may be redistributed, and under what conditions;
    `docs/HISTORY.md` records where some of them came from.
  - **Then**: add `LICENSE` (and SPDX headers if wanted), a `THIRD_PARTY.md` (or a section in the
    README) listing the bundled code with its license and origin, and a "License" section in the
    README (step 8).

## 9. Decouple from ROOT where possible

Goal: the libraries (Tools/MePhys, Glauber, Knockout, DIS) build and run without ROOT; ROOT remains
an optional dependency, only for the programs that need it (fits with Minuit, ROOT files and
histograms). This removes a large dependency (ROOT 6 with mathmore and minuit2, built with its own C++
standard) from everything that only computes cross sections.

**Inventory** (2026-11-22, libraries; programs in parentheses):

| ROOT feature | Where | Uses |
|---|---|---|
| `TVector3` | Tools (kinematics, spinors, deuteron), Glauber, Knockout, DIS | ~390 (17 programs) |
| `TLorentzVector` | mostly the Tools kinematics and spinors | ~166 |
| `TLorentzRotation`, `TRotation`, `TVector2` | kinematics, spinor boosts, `DistMomDistrGrid`, `Rho/PionTCross` | ~85, ~31, ~6 |
| `TString` | the format-string parsing of `TKinematics*`, `TDeuteron` names | ~114 |
| `TObject`/`TNamed`, `ClassDef`/`ClassImp`, LinkDef dictionaries, `TRootIOCtor`, `Streamer` | 16 classes (`T*` in Tools, `TMFSpinor`) | ROOT I/O is never used in the libraries |
| `TMath` | constants and small functions | a few dozen |
| mathmore: `ROOT::Math::sph_legendre`, `ROOT::Math::Polynomial` | `TDeuteron` (spherical harmonics), `TKinematics2to3(WithLabAngles)` (quadratic roots) | 6 |
| `TMatrixT`, `TVectorT`, `TDecompLU` | `TYukawaPWF` (a linear solve) | 1 file |
| Minuit2 / `TF1` / fitting | programs only | ~15 programs (`F2_largex_fit`, `deepsdis` fits, `GPD/diDVCS*`) |
| `TFile`, `TTree`, histograms, `TCanvas` | programs only | ~5 programs |

**Phases** (each its own PR, full test suite including `slow` after each; values must not change):

1. [x] **Drop the ROOT dictionaries and base classes** (done 2026-11-22). No program writes these
   objects to ROOT files (the `TTree`s of `progs/GPD` hold plain numbers). Removed: `ClassDef`/
   `ClassImp`, the 16 LinkDef files and the `rootcling` step (`nucqcd_root_dictionary` in
   `cmake/NucQCD.cmake`), the `.pcm` files, the `TRootIOCtor` constructors, and the (already empty)
   `Streamer`s with their commented-out old bodies (~1000 lines in all). The `TObject` bases are gone;
   the three `TNamed` classes (`TKinematics`, `TKinematics2to3`, `TDeuteron::Wavefunction`) keep a name
   and title as `std::string` members with `GetName`/`GetTitle`/`SetName`/`SetTitle`, so the
   constructors that programs call are unchanged. `TKinematics::Print` takes `const char*` instead of
   ROOT's `Option_t*`. The derived wave functions constructed their base through the ROOT I/O
   constructor (`TWavefunctionImplementation(0)`); it now has a default constructor. The Doxyfile no
   longer needs the `ClassDef` predefinitions. All tests pass with unchanged references.
2. [x] **Small math** (done 2026-11-22). `TMath` → `<cmath>` (`std::sqrt`, `std::pow`, `std::copysign`
   for `Sign`, ...) and the `nucqcd::pi`/`nucqcd::rad_to_deg` constants (the same doubles); the
   unnormalized `TMath::Gaus` of the test wave function written out. `ROOT::Math::sph_legendre(l,m,θ)`
   → `gsl_sf_legendre_sphPlm(l,m,cos θ)`, as in `SphericalHarmonicCos`. The new `Tools/MathUtil.hpp`
   has `nucqcd::quadratic_real_roots` (the `Polynomial` roots in `TKinematics2to3(WithLabAngles)`
   are quadratics; ROOT solves them with `gsl_poly_complex_solve_quadratic` and keeps the roots with
   imaginary part 0, and so does this function, whereas `gsl_poly_solve_quadratic` differs in the last
   bit when b = 0) and `nucqcd::solve_3x3` for the `TYukawaPWF` boundary conditions, a copy of
   `TDecompLU`'s algorithm (Crout with implicit pivoting) rather than `gsl_linalg_LU`, which differs
   by up to 1e-9 relative. `MathUtil.cpp` is compiled with `-ffp-contract=off`: with the fused
   multiply-adds of `-march=native` the D-wave coefficients changed in the last bits, and w(r) near
   r = 0, where the terms cancel, by up to 2e-4 relative. Both functions were compared with ROOT on
   millions of random inputs without a difference, and the deuteron wave functions (all
   parametrizations, r and p space, off-shell parts), spherical harmonics and 2→3 kinematics were
   bitwise the same before and after; `unit_mathutil` keeps ROOT reference values. MathMore is no
   longer a required ROOT component. The programs still use `TMath` (phase 6).
3. [x] **`TString` → `std::string`** (done 2026-11-22). The format-string parsing of `TKinematics` and
   `TKinematics2to3`, the `GetVar`/`GetVarArray`/`GetPhysicalVarArray`/`GetVarByGetter`/`GetVarName`
   interfaces (now `const std::string&` and `std::string`; callers pass literals or `std::string`, so
   no program changed) and the wave-function names in `TDeuteron`. New `Tools/StringUtil.hpp`
   (`to_lower`, `iequals`, `icontains`, `split`) reproduces the TString behaviour: case-insensitive
   comparisons where `kIgnoreCase` was used (`TKinematics::GetVarArray` stays case sensitive, as
   before), and `split` leaves out empty tokens, as `TString::Tokenize` does (checked with ROOT).
   Also fixed: `TKinematics2to3::ReadFormatString` copied each token into a `char[11]` with `strcpy`
   (overflow for tokens longer than 10 characters), and the token lists leaked when a format string
   was rejected. The dead, commented-out `HashName` functions are removed. Tests: `unit_stringutil`,
   and mixed-case and empty-token format strings in `unit_kinematics`.
4. [ ] **Three- and four-vectors and Lorentz transformations.** In two parts:
   - [x] **4a: the types and their tests** (done 2026-11-23). `Tools/Core` has `nucqcd::Vector2`,
     `Vector3` (`Vector3.hpp`), `LorentzVector` (`LorentzVector.hpp`), `Rotation` and
     `LorentzRotation` (`Rotation.hpp`), with the subset of ROOT's interface that nucqcd uses (an
     inventory of the method calls in libraries and programs) under the same names, so that the
     switch is a rename. The formulas and their order of operations are ROOT 6.40's
     (`math/physics`, checked against the installed headers), and so are the details that matter
     for the results: `TMath::ATan2` (0 or ±π/2 for x = 0, where `std::atan2` can give ±π),
     `TVector2::Phi` in [0, 2π) versus `TVector3::Phi` in [-π, π], the implicit conversions
     `TRotation` → `TLorentzRotation` and `TVector3` → `TLorentzRotation` (a boost), active
     transformations with `a.Transform(b)` = b·a. What ROOT has inline in its headers is inline here
     too; what ROOT compiles into `libPhysics` is in `Vector3.cpp`, `LorentzVector.cpp` and
     `Rotation.cpp`, compiled with `-ffp-contract=off`, because ROOT's `libPhysics` contains no fused
     multiply-adds (checked with `objdump`). Differences: no `TObject` base, a bad index throws
     `std::out_of_range`, no warnings printed for zero axes or vectors (the behaviour is ROOT's), and
     `SetBoost` is private as in ROOT. `tests/unit/unit_vectors.cpp` compares every method with ROOT
     on 20000 random cases per class (zeros of both signs, ±1, magnitudes 1e-12 to 1e12, angles
     including multiples of π/2 and up to ±100, velocities up to and beyond 1, random chains of
     rotations and boosts) and requires bitwise equality (several million comparisons, all equal; a
     deliberately wrong `atan2` is caught). The test file is compiled with `-ffp-contract=off`:
     with `-march=native` the compiler chooses fused multiply-adds per call site, so two identical
     inline functions can differ in the last bits (1–4 ulp were seen for `VectorMultiplication`).
     That also means the converted code cannot be guaranteed bitwise equal at every call site; the
     regression references are the check there (expected bitwise or within ~1e-14). The test needs
     ROOT and stays as long as ROOT is available.
   - [ ] **4b: convert the code**, library by library (Tools kinematics and spinors, then Glauber,
     Knockout, DIS, GPD_Partons; then the programs in phase 6), with the regression and slow tests
     after each. `FourVector.hpp`'s `ToFourVector`/`ToLorentzVector` switch to `nucqcd::LorentzVector`
     with the Tools conversion. Add `namespace nucqcd` to the classes touched (code layout item 4)
     where it does not double the size of the change.
5. [ ] **Build system.** `find_package(ROOT)` becomes optional; the `root_libs` interface target is
   linked only by programs that need it; the libraries and tests no longer link ROOT. The
   `ROOTConfig`/`thisroot.sh` paragraph in CLAUDE.md and the README change accordingly.
6. [ ] **Programs.** Fits: keep ROOT for them, or switch to standalone Minuit2 (available as its own
   CMake project) if a ROOT-free build of those programs is wanted. Programs that write ROOT files or
   histograms keep ROOT (or move to text output, or to the attic if unused). Mark in `progs/` which
   programs need ROOT, and build them only when it is found.

Phases 1 and 2 are independent and small; 3 and 4 can go in either order; 5 and 6 come last.

## 10. LHAPDF interface for the nucleon structure functions

Add LHAPDF6 sets as a parametrization in `NucleonStructure` (and so in `NuclStructure`,
`DeuteronStructure` and the programs that take a structure-function name), the C++ counterpart of
`LHAPDFStructure` in `taggedpd/nucleon.py` of the poldeut repository.

- [x] **Physics** (done 2026-11-22, `Tools/LHAPDFStructure`): leading-order parton model,
  F₂ = Σ_q e_q² x(q + q̄) (nf = 5 by default) from an unpolarized set and g₁ = ½ Σ_q e_q² (Δq + Δq̄)
  (nf = 3 by default) from a polarized set; the neutron by isospin (u <-> d). F1 = F2 (1 + 4M²x²/Q²) /
  (2x(1+R)) with R1998 by default, as for MSTW and CTEQ in `NucleonStructure` (so that
  `LHAPDF:MSTW2008lo68cl` and `MSTW` can be compared); a constant R and no target-mass factor
  (poldeut's default) are options (`nucqcd::LHAPDFOptions`). g1+g2 from the Wandzura-Wilczek relation,
  tabulated on (x, Q²) nodes the first time it is asked for (as in poldeut's `ww_grid`). 0 outside
  0 < x < 1; outside a set's range LHAPDF's own extrapolation is used (documented).
- [x] **Interface:** `NucleonStructure("LHAPDF:CT18LO")`, `"LHAPDF:CT18LO:NNPDFpol11_100"`,
  `"LHAPDF::NNPDFpol11_100"` (polarized only), members as `"CT18LO/3"`; or
  `NucleonStructure(nucqcd::LHAPDFOptions)` for the options. `getF_xQ`/`getF_xQW` (x is used) and the
  new `getG1`/`getG1plusg2` members (they throw for parametrizations without a polarized LHAPDF set;
  GRSV stays in the static `getG1_grsv2000`). Loading a set is slow and `DeuteronStructure` makes a
  `NuclStructure` per kinematic point, so `LHAPDFStructure::get` keeps one object per set of options
  for the whole program, shared by all `NucleonStructure`s.
- [x] **Build:** optional (`LHAPDF_ENABLED`, ON by default: used when found), found through LHAPDF's
  pkg-config file, also in `~/.local`. Without LHAPDF the library builds, `LHAPDFStructure::available()`
  is false and the "LHAPDF:" names throw. `scripts/install_partons.sh` installs the sets CT18LO and
  NNPDFpol11_100.
- [x] **Tests** (`unit_lhapdf`, as poldeut's `tests/test_lhapdf.py`): parton-model sums and isospin
  against direct LHAPDF calls, F1 from F2 and R (both options), zero outside 0 < x < 1, one set only,
  malformed names, copies, the WW table against a direct integration (GSL QAG) within poldeut's
  tolerance, a rough comparison with SLAC F2 and GRSV g1. Skipped without LHAPDF or the sets
  (checked in a build with `-DLHAPDF_ENABLED=OFF`).
- [x] **Explicit combinations of parametrizations** (done 2026-11-22). `NucleonStructure("<low>+<high>")`
  or `"<low>+<high>@<W²>"` (W² boundary in GeV², default 9), or the constructor
  `NucleonStructure(low, high, Wsq_switch)`: low for W² <= boundary, high above, e.g. `"CB+SLAC"` or
  `"CB+LHAPDF:CT18LO@4"`. `NuclStructure("CB")` is now the explicit combination `"CB+SLAC@9"` instead
  of a rename in its constructor: the values are unchanged (all regression references pass), but
  `getName()` now returns "CB" also above 9 GeV² (until 2026-11 "SLAC"). `DeuteronStructure` asks
  the new `isF2Fit()` (true where SLAC is used, so F1 is built from F2 with R = 0.18 and the off-shell
  kinematics) instead of comparing names. `NucleonStructure("CB")` itself still means CB up to
  25 GeV² and throws above.
- [x] **MSTW through LHAPDF** (done 2026-11-22). The bundled MSTW grid code (`Tools/mstwpdf.h/.cpp`)
  and grids (`share/mstw_grids`, LO, NLO and NNLO; only LO was used) are removed. `"MSTW"` in
  `NucleonStructure` (and so `NuclStructure`, `DeuteronStructure`) is now an alias of
  `"LHAPDF:MSTW2008lo68cl"` (central member), with the same F1 convention; `GPD_T_Nucl_grid` (the
  forward limit of the transversity GPDs, GK model) takes u, d, ū, d̄ from the same set, and also
  accepts `"LHAPDF:<set>"`. F2 compared on 14 x values in [1e-4, 0.95] and 9 Q² values in
  [1, 1e4] GeV²: equal at the grid nodes (1e-16), up to 1.5e-3 relative between them (largest at
  x = 1e-4, Q² = 2 GeV², just above the charm threshold), ~0.1% below Q² = 1 GeV²
  (extrapolation): both interpolate the same grid, differently. The MSTW references in
  `reg_nucleon_structure` were re-baselined (changes up to 1.2e-4 at Q² = 2 GeV², ~1e-9 at 5 GeV²);
  `unit_lhapdf` checks two old values at a grid node. "MSTW" and `PARTONS_ENABLED` now need LHAPDF
  (CMake stops for PARTONS without it); the install script installs the set.

## 11. Run bookkeeping (provenance) for all programs

Some programs call `Bookkeep(argc, argv, arg_names)` (`Tools/Utilfunctions`) first thing in `main`, to
record how a result was produced. Make this better, then use it in every program.

**Current state** (2026-11-22): 8 of the 77 programs outside `progs/attic` call it (`poldeut`,
`b1_calc`, `pionT_calc12`, `Ca_dens`, `PxPz_exp_sample`, ...). It prints the date and time, each
argument next to a name from a string array, and the output of `git rev-parse --short HEAD`, run with
`popen` in `nucqcd::share_dir` at run time. Shortcomings:
- The git hash is that of the checkout when the program runs, not of the build. A binary built from
  an older commit, or with uncommitted changes, reports the wrong code; "dirty" is never recorded; it
  fails silently without git or when the share directory is not a checkout.
- `arg_names` must have `argc` entries: with fewer (a program called with extra arguments) it reads
  past the array (undefined behaviour).
- The lines go to stdout without a comment marker, mixed with the data that most programs print to
  stdout, so readers (gnuplot, numpy) have to skip them by hand.
- Not recorded: host, working directory, build type and compiler flags, versions of ROOT, GSL, LHAPDF
  (and the PDF sets used), relevant environment variables (`NUCQCD_GRID_DIR`), number of MPI ranks,
  run time and exit status.

1. [x] **Design and implementation** (done 2026-11-22). Decided with the author: the log goes to
   **stdout** (where most programs write their data, so it stays in the data file), every line starts
   with **`# `** (skipped by gnuplot, numpy.loadtxt, pandas with comment='#'; some programs' own headers
   use `// `, to be aligned in the roll-out), and **argument names are an optional list** checked
   against argc (an argument parser would be a separate item). Implemented as `nucqcd::RunLog`
   (`Tools/RunLog`):
   `nucqcd::RunLog runlog(argc, argv, {"wave function", "Q2 [GeV^2]"});` as the first statement of
   `main` writes the quoted command line (pastable into a shell), the named arguments (extra arguments
   as "argument i", missing ones as "(missing)"), date and time with time zone, host, working directory,
   `NUCQCD_GRID_DIR` if set, and the build; at the end of `main` or at `exit()` (an atexit handler) the
   end time and wall time. Only MPI rank 0 writes (rank from the launcher's environment, so no MPI
   dependency in libMePhys). Build information (`nucqcd::build_info`): `git describe --always --dirty`
   and the commit date, written into the build tree at every build by `cmake/BuildInfo.cmake` (the
   file only changes with the git state, so nothing is recompiled otherwise), plus the build type,
   compiler, flags and the ROOT, GSL and LHAPDF versions from the configuration. Tests: `unit_runlog`
   (marker, quoting, names, end written once, at `exit()` in a child process, MPI rank, build info).
2. [ ] **Roll-out.** Done (2026-11-22): `Bookkeep` replaced by `RunLog` in the 8 programs that used it
   (`poldeut`, `EMT_deut`, `b1_calc`, `pionT_calc12`, `Ca_dens`, `Nucl_dens_profile`,
   `PxPz_exp_sample`, `conv_V_test`; their "exec name" entries dropped, since `RunLog`'s names start at
   argv[1]; `poldeut` and `EMT_deut` had too few names and read past the array), and `Bookkeep`
   removed from `Utilfunctions`. To do: add `RunLog` to all other programs (not `progs/attic`). Check per program that its data output still parses
   (commented header lines); change the programs' own `// ` header lines to `# `. Can be combined with
   ROOT phase 6, which touches the programs too.
