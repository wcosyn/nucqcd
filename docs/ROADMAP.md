# Roadmap: cleanup and modernization of nucqcd

This roadmap follows from an assessment of the repository done in November 2026 (added in #3; status
updates also in #40 and #43, otherwise with the work; `docs/HISTORY.md` came with #6). The steps are
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
| 4 | Bug fixes | ✅ done: listed bugs (PRs #8, #9, #12); found later: stale FSI+CT grid (#18), kinematics (#19, #20, #33), DoubleNModel/Pair (#21, then moved to the attic in #22), π⁻ Glauber parameters (#28), TDeuteron errors (#31), He3 spin sum (#47) and factor A in `He3Cross` (#48), `TMPI` not compiling with MPI-3 (fixed in the attic with step 6, "Retire MPI", #54), `RunLog` reading freed memory at `exit()` (#57), `cube_romb` overflow beyond 10 levels in 3D (#60). New bugs found later are added here |
| 5 | Low-risk modernization | ✅ done: `-fopenmp`, Fortran interop (#11), constants (#13, #14), grid cache (#15–#17), ownership (#23), exceptions (#24), Doxygen (#10, #25–#30, #32, #34–#39: no warnings left) |
| 6 | Structural refactors | 🔶 done: constants macros (#41), `NuclStructure` merge (#42), code layout (`Tools` split #52, prefixed includes and `.hpp` #53, namespace `nucqcd` #65, the `T*` classes renamed #76, the `TNamed` leftovers #78), `Model`/`Cross` renamed (#53), MPI retired (#54), integrator benchmark (#61), precision reporting (#62), integrator option 0 removed (#69), CRLF to LF (#66), long files split (#70), `GetVar*` table (#71), model names as enumerations (#72), the long `TKinematics2to3` functions split (#74) and the crash they showed fixed (#75), the thick Glauber grids deduplicated (#77); `rombergerN` removed from the libraries (#82); copy-paste comments and a duplicate header guard (#86); the integration method as an enumeration with one dispatch function (#87, #88); the Glauber abstract classes reduced from seven to two (#92); the raw grid arrays replaced by `Grid3D` and `std::vector` (#93, #94); the long Glauber files split (#96); the off-axis thick grids in the frame of the particle (#97, 20 to 100 times faster, values change by up to a few 1e-3); to do: parallel grid points, the hard-coded paths, optional other integration algorithms |
| 7 | Unbundle Cuba and cubature | ✅ done: Cuba (#67), cubature (cubature 1.0.4 downloaded and patched; results of integrands with steps change at the 1e-6 to 1e-2 level, `-DCUBATURE_LEGACY_SPLIT=yes` reproduces the old ones) |
| 8 | Repository README | ⬜ (last, after 9, 10 and 11) |
| 9 | Decouple from ROOT where possible | ✅ done: phases 1–5 (#44, #45, #50, #55, #56, #58), phase 6: GPD generators (#59) and the fits on the standalone Minuit2 (#73). ROOT is only used by the test `unit_vectors` |
| 10 | LHAPDF interface for the nucleon structure functions | ✅ done (#46), including MSTW through LHAPDF |
| 11 | Run bookkeeping (provenance) for all programs | ✅ done: design and `RunLog` (#49), `Bookkeep` replaced in its 8 programs (#51), `RunLog` in all other non-attic programs (#64) |

## Where we are / next up

*Updated 2026-11-26, night (after PR #97; no PR is open).*

**Done** (PRs #44-#97, details in the steps below): ROOT decoupling (ROOT is now only used by the test
`unit_vectors`; the fits run on the standalone Minuit2, #73); LHAPDF structure functions (step 10); He3 fixes; run
bookkeeping in all programs (step 11, #49, #51, #64); Cuba unbundled (#67); the code layout (`Tools` split, includes
from the repository root, `.hpp`, namespace `nucqcd`, the `T*` classes renamed, the `TNamed` leftovers removed:
#52, #53, #65, #76, #78); MPI retired (#54); the integrator benchmark (#61), `cube_romb` overflow (#60), precision
reporting (#62), the removal of the dead integrator option 0 (#69) and of `rombergerN` from the libraries (#80-#82:
dead chains removed, `RhoDeuteron` and `InclusiveCrossRes` archived, the plane-wave grid on `cube_adaptive`, the
routine itself in `progs/attic/rombergerN.hpp`); LF line endings (#66); the two longest files and two long functions
split (#70, #74) with a crash fix (#75); the `GetVar*` table (#71); model names as enumerations (#72); the safety net
for all Glauber grid variants (#68) and the deduplication of the two thick grids (#77). Since then (2026-11-26): the
`InclusiveCrossRes` "FSI returns 0" lead investigated, not a bug (#84); the rho t- versus z-integral bug found and fixed,
with the derivation in `docs/rho_jacobian.md` (#85); copy-paste comments and a duplicate header guard cleaned up (#86);
the integration method as an enumeration with one dispatch function, `numint::cube` (#87, #88); the W > 2 GeV cut checked
on the input of the rho programs (#90); a note on the vintage of the physical constants, which are not updated (#91); the
seven abstract Glauber grid classes reduced to two (`GridVariant`, #92); the raw pointer arrays of the Glauber grids,
`DistMomDistrGrid`, `FsiCorrelator` and `setPionGlauberData` replaced by `Grid3D` and `std::vector` (#93, #94); the four longest Glauber files split (#96); the thick Glauber grids off the z axis integrated in the frame
of the particle, 20 to 100 times faster (#97).

**Options for the next step** (my recommended order; the author decides):
1. **Open leads** (step 4, "Open leads (to investigate)"):
   - ~~the rho transparency t- versus z-integrals~~ (2026-11-26: the z version had the residual nucleus mass instead of
     the Jacobian |dt/dz|/(2 nu); fixed in #85, checked with the Glauber distributions for helium, z/t within 2.2%; a carbon
     or finer-grid check and the effect on the z-version transparencies are open; the W > 2 GeV cut is checked on the input of
     `rhoset` and `rhosingle`, see the item);
   - ~~the `InclusiveCrossRes` FSI part that returns exactly 0~~ (investigated 2026-11-26: not a bug, the
     settings of `increstest` cancel by construction; with `fixprop` 1, `symm` 0 the FSI is ~6% of the plane wave;
     reviving the class is still open, see the item);
   - ~~the **PARTONS setup** on the second computer~~ (2026-11-27: the VGG99 PDF set was never given to PARTONS v5, fixed; the
     NaN at the corner of the nucleon grid fixed, see the item).
2. **Glauber grids:** the off-axis speed is done for the thick grids (#97: the integrals in the frame of the particle; values
   changed by up to a few 1e-3, towards the converged value). Open: filling grid points in parallel (needs the hit point out of the
   members), the thin `GlauberGrid` (not converted, see the item), the missing SRC normalization term on the z axis (see the item).
   The hierarchy dedup is done (#77, #92, #93, #94);
   only the named `printFsi_*_grid` functions are left, low priority.
3. **Smaller items:** other integration algorithms in `numint::cube` (step 6, optional part 3), the hard-coded absolute paths in some programs (step 6, copy-paste leftovers), the longest files of the libraries other than `Glauber/` (see the item on long files, step 6),
   `SolvePhotonAndParticle1` (134 lines), reporting the Suave bug to the Cuba authors, a format string that gives
   no scenario at all is still accepted (step 6, long functions), `Resonance` stays in the DIS library although
   only the archived class used it.
4. **Last:** the license decision (step 8), then the README (currently a stub; step 8 lists what it must cover).

**Open choices for the author:** the license (step 8, see the item there); which of the open leads to take up first; the README (step 8)
comes last.

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

- [x] **`RunLog`: invalid memory access at `exit()`** (#57; found 2026-11-23 while making ROOT optional,
  step 9 phase 5). The `atexit` handler that writes the end of a `RunLog` when a program calls
  `exit()` iterated a function-local static `std::set` that was already destroyed: the handler was
  registered before the set was first constructed, so the set's destructor ran first. Undefined
  behaviour in every program that calls `exit()` with a live `RunLog` (valgrind: invalid reads in
  `finish_at_exit`); with ROOT linked it happened to work, without ROOT `unit_runlog` crashed with
  SIGSEGV. Fixed: the registry and its mutex are a never-destroyed object. No values affected (only
  the "# finished" line of the run log).

- [x] **`numint::cube_romb`: 32-bit overflow beyond 10 refinement levels in 3D** (#60; found 2026-11-23
  with the integrator benchmark, step 6). `_cube_romb_step` counted the (2n+1)^N grid points of a
  refinement level in an `unsigned`, which overflows from level 11 on in 3D (2049³ > 2³²; from level
  8 in 4D), so the loop summed the wrong points: the 3D step function (unit sphere in the unit cube)
  at requested precision 1e-7 gave 0.0087 instead of π/6 (error flag `EMAXITER`). Fixed: 64-bit
  count and index. Below the overflow the results are bitwise unchanged (same case at 1e-5:
  identical value); with the fix the 1e-7 case gives π/6 to 1.7e-7 (still `EMAXITER`: 12 levels do not
  reach 1e-7 on a step). Affects only calls that need more than 10 levels in 3D (more than 4·10⁹
  evaluations per level), i.e. very high precision or hard integrands; such calls ran for minutes.
  The evaluation counter that `cube_romb` returns (`unsigned &neval`) still wraps (it is part of the
  interface; only informational).

### Open leads (to investigate)

- [x] **PARTONS setup on the author's second computer (added 2026-11-23, fixed 2026-11-27).** The GPD programs stopped
  with "Info file not found for PDF set 'UNDEFINED'". Two causes, neither in the template of
  `partons.properties` (PARTONS v5 reads no PDF set from it; the template starts PARTONS fine):
  - since 2020 PARTONS has no built-in MSTW grid: `GPDVGG99` takes its LHAPDF set as the module parameter `setName`
    (default `UNDEFINED`), which no program set. All nine places that create a VGG99 module now call
    `setSetName("MSTW2008nlo68cl")`, the set of the old built-in grid (`mstw2008nlo.00.dat`);
  - `GPDMMS13` hard-codes `MSTW2008nlo68cl`, which the install script did not fetch. It does now (next to the
    sets of #46).
  The script also warns when an existing `partons.properties` has v4 keys (`environment.configuration.file.path`,
  `grid.directory`), but still never overwrites it. Cached VGG99 grids in `share/gpd_deutgrids/` were made with the
  old built-in set; delete them to recompute.
  **The NaN ("Unphysical x given: nan"), found 2026-11-27:** it came from the degenerate corner |x| = xi = 1 of the
  nucleon grid in `GPD_V_Nucl_grid`, where the integration limits of VGG99 and MMS13 are 0/0. The old PARTONS passed
  the NaN on and the grid zeroed it (`isnan` check); LHAPDF 6 throws instead and aborted the whole grid. The corner is
  now skipped (set to 0, as before). No other grid point throws. Note that `conv_V_test` does not build the nucleon grid
  at all when xi exceeds the deuteron limit for the given t, so such runs are trivially zero (e.g. xi 0.2, -t 0.5).
  `GPDMMS13` provides only H and E in PARTONS, so `conv_V_test` with MMS13 stops at the axial (Ht, Et) grid: use MMS13
  only for the vector part.
- [x] **`InclusiveCrossRes`: the FSI part returns exactly 0: not a bug (added 2026-11-25, investigated 2026-11-26).** The resonance model for inclusive
  DIS on the deuteron (W. Cosyn, W. Melnitchouk and M. Sargsian, Phys. Rev. C 89, 014612 (2014)) was archived in
  #82 (`progs/attic/threeres`, not built). It compiles once the stale constructor call
  `structure(elec, proton, strucname)` is changed to `structure(proton, strucname)` (`DeuteronStructure` lost its
  `ElectronKinematics` argument). The plane-wave part works (`calc_F2Dinc` 0.30 for the proton at x = 0.3,
  Q^2 = 5 GeV^2, Paris, Alekhin, resonances 1.4 and 2.4 GeV with coefficients +-1/sqrt(2), beam 5 GeV).
  **Finding:** the four resonance-pair integrals `(res1, res2)` of `calc_F2DincFSI` are *not* zero (about 6.7e-11 each
  at x = 0.3, before the factor chi2 and the prefactors) but have the signs +, -, -, +, because the pair weight is
  `c_1 c_2*` with `c_1 = -c_2 = 1/sqrt(2)`. The sum vanishes exactly whenever the amplitude of a pair does not depend on
  the two resonances apart from that coefficient. In the code the resonance masses enter the FSI amplitude only through
  (a) `sigma(Q^2)`, via `sigma_slope` (the example command has slope 0), (b) the shift of the spectator
  longitudinal momentum for `fixprop` = 1 with `symm` = 0 or -1, and (c) the `kin_choice` cut. With `fixprop` = 0 the
  integrands use the integrated `Wxprime2`, not the resonance masses. Checks at x = 0.3 (proton, same settings as above, off-shell type 3):
  sigma slope 0, `symm` 1 or 0, `fixprop` 0 or 1 with `symm` = 1: exactly 0; sigma slope 1 GeV, `fixprop` 0:
  -4.3e-6 (a 1e-5 residue of the cancellation); **`fixprop` 1, `symm` 0, slope 0: 0.0181 and 0.0183 (the two outputs), about 6% of the
  plane-wave 0.30**, a physically plausible size. So `increstest 1 3 0 5. paris Alekhin 40. 0. 1.4 2.4` (`symm` 1, `fixprop` 0,
  slope 0) is a configuration where the model gives zero by construction, not a broken integral (the 280 s without an output line
  was the slow first point at x = 0.025). The paper's model avoids the cancellation through the mass-dependent resonance
  parameters and the `fixprop`/`symm` prescriptions; I could not confirm which settings the paper's curves used
  (**to check by the author**). Minor findings in the class: `abs(structfactor1<1E-09)` has a misplaced bracket (works by
  accident: abs of a bool), `calc_F2Dinc` prints "blaaa" on stdout, the convergence loop of `get_prz` runs once so it never
  converges (`otherWx2` is then the first guess, intended or not is unknown). **Still to do if the class is revived:**
  (1) compare with the curves of the paper using the settings above, (2) decide the integrator (the nested `rombergerN`
  of `progs/attic/rombergerN.hpp`, or `cube_adaptive`, see the item on `rombergerN`), (3) write a regression test (use
  `fixprop` 1, `symm` 0) before moving it back into the DIS library. Not urgent: the class is not used by anything.
- [x] **Rho transparency: t- and z-integrated results disagree: cause found and fixed (added 2026-11-23, fixed 2026-11-26, from the author's
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
  - **Cause (found 2026-11-26):** the z-integrated `getCrossz` multiplied by `massX`, the invariant mass of the
    whole residual nucleus system (about M_A), in `RhoTCross::getfrontfactor`. The Jacobian between the two ways of
    eliminating the on-shell delta function of the final nucleon (at fixed missing momentum) is |dt/dz| = 2 nu K, with
    K = `massi` = M_A - M_(A-1) - E_exc, the off-shell energy of the struck nucleon (0.922 GeV for carbon, not the
    on-shell m_N), for pm = 0 and a pm- and angle-dependent mass otherwise (K ~ 0.7-1.7 GeV in a random sample), not
    M_A. The derivation is in `docs/rho_jacobian.md`. The z integrals were therefore a factor ~M_A/m_N ~ A too large, which is the "goes like A" of the 2012 commit
    messages, and it cancels in the Glauber/plane-wave ratios. The archived `RhoDeuteron` (where the Jacobians were
    "finally figured out") divides the t version by the off-shell nucleon energy `Eoff`, a nucleon-sized quantity and not
    M_A, which points the same way (not checked in detail).
  - **Which results were affected:** only the z branch of `getfrontfactor` changed, so the t-integrated results, and
    with them the published transparencies (the last calculations integrated over t, per the README of the paper's
    directory), are unchanged; the absolute z-integrated cross sections were too large by ~M_A/massi. The
    transparencies of the z version barely noticed (a ratio of RMSGA to PWA over the same range, and `massX` was nearly
    constant); how much they move with the corrected factor was not quantified.
  - **Check:** a Monte Carlo that copies the kinematics, cuts and prefactors of `getCrosst`/`getCrossz` with a Gaussian
    in place of the Glauber momentum distributions (carbon-like M_A = 11.18 GeV, Q^2 = 1.5 GeV^2, nu = 2.5-3.5 GeV,
    t in -0.4..-0.1, z in 0.9..1): z/t = 11.66 before; with the exact Jacobian 0.997 (statistical error 0.05%), and
    0.988 at Q^2 = 2.5, nu = 3.0-3.8. The remaining 0.3-1.2% was not traced (edge singularities of the root
    solutions in the Monte Carlo, or a small leftover; **not checked**). With massX removed the ratio was 1.03, 1.02,
    0.98 for three kinematic settings, so a constant mass is not the fix. For a point nucleus (pm -> 0) the author's
    commented analytic pair in `rhoset.cpp` agrees to 4 digits (0.9996), and its z/t ratio 2 nu M_A is this Jacobian.
  - **Fix:** `getfrontfactor` takes C0, Cz and Cx and returns `|dt/dz|/(2 nu)` times the elementary factor for z,
    `q |C0 px - Cz nu px/q - Cx nu (z - pz/q)| / |Cz px - Cx pz|` (derived in the code comment). The t version is
    unchanged. **Changes values:** every z-integrated rho result (`rhoset` with torz 0) gets smaller by roughly
    M_A/m_N (about 12 for carbon); the t results and the transparencies do not move. No test covers `RhoTCross`.
  - **Normalisation (checked against the paper, 2026-11-26, `docs/rho_jacobian.md` section 4):** the prefactor and front
    factor of `getCrosst` are Eq. (13) of Cosyn and Ryckebusch, Phys. Rev. C 87, 064608 (2013), arXiv:1301.1904 (formulas
    checked against the author's TeX source), so the t version is the reference; the paper
    has no z integral or t-z Jacobian. The code gives each root of the on-shell condition the weight 1, i.e. treats
    `dphi* = dphi_rho`; the exact weight is `sqrt(lambda)/(2 q kappa)` (mean 1.015, large excursions at the edges, exactly 1
    when integrated over the azimuth). Not changed; whether it matters with the azimuth-dependent Glauber distributions
    was not studied.
  - **Full comparison with the Glauber distributions (2026-11-26, `docs/rho_jacobian.md` section 6):** helium, Q² = 2 GeV²,
    nu = 3.0-3.2 GeV, 2 x 3 Gauss-Legendre points per version, real `RhoTCross`: z/t = 0.999-1.022 for all nine
    cases (plane wave 1.014) and transparencies from t and z within 0.011. A coarse test (cut edges, quadrature), helium
    only: carbon and a finer grid are still to do if wanted. The "before" value of the unfixed code was not run (the
    Monte Carlo gives the factor).
  - **W > 2 GeV cut (2026-11-26):** the analysis of the paper also cuts W > 2 GeV, which `getCrosst`/`getCrossz` do not
    (it is a limit on the nu range at fixed Q^2). `RhoTCross::checkWcut(Q2, nu_min)` (with `getW` and `Wmin`) now throws
    `std::invalid_argument` when W(nu_min) < 2 GeV, W^2 = m_p^2 + 2 m_p nu - Q^2 (nucleon at rest, proton mass:
    inferred, to confirm), and `rhoset` and `rhosingle` call it unless `nocuts` is set; unit test in `unit_errors`.
    Existing calculations whose input violated the cut would now stop: none of the helium test runs did (W = 2.12 GeV at
    nu = 3.0, Q^2 = 2). `pionT_calc12` (the pion program) has no such check (the pion analysis has other cuts, not looked at).
  - **Not verified:** the absolute normalisation of the t version was compared with the paper's Eq. (13) (it
    matches), not with an independent calculation.

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
  message (`WeakQEHadronCurrent.cpp` had CRLF line endings then; converted to LF in #66).
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

- [x] **Constants macros** (#41)**:** all uses of the macro aliases (about 1700, in 143 files: libraries,
  programs, tests) are replaced by the `nucqcd::` names, and the aliases are removed from
  `constants.hpp.in`, including `HOMEDIR` (now `nucqcd::share_dir`). Only `PI` is left, a literal,
  until `numint` drops its own. Done with a script that skips string literals and keeps line
  endings; `progs/GPD/diDVCS_gen_gridonly.cpp` keeps its own `#define ALPHA` (same value).
  `ClassGridThick_SCX` takes its share directory as `const char*` (it took `char*`, which only
  accepted the string literal `HOMEDIR` by a GCC extension). Values are unchanged: all tests,
  including `slow`, pass with the same references. `unit_constants` now fails to compile if an old
  alias comes back.
- [x] **Update the constants to current PDG values: decided not to (author, 2026-11-26).** The constants
  are not of one PDG edition but a mix of vintages, as they were when the code was written. Checked against the sources:
  the proton mass 938.27231 MeV is the 1986 CODATA value that the PDG listings of the 1990s quote (current PDG: 938.27209
  MeV, from CODATA 2018), `alpha_em` 0.00729735253 is the CODATA 1998 value, and the neutron mass 939.56536 MeV matches no
  edition that was found (the 1986 CODATA/PDG 1998 value is 939.56563 MeV, so the digits may have been swapped when it
  was typed; current PDG: 939.56542 MeV); the boson masses, `G_F`, `V_ud` and `sin2_theta_w` look like about PDG 2008-2012
  and the K0 and Sigma masses like about 2000 (judged from the numbers only). The changes to the current values are below
  1e-6 relative for the nucleon masses, `alpha_em` and `hbarc` and of the order of 1e-4 to 1e-3 for the boson masses,
  `V_ud`, `sin2_theta_w` and a few hadron masses (rough estimates, not re-checked), far below the accuracy of the
  calculations (numerical precisions of 1e-3 to 1e-5, model uncertainties of the order of percent), so the values stay.
  The note is in `Tools/Core/constants.hpp.in`. The rho mass was already moved to PDG 2022-2026 (step 5).
- [x] **Code layout** (sharpened 2026-11-22 with the author; all four parts done: #52, #53, #65, and the rename of the `T*` classes #76). Headers next to their sources is fine and
  stays (no general `include/`/`src/` split: that pays off for libraries installed for outside
  projects, which nucqcd's are not, and moving every file would cost the `git blame` history). What
  does not follow good practice, in this order:
  1. ✅ **Split `Tools` (library `MePhys`)** (done 2026-11-23, #52): the 36 headers plus Fortran are now
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
  2. ✅ **Prefixed includes, with the repository root as the include root** (done 2026-11-23, #53).
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
  3. ✅ **One header extension: `.hpp`** (done 2026-11-23, same PR, #53): the 19 `.h` headers in Tools
     were renamed with `git mv`, and their `\file` tags and mentions updated. Only the third-party
     `Integrator/adaptive/cubature.h` and Cuba keep `.h`.
  4. [x] **Namespace `nucqcd`** (done 2026-11-25, #65). Instead of class by class, one scripted change:
     the declarations and definitions in all headers and sources of `Tools/{Core,Kinematics,WaveFunctions,
     StructureFunctions}`, `Glauber`, `Knockout`, `DIS` and `GPD_Partons` (147 files; the files that already
     had a `namespace nucqcd` block were left alone) are wrapped in `namespace nucqcd { ... }` after their
     includes, so classes, free functions (`interpolate`, `rombergerN`, ...), enums and typedefs are
     `nucqcd::`-qualified. The 130 programs, tests and program-local headers that include a library header
     got `using namespace nucqcd;` after their includes (the "compatibility" step: no program had to be
     rewritten; new programs can qualify names or do the same). Forward declarations of library classes
     are written `namespace nucqcd { class X; }`. Left alone: `Integrator` (own namespace `numint`), the
     Fortran, the `extern "C"` declarations (C linkage has no namespaces), the macros, and the existing
     explicit `nucqcd::` qualifications in the libraries (redundant now, harmless). No values change: all
     tests including `slow` pass with unchanged references; `make doc` gives no warnings (the Doxygen groups
     of `EMQECross.hpp` and `WeakQECross.hpp` had to be closed outside the namespace).
- [~] **Glauber grid hierarchy:** deduplicate it. **Done for the two thick grids (2026-11-25, #77), the abstract classes are not touched.** `GlauberGridThick` and `GlauberDecayGridThick` are
  largely copy-pasted, as are the {plain, Thick} × {CT, Decay} abstract classes.
  **Safety net done (2026-11-25, #68):** `reg_glauber_variants` [slow, ~90 s] covers what the other tests
  did not: `GlauberGrid` (the thin-nucleus grid only works off axis, it throws "use OneGlauberGrid" for a
  particle along z), `GlauberGridThick` with the SRC grids, `GlauberDecayGridThick` (all 8 grids; rho,
  rho→ππ pair, proton; along z and off axis), `OneGlauberGrid` and `ROMEAGrid`; tiny grids, 4 points each,
  references captured before any refactor. Until then `GlauberDecayGridThick` (used by `RhoTCross` and
  `RhoDeuteron`) was covered only by a compile check.
  **Measured scope** (same date): `GlauberDecayGridThick.cpp` has 1005 non-blank lines, 738 of them
  identical to `GlauberGridThick.cpp` after renaming; the headers share 146 of 237. The abstract classes
  are one interface layer per named grid (`Src`, `Ct`, `Decay`, ...) with four interpolation wrappers each.
  The two concrete classes cannot simply be merged: the decay class computes 8 grids at once (the 4
  extra ones with decay), so using it for plain grids would double their cost, change the cache files and
  the adaptive-integration sampling (values within the precision, not bitwise). The dedup therefore needs
  a shared implementation parameterized on the number of variants, with the grids and file names left as
  they are; the test above must stay bitwise equal. Not started.
  **Dedup of `GlauberGridThick` and `GlauberDecayGridThick` (#77):** the two classes differed only in the
  number of variants in the first index of the grid arrays (2: RMSGA, +SRC; 4: also +decay) and in the integral
  of one grid point; everything else was copied. New `Glauber/GlauberThickGridCommon` (a class that derives
  virtually from `AbstractFsiCTGrid`, so that it can override its pure virtual functions) holds, for a number
  of variants given to its constructor, the two arrays (allocation, destruction), `constructAllGrids` and
  `constructCtGrid` (one loop `constructGrids(ctOnly)`), `readinFsiGrid`/`readinFsiCtGrid`/`writeoutFsiGrid`/
  `writeoutFsiCtGrid`, `getFsiGridN_interp`, `print_grid` and `interpolateVariant` (the proton and neutron
  factors of one variant). The two classes keep what differs: `calcGlauberphases*`, the integrands
  `klaas_int_*`, the file names, and the named interpolation and print functions, now one line each. To
  make it possible `AbstractFsiCTGridThick` inherits `AbstractFsiCTGrid` virtually (as
  `AbstractFsiCTDecayGrid` already did). The loops are the old ones, including their quirk (the r = 0 shortcut
  reuses the loop variables j and k, which ends the loops). Also removed: the integrand functions that only
  the old `rombergerN` branches called (declared in `GlauberGrid`, `GlauberGridThick`, `GlauberDecayGridThick`,
  `OneGlauberGrid`; in `OneGlauberGrid` they were live code, unreachable since option 0 is gone, which also
  removes its hard-coded `rombergerN` uses), the commented-out copies of them, and the commented
  `rombergerN` calls. Lines: `GlauberGridThick.cpp` 1104 → 519, `GlauberDecayGridThick.cpp` 1008 → 443,
  `OneGlauberGrid.cpp` 278 → 185, `GlauberGrid.cpp` 754 → 664, plus the headers (about 250 lines) and the new
  class (about 250). **Checked:** all 23 cache files that `reg_glauber_variants` and `reg_glauber_grid`
  write (thin, thick, decay: along z and off axis) are byte-identical before and after (the cache file format
  is unchanged), the references of both tests pass, all other tests unchanged. One output change: the
  decay grid's `constructCtGrid` printed "Mem reservation for treshold array" to stdout; it no longer does.
  **The abstract classes (2026-11-26):** the seven classes `AbstractFsiGrid`, `AbstractFsiCTGrid` and the five between
  them and the Glauber classes (`AbstractFsiGridThick`, `AbstractFsiDecayGridThick`, `AbstractFsiCTGridThick`,
  `AbstractFsiCTDecayGrid`, `AbstractFsiCTDecayGridThick`: 10 files, one named interface per variant, each with four
  interpolation wrappers, repeated, in a diamond of virtual bases; `GlauberDecayGridThick` could not even call
  `getFsiDecayGridFull_interp3` without a cast, it was ambiguous) are now two. `AbstractFsiGrid` has the enumeration
  `GridVariant` (bit 0 SRC, bit 1 decay, bit 2 CT; the value is the grid number of `GlauberDecayGridThick`), a virtual
  `getGridNumber(variant)` (-1 if the class has no such grid; default: only `Fsi`), `hasVariant`, and the 40 named
  functions (`getFsi[Ct|Src|SrcCt|Decay|...]GridFull_interp[|vec|3|2|1]`) as inline one-liners that call
  `getFsiGridN_interp*` with the grid number of the variant (a variant that the class does not have throws
  `std::logic_error`; before, it did not compile). The call sites (about 60) are unchanged. `GlauberThickGridCommon`
  took over the `FsiCorrelator`, the thick nucleus, `number_of_grids` (2 x variants), the grid numbers and the `Thick.`
  prefix of the file names (`GlauberDecayGridThick` adds its decay string itself); `GlauberGrid` and `ROMEAGrid` moved
  the bodies of their named interpolations into `getFsiGridN_interp`. The inheritance is a single chain
  (`AbstractFsiGrid` <- `AbstractFsiCTGrid` <- `GlauberThickGridCommon` <- the thick classes), no virtual bases. The
  `dynamic_cast`s to the thick abstract classes in `EMQEHadronCurrent` and `WeakQEHadronCurrent` are gone (the grid
  pointer has the functions now). 1022 lines removed, 197 added. **Checked:** a probe (before and after, scratch
  program) that builds each concrete class with a tiny grid and prints `number_of_grids`, the cache file names, all named
  interpolations and `getFsiGridN_interp` at a point (12 digits) is identical (71 lines), and so are the 11 cache files
  it writes (names and md5); `ctest -L slow` with the 23 cache files of the two Glauber regression tests: identical names
  and md5 to master; `ctest -LE slow` 136/136.
  **The grid arrays of the Glauber classes (2026-11-26, part A of the raw pointer arrays):** the `std::complex<double> *****`
  arrays (`fsi_grid`, `fsi_ct_grid` of `GlauberThickGridCommon`: `[variant][nucleon][r][cos theta][phi]`, and of `GlauberGrid`:
  `[shell][m][r][cos theta][phi]`, ragged in `m`), the 3D `fsi_grid` of `ROMEAGrid`, the `int **` threshold arrays and the
  `complex<double> **` optical potential were allocated with a `new` per row and freed with matching loops. New
  `Tools/Core/Grid3D.hpp`: `Grid3D<T>`, an owning 3D array in one contiguous `std::vector` with the `grid[i][j][k]` syntax
  (value-initialized, copyable, movable; iterators in the order of the cache files). The 5D arrays are
  `std::vector<std::vector<Grid3D<std::complex<double>>>>` (`[variant][nucleon]` or `[shell][m]` outside, a Grid3D for
  r, cos theta, phi), the threshold arrays are local `std::vector<std::vector<int>>` in the loops that use them (they were
  members that were allocated and freed inside one function; the check "Treshold should only be initialized once!" is
  gone with them), the optical potential is a `std::vector<std::vector<std::complex<double>>>` that is cleared at the end of
  `constructAllGrids` (the getter returns a const reference). `Interp3d` and `getInterp` take a `Grid3D`. The manual
  new/delete code of 4 files is gone; the index uses (about 250) are unchanged. **Checked:** `ctest -LE slow` 139/139 (3 new
  tests in `unit_grid3d`: indexing and order, copy and move, `Interp3d`); `ctest -L slow` 11/11; the before/after probe
  of #92 (all named interpolations to 12 digits, file names, `number_of_grids`) identical (71 lines); of the 23 cache
  files of the Glauber regression tests 19 are byte-identical to master and **4 differ in the last bit**: the thin
  `GlauberGrid` and `OneGlauberGrid` files (FSI and CT each), in the block of the product over all nucleons (64 of the 700
  numbers of the probe file differ by 1 ulp, 2e-16 relative). **Cause (tested):** FMA contraction (`-O3 -march=native`)
  chose differently for the changed array access (`GlauberGrid.o`: 38 fused multiply-add instructions, 37 on master); with
  master and the branch both built with `-ffp-contract=off` (scratch builds of all Glauber sources) all 11 probe cache files and
  the probe output are byte-identical. So the logic is the same; cached grids from before stay valid (the file format is
  unchanged, the values differ at 1e-16).
  **Part B (2026-11-26): `DistMomDistrGrid`, `FsiCorrelator`, `FastParticle`.** `DistMomDistrGrid`: the 4D `rhogrid` and
  `rhoctgrid` (`double ****`, `[grid][p][cos theta][phi]`) are `std::vector<Grid3D<double>>`, `rhopwgrid` a `std::vector<double>`
  (the getters return const references; `Interp3d` for `double` takes a `Grid3D`), and the destructor, the copy constructor and
  the assignment operator (about 160 lines that copied every member by hand; the assignment also did not free the arrays it
  overwrote) are gone: the class is copied, moved and destroyed by the compiler. `FsiCorrelator`: the Newton solver for the
  correlation grid (`functionmatrix`, `x`, and the Jacobian, gradient, step vectors of `newt`/`lnsrch`, with their manual
  `corrmaintenance` frees on every exit) works on `std::vector`s; `LUsolve` got an overload for `std::vector` (it builds the
  row pointers and calls the pointer version, so the arithmetic is the same; `fitdata` keeps its raw arrays and its `double *`
  result). `FastParticle::setPionGlauberData` (used by `progs/pionprog` to make the pion table): the local data tables are
  `std::vector<DataTable>` (a table of three columns; `table[j]` is the pointer to column j, so the 50 uses and the calls of
  `fitdata` are unchanged); the old arrays were freed on the normal path but leaked when a data file could not be opened.
  **Checked:** `ctest -LE slow` 140/140 and `ctest -L slow` 11/11 with the 23 cache files identical to master; and, because
  none of these three has a test of its own, scratch comparisons master against the branch (the probes are not in the
  repository): `setPionGlauberData` for 30 momenta (all six outputs, 17 digits): identical; `DistMomDistrGrid` with FSI for a
  rho in carbon (the flow of `RhoTCross`: fill, interpolate for 5 grids at 3 momenta, copy constructor, assignment, a
  `std::map` entry, assignment over a filled object; 80 values at 17 digits, and the 6 cache files): identical;
  `FsiCorrelator` for carbon on 6 x 6, 12 x 12 and 20 x 20 grids (the solver only runs when `share/input/CorrGrid.*` is
  missing, which is never the case for the shipped nuclei, so no test runs it): printed interpolated values, the written cache file
  (md5) and the log (iteration counts included) identical. New test: `LUsolve` for vectors and for pointers.
  Still to do for this item: only the named print functions (`printFsi_*_grid`, declared in the concrete classes only now, with their
  per-class header text, not in the abstract ones; low priority).
  **Found on the way, fixed (2026-11-25, #69):** in the Glauber grids `integrator == 0` (the old
  `rombergerN`) was dead: the `rombergerN` calls were commented out in `GlauberGrid`, `GlauberGridThick`,
  `GlauberDecayGridThick` (and in `DistMomDistrGrid`, `EMQEHadronCurrent`, `RhoDeuteron`), so the results
  stayed 0 and the grids silently came out as "no FSI" (1 − σ·0). Nothing known was affected: the tests and
  the published runs use 2. **The option is removed** (author's decision): the `integrator == 0` branches
  are gone from those classes and from the live ones in `OneGlauberGrid` and `RhoTCross`, and every
  constructor with an `integrator` option (`AbstractFsiGrid` and so all grids, `DistMomDistrGrid`,
  `EMQEHadronCurrent`, `WeakQEHadronCurrent`, `EMQECross`, `WeakQECross`, `RhoTCross`, `RhoDeuteron`,
  `PionTCross`) now rejects anything but 1 (`numint::cube_romb`) and 2 (`numint::cube_adaptive`) with
  `std::invalid_argument` (`nucqcd::checked_integrator_flag`, `Tools/Core/IntegratorFlag.hpp`). A program
  that still passes 0 now fails at construction instead of giving wrong results or (in `OneGlauberGrid`,
  `RhoTCross`) using `rombergerN`. Not touched: `DQEinclusive`, whose `integrator` argument is a different
  switch (its case 0 is `cube_adaptive`), and the hard-coded `rombergerN` uses (next item). Tested in
  `unit_errors`; all tests including `slow` pass with unchanged references.
- [~] **Speed of the Glauber grids (added 2026-11-23; off-axis thick grids done 2026-11-26, #97; parallel grid points open).** From the author's experience,
  grids for particles that do not all move along the z axis are slow to compute (cached afterwards,
  but still). The `numint` integrators used there (`cube_adaptive`, `cube_romb`) are serial; only the
  explicit `numint::cuhre/vegas/...` wrappers use Cuba, which forks workers on one node
  (`CUBACORES`). First measure where the time goes, then compare: (a) filling the grid points in
  parallel threads (each point is an independent integral, but the integrand reads the point from
  members such as `r_hit`, so it has to be passed in explicitly first); (b) the vectorized
  `hcubature_v` (see the cubature item under step 7), which hands over batches of points that can be
  evaluated in threads; (c) Cuhre with `CUBACORES`. Best done together with the deduplication above.
  (b) and (c) change values within the precision, so they need the author's OK.
  **Measured (2026-11-23, `docs/integrator_benchmark.md`):** a 12×6×3 carbon grid for a 1500 MeV proton
  takes 0.4 s along z but 76 s off axis (θ = 0.5 rad) with `cube_adaptive` at 1e-3; off axis the
  integrals run into the hard-coded cap of 2·10⁶ (5·10⁶ with CT) evaluations per grid point, so asking for
  1e-5 or 1e-6 gives the same grid (accuracy ~1e-3) and a 100× larger cap did not finish in 2 hours.
  So the off-axis cost is "every point hits the cap"; `cube_romb` at 1e-3 took 16 s there with similar
  accuracy. Next: find out which part of the off-axis integrand is hard (the cut at the hit point, the
  profile function's oscillations?), consider a better variable choice, `cube_romb` for these grids, or
  parallel grid points.
  **Done (2026-11-26, #97): the thick grids off axis in the frame of the particle.** The cause was the integrand, not the
  integrator: in the spherical coordinates of the lab the Glauber phase of one particle is a step along a plane
  (z along the momentum > z of the hit point) times a Gaussian in the transverse distance to the particle's line, a narrow tilted tube;
  neither lines up with the integration variables, so no integrator converged within the cap. For a particle along z the code
  already used cylindrical coordinates around the momentum (`klaas_int_bound`), which turn the step into a limit of the
  z integral and the tube into a Gaussian in b. `GlauberThickGridCommon::integrateTilted` does the same for one particle in any direction
  (variables b, z, phi' around the momentum, phi' from the direction of the b vector of the hit point; the lab coordinates of the point,
  needed for the SRC factors, are rebuilt from the frame; an incoming particle integrates the other half of z), for
  `GlauberGridThick` and `GlauberDecayGridThick` (several particles still use the lab coordinates; `setUseTiltedFrame(false)` selects them
  for one particle, as before, for comparisons). **Measured** (12 x 6 x 3 carbon grid, 1500 MeV proton, theta 0.5 rad, `cube_adaptive` 1e-3):
  61 s before, 0.6 s with the plain FSI/CT variants and 3.2 s with the SRC grids (the normalization integral, see below);
  the 4 x 4 x 3 rho (type 4) grid: 0.6 s. **Values:** the old values were cap-limited, with errors of up to 2.3e-3 (relative to the converged value)
  although the integrator reported 1e-3: new at 1e-3 versus new at 1e-7 differ by 7e-5, old at 1e-3 versus the old integrand with `cube_romb` at 1e-4 by 2e-3,
  the new integrand with `cube_romb` 1e-4 and `cube_adaptive` 1e-6 agree to 1e-7 at the point checked, the two old ones to 6e-5; new versus old, both tight: 8e-5 (adaptive 1e-6)
  and 1.5e-4 (romb 1e-4) for all eight decay variants including SRC. **A trap on the way:** the SRC variants of the lab
  integrand contain the term N = integral of density * correlation over the whole nucleus (1 - the FSI part is then N - amplitude * integral);
  `klaas_int_bound` (particle along z) takes N = 1. First version of the new path did too and differed from the old by up to 8.5e-3 in the SRC variants;
  the normalization is now computed (an extra integral per grid point, at precision prec/100 because the correlation peak at the hit point makes
  the integrator stop too early at prec: errors of 3e-3 at prec 1e-3). The test `reg_glauber_tilted` compares the two coordinate systems (proton, rho, rho pion pair, incoming proton);
  the references of the off-axis entries of `reg_glauber_grid` and `reg_glauber_variants` were regenerated (largest differences 1e-3 absolute in the former, 9e-5 in the latter; the entries on the z axis and
  of the thin grid did not move).
  **Grid files in `share/grids` computed off axis before this are not recomputed** (the file name does not record the physics): delete them to get the new values.
  **Open leads:** (1) filling grid points in parallel (the integrands still read the hit point from members of the grid and of the particle; it would be
  all grid points independent, about 8x on 8 cores, values unchanged); (2) N does not depend on phi of the hit point or the particle, so one integral per (r, cos theta) and proton
  would do instead of one per grid point (5x fewer in the 120 x 36 x 5 grid; most of the 3.2 s above); (3) **the SRC variants for a particle along z omit N**:
  they use `1 - amplitude * integral` where the off-axis (old and new) formula has `N - amplitude * integral`; N * (correlation grid at the hit point)
  differs from 1 by about 1e-2 at some points (the 8.5e-3 above), so the on-axis SRC grids may be off by that much: to quantify and, with the author's OK, fix
  (changes values of A(e,e'p) with SRC); (4) the thin `GlauberGrid` has the same lab-coordinate integrand with the single-particle wavefunction as weight
  and was not converted (its test takes 20 s).
- [x] **Remove the integrators `rombergerN` and `cube_romb`?** (`rombergerN` removed 2026-11-25, #82; `cube_romb` stays) (added 2026-11-23). The author
  remembers the adaptive integrators doing better, especially for multidimensional integrands with
  non-trivial boundaries or cuts inside the integrand. **Benchmark done 2026-11-23
  (`docs/integrator_benchmark.md`, `progs/benchmark/integrators.cpp`, #61): the memory is not
  confirmed by these tests.** On integrands with a step or a kink inside the domain, `cube_adaptive`
  (and upstream `hcubature` 1.0.4) is the least reliable (stalls at `maxEval` with 1e-5 to 1.5e-3 and
  reports success; 1.0.4 even claims convergence on a kink), while the Romberg integrators converge,
  at a high price. `cube_adaptive` is much faster on smooth integrands and peaks (1400× on the
  along-z Glauber grid); `cube_romb` was faster on the off-axis grid at 1e-3 and is comparable on the
  A(e,e'p) amplitudes. Cuhre misses narrow peaks at loose tolerance. **Recommendation: keep
  `cube_romb` and `cube_adaptive`** (author's decision); `rombergerN` could still go case by case. The
  benchmark plan was: compare `rombergerN`, `cube_romb`, `cube_adaptive` (and Cuhre) on accuracy versus number of
  evaluations and time: (1) test integrands with known values in 2 and 3 dimensions: smooth ones, ones
  with a step or cut inside the domain, a non-rectangular domain mapped onto a box (as with the
  impact-parameter and z integrations in the Glauber phases), and a peaked one; (2) real integrands:
  a Glauber grid point for a particle along z and one off-axis, an A(e,e'p) amplitude from `Model`,
  and an `InclusiveCrossRes` integral (not done: it uses `rombergerN` only, without a flag), each
  against a high-precision reference. If they are removed: the `integrator` flag of the grid
  and cross-section classes no longer has the `rombergerN` option 0 (removed in #69, see "Glauber grid
  hierarchy"; 1: `numint::cube_romb`, 2: `cube_adaptive`; the tests use 2). Hard-coded uses have to be replaced, which changes values
  (author's OK per class). **Survey and first step (2026-11-25, #80):** of the 25 calls in the libraries and the programs, 6 (in 12 functions) were dead: the
  nested chains `intPm*` → `intCosTheta*` → `intPhi*` of `RhoTCross` (6 functions) and of `RhoDeuteron`
  (t and z, 6 functions) were only called by the `integrator == 0` branches removed in #69 (and
  `OneGlauberGrid` in #77); they are deleted (no value changes; `RhoTCross` has no `rombergerN` left).
  **Outcome (2026-11-25, #81, #82; the author chose: archive what is abandoned, replace the rest):**
  - `RhoDeuteron` (rho production on the deuteron, an abandoned calculation) moved from `Knockout` to
    `progs/attic/rhodeuteron` with `rhodint` (from `progs/RhoTransparency`); both still build. `rhod`
    (`progs/attic/mpi`, built only with MPI, not compiled on this computer) includes it from there.
  - `InclusiveCrossRes` (6 `rombergerN` calls; the resonance model for inclusive DIS, W. Cosyn, W. Melnitchouk and
    M. Sargsian, Phys. Rev. C 89, 014612 (2014); already not built and no longer compiling: `DIS/CMakeLists.txt` had it
    commented out, `DeuteronStructure` has a different constructor now) moved with `threeres` and `increstest` to
    `progs/attic/threeres`, not built (its `CMakeLists.txt` says so). `DIS/Resonance` stays in the DIS library.
  - `DistMomDistrGrid::constructpwGrid` (the plane-wave momentum distribution grid, the only live library use):
    the 1D radial integral is now a 1D `numint::cube_adaptive` call (1e-7 relative). The regression test
    `reg_distmomdistr_pw` (carbon s and p shells, 9 momenta each, recorded with `rombergerN` in #81) passes at
    1e-6 relative, so the values moved by less than that.
  - `rhoset`: the `integrator == 0` branch (its own option, the nested `rombergerN`) and its helper functions
    `nuint`/`tzint` are removed, as option 0 was in the classes (#69); the program checks the flag
    (`checked_integrator_flag`). 1 and 2 (`cube_romb`, `cube_adaptive`) were already there and are unchanged.
  - `rombergerN` itself left `Tools/Core/Utilfunctions.hpp` for `progs/attic/rombergerN.hpp`, which the archived
    `RhoDeuteron`, `rhodint` and `InclusiveCrossRes` include; so does the benchmark (`progs/benchmark/integrators.cpp`,
    method 0, so that `docs/integrator_benchmark.md` stays reproducible). Nothing in the libraries uses it.
  - Tests: `reg_rombergern_users` became `reg_distmomdistr_pw` (the `RhoDeuteron` case goes with the archive);
    the `RhoDeuteron` check was removed from `unit_errors`, the `InclusiveCrossRes` assertion from `unit_ownership`.
  **Alternatives that were considered** for the live uses (they differ in how much of the code
  changes and in what the values do): (a) *keep the nested structure, change the driver* to `numint`'s 1D
  Romberg: the least change to the physics code (the integrands stay), values move at the level of the
  stopping rule; (b) *one multi-dimensional call* of `cube_adaptive` (or `cube_romb`) over the nested
  variables, as in the Glauber grids and the A(e,e'p) amplitudes: fewer evaluations for smooth
  integrands (benchmark: up to 1000x), but the limits of `InclusiveCrossRes` depend on pr (a
  non-rectangular domain: needs a change of variables) and cuts inside the integrand are where
  `cube_adaptive` was the least reliable in the benchmark; (c) *GSL's adaptive 1D routines*
  (`gsl_integration_qag`/`cquad`), nested: robust against cuts, an error estimate per level, a new
  dependency use in `Tools` (GSL is already required); (d) *keep `rombergerN`* for the nested 1D chains and
  only move it into `Integrator/numint` as a documented legacy routine, so that the question "which
  integrator" stays open for the old programs. A regression test with the current values (a few
  kinematic points of each of the four cases) should be written first, so that the effect of any
  replacement is measured, not guessed: none of the four has one today.
  `cube_romb` in `Model.cpp:744` and `progs/momdistr/momdistr.cpp`. Then remove
  `Integrator/numint/cube_romb.cpp` (and `cube_tr`/`fubini`/`vfubini` if unused) and its tests in
  `unit_numint`, and simplify the `integrator` flag (an enum, see "Model selection" below).
- [x] **Better names for `Model` and `Cross`** (done 2026-11-23, in the include PR #53): `Knockout/Model`
  and `Knockout/Cross` are now `EMQEHadronCurrent` and `EMQECross` (files and classes, `git mv`),
  following `WeakQEHadronCurrent`/`WeakQECross` and `NucleonEMOperator`/`NucleonWeakOperator`; the
  class names were replaced in the 19 files that include them (library, programs including the
  attic, tests). Also fixed: the `\file` tag of `DQEinclusive.hpp` said `Cross.hpp`.
- [x] **The requested precision often does not control the accuracy** (reporting done 2026-11-25, #62; the caps themselves are unchanged) (found with the integrator
  benchmark, 2026-11-23). `cube_adaptive` returns 0 (success) when it stops at `maxEval`, and the
  `numint` wrapper (`numint2adapt.hpp`) discards the error estimate, so callers never learn that the
  precision was not reached. In practice the evaluation caps decide: `maxEval` of
  `EMQECross::getDiffCross` (the 24-component amplitude integral always runs to it; the programs mostly
  pass 2·10⁴, which gives 1e-5 to 1e-4 relative, see the benchmark), and the hard-coded 2·10⁶/5·10⁶ in
  the Glauber grids (off-axis grids stop there at ~1e-3). To do: return or log the error estimate
  and a "cap reached" flag from `cube_adaptive` (and check it, or at least warn, in the classes);
  document what `prec` and `maxEval` mean for each class. No values change by reporting. **Done:** `cube_adaptive` checks the cubature error estimate against the request (the cubature
  code's own test) and counts the calls that stopped without meeting it (`numint::adaptive_stats()`,
  `Integrator/numint/adaptive_stats.hpp`); the return value is unchanged. The first such call prints a
  warning to stderr (`adaptive_stats().warn = false` silences it), and `nucqcd::RunLog` writes
  `# cube_adaptive: N integrals, M stopped at maxEval ... (worst estimated relative error ...)` before its
  "finished" line when `cube_adaptive` was used. Tested in `unit_numint`. Still to do: per class,
  document what `prec` and `maxEval` mean, and decide case by case whether to raise caps.
- [x] **Model selection** (done 2026-11-25, #72). The "about 700 string comparisons" were mostly the
  variable names of `TKinematics*` (now one table, #71, and format strings that are parsed by design);
  the comparisons of model names in the libraries were about 25, in four places, and are now
  enumerations parsed in one place each:
  - `DeuteronModel` (`Tools/WaveFunctions/WaveFunctionModels.hpp`): the factory
    `TDeuteron::Wavefunction::CreateWavefunction` and `OldDeuteron`; `He3Model`: `He3wf`. The name
    constructors and factories stay (programs read names from the command line) and call
    `parse_deuteron_model` / `parse_he3_model`; new overloads take the enumeration, and the classes switch on
    it. The error message of an unknown deuteron name is the same as before.
  - `StructureFunctionModel` (`Tools/StructureFunctions/StructureFunctionModel.hpp`): the plain names of
    `NucleonStructure` ("SLAC", "Alekhin", "CB", "CTEQ", "HMRS", "MSTW"); the grammar (`LHAPDF:<set>`,
    `<low>+<high>@<W2>`) stays a string and is parsed as before. New constructor
    `NucleonStructure(StructureFunctionModel)`.
  - Names stay case sensitive for He3 and the structure functions and case insensitive for deuterons, as
    they were. **Behaviour changes** (found while doing it, all turning silent or late failures into early
    ones): (1) a misspelled `NucleonStructure` name used to construct and only threw when a structure
    function was asked for, it now throws `std::invalid_argument` at construction, listing the valid names
    (nothing in the repository or the tests relied on it); (2) `OldDeuteron` with an unknown name used to
    print a message and `assert(1==0)`, which does nothing in a release build, leaving an empty
    parametrization (zero wave function); it throws now; (3) the He3 error message now also lists the
    valid names. Tested in the new `unit_models`; all tests including `slow` pass with unchanged references.
  - Still strings (not model names in this sense, or external): the structure-function grammar, the GPD
    model names handed to PARTONS (`GPD_Partons`), `pdf_name` of `GPD_T_Nucl_grid` ("MSTW" / "LHAPDF:...",
    parsed by `LHAPDFStructure`), and the kinematics format strings.
- [x] **Merge `NuclStructure` into `NucleonStructure`** (#42). The parametrizations are now implemented
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
- [x] **Retire MPI (`TMPI`)** (done 2026-11-23, #54, author's decision). `TMPI` (from strangecalc) only
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
- [x] **Rename the files and classes that start with `T`** (done 2026-11-25, #76, the author's request and
  choice of names): the prefix was a ROOT naming convention (the classes came from the strangecalc code,
  built on ROOT), and ROOT is no longer used by the libraries (step 9). Done with `git mv` for the files (so
  `git blame` follows) and one scripted replacement, with word boundaries, of the names in the libraries,
  programs, tests, CMake lists and CLAUDE.md, in a commit of its own that is in `.git-blame-ignore-revs`;
  no code changed:
  | old | new |
  |---|---|
  | `TKinematics`, `TKinematics2to2`, `TKinematics2to3`, `TKinematics2to3WithLabAngles` | `Kinematics`, `Kinematics2to2`, `Kinematics2to3`, `Kinematics2to3WithLabAngles` (also the files `Kinematics2to3_*.cpp`) |
  | `TElectronKinematics`, `TLeptonKinematics` | `ElectronKinematics`, `LeptonKinematics` |
  | `TSpinor`, `TVSpinor`, `THelicitySpinor`, `THelicityVSpinor` | `Spinor`, `VSpinor`, `HelicitySpinor`, `HelicityVSpinor` |
  | `TLorentzQuaternion` | `LorentzQuaternion` |
  | `TDeuteron` (and its nested `Wavefunction`, `Polarization`) | `Deuteron` |
  | `TMFSpinor` | `MeanFieldSpinor` |
  | `TYukawaPWF`, `TGrossPWF`, `TCstPWF` | `YukawaWavefunction`, `GrossWavefunction`, `CstWavefunction` |
  | `TInterpolatingWavefunction`, `TWavefunctionImplementation` | `InterpolatingWavefunction`, `WavefunctionImplementation` |
  The include guards got the new names too. The local struct `Kinematics` of `NuclStructure.cpp`, which
  would have collided with the new class, is `OnShellPoint`. **Not renamed:** `TMPI` (the attic copy in
  `progs/attic/mpi`), and the ROOT classes in the comments. **This ROADMAP and `docs/HISTORY.md` still use
  the old names where they describe what was done before this rename** (PR numbers #1-#75 refer to the old
  names). **The `TNamed` leftovers are gone too (2026-11-25, #78):** `Kinematics`, `Kinematics2to2`,
  `Kinematics2to3` and `Kinematics2to3WithLabAngles` had a name and a title (members, `GetName`, `GetTitle`,
  `SetName`, `SetTitle`, and the first two arguments of every constructor) that nothing used: 134 of the
  constructor calls (libraries, programs, tests) passed `"", ""`, three more were the names of the KY, KN and
  YN subsystems and `transp_eli` ("kinematics_name"), and the only readers were the first line of the two
  `Print`/`VariableInfo` outputs. The arguments and members are removed (a script took the two leading string
  literals out of the 134 calls, the compiler checked the rest); the output of `Kinematics::Print` and
  `Kinematics2to3::VariableInfo` now starts with the class name instead of "name: title". The virtual
  `Clone(const char *newname)` of ROOT's `TObject` (17 declarations: the kinematics, `Spinor` and its
  subclasses, `MeanFieldSpinor`, `ElectronKinematics`, `LeptonKinematics` and the wave function
  implementations) is `Clone()`: nobody passed a name. `Deuteron::Wavefunction` keeps its name
  (`GetName`: "CDBonn", "Paris", ..., which names the parametrization in an error message and is tested in
  `unit_stringutil` and `unit_models`) and loses the title, `SetName` and `SetTitle`.
  Checked: full build with `BUILD_PROGS` and PARTONS, `ctest -LE slow` 133/133, `ctest -L slow` 11/11,
  `make doc` no warnings, all with unchanged references.
- [~] **Very long files and functions** (files split 2026-11-25, #70; functions not): the files were split
  without touching any code (the functions were moved verbatim; all tests pass, references unchanged).
  `Tools/Kinematics/TKinematics2to3.cpp` (6.5k lines) is now six files: `TKinematics2to3.cpp` (class
  description, constructors, assignment, 580 lines), `_ReadFormatString`, `_Update`, `_Variables` (setting,
  stepping and the `GetVar*` by name), `_Getters` (vectors, angles, CM quantities, `operator==`,
  `GetVarByGetter`) and `_Frames` (boosts, rotations, bounds). **Update (#71):** the `GetVar*` by name are done, see below. `Glauber/NucleonEMOperator.cpp` (4.2k
  lines) was 3.9k lines of in-medium form-factor tables: they moved to `NucleonEMOperatorTables.cpp`,
  leaving 250 lines of code. **`GetVar*` by name done (2026-11-25, #71):** `GetVar`, `GetVarArray`, `GetPhysicalVarArray` and
  `GetVarByGetter` each had the same chain of 89 `if (varName == ...)` branches (152 names, with aliases
  such as `ynphincm`/`nyphincm`). They now call `TKinematics2to3::FindVarGetter`, one table in
  `TKinematics2to3_VarNames.cpp` (a private static member, since some getters are private); about 1330
  lines less. Checked against the old code with a scratch program on all 152 names (`GetVar` in lower and
  capitalized form, `GetVarArray`, `GetPhysicalVarArray` on a 5-step range): identical to the last digit;
  new unit test in `unit_kinematics`. Two small behaviour changes: (1) **a typo in the table:** the
  documented example `kyphikcm_deg` was spelled `kyhpikcm_deg`, so `GetVar("kyphikcm_deg")` threw
  "unknown"; both spellings are accepted now (the old one as an alias); (2) the error message of
  `GetPhysicalVarArray` named `GetVarArray` and now names itself. **The two long functions are done (2026-11-25, #74).**
  - `TKinematics2to3::ReadFormatString` (1.2k lines, one 15-25-line branch per variable name) is about 270
    lines: local helpers for the groups of variables (`energyVar`, `particleEnergy`, `labTheta`, `labPhi`,
    `cmAngle`, `pairS`, `tg`, `td`) hold the checks, a chain of one line per name calls them. The checks
    of each group were compared with the old branches (the 18 CM-angle ones by a script).
  - `UpdateKinematics` (0.9k lines) is a dispatcher (about 100 lines with the setup and the final
    assignment) and steps: `SolveEnergyScenario` (`SolvePhotonAndParticle1`, `RotateTo23Frame`,
    `SolveParticle2InCM`, `SolveParticle2InLab`) and `SolveMandelstamScenario` (`SelectMandelstamVariables`,
    `SolveMandelstamParticle1`, `RotateTo23Frame`, `SolveMandelstamParticles23`), with the intermediate
    results in a private struct `UpdateWork`; the four copies of "cosine just outside [-1,1]: clamp or
    unphysical" and of "sine from cosine" are `ClampCosine` and `SineFromCosine`; the 165 lines of
    commented-out code of the (refused) scenario 2 are gone. Private member functions, declared in the header.
  - **Checked against the old code** with a scratch program that builds a `TKinematics2to3` for 900k
    format strings (random strings with random values, and strings made of the three scenarios' variables
    with the values of 4000 physical events, all three axis conventions, isospin 1-6) and prints the
    exception message, or the physical flag, the variable names and 34 variables: the output is identical,
    including the messages and the 63 inputs that crash (below), except for 14 of the Mandelstam cases (angles
    that are 0 or pi up to rounding: `acos` of 1 - 1e-16 is 1.5e-8 in one build and 0 in the other; the
    floating-point contraction differs between the restructured functions) and one value that differs in
    the last digit. New tests in `unit_kinematics`: events rebuilt from a CM-angle format string and from a
    Mandelstam one, and the refusals (duplicates, conflicts, unknown names, no Q^2, the azimuth that is zero by
    definition).
  - **Found on the way (1 fixed in #75, 2 not):** (1) a format string that does not give
    scenario 1 all of its variables (e.g. `qsquared` twice, or two energies and no angle) is accepted and
    `UpdateKinematics` then dereferences a null `fE1Var`/`fTh1Var`/`fTh2Var`/`fPhi2Var`: a segmentation
    fault (63 of 900k fuzzed strings, e.g. `qsquared:ey:ykcosthycm:wlab:ykphiycm:en` for kN); **fixed (2026-11-25, #75):** `ReadFormatString` now throws "Format string is not sane!" when the
    variables of scenario 1 are missing (after the loop, so the other messages are unchanged); tried on
    the same 900k strings: the 63 crashes became exceptions and nothing else changed; tested in
    `unit_kinematics`. (A string that gives no scenario at all, `fKinSolution[0] == '-'`, is still accepted and
    leaves the vectors at zero: no crash, not changed.) (2) Quirks kept as they were: `pn` with an axis convention other than kN
    stores `fEn` (not `fPn`) with the name `P` (`pk` and `py` store their own); the checks in the CM-angle
    branches compare `fKinSolution[0] == 2` (the number, never true) instead of `'2'`.
  Still to do: the 1.0-1.2k-line `GlauberGridThick.cpp`, `GlauberDecayGridThick.cpp`, `GlauberGrid.cpp`,
  `DistMomDistrGrid.cpp` (with the grid dedup above), and `SolvePhotonAndParticle1` (134 lines) could be
  split further.

  **Glauber files (2026-11-26, #96):** the four longest files of `Glauber/` were split the same way, the functions moved verbatim into a new file
  of the class: `FastParticle.cpp` 986 lines -> 551 + `FastParticle_ScatteringData.cpp` 465 (`setGlauberParameters`, `setPionGlauberData`,
  `interpPionGlauberData` and the table struct of the latter), `DistMomDistrGrid.cpp` 636 -> 333 + `DistMomDistrGrid_Construct.cpp` 340
  (`constructAllGrids`, `constructCtGrid`, `constructpwGrid`, `intRhoRpw` and the two integrands), `FsiCorrelator.cpp` 585 -> 284 +
  `FsiCorrelator_Solver.cpp` 328 (the matrix, the correlation function and the Newton solver), `GlauberGrid.cpp` 540 -> 184 +
  `GlauberGrid_Construct.cpp` 380 (`constructAllGrids`, `constructCtGrid`, the Glauber phases of one grid point, the integrands). A splitter script
  checked that every line of the original is in exactly one of the two files. **Checked:** `ctest -LE slow` 141/141, `ctest -L slow` 11/11 and the 23 grid cache
  files **byte-identical** to master (moving functions between files can change the last bits through the compiler's choice of fused multiply-adds, as in
  #93; here it did not); `setPionGlauberData` for 30 momenta and `FsiCorrelator` on 6 x 6, 12 x 12 and 20 x 20 grids (log and cache file): identical
  to master. The doxygen limit `DOT_GRAPH_MAX_NODES` was raised from 50 to 100 (the include graph of `ErrorMessage.hpp` got too large with the new files).
  The longest files left in `Glauber/`: `FastParticle.cpp` 551, `GlauberGridThick.cpp` 451 and `GlauberDecayGridThick.cpp` 381 lines (and the data
  tables `NucleonEMOperatorTables.cpp` and `global2009.f`); elsewhere in the libraries `Tools/WaveFunctions/Deuteron.cpp` (2169), `Tools/Kinematics/Kinematics.cpp`
  (1745), `Knockout/DQEinclusive.cpp` (1275), `DIS/InclusiveCross.cpp` (1173), `Tools/Core/Spinor.cpp` (1165) and `GammaStructure.cpp` (1072) are
  longer than 1000 lines (not looked at).
- [x] **Comments that belong to another file (copy-paste leftovers)** (crawl of 2026-11-26, fixed the same day, mostly
  comments; **partly done**, see the rest below). A scripted crawl of the tracked C++ files outside `progs/attic` and
  `tests/` (file tags, header guards, the first comment block of every file against the file name, usage lines against
  the program name), reading the hits by hand; what had been cleaned in earlier items is not repeated. Fixed:
  - **Header guard collision:** `Knockout/EMQEHadronCurrent.hpp` and `Knockout/WeakQEHadronCurrent.hpp` both used
    `MODEL_HPP` (left over from the `Model` rename), so a translation unit that included both would silently have
    lost the second one (none did). Now `EMQEHADRONCURRENT_HPP` and `WEAKQEHADRONCURRENT_HPP`; `EMQECross.hpp`
    `CROSS_HPP` -> `EMQECROSS_HPP`. A scan of all headers finds no duplicate guard left.
  - **Wrong `\file` tag or brief:** `Glauber/NucleonWeakOperator.hpp`, `DIS/InclusiveCross.hpp`, the archived
    `progs/attic/threeres/InclusiveCrossRes.hpp`, the first line of `Tools/WaveFunctions/YukawaWavefunction.cpp`, and
    the brief of `GPD_Partons/Deut_GPD_T_set.hpp` ("vector helicity amplitudes", a copy of the `V` set: now
    "transversity", **inferred from the class name and `Deut_Conv_GPD_T`, to confirm**).
  - **Usage lines naming another program or giving its arguments:** `progs/neutrinoQE/ddiff.cpp`, `ddiff_minerva.cpp`,
    `ddiff_minerva_N.cpp` (now the full argument list of the `RunLog` names: they also gave `./neutest` and left out
    the ROMEA flag and the experiment), `neutestint.cpp`, `progs/Obstest/Ca_dens.cpp`, `Nucl_dens_profile.cpp`,
    `PxPz_exp_sample.cpp` (all said `./observables [Q2] [omega] [pm]`), `progs/onenucl/A_QE_cross.cpp`,
    `progs/pionestimate/pionestimate.cpp`, `progs/F2_largex_fit/F2_bootstrap.cpp` (it said `F2_fit` and two
    parameters, the code reads one). `Nucl_dens_profile.cpp` no longer says it is for Calcium (it takes the
    nucleus as an argument).
  - `CLAUDE.md`: the `gh pr edit` note (the Ubuntu `gh` 2.46 is too old; 2.106 from GitHub's apt repository works).

  **Still to do (not comments, or the author's call):**
  - **Hard-coded absolute paths of the author's machine:** `progs/deepsdis/fitscript.py` (`/home/wim/Code/trunk/bin/...`,
    `/home/wim/Calculations/Bonus/results/`), `progs/deepsdis/bonus_F2nfit.cpp`, `bonus_F2nfit_allangles.cpp`,
    `bonus_F2nfit_sys.cpp` (`outdir="/home/wim/Calculations/Bonus/results/parabfits..."`), `progs/onenucl/transp_eli.cpp`
    (`$HOME/Code/trunk/share`, a code path that does not exist: use `nucqcd::share_dir`); these programs fail or write
    elsewhere on another machine. Changes behaviour (paths), so separate.
  - **Library-wide Doxygen blocks that sit in one class header** (low priority): the `\defgroup DIS` text in
    `DIS/DeuteronCross.hpp` and the `\mainpage`/`\defgroup Glauber` text in `Glauber/AbstractFsiGrid.hpp`.
  - **Not covered by the crawl:** the Fortran files, `progs/attic`, `tests/`, the Python scripts apart from the usage
    check, and comments inside function bodies; a reader pass of the remaining `progs/deepsdis` files (many are variants
    of each other) is worth doing. Only the `\brief` lines of the `GPD_Partons` headers were checked (one wrong: `Deut_GPD_T_set.hpp`), not the rest of their comments.

- [x] **Integration method as an enumeration with one dispatch function** (added and done 2026-11-26, parts 1 and 2:
  #87 and the next PR; an optional part 3 below). The choice between `numint::cube_romb` and `numint::cube_adaptive` was
  an `int integrator` (1 or 2) in the constructors of the Glauber grids, `RhoTCross`, `PionTCross`, the EM/Weak QE
  classes and some programs, validated by `checked_integrator_flag` and then dispatched by `if (integrator == 1)
  cube_romb(...) else cube_adaptive(...)` at about 80 places in 15 files. `cube_adaptive` takes two more arguments
  (minimum and maximum number of evaluations); the call sites hard-code about 14 different pairs (`1E03,maxEval`,
  `5E02,maxEval`, `2E03,2E06`, `1E04,1E07`, ...); the integer is part of the grid cache file names (`integr2`,
  `intgr2`). Now, in the style of "Model selection" (#72):
  - an enumeration, `numint::IntegrationMethod { Romberg = 1, Adaptive = 2 }`
    (`Integrator/numint/integration_method.hpp`), with parse functions `nucqcd::integration_method(int)` and
    `(std::string)` ("romberg", "adaptive", "1", "2"; `Tools/Core/IntegratorFlag.hpp`), which replace
    `checked_integrator_flag` (removed);
  - a small options struct `numint::IntegrationOptions` (method, `epsabs`, `epsrel`, minimum and maximum evaluations) and
    one `numint::cube(options, f, lower, upper, result, neval[, dirty])` (`Integrator/numint/cube.hpp`) that hides the
    signature difference (`cube_romb` has no `dirty`/`compare` arguments for `double`/`dcomplex` results, so there are two
    forms); the evaluation caps of each call site are unchanged, only named in the options;
  - the integer stays in the cache file names (the 23 file names of the sandbox grids are identical before and after).
  - **Part 1 (#87):** the enumeration, options, `numint::cube`, the 13 Glauber grid classes and `DistMomDistrGrid` (19
    dispatch sites; the outer `if (integrator == 1 || integrator == 2) {...} else throw "integrator type not
    implemented"` blocks are gone, a closed enumeration makes them unreachable). A unit test checks that `numint::cube`
    returns bit-identical results and counts to calling the two routines.
  - **Part 2:** `RhoTCross`, `PionTCross`, `EMQECross`, `EMQEHadronCurrent`, `WeakQECross`, `WeakQEHadronCurrent` (15
    dispatch sites, 13 unwrapped `if/else throw` blocks; the `integrator` members and constructor arguments are an
    `IntegrationMethod`; `EMQECross` and `WeakQECross` had declared the member as a `double`), the programs and tests
    that construct them (a literal `2` became `IntegrationMethod::Adaptive`, the sdiff programs' `integrator` struct
    member and parameters changed type, `rhoset` converts its integer argument with `integration_method` and uses
    `numint::cube`), the archived programs that still build. `DQEinclusive` is untouched (its flag has five values,
    see part 3). Not changed on purpose: `RhoTCross`/`PionTCross` ignore their own `integrator` argument when they build
    their Glauber grids (a literal `Adaptive`, as before; and `rhoset` passes `Adaptive` to `RhoTCross` whatever its
    argument 9 says).
  - **Bug found on the way (part 1):** `progs/RhoTransparency/rhotest.cpp` passed the precision and the integrator to
    `GlauberDecayGridThick` in the wrong order (precision 2, integrator 0, which throws since #69); fixed by swapping
    (the program was not run, it needs a long Glauber computation).
  - **Part 3 (optional, to decide; the author asked how much work it is):** other algorithms in the enumeration. The
    library has Cuba's `cuhre`, `vegas`, `suave` and `divonne` (`Integrator/numint/numint2Cuba.hpp`: `void`, with
    `nvec`, `flags`, a seed, algorithm-specific parameters (key, `nnew`/`flatness`, `nstart`/`nincrease`/`nbatch`,
    `key1..3`, `maxpass`, `border`, ..., and a state file), and `nregions`, `neval`, `fail`, the error and `prob` as
    output arguments) and `cube_tr` (a trapezoid rule, used by nothing outside the library). They are called directly in
    `progs/neutrinoQE/sdiff_Q2*.cpp` (the "flux integrator": 0 vegas, 1 cuhre, 4 divonne), `progs/2bodymom`
    (`cuhre`, `suave`), `progs/benchmark/integrators.cpp` (`cuhre`) and `Knockout/DQEinclusive` (a flag 0 adaptive
    cubature, 1-4 Cuhre, Suave, Vegas, Divonne: a different numbering than 1 and 2 of the grids, so it needs a mapping).
    Work: new enumerators (`Cuhre`, `Vegas`, `Suave`, `Divonne`) with cases in `numint::cube` (about 150 lines in
    `cube.hpp`), the mapping of `fail` onto the return value, the algorithm-specific parameters with the defaults that
    the programs use today (key 11, flags 0, seeds) in an optional struct, a unit test per method with fixed seeds
    and loose tolerances (Vegas and Suave are stochastic), and moving the direct calls of the programs onto it (separate
    PR: `DQEinclusive` and the sdiff programs). Open choices: which algorithms are worth offering for the Glauber and
    Knockout integrands (the benchmark program of #61 can be extended to tell), whether the error estimate and `prob`
    should be returned, and the interplay of Cuba's forked workers (`CUBACORES`) with grids that are computed point by
    point inside loops. Not started.

## 7. Unbundle Cuba and cubature

This is a build-system change only, and the `unit_cuba` tests verify it.

- [x] **Cuba unbundled** (done 2026-11-25, #67). `Integrator/Cuba-4.2.1` and the `Integrator/Cuba`
  symlink (74 files) are removed. `Integrator/CMakeLists.txt`:
  - uses an installed Cuba if it finds one (`-DCUBA_ROOT`, `$CUBA_ROOT`, `~/.local`, the system), with
    `<Cuba/cuba.h>` provided from `<build>/Integrator/cuba-include`; it says when that Cuba is not in
    `~/.local` (a system package cannot carry our patch);
  - otherwise downloads Cuba **4.2.1** (SHA256 pinned; `-DCUBA_TARBALL=<file>` for machines without
    network, `-DCUBA_DOWNLOAD=yes` to ignore an installed one), applies the patch and builds it in the
    build tree with the same flags as before (`-std=gnu17`, serial `make`);
  - keeps the `cuba` INTERFACE target, so nothing that links it changed.
  4.2.1 was kept, not 4.2.2, so the code is the bundled one except for our patch (no values change); the
  numint wrappers need the 4.2 API.
- [x] **The patch** is `cmake/patches/cuba-4.2.1-nucqcd.patch` (`MaxL` → `MinL` in `src/suave/Fluct.c`,
  step 4, still in 4.2.2 (checked 2026-11-21), and the `lib:` target of `makefile.in` without the headers
  dependency, a local change the bundle also had). It was made by diffing the bundle against the
  pristine tarball: those two were the only differences (apart from the demo output and objects).
  `unit_cuba` tests Suave with a small flatness. Reporting the bug upstream is still to do.
- [x] **`scripts/install_cuba.sh [prefix] [tarball]`** installs the patched 4.2.1 (with `-fPIC`) into
  `~/.local` or another prefix, as `install_partons.sh` does for PARTONS.
- Tested: configure and build from scratch with the download, and with `scripts/install_cuba.sh` into a
  scratch prefix plus `-DCUBA_ROOT`: `unit_cuba` and `unit_numint` pass in both; the main build with the
  download passes all tests including `slow` (unchanged references). Not tested: a system Cuba package
  (none installed here), and a build without network access other than through `CUBA_TARBALL`.

- [x] **cubature unbundled (added 2026-11-23, done 2026-11-26): the adaptive integrator is now cubature 1.0.4, values change at the 1e-6 to 1e-2
  level for integrands with steps (author's decision), `CUBATURE_LEGACY_SPLIT` reproduces the old results bit for bit.** The bundled
  `Integrator/adaptive/cubature.{c,h}` (1435 lines) was Steven G. Johnson's pre-1.0 `adapt_integrate` from October 2010 with a `minEval`
  argument that upstream does not have (no convergence before `minEval` evaluations; the code passes 1e2 to 1e4). It is replaced by **cubature 1.0.4**
  (github.com/stevengj/cubature; **GPL-2-or-later, as the bundled copy was: no license change**), downloaded at configure time (SHA256 pinned,
  `-DCUBATURE_TARBALL=<file>` for machines without network, a wrong tarball is refused), patched with `cmake/patches/cubature-1.0.4-nucqcd.patch`
  and its one source file `hcubature.c` compiled into `libIntegrator` (and `IntegratorST` with `STATIC_BUILD`); `<cubature/cubature.h>` comes from
  `<build>/Integrator/cubature-include`. The unpacked and patched source is redone when the patch changes. A system `libcubature-dev` is not used:
  it has no `minEval`. `numint2adapt.hpp` (`cube_adaptive`) calls `hcubature_minEval` with `ERROR_INDIVIDUAL` (the same convergence test: every
  component below the absolute or the relative tolerance); the integrand returns `int` now. The old directory is gone; the unit test and the
  archived `inttest` use the new interface.
  **The patch.** (1) `hcubature_minEval`: `hcubature` with the minimum number of evaluations (`hcubature` itself is unchanged, `minEval` = 0).
  (2) The macro `NUCQCD_LEGACY_SPLIT` (CMake option `CUBATURE_LEGACY_SPLIT`, off) selects the rule of the 2010 version to choose the dimension to
  split a region; by default it is the rule of 1.0.4. Comparing the two sources, this rule is the only difference in the algorithm: 1.0.4 breaks
  near-ties of the fourth differences in favour of the longest side of the region (the fix of a hang with unlimited `maxEval`, upstream #14, which nucqcd
  never uses: it always passes a limit). For smooth integrands this changes nothing but the last digits; for integrands with steps (the Glauber integrals
  have cuts) it changes the adaptive subdivision, and so the results at the level of the integration error. **The first version of this item kept
  the 2010 rule so that nothing changed; the author then asked for upstream's behaviour and regenerated references (2026-11-26).**
  **What changed (default build vs master).** The 23 grid cache files of the slow tests: 1 identical (the ROMEA file), in the others the
  value changes per file (absolute, maximum over the file) are between 2.8e-6 and 5.4e-2, the median 4.8e-5: the thin off-axis `GlauberGrid` (tiny 4 x 4 x 3
  grid, 1e-3 requested) has entries off by up to 6e-2 and the 120 x 36 x 5 thick grid by up to 6e-3 per grid value, the other tiny grids by 1e-6 to 5e-4. In
  the references of the regression tests (**343 entries re-baselined**; 65 entries that moved by less than 1e-9 relative, last-bit noise, were restored to
  master's values, among them the 4 changed entries of `reg_nucleon_structure`, which does not use the integrator): `reg_glauber_grid` (A(e,e'p) of carbon)
  20 of 20 entries, up to 3.3e-3; `reg_glauber_variants` 254 of 392, up to 6.6e-3 on the real parts of values of order one (the largest relative changes,
  up to 4e-1, are on imaginary parts of 2e-3); `reg_glauber_points` 22 of 38, up to 9e-5; `reg_inclusive_dis` 21 of 28, F2D FSI pieces by about 4e-4 and the
  neutron Azz at x = 0.8 (a number of 7e-7) by 8%; `reg_deuteron_tagged` 6 of 23, up to 7e-7; `reg_distmomdistr_pw` 17 of 18, up to 2e-7; `reg_eep_planewave` 2 of 90, up to
  3e-6; `reg_weakqe` 1 of 22, 1.4e-9. The size of the grid changes is above the 1e-3 precision that the grids request, so on integrands with steps at least one of the two versions is not
  converged where it claims to be (in a scratch test of a step and a kink: the new rule was better on the kink, worse on the step: 6e-3 against 3e-3 at 2e4
  evaluations). Results of published calculations were computed with the 2010 rule: reproduce them with `-DCUBATURE_LEGACY_SPLIT=yes`.
  **Checked:** with the default build, `ctest -LE slow` 141/141 and `ctest -L slow` 11/11 against the regenerated references. With
  `-DCUBATURE_LEGACY_SPLIT=yes` and **master's original references**: 141/141 and 11/11, and the 23 grid cache files **byte-identical to master**; a
  scratch comparison of the old bundled `adapt_integrate` and the patched `hcubature_minEval` on 108 cases (6 integrands, 3 tolerances, 3 values of `minEval`,
  2 of `maxEval`, scalar and 6-component vector integrands): with the legacy rule **bitwise identical in the value, the error estimate and the number of
  evaluations in all 108**, with the default rule in 26. New unit tests: `hcubature` on the inttest integrand, and `hcubature_minEval` (minEval = 0 is
  `hcubature`, the number of evaluations exceeds `minEval`). The offline tarball path and the `STATIC_BUILD` configure were tested with a scratch configure.
  **Not tested:** the `STATIC_BUILD` build itself (only its configure); a build without network from a clean tree beyond the tarball configure.
  Still open: `pcubature` (excellent on smooth integrands, misses narrow peaks) and `hcubature_v` for parallel evaluation (see "Speed of the Glauber
  grids", step 6) are in the downloaded source but not used. **Benchmark (2026-11-23, `docs/integrator_benchmark.md`):** 1.0.4's `hcubature` gives practically the
  same results as the old code on smooth integrands, is not better on steps, and on a kink reports convergence while wrong (estimated 1e-5, true 6e-4).

## 8. Repository README

Write this last, after the refactors, so it describes the final layout and uses the new name.

- [ ] Cover:
  - what the code computes, and the libraries and how they depend on each other
  - dependencies and how to install them: GSL, a Fortran compiler, CMake; Cuba and Minuit2 are downloaded and
    built by CMake (or installed with `scripts/install_cuba.sh`, `scripts/install_minuit2.sh`); optionally LHAPDF
    with its PDF sets and PARTONS (`scripts/install_partons.sh`), Doxygen; ROOT only for the test `unit_vectors`;
    MPI only for two archived programs
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
    the cubature code that is downloaded and compiled into `libIntegrator` (`hcubature.c` of cubature 1.0.4) is GPL-2-or-later, so binaries that link them are covered by
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
   - [x] **4a: the types and their tests** (done 2026-11-23, #55). `Tools/Core` has `nucqcd::Vector2`,
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
   - [x] **4b: convert the code** (done 2026-11-23, #56). One scripted rename over the whole tree
     (libraries, programs, tests; line endings kept): `TVector3`, `TVector2`, `TLorentzVector`,
     `TRotation`, `TLorentzRotation` → `nucqcd::Vector3`, `Vector2`, `LorentzVector`, `Rotation`,
     `LorentzRotation`, their ROOT includes → `Tools/Core/{Vector3,LorentzVector,Rotation}.hpp`,
     and forward declarations → `namespace nucqcd { class ...; }` (100 files, ~980 lines). Done in one
     go rather than library by library because the Tools classes return these types to everything
     else. Besides the rename: `TSpinor.cpp` used `Double_t` (now `double`); two files had relied on
     the ROOT headers for `<cstring>`/`<iterator>`; the dictionary-era `#ifdef ROOT_Rtypes` blocks in
     `TSpinor.hpp`/`TLorentzQuaternion.hpp` (full includes when a ROOT header came first, forward
     declarations otherwise) now always forward-declare. `progs/GPD/diDVCS_gen_gridonly.cpp` links
     only ROOT and stays with ROOT's types (phase 6). No other method was missing. **The libraries
     now include no ROOT header.** Results: all tests pass with unchanged references. Re-baselining
     to compare all digits: every value computed without integration is bitwise the same; the
     adaptive 3D integrals of `reg_eep_planewave` (A(e,e'p) plane wave, precision 1e-5) moved by
     1e-15 to 3e-11 relative (54 of the 72 `_pw` values; the 18 `_free` values, without integration, are unchanged) and three
     `reg_glauber_grid` values by ≤ 7e-16: last-bit differences where the compiler contracts
     differently into fused multiply-adds at inline call sites (see 4a), slightly amplified by the
     adaptive subdivision. (`reg_inclusive_dis` and `reg_nucleon_structure` re-baseline differently
     from their references on this machine on master too, at 1e-15: not from this change.) Cached
     Glauber grids stay valid (last-bit differences only). `namespace nucqcd` for the other classes
     (code layout item 4) is still to do.
5. [x] **Build system** (done 2026-11-23, #58). ROOT is optional: `ROOT_ENABLED` (ON by default, like
   `LHAPDF_ENABLED`: used when found), `find_package(ROOT ... QUIET COMPONENTS Minuit)` (Minuit2 is
   not used by anything). The libraries no longer link ROOT (`ToolsCore` linked `root_libs`, and
   through it everything else). Programs and tests that need ROOT list `root_libs` in their `LINK`;
   `nucqcd_add_program(s)` and `nucqcd_add_test` skip them when ROOT is not found, and the
   configuration lists them: the 9 fit programs (`TMinuit`/`TFitter`), the 5 `diDVCS*` event
   generators (phase 6 removes ROOT from them) and `unit_vectors`. `wftest` included ROOT's
   `TVector.h` without using it. `RunLog` reports the ROOT version only when ROOT was found.
   `scripts/install_partons.sh` no longer needs ROOT for its test build. Checked: a build with
   `-DROOT_ENABLED=no`, `BUILD_PROGS` and `PARTONS_ENABLED` builds everything else, no library
   links ROOT (`readelf`), and the tests pass there too (after the `RunLog` fix, step 4).
6. [ ] **Programs.** What is left of ROOT (inventory 2026-11-23):
   - [x] the `diDVCS*` event generators in `progs/GPD` (done 2026-11-23, #59; author's
     decision): no ROOT left. `progs/GPD/event_random.hpp` replaces `TRandom3` and
     `TF1::GetRandom`: one `std::mt19937_64`, seeded differently for every run as ROOT's
     `SetSeed(0)` did (from `std::random_device`), or from the environment variable `NUCQCD_SEED`
     to repeat a run; the seed goes to stderr (stdout, with the events, is unchanged).
     `uniform(a, b)` and `breit_wigner(mean, gamma)` use ROOT's formulas (`TRandom::Uniform`,
     `TRandom::BreitWigner`, with a uniform number in (0, 1] as `TRandom3::Rndm`); the sin²θ and
     cos²θ decay angles are sampled exactly by accept-reject (`TF1::GetRandom` inverted a binned
     CDF with interpolation). Checked on 10⁶ samples: ⟨sin²θ⟩ = 0.7499 and ⟨cos²θ⟩ = 0.7501
     (exact 0.75), the histogram matches the exact bin fractions, Breit-Wigner median and quartiles
     as expected, `NUCQCD_SEED` reproducible. The event samples therefore differ from ROOT's event by
     event but not in distribution. `TMath` → `<cmath>`/`std::min`/`std::max`, `TMath::Pi()` and
     `DegToRad()` → `nucqcd::pi`, `nucqcd::deg_to_rad`; `Double_t` → `double`; in
     `diDVCS_gen_gridonly` ROOT's vectors → `nucqcd::` ones (it now links `ToolsCore`). `TFile`,
     `TTree`, `TF2` were included but unused (the events go to `cout`); `diDVCS` and `diDVCS_Q2`
     only had the includes. Not run end to end: on this computer PARTONS stops with "Info file not
     found for PDF set 'UNDEFINED'" (see the PARTONS setup lead in step 4).
   - [x] **the fits** (done 2026-11-25, #73; `progs/F2_largex_fit`: `F2_fit`, `F2_bootstrap`, `EMC_fit`;
     `progs/deepsdis`: `deepsfit`, `bonusscatterfit`, `bonusfit_rig`, `bonus_F2nfit(_sys|_allangles)`): they
     used ROOT's `TMinuit` through `TFitter` (and peeked into TMinuit's status strings), plus `TMatrixT`
     for a covariance matrix that was only printed. They now use the **standalone Minuit2** (the ROOT sources
     of Minuit2 with their own CMake project, github.com/GooFit/Minuit2 tag v6-40-02 = ROOT 6.40.02,
     LGPL-2.1, 6 s to build, static libraries): `cmake/Minuit2.cmake` uses an installed one
     (`-DMINUIT2_ROOT`, `$MINUIT2_ROOT`, `~/.local`, the system; not ROOT's own libMinuit2, which needs
     ROOT; `scripts/install_minuit2.sh`) or downloads and builds the pinned release (SHA256;
     `-DMINUIT2_TARBALL`, `-DMINUIT2_DOWNLOAD`, `-DMINUIT2_ENABLED=no` skips the fits). Instead of
     porting each program to Minuit2's interface (ten copies of the same ~150 lines), `Fitting/MinuitFitter.hpp`
     (`nucqcd::MinuitFitter`) keeps the part of the `TFitter` interface they used (external parameter
     numbers, limits `0,0` = none, the old `fcn(npar, gin, f, par, iflag)`, the commands `MIGRAD`, `SIMPLEX`,
     `MINOS`, `SET PRINT`, `SET STR`, `SHO COV`, ...); the programs changed by a script (includes, the
     creation of the fitter, `static_cast<TFitter*>(...)->GetMinuit()->fCstatu.CompareTo("CONVERGED ")` →
     `Converged()`, the `fCovmes` tests → `CovarianceAccurate()` / `CovarianceForcedPosDef()`,
     `TMatrixT` → `std::vector`, `TMath::Abs` → `std::abs`, `Int_t`/`Double_t`). No program links ROOT any
     more. **Differences** (Minuit2 is a rewrite, not the same algorithm; see the header): `IMPROVE` does
     nothing (no equivalent), `SET RAN` and `SET LIM` are ignored, a fixed parameter reports error 0 as in
     TMinuit (Minuit2 reports its step size, which made the "fix parameters at their limits" loop of the
     programs endless: found when the first test run wrote a 6 GB log). **Checked** against the ROOT binaries
     on the same data: `EMC_fit` gives the same fit to every printed digit (FCN 88.1853, m1 = 0.103422,
     error 1.92254e-3, identical data/fit table). `F2_fit stat noSLAC wHT xcut` finds the same parameter values
     for the four parameters that stay free (1.00597, 1.8468, -6.96952, -0.08836), with parabolic errors
     about 1.4 times larger; it differs in the status of the error matrix: Minuit2 calls it accurate, so
     the program fixes the two parameters that sit at their limits (5 and 10), where the TMinuit run did not
     (and with fixed parameters the summary lines at the end of these programs, which pass `nFreePars` values
     to the FCN, evaluate garbage: the chi-square printed there, 173787, is not the fit's; a pre-existing flaw,
     the fixed parameters are not in the array, not fixed here). **Not checked:** `deepsfit`,
     `bonusscatterfit`, `bonusfit_rig`, `bonus_F2nfit*` and `F2_bootstrap` (they need the DEEPS and BONuS
     data and take long; they compile and link) and `parabfit`, which was never in the CMake lists (it
     compiles with the port). Not covered by tests: the fits are programs, `unit_minuit` is not written.
   - `RunLog` in the other programs and their `// ` headers (step 11 part 2) touch the same files.
   With the generators and the fits done, ROOT is only needed for `unit_vectors` (its comparison with ROOT's vector classes).

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

1. [x] **Design and implementation** (done 2026-11-22, #49). Decided with the author: the log goes to
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
2. [ ] **Roll-out.** Done (2026-11-22, #51): `Bookkeep` replaced by `RunLog` in the 8 programs that used it
   (`poldeut`, `EMT_deut`, `b1_calc`, `pionT_calc12`, `Ca_dens`, `Nucl_dens_profile`,
   `PxPz_exp_sample`, `conv_V_test`; their "exec name" entries dropped, since `RunLog`'s names start at
   argv[1]; `poldeut` and `EMT_deut` had too few names and read past the array), and `Bookkeep`
   removed from `Utilfunctions`. **Done (2026-11-25, #64):** `RunLog` is the first statement of `main` in the other 67 programs outside
   `progs/attic` (77 programs in all; the argument names were read from how each `argv[i]` is used, units only
   where the code shows them; programs without arguments pass none, and four `int main()` programs now take
   `argc, argv`). None of them printed `//` header lines, so no output text changed, and none writes the
   log to a data file: it goes to stdout, before and after the program's own output. Checked that the
   consumers in the repository skip `#` (`np.loadtxt` in `bonusreadin.py`, `gridplot.py`, the event
   parsers of `2bodymom`); `progs/deepsdis/fitscript.py` (Python 2, hard-coded paths, parses the stdout of
   `bonusextrapolate` and `parabfit` by token position) now drops `#` lines first. Not checked: scripts
   outside the repository that read the stdout of `transp_eli` (it writes a `#__PYHDR_START__` JSON header; the `# `
   lines now come before it), `transp`, `JlabPAC16_calc` and the other `onenucl` programs. Built with
   `PARTONS_ENABLED`; the GPD programs were compiled but not run.
