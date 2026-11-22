# History of the repository

This file records where the code came from and what happened to it, so that old commits, branches and
results can still be traced. The current cleanup plan is in [ROADMAP.md](ROADMAP.md).

## Hosting

| Period | Where | Notes |
|---|---|---|
| 2012-03 – 2017 | Subversion (SVN) | First commit `bf133f9`, 2012-03-06. Working copies were checked out as `~/Code/trunk`, and some old paths in `progs/` still refer to it. |
| 2017 – 2022 | `github.ugent.be/wcosyn/physics-code` (Ghent University GitHub Enterprise) | The SVN repository was converted to git with its full history. The first git merge is UGent PR #1 (`cmake_prog`, 2017-11-10). |
| 2022-09 – now | [github.com/wcosyn/nucqcd](https://github.com/wcosyn/nucqcd) | The UGent account ended in 2022 after a change of institution, so everything was migrated here, still as `physics-code`. Renamed to `nucqcd` ("QCD in nuclei") on 2026-11-20; GitHub redirects the old URLs. |

Traces of the earlier hosts in the git history:

- **SVN:** there are no merge commits before 2017-11 (SVN merges became ordinary commits), except for
  `b7be145` (2014-10-13), a converted merge of the SVN branch `LDRD-branch-0.01`. The lightweight tag
  `LDRD-0.01` on that commit is the converted SVN tag. The conversion kept no `git-svn-id` lines,
  so SVN revision numbers cannot be recovered from the commits.
- **UGent GitHub:** twelve merge commits from 2017-11-09 to 2022-06-16 have messages like
  `Merge branch 'master' of https://github.ugent.be/wcosyn/physics-code`. They come from pulling
  between machines. PR and issue numbers mentioned before 2022 (e.g. "pull request #1") refer to the
  UGent server, which is no longer accessible.
- The git history was pushed to github.com as it was (the UGent merge commits are still there), so
  commit hashes from the UGent years are unchanged, and hashes cited in papers or notes still resolve.

## Contributors

- **Wim Cosyn** (2012 – now): the libraries and most programs.
- **Camille Colle** (2013 – 2016): semiclassical single-charge-exchange (SCX) grids and the
  2-nucleon knockout / pair c.m. momentum programs (`progs/2bodymom`, the `progs/camille` programs,
  `ClassGridThick_SCX` in Glauber). Camille also updated the bundled Cuba (then 3.3) and `numint`.
- **Philip Hacker** (2016-11 – 2017-01): neutrino QE scattering, including the weak knockout
  currents in `Knockout/WeakQE*`, the hydrogen cross section, and comparisons with MINERvA and MiniBooNE
  (`progs/Minerva`, later folded into `progs/neutrinoQE`).
- **Alan Sosa** (2024-07 – 2024-09): b1–b4 tensor structure functions in `progs/deepsdis/b1_calc.cpp`.
  This work was not merged; it is kept in the tag `archive/Alan` (see below).

## When the main parts were added

| Part | First commit |
|---|---|
| Tools (MePhys), Glauber, Knockout, DIS, `progs/deepsdis` | 2012-03 (initial import) |
| `progs/2bodymom` (2N knockout) | 2014-08 |
| `progs/poldeut` (light-front deuteron, spectator tagging) | 2016-03 |
| CMake build for all programs, `progs/neutrinoQE`, `progs/RhoTransparency` | 2017-11 |
| `progs/GPD` (deuteron GPDs, diDVCS) | 2018-01 |
| `progs/F2_largex_fit` and `DIS/evolution` (QCD evolution, Fortran) | 2018-06 |
| `GPD_Partons` (PARTONS interface) | 2019-03 |
| Test suite (`tests/`), cleanup roadmap | 2026-11 |

## Archived branches

The unmerged branches were archived as tags on 2026-11-20 and then deleted
([ROADMAP.md](ROADMAP.md), step 2). To work on one again, run `git switch -c <name> archive/<name>`.

| Tag | Branch content |
|---|---|
| `archive/LDRD-branch-0.01` | 2014 SVN branch for the LDRD project (He3 tagging): copy constructors and kinematics changes. It also removed most of the repository at the time. |
| `archive/evolution_integration` | 2018 branch for the QCD evolution. Everything in it is also on master: `std::isnan` in `F2_bootstrap.cpp` (same as `f246735`) and `evolution_speedy.f90` (from `fb67bb0`). |
| `archive/Alan` | 2024: extraction of b1–b4 and new structure-function outputs in `b1_calc.cpp`. Its `cos 2θq` fix of the Azz angular factor was applied to master separately. |

The merged branches `fix-ww-integral`, `pionT`, `brandon` and `partons_integration` were deleted.

## Programs behind published or finished calculations

This is not a complete list. Add to it when you trace a result back to the program that produced it.

- **`progs/F2_largex_fit`** (`F2_fit`, `F2_bootstrap`), 2018: fit of the x > 1 F2 of carbon
  (JLab data from Fomin et al., optionally with SLAC data), with and without higher twist and with
  x cuts, plus bootstrap error estimates. The leading-twist part is evolved from Q² = √18 GeV²
  with `DIS/evolution`: `evolve_ns_adam` in `wrappers.f90` calls the non-singlet `evolve_ns_speedy`
  from `evolution_speedy.f90`, optimized for large x. These two programs are the only users of
  `libEvolution`. Its GPD and singlet evolution routines are not called from C++.
- **`progs/deepsdis`**: DEEPS and BONuS spectator-tagging analyses, including the BONuS paper plots
  (`bonustest`), and b1 calculations compared with Kumano's model (`b1_calc`, `b1_calc_comp`).
  The experiment data are in `share/expdata/`.

## Unfinished code

- **A(e,e'NN) with correlated pairs** (`Knockout/DoubleNModel`, `Knockout/DoubleNCross`,
  `Tools/pair`, `progs/twonucl`), in the initial import of 2012-03 and not developed further. The
  Vegas logs committed with `progs/twonucl` (`corrpn`, `uncorrpn`, `corrtest`, `uncorrtest`,
  removed in 2017) show that the 6D integral did not converge at its ¹²C(e,e'pn) kinematics, and no
  results came from it. Several bugs were fixed in 2026-11 ([ROADMAP.md](ROADMAP.md), step 4).
