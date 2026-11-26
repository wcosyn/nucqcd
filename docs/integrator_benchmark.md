# Integrator benchmark (2026-11-23)

Roadmap step 6 asks whether `rombergerN` and `numint::cube_romb` can be removed, after a benchmark
that checks the author's memory: that the adaptive integrators did better, especially on
multidimensional integrands with non-trivial boundaries or cuts. This file records the benchmark.
Program: `progs/benchmark/integrators.cpp` (`bin/integrators A`, `B`, `C`), built with `BUILD_PROGS`.
Machine: 12 cores, Release build (`-O3 -march=native`), one core per run (Cuba with `CUBACORES=0`).

The integrators:

| Name | Where | Method | Stops when |
|---|---|---|---|
| `rombergerN` | `Tools/Core/Utilfunctions.hpp` | 1D Romberg, nested per dimension | relative change < `acc`, or absolute change < `acc`·1e-7, or the result is < 1e-5 of a running `estimate`; at most `max` steps (code: 3 to 8); returns 0 for intervals shorter than 1e-6 |
| `cube_romb` | `Integrator/numint/cube_romb.cpp` | Romberg on the full N-dim grid, all dimensions refined together | converged, or 12 refinement levels (`EMAXITER`); no limit on evaluations |
| `cube_adaptive` | `Integrator/numint/numint2adapt.hpp` → `hcubature` of cubature (until 2026-11 the bundled `Integrator/adaptive/cubature.c`, S. G. Johnson, 2010, with a local `minEval`; now the downloaded 1.0.4 patched to behave the same, see ROADMAP step 7) | h-adaptive Genz-Malik cubature | all components below the absolute *or* relative tolerance, or `maxEval` reached; **returns 0 (success) in both cases**, and the wrapper discards the error estimate |
| Cuhre | Cuba 4.2.1 via `numint::cuhre` | globally adaptive cubature | converged, or `maxEval` (`fail` = 1) |

## Part A: test integrands with known values on [0,1]³

Requested relative precisions 1e-3, 1e-5, 1e-7; `maxEval` 2·10⁷ for `cube_adaptive` and Cuhre;
`rombergerN` with 3 to 8 steps per dimension, as in the code. Relative error with respect to the exact
value (for the Glauber-like case: nested GSL Gauss-Kronrod, 1e-11). "flag": the integrator's own
failure report. Condensed (full table: run `bin/integrators A`):

| Integrand | Prec | `rombergerN` | `cube_romb` | `cube_adaptive` | Cuhre |
|---|---|---|---|---|---|
| smooth exp(−r²) | 1e-3 | 5e-7 (729) | 8e-10 (4913) | 2e-7 (166) | 2e-11 (381) |
| | 1e-7 | 7e-13 (4·10⁴) | 5e-13 (4·10⁴) | 2e-11 (3664) | 2e-11 (381) |
| step: inside the unit sphere | 1e-3 | 4e-4 (1.4·10⁵) | 3e-4 (4·10⁴) | **1.5e-3 (2·10⁷, flag 0)** | 1.3e-4 (2.6·10⁵) |
| | 1e-5 | 2.3e-5 (1.3·10⁷) | 2.6e-8 (1.1·10⁹, 16 s) | **1.5e-3 (2·10⁷, flag 0)** | 1.6e-6 (2·10⁷, flag 1) |
| | 1e-7 | 2.8e-5 (1.3·10⁷) | **0.98 (overflow bug, flag `EMAXITER`)** | **1.5e-3 (flag 0)** | 1.6e-6 (flag 1) |
| kink: (1−r²) inside the sphere | 1e-3 | 1.4e-5 | 1e-4 | 2e-5 (2·10⁷, flag 0) | 1.8e-4 |
| | 1e-7 | 1.3e-7 (4·10⁵) | 1.9e-8 (1.7·10⁷) | **2e-5 (2·10⁷, flag 0)** | 1e-5 (flag 1) |
| narrow peak (σ = 0.02) | 1e-3 | 1.3e-4 | 6e-15 | 1.6e-4 (8284) | **0.99 (misses the peak, flag 0)** |
| | 1e-7 | 3.1e-7 | 5e-14 | 5e-10 (8·10⁵) | 6e-5 |
| simplex mapped onto the cube | 1e-7 | 2e-11 | 2e-11 | 2e-10 | 1e-11 |
| Glauber-like: Woods-Saxon over b, φ, z with a cut | 1e-3 | 8e-5 | 8e-7 | 1.5e-4 (4522) | 2e-5 |
| | 1e-5 | 7e-8 (1.4·10⁶) | 1.5e-6 (2.7·10⁵) | 1.1e-6 (2·10⁷, flag 0) | 1.3e-6 (10⁵) |
| | 1e-7 | 3e-8 (6.7·10⁶) | 1.1e-10 (1.4·10⁸, 5 s) | **1.1e-6 (2·10⁷, flag 0)** | 9e-8 (7.8·10⁶) |

(numbers in brackets: evaluations)

Upstream cubature 1.0.4 (`hcubature`, `pcubature`; compiled separately for this comparison, not part of
nucqcd) on the same integrands, `maxEval` 2·10⁷:

- `hcubature` is practically the same as the bundled 2010 code on smooth integrands (identical results
  for the peak and the simplex). On the step it stalls at 3.5e-3 (its own estimate: 3.6e-4); on the
  kink it **reports convergence while being wrong** (at 1e-5: estimated 1e-5, true 6.3e-4, after
  3.2·10⁵ evaluations). It also returns 0 when `maxEval` is reached.
- `pcubature` (p-adaptive Clenshaw-Curtis) is excellent on smooth integrands (1e-15 on the simplex) and
  better than `hcubature` on the kink (3e-9 at 1e-7), but **misses the narrow peak completely** (27
  evaluations, error 1.0).

## Part B: the real integrands, through the classes' `integrator` flag

### Glauber grid (`GlauberGridThick`, carbon, 12×6×3 grid, 1500 MeV proton)

Max. relative deviation of 4 grid points (FSI factor with and without CT) from `cube_adaptive` at
1e-6; every grid computed from scratch (empty `NUCQCD_GRID_DIR`).

| Proton | Method | Prec | Deviation | Time |
|---|---|---|---|---|
| along z | `cube_adaptive` | 1e-3 | 6.6e-5 | 0.4 s |
| | `cube_romb` | 1e-3 | 1.1e-4 | 580 s |
| | `cube_adaptive` | 1e-5 | 3e-7 | 12 s |
| | `cube_romb` | 1e-5 | — | > 2 h (stopped) |
| θ = 0.5 rad | `cube_adaptive` | 1e-3 | 1.1e-3 | 76 s |
| | `cube_romb` | 1e-3 | 1.3e-3 | 16 s |
| | `cube_adaptive` | 1e-5 | 0 (same as the reference) | 123 s |
| | `cube_romb` | 1e-5 | — | > 600 s (stopped) |

The grid code caps `cube_adaptive` at 2·10⁶ (5·10⁶ with CT) evaluations per grid point. Raising the cap
100×: the along-z reference at 1e-6 takes 2850 s instead of 81 s (the deviations above stay the same,
so the capped reference is good to ~1e-7), and the off-axis reference does not finish in 2 hours. So
**off-axis grids are limited by the evaluation cap, not by the requested precision**: requesting 1e-5
or 1e-6 gives the same grid, whose accuracy is about 1e-3 (the deviation of both methods at 1e-3). This
matches the author's memory that off-axis grids are slow (see "Speed of the Glauber grids" in the
ROADMAP).

### A(e,e'p) cross section (`EMQECross` → `EMQEHadronCurrent`, carbon, Q² = 1.8 GeV², p_m ≈ 150 MeV/c)

The amplitude integral is a 3D integral of 12 complex amplitudes (24 components). `cube_adaptive`
needs *all* components converged, so it runs until `maxEval` (the argument of `getDiffCross`) whatever
the requested precision: the results at 1e-3, 1e-5 and 1e-6 are identical and take the same time; the
absolute tolerance 1e-8 is not the reason (removing it changes nothing). **The accuracy is set by
`maxEval`**: plane wave and RMSGA, relative deviation from `maxEval` = 2·10⁷:

| Shell | FSI | 2·10³ | 2·10⁴ | 2·10⁵ | 2·10⁶ | time at 2·10⁴ |
|---|---|---|---|---|---|---|
| 0 (1s½) | plane wave | 7e-4 | 1.4e-5 | 3.6e-7 | 7e-8 | 0.04 s |
| 0 | RMSGA | 1.1e-3 | 1.9e-5 | 9.3e-6 | 2.4e-6 | 0.05 s |
| 1 (1p3/2) | plane wave | 2.3e-3 | 3.0e-5 | 1.6e-7 | 8.5e-8 | 0.09 s |
| 1 | RMSGA | 1.1e-3 | 1.0e-4 | 1.7e-5 | 3.4e-6 | 0.10 s |

The programs mostly pass 2·10⁴ (`progs/onenucl`, `progs/Obstest`), i.e. 1e-5 to 1e-4.
`cube_romb` on the same cross section (plane wave, relative to `cube_adaptive` with 2·10⁶): at 1e-3
1.6e-6 / 1.7e-5 (0.4 s / 1.2 s), at 1e-5 6.8e-7 / 8.7e-8 (3 s / 40 s) for shells 0 / 1.

## Conclusions

1. **The memory "adaptive is better for cuts" is not confirmed by these tests.** For integrands with a
   step or a kink inside the domain, the h-adaptive cubature (`cube_adaptive`, and upstream `hcubature`
   too) is the *least* reliable: it stalls at `maxEval` with errors of 1e-5 to 1.5e-3 and reports
   success, or (1.0.4) claims a convergence it has not reached. Romberg (`cube_romb`, `rombergerN`)
   converges there, at a high price. It may be that older versions of the code or other integrands
   behaved differently.
2. **Speed:** for smooth integrands and narrow peaks, `cube_adaptive` is by far the cheapest; for the
   along-z Glauber grid it is 1400× faster than `cube_romb` at the same accuracy. `cube_romb`'s cost
   explodes with precision in 3D (all dimensions refined together), and it has no evaluation limit.
   For the off-axis grid at 1e-3, `cube_romb` was faster (16 s against 76 s, similar accuracy). For the
   A(e,e'p) amplitudes the two are comparable at equal accuracy (`cube_romb` at 1e-3: 1.6e-6 in 0.4 s;
   `cube_adaptive` with 2·10⁵ evaluations: 3.6e-7 in 0.3 s; shell 0, plane wave).
3. **Cuhre** is good on smooth integrands and steps, but misses narrow peaks at loose tolerance (as does
   `pcubature`) without reporting it.
4. **Removing `cube_romb` is not recommended** on this evidence: it is the robust choice for
   discontinuous integrands, and competitive for the A(e,e'p) amplitudes and off-axis grids at modest
   precision (it needs the overflow fix #60 for high precision in 3D). `rombergerN` (flag 0) is
   mostly commented out in the classes; its ad hoc stopping rules (absolute 1e-7·acc, 1e-5 of a running
   estimate, zero for intervals < 1e-6) make its precision hard to state; the hard-coded uses
   (`InclusiveCrossRes`, `RhoDeuteron`, `RhoTCross`, `OneGlauberGrid`, `DistMomDistrGrid`) could move to
   `cube_romb` or `cube_adaptive` case by case, with value changes to be checked.
5. **The requested precision often does not control the accuracy.** The evaluation caps do: `maxEval`
   in `getDiffCross` for A(e,e'p), and the hard-coded 2·10⁶/5·10⁶ in the Glauber grids. Callers get no
   warning, because the `numint` wrapper of `cube_adaptive` returns 0 and drops the error estimate.
   Worth fixing independently of the integrator choice (see the ROADMAP).
6. **Unbundling cubature (step 7):** moving to upstream 1.0.4 does not improve accuracy (worse on kinks);
   it is a maintenance question only. Done 2026-11-26: the build downloads 1.0.4 (its rule to choose the dimension to split is the difference measured above,
   it changes results of integrands with steps by 1e-6 to 1e-2; `-DCUBATURE_LEGACY_SPLIT=yes` selects the 2010 rule and reproduces the old results bit for bit). `pcubature` could be considered for smooth integrands.

Bugs found: `cube_romb` 32-bit overflow beyond 10 refinement levels in 3D (#60).
