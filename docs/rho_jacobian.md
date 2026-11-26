# The t and z integrations of the rho cross section (`RhoTCross`)

This note derives the Jacobian that connects the two ways of computing the integrated incoherent rho
electroproduction cross section in `Knockout/RhoTCross.cpp`: `getCrosst` (integrand per `nu` and `t`, used by
`rhoset` with `torz` 1) and `getCrossz` (per `nu` and `z`, `torz` 0). Both must give the same integral over the same
region (`t` in -0.4..-0.1 GeV², `z` in 0.9..1). They did not until 2026-11 (ROADMAP step 4, "Rho transparency: t- and
z-integrated results disagree"). Everything is in GeV. The derivation is mine (inferred from the code, not from a
paper): **marked for the author to confirm**, in particular the reading of what `C0` and `massi` mean.

## Summary: what was wrong and what was changed

- **Before (code up to 2026-11):** `getCrosst` and `getCrossz` differ in their prefactors (`1/(2q)` and `nu/q`, ratio
  `2 nu`) and in `getfrontfactor`, which multiplied the t version by 1 and the z version by `massX`,
  `sqrt(E_X² - px² - (q-pz)²)`, the invariant mass of the whole residual nucleus system (about `M_A`). For the two
  integrals to agree, the ratio of the z and t weights has to be `|dt/dz|`, so `massX` had to equal `|dt/dz|/(2 nu)`.
- **The mistake:** `|dt/dz|/(2 nu)` is about 0.7-1.7 GeV (`massi`, the off-shell energy of the struck nucleon, at pm = 0),
  not `M_A`. `massX` is right only if the whole nucleus is the struck object. So the z integrals were too large by about
  `M_A/massi` (12 for carbon), the "goes like A" of the 2012 commit messages. It cancels in the Glauber/plane-wave ratios.
- **The correction:** in the z case, `getfrontfactor` returns the exact `|dt/dz|/(2 nu)` of equation (4) below instead of
  `massX`. The t version is unchanged.

## Which results were affected

- **Published transparencies: not affected.** The correction changes only the z branch of `getfrontfactor`
  (`torz` = 0, used by `getCrossz`). For the t branch (`torz` = 1, used by `getCrosst`) the function returns the same
  value as before, `λ/m_N² e^{6t}`. According to the README of the paper's calculation directory, the last
  calculations (`rhoset`, with the cuts) integrated over t, so they did not pass through the corrected factor.
- **Why the transparencies looked fine even in the z version:** the transparency is the ratio of the RMSGA to the PWA
  cross section, integrated over the same ranges, and the wrong factor multiplied both. The factor `massX` was
  also nearly constant over the integration region, so it cancelled to a good approximation. How much the z-version
  transparencies move with the corrected, pm-dependent factor was **not quantified** (it needs a full run with the
  Glauber distributions in both versions, see the ROADMAP item).
- **Affected:** absolute z-integrated cross sections (`rhoset` with `torz` 0, `rhosingle` with the z option), which
  were too large by about `M_A/massi`.

## 1. Setup

- Photon: energy `nu`, momentum `q = sqrt(Q² + nu²)` along z, `t = (q_γ - p_ρ)²`.
- Struck nucleon with missing momentum `pm` (polar angle θ, azimuth φ with respect to q). Its energy is the
  *off-shell* energy `massi = M_A - M_(A-1) - E_exc` (the residual A-1 system is on shell), so `massi` is about
  0.922 GeV for carbon, not the on-shell `m_N`. This is the `massi` of `klaas_rho_t`/`klaas_rho_z`.
- `C = P_N + q_γ` with components `C0 = nu + massi`, `Cz = q + pm cosθ`, `Cx = pm sinθ cosφ`, `Cy = pm sinθ sinφ`;
  `s = C²` (the invariant mass squared that enters the elementary cross section).
- Rho: `E_ρ = z nu`, momentum `(px, 0, pz)`; the azimuth of the rho is fixed to define the x axis, the azimuth φ of
  `pm` takes its place.
- The final nucleon is on shell with the vacuum mass `m_N`:

      f(z, t) = (C - p_ρ)² - m_N² = s + m_ρ² - m_N² - 2 (C0 E_ρ - Cz pz - Cx px) = 0.        (1)

  This is the delta function `delta(f)` that the code removes by solving a quadratic for `pz` (two roots, hence
  the two branches in `klaas_rho_t`/`klaas_rho_z`).
- `t = -Q² + m_ρ² - 2 E_ρ nu + 2 q pz` gives `pz(z, t) = (A + z nu²)/q` with `A = (t + Q² - m_ρ²)/2`, and
  `px² = E_ρ² - pz² - m_ρ²`.

## 2. Phase space and the two ways to remove the delta function

The rho phase space in these variables is

    d³p_ρ / (2 E_ρ) = nu/(4 q) dz dt dφ_ρ                                                   (2)

(use `d³p/E = p dE dΩ` and, at fixed `E_ρ`, `dt = 2 q p_ρ dcosθ_ρ`). The cross section contains
`d³p_ρ/(2E_ρ) delta(f) G`, with `G` everything else. Integrating over `z` at fixed `t` (the t version) or over `t` at
fixed `z` (the z version) removes the delta function with a different Jacobian:

    W_t = nu/(4q) 1/|df/dz|,        W_z = nu/(4q) 1/|df/dt|,        W_z/W_t = |df/dz|/|df/dt| = |dt/dz|_{f=0}.      (3)

So for the two integrals to be equal, the weights that multiply `G` in the code must have the ratio `|dt/dz|`.

## 3. The Jacobian

With `d pz/dt = 1/(2q)`, `d pz/dz = nu²/q`, `d px/dt = -pz/(2 q px)`, `d px/dz = nu² (z - pz/q)/px`, differentiating (1)
gives `|df/dz| = 2 nu kappa` and `|df/dt| = tau/q` with the two quantities (kappa has the dimension of a mass)

    kappa = |C0 - Cz nu/q - Cx nu (z - pz/q)/px|,      tau = |Cz - Cx pz/px|,

and therefore, with `K := |dt/dz| / (2 nu)` (also a mass),

    K = q kappa / tau = q |C0 px - Cz nu px/q - Cx nu (z - pz/q)| / |Cz px - Cx pz|.         (4)

This is what `RhoTCross::getfrontfactor` returns, times the elementary-cross-section factor, in the z case.

**Special case pm = 0** (`Cz = q`, `Cx = 0`): `kappa = C0 - nu`, `tau = q`, so `K = kappa = massi`, the off-shell energy of the struck nucleon (0.922
GeV for carbon). **Not** the on-shell `m_N` = 0.938 GeV, and not `M_A`. For pm > 0 it depends on pm and on the angles:
in a random sample (pm up to 0.3 GeV, nu 2.5-3.5 GeV, Q² 1.5 GeV²) it ranged from about 0.7 to 1.7 GeV.

Before 2026-11 the code used `massX = sqrt(E_X² - px² - (q-pz)²)`, the invariant mass of the whole residual nucleus
system (about `M_A`), in this place. That made the z integrals too large by about `M_A/K`, about 12 for carbon (the
"goes like A" of the 2012 commit messages). For a point nucleus (pm = 0 with the nucleus as the struck object,
`massi = M_A`) the old factor is right, which is probably how it came about; the analytic pair that is commented out in
`progs/RhoTransparency/rhoset.cpp` (`0.5 Erho/q` for t, `nu M_A Erho/q` for z) is this case.

## 4. Comparison with the paper, and what the "4 kappa" was

The cross section implemented in `RhoTCross` is Eq. (13) of W. Cosyn and J. Ryckebusch, Phys. Rev. C 87, 064608
(2013), arXiv:1301.1904, in the equation numbering of the arXiv version (formulas checked against the author's TeX
source `rhopaper_final.tex`, equations `eq:A`, `eq:factorize`, `eq:distormomentum`, `eq:elem`, `eq:cross`):

    d sigma^{eA}/(dE_e' dOmega_e' dt dphi_rho) = sum_alpha ∫ d³pm  alpha_EM/(4 (2 pi)²) · E_e' rho^D_alpha(pm)
        · [s² - 2 s (m_N² - Q²) + (m_N² + Q²)²] / [q E_e Q² (1 - eps) m_N²] · (d sigma/d|t|)_0 e^{-beta t},   beta = 6 GeV^-2.

- The prefactor of `getCrosst` is this one: `alpha (Ebeam-nu)/(2 (2 pi)² Ebeam Q² (1-eps)) · 1/(2q)` equals
  `alpha/(4 (2 pi)²) · E_e'/(q E_e Q² (1-eps))`, and `getfrontfactor` is the bracket with `m_N²` and `e^{6t}`. So the
  **t version is the paper's formula**. The paper has no z integral and no t-z Jacobian (z only appears as the cut
  z > 0.9), so the z version needs the Jacobian of section 3, which is what the fix supplies.
- Eq. (12) of the paper, `d sigma^{gamma* N}/(d|t| dphi*) = m_N²/(8 pi² λ) · |M|²` (spinors normalised to
  `ubar u = 1`), is the standard two-body cross section; the `1/λ` is the flux/phase space (`λ` the polynomial above),
  and `m_N²` comes from the spinor normalisation. That is why Eq. (13) contains `λ/m_N²` and no separate `1/sqrt(λ)`
  flux factor, and why the weights `W_t`, `W_z` of section 2, which are only the phase-space part, differ from the
  code's weights by a factor that includes masses (the `4 kappa` of an earlier version of this note, which compared
  only the phase-space part; it is not an error).
- What remains: Eq. (12) is differential in the CM azimuth phi* (about the gamma* N axis), Eq. (13) in the lab azimuth
  phi_rho (about q), and the code gives each of the two roots the weight 1, i.e. it treats `dphi* = dphi_rho`.
  From (5), the exact weight per root is `|dphi*/dphi_rho| = 4 sqrt(λ) W_t = sqrt(λ)/(2 q kappa)` (my derivation;
  the identification of the integrands root by root follows from (5) but was not checked separately). Numerically this
  has mean 1.015 over roots with z > 0.9 (pm uniform in 0-0.33 GeV, nu 2.5-3.5 GeV, Q² 1.5 GeV², t in -0.4..-0.1) with
  large excursions at the edges of the root structure (min 0, max 42), and it integrates to exactly 1 over a full
  azimuth, so for a distribution that does not depend on the azimuth the weight 1 is exact on average. Whether it
  matters for the Glauber distributions (which depend on the angle between pm and the rho) was not studied.

The same equations are in `rhopaper_print.tex` (the latest version, after the journal review) with only cosmetic
changes (the label alpha on the amplitude, square brackets around the polynomial). The paper's cuts are
z > 0.9, W > 2 GeV and -0.4 < t < -0.1 GeV²; `getCrosst`/`getCrossz` apply the z and t cuts, but not the W cut,
which is a limit on the nu range at fixed Q²: since 2026-11-26 `rhoset` and `rhosingle` check it on the input
(`RhoTCross::checkWcut(Q2, nu_min)`, W² = m_p² + 2 m_p nu - Q² for a nucleon at rest, W growing with nu, so the lowest nu
decides) and throw `std::invalid_argument` if W(nu_min) < 2 GeV, unless the cuts are switched off (`nocuts`).
The use of the proton mass (and not the average nucleon mass) is inferred (the difference is below 1e-3 GeV²): to confirm.

## 5. Numerical checks (scratch programs, 2026-11-26; not in the repository)

- (4) against a finite difference of the t-version root `z(t)`: 2000 random points, maximum relative difference
  0.7% (finite-difference noise near the edges of the root structure).
- (5): the two sides agree to 0.3% in three kinematic points with an allowed `t` (4e6 azimuth points; the integrable
  `1/sqrt` edges make the sum converge slowly from below); outside the allowed range both root sets are empty.
- Code t weight against (5): integrates over the azimuth to `4 kappa_avg` times `pi/(2 sqrt(λ))`, `4 kappa_avg` = 3.3-4.2 GeV
  (3.69 GeV = `4 massi` at pm -> 0). Explained in section 4: the code weights are not the bare phase-space weights.
- Monte Carlo with the code's kinematics and cuts and a Gaussian momentum distribution: z/t = 11.66 with `massX`,
  0.997 with (4) at Q² = 1.5 GeV², nu 2.5-3.5 GeV; 0.988 at Q² = 2.5 GeV², nu 3.0-3.8 GeV.

## 6. Full comparison with the Glauber distributions (2026-11-26)

A small scratch driver (not in the repository: a Gauss-Legendre loop that calls `RhoTCross::getCrosst` over (nu, t)
and `RhoTCross::getCrossz` over (nu, z)) with the real `RhoTCross`: nucleus 0 (He), Q² = 2 GeV², nu = 3.0-3.2 GeV,
Ebeam = 4.023 GeV, t in -0.4..-0.1 GeV², z in 0.9..1, `usersigma` = 20 mb, rho particle type 4, `p_dil` = 0, cuts on,
inner `cube_adaptive` with maxEval = 3000 at precision 1e-3, **2 x 3 Gauss-Legendre points** per version (about 32
minutes each on one core, the cost is the `DistMomDistrGrid` per rho-angle bin). Indices as in the README of the
paper's calculation directory (0 RMSGA, 1 +SRC, 2 +decay, 3 +SRC+decay, 4 +CT, 5 +CT+SRC, 6 +CT+decay, 7 the works,
8 plane wave). Integrated results (arbitrary units, common factors as in the code) and ratios:

| Index | Case | t integral | z integral | z/t | T from t | T from z |
|---|---|---|---|---|---|---|
| 0 | RMSGA | 1.5224e-07 | 1.5559e-07 | 1.022 | 0.802 | 0.808 |
| 1 | RMSGA+SRC | 1.5296e-07 | 1.5621e-07 | 1.021 | 0.806 | 0.812 |
| 2 | RMSGA+Decay | 1.2764e-07 | 1.2889e-07 | 1.010 | 0.672 | 0.670 |
| 3 | RMSGA+SRC+Decay | 1.2811e-07 | 1.2927e-07 | 1.009 | 0.675 | 0.672 |
| 4 | RMSGA+CT | 1.6422e-07 | 1.6645e-07 | 1.014 | 0.865 | 0.865 |
| 5 | RMSGA+CT+SRC | 1.6427e-07 | 1.6647e-07 | 1.013 | 0.865 | 0.865 |
| 6 | RMSGA+CT+Decay | 1.4312e-07 | 1.4299e-07 | 0.999 | 0.754 | 0.743 |
| 7 | RMSGA+CT+SRC+Decay | 1.4294e-07 | 1.4281e-07 | 0.999 | 0.753 | 0.742 |
| 8 | plane wave | 1.8986e-07 | 1.9248e-07 | 1.014 | 1 | 1 |

With the old `massX` the z column would be larger by about `M_A/massi` (about 4 for helium; inferred from the Monte
Carlo, the unfixed code was not run here). Caveats: the 2 x 3 grid is very coarse and some points are exactly zero
because of the cuts (one t point and one z point in the run), so a 1-2% difference can be quadrature and cut-edge
effects: this shows that the factor ~A is gone, not that the remaining 1-2% is zero. Helium only. The nu range was
chosen so that W² = m_N² + 2 m_N nu - Q² stays above 4 GeV² (W > 2 GeV; the programs now check this, see the remark on the cuts above).
