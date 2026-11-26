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

## 4. What the weights are, and a remaining question

Equations (2) and (3) give the weights of the delta-function-derived integrand, `W_t = 1/(8 q kappa)` and
`W_z = nu/(4 tau)`. The code does not use these. It uses `1/(2q)` (t) and `(nu/q) K = nu kappa/tau` (z, after the fix),
which are `4 kappa` times `W_t` and `W_z`. The two code weights are consistent with each other, so t and z agree,
but both differ from the delta-function weights by the same factor `4 kappa`, which is pm dependent.

There is an independent check of the delta-function weights. Integrating over the azimuth φ at fixed `pm`, `theta`,
`nu`, `t` must give the Lorentz-invariant two-body phase space of the photon-nucleon system (go to its CM frame, where
`d³p_ρ/(2E_ρ) delta(f) = |p*|/(4 sqrt(s)) dΩ*`, `dt = 2 |q*||p*| dcosθ*`, `|q*| = sqrt(λ)/(2 sqrt(s))`):

    ∮ dφ Σ_roots W_t = pi / (2 sqrt(λ)),   λ = s² - 2 s (m_N² - Q²) + (m_N² + Q²)²,     (5)

for `t` inside the physical range of that `s` (zero outside). `λ` is exactly the polynomial in the front factor of
`getfrontfactor`. Equation (5) tests the derivation of `df/dz` (both the formula and the root structure).

## 5. Numerical checks (scratch programs, 2026-11-26; not in the repository)

- (4) against a finite difference of the t-version root `z(t)`: 2000 random points, maximum relative difference
  0.7% (finite-difference noise near the edges of the root structure).
- (5): the two sides agree to 0.3% in three kinematic points with an allowed `t` (4e6 azimuth points; the integrable
  `1/sqrt` edges make the sum converge slowly from below); outside the allowed range both root sets are empty.
- Code weight against (5): the t weight of the code integrates over the azimuth to `4 kappa_avg` times `pi/(2 sqrt(λ))`,
  with `4 kappa_avg` (the azimuthal harmonic mean of kappa) = 3.65-3.70 GeV at pm ~ 0.01-0.02 GeV (`4 massi` = 3.69 GeV)
  and 3.3-4.2 GeV for pm = 0.09-0.2 GeV. So the t weight (and with the fix the z weight) is the invariant weight times about `4 massi` (an overall constant that
  might be a convention absorbed elsewhere, **not checked**) times a pm-dependent factor of up to about ±12%
  (**the author should decide** whether to replace both weights by the delta-function weights; this would change the
  t results too, by that factor).
- Monte Carlo with the code's kinematics and cuts and a Gaussian momentum distribution: z/t = 11.66 with `massX`,
  0.997 with (4) at Q² = 1.5 GeV², nu 2.5-3.5 GeV; 0.988 at Q² = 2.5 GeV², nu 3.0-3.8 GeV.
