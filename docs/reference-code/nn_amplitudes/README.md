# Nucleon-nucleon amplitudes (Fortran reference code)

These files are not compiled or used by any library or program in this repository. They were added in
2013 (commit `9b1784a`) as reference material for charge-exchange and NN rescattering, and were kept in
`share/input/` until 2026.

- `bystricky.f`: NN elastic scattering helicity amplitudes `f_M(s,t,...)` in the Bystricky a–e
  parametrization (Rev. Mod. Phys. 65, 47 (1993)). `Bystricky_dat` reads the tabulated amplitudes below.
- `fpn_chex.f`: the pn charge-exchange amplitude `f_pn_chex`.
- `B_ampl_pp.dat`, `B_ampl_pn.dat`, `B_ampl_pp_nocoul.dat`: Bystricky amplitude tables (pp with and
  without Coulomb effects, pn). `bystricky.f` opens them by relative path, so run it from this directory.
