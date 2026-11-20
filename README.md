# nucqcd

Code for high-energy reactions on nuclei (JLab to EIC energies), used to study QCD and the strong
force in nuclei. Formerly known as `physics-code`; see [docs/HISTORY.md](docs/HISTORY.md).

Short version for now...

The manual subfolder has a tex file manual.tex with some info.

This repository contains 5 libraries, each in its subfolder
- Integrator has several integration routines (3rd party)
- MePhys contains a lot of standard physics and math things
- Glauber has all the glauber code and mean-field nucleus stuff
- Knockout has quasi-elastic models
- DIS has deep-inelastic models

progs subfolder contains a lot of executables, both for testing and used for calculations (documenting in progress)
