# libEvolution: QCD evolution (Fortran 90).
# Built with CMake's default flags only (-O3 for Release), as it always was: it is not linked to
# MePhys/Integrator, through which the project flags used to reach the other libraries.
# (-march=native could change its floating-point results slightly; the !$omp directives in
# evolvers.f90 and wrappers.f90 are therefore inactive.)
set_property(DIRECTORY PROPERTY COMPILE_OPTIONS "")

physcode_add_library(Evolution
  SOURCES
    kinds.f90
    utils.f90
    decouple.f90
    evolvers.f90
    msbar.f90
    vector_nonsinglet.f90
    vector_singlet.f90
    wrappers.f90
    evolution_speedy.f90)
