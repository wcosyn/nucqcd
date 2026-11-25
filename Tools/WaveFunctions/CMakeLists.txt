# libToolsWaveFunctions: deuteron wave functions (Deuteron and its implementations YukawaWavefunction,
# GrossWavefunction, CstWavefunction, InterpolatingWavefunction; OldDeuteron) and the He3 wave function (He3wf, a
# wrapper around the Fortran getwavemod_tm.f with its helpers clebsch.f and spline.f).
# Deuteron.hpp uses Kinematics2to3, hence the link to ToolsKinematics.

nucqcd_add_library(ToolsWaveFunctions
  SOURCES
    Deuteron.cpp
    WavefunctionImplementation.cpp
    YukawaWavefunction.cpp
    GrossWavefunction.cpp
    CstWavefunction.cpp
    InterpolatingWavefunction.cpp
    OldDeuteron.cpp
    WaveFunctionModels.cpp
    He3wf.cpp
    # Fortran
    getwavemod_tm.f
    clebsch.f
    spline.f
  LINK ToolsKinematics)
