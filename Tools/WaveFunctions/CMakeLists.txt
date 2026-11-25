# libToolsWaveFunctions: deuteron wave functions (TDeuteron and its implementations TYukawaPWF,
# TGrossPWF, TCstPWF, TInterpolatingWavefunction; OldDeuteron) and the He3 wave function (He3wf, a
# wrapper around the Fortran getwavemod_tm.f with its helpers clebsch.f and spline.f).
# TDeuteron.hpp uses TKinematics2to3, hence the link to ToolsKinematics.

nucqcd_add_library(ToolsWaveFunctions
  SOURCES
    TDeuteron.cpp
    TWavefunctionImplementation.cpp
    TYukawaPWF.cpp
    TGrossPWF.cpp
    TCstPWF.cpp
    TInterpolatingWavefunction.cpp
    OldDeuteron.cpp
    WaveFunctionModels.cpp
    He3wf.cpp
    # Fortran
    getwavemod_tm.f
    clebsch.f
    spline.f
  LINK ToolsKinematics)
