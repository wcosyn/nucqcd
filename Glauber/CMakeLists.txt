# libGlauber: mean-field nuclei, nucleon current operators, Glauber/RMSGA FSI grids.
nucqcd_add_library(Glauber
  SOURCES
    MeanFieldNucleus.cpp
    MeanFieldNucleusThick.cpp
    NucleonEMOperator.cpp
    NucleonEMOperatorTables.cpp
    NucleonWeakOperator.cpp
    FsiCorrelator.cpp
    FastParticle.cpp
    AbstractFsiGrid.cpp
    AbstractFsiCTGrid.cpp
    AbstractFsiCTDecayGrid.cpp
    AbstractFsiGridThick.cpp
    AbstractFsiDecayGridThick.cpp
    AbstractFsiCTGridThick.cpp
    AbstractFsiCTDecayGridThick.cpp
    GlauberGrid.cpp
    GlauberGridThick.cpp
    GlauberDecayGridThick.cpp
    OneGlauberGrid.cpp
    TMFSpinor.cpp
    DistMomDistrGrid.cpp
    ROMEAGrid.cpp
    ClassGridThick_SCX.cpp
    global2009.f          # Fortran
  LINK ToolsWaveFunctions Integrator)
