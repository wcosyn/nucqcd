# libGlauber: mean-field nuclei, nucleon current operators, Glauber/RMSGA FSI grids.
nucqcd_add_library(Glauber
  SOURCES
    MeanFieldNucleus.cpp
    MeanFieldNucleusThick.cpp
    NucleonEMOperator.cpp
    NucleonEMOperatorTables.cpp
    NucleonWeakOperator.cpp
    FsiCorrelator.cpp
    FsiCorrelator_Solver.cpp
    FastParticle.cpp
    FastParticle_ScatteringData.cpp
    AbstractFsiGrid.cpp
    AbstractFsiCTGrid.cpp
    GlauberGrid.cpp
    GlauberGrid_Construct.cpp
    GlauberGridThick.cpp
    GlauberDecayGridThick.cpp
    GlauberThickGridCommon.cpp
    OneGlauberGrid.cpp
    MeanFieldSpinor.cpp
    DistMomDistrGrid.cpp
    DistMomDistrGrid_Construct.cpp
    ROMEAGrid.cpp
    ClassGridThick_SCX.cpp
    global2009.f          # Fortran
  LINK ToolsWaveFunctions Integrator)
