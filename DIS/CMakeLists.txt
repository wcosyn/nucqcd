# libDIS: inclusive and tagged DIS on the deuteron and He3, polarized deuteron, resonances.
nucqcd_add_library(DIS
  SOURCES
    InclusiveCross.cpp
    # InclusiveCrossRes.cpp   # only used the Romberg routines, no longer used
    DeuteronCross.cpp
    Resonance.cpp
    He3Cross.cpp
    Poldeut.cpp
    Deut_Conv_EMTLF.cpp
  LINK Knockout ToolsStructureFunctions Integrator cuba)

add_subdirectory(evolution)
