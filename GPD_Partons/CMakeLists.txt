# libGPD: deuteron vector and transversity GPD convolutions, built on PARTONS.
# Only included with -DPARTONS_ENABLED=yes (the PARTONS libraries are found in the top level).
nucqcd_add_library(GPD
  SOURCES
    Deut_Conv_GPD_T.cpp
    Deut_Conv_GPD_V.cpp
    Deut_GPD_T_set.cpp
    Deut_GPD_V_set.cpp
    DiDVCS.cpp
    GPD_T_Nucl_grid.cpp
    GPD_V_Nucl_grid.cpp
    TransGPD_set.cpp
    TwoVector_Deut.cpp
    TwoVector_Nucl.cpp
  LINK partons_libs ToolsStructureFunctions Integrator Glauber cuba)
