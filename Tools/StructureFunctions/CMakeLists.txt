# libToolsStructureFunctions: nucleon structure functions and PDFs (NucleonStructure, its kinematics
# wrapper NuclStructure, LHAPDFStructure) and the off-shell deuteron structure functions
# (DeuteronStructure), with the Fortran parametrizations (Alekhin, Christy-Bosted F1F209, CTEQ, GRSV,
# HMRS, R1990/R1998). DeuteronStructure uses the kinematics classes, hence the link to ToolsKinematics.

nucqcd_add_library(ToolsStructureFunctions
  SOURCES
    NucleonStructure.cpp
    NuclStructure.cpp
    LHAPDFStructure.cpp
    DeuteronStructure.cpp
    # Fortran
    alekhin.f
    F1F209_test.f
    cteq.f
    grsv2000pdf_g1.f
    R199x.f
    hmrs-b.f
  LINK ToolsKinematics)
if(TARGET PkgConfig::LHAPDF)
  # LHAPDFStructure.cpp uses LHAPDF only when NUCQCD_HAVE_LHAPDF is defined (it throws otherwise)
  set_property(SOURCE LHAPDFStructure.cpp APPEND PROPERTY COMPILE_DEFINITIONS NUCQCD_HAVE_LHAPDF)
  target_link_libraries(ToolsStructureFunctions PUBLIC PkgConfig::LHAPDF)
  if(STATIC_BUILD)
    target_link_libraries(ToolsStructureFunctionsST PUBLIC PkgConfig::LHAPDF)
  endif()
endif()
