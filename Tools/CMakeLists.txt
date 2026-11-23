# Tools: general physics and math, in four libraries (roadmap step 6, "Code layout"), each linking
# the ones before it:
#   ToolsCore               constants, errors, string/math utilities, grid cache, RunLog, the Dirac
#                           and Lorentz algebra (FourVector, Matrix, GammaStructure, TSpinor, ...)
#   ToolsKinematics         TKinematics*, TElectronKinematics, TLeptonKinematics, LightConeKin*
#   ToolsWaveFunctions      deuteron (TDeuteron and its parametrizations, OldDeuteron) and He3 (He3wf)
#   ToolsStructureFunctions nucleon/deuteron structure functions and PDFs (NucleonStructure, ...)
# (StructureFunctions needs Kinematics, not WaveFunctions.) The names carry a "Tools" prefix
# because a libCore.so would collide with ROOT's.
#
# MePhys, the former single library, is an INTERFACE target that links all four, so programs that
# link MePhys need no change.

add_subdirectory(Core)
add_subdirectory(Kinematics)
add_subdirectory(WaveFunctions)
add_subdirectory(StructureFunctions)

add_library(MePhys INTERFACE)
target_link_libraries(MePhys INTERFACE ToolsCore ToolsKinematics ToolsWaveFunctions ToolsStructureFunctions)
if(STATIC_BUILD)
  add_library(MePhysST INTERFACE)
  target_link_libraries(MePhysST INTERFACE
    ToolsCoreST ToolsKinematicsST ToolsWaveFunctionsST ToolsStructureFunctionsST)
endif()
