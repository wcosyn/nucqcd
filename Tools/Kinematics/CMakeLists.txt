# libToolsKinematics: kinematics of 2->2 and 2->3 reactions (TKinematics, TKinematics2to2,
# TKinematics2to3(WithLabAngles)), the lepton side (TElectronKinematics, TLeptonKinematics), and
# light-cone kinematics (LightConeKin2to2/2to3).

nucqcd_add_library(ToolsKinematics
  SOURCES
    TKinematics.cpp
    TKinematics2to2.cpp
    TKinematics2to3.cpp
    TKinematics2to3_ReadFormatString.cpp
    TKinematics2to3_Update.cpp
    TKinematics2to3_Variables.cpp
    TKinematics2to3_Getters.cpp
    TKinematics2to3_Frames.cpp
    TKinematics2to3WithLabAngles.cpp
    TElectronKinematics.cpp
    TLeptonKinematics.cpp
    LightConeKin2to2.cpp
    LightConeKin2to3.cpp
  LINK ToolsCore)
