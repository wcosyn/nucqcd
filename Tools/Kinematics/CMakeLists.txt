# libToolsKinematics: kinematics of 2->2 and 2->3 reactions (Kinematics, Kinematics2to2,
# Kinematics2to3(WithLabAngles)), the lepton side (ElectronKinematics, LeptonKinematics), and
# light-cone kinematics (LightConeKin2to2/2to3).

nucqcd_add_library(ToolsKinematics
  SOURCES
    Kinematics.cpp
    Kinematics2to2.cpp
    Kinematics2to3.cpp
    Kinematics2to3_ReadFormatString.cpp
    Kinematics2to3_Update.cpp
    Kinematics2to3_Variables.cpp
    Kinematics2to3_Getters.cpp
    Kinematics2to3_Frames.cpp
    Kinematics2to3_VarNames.cpp
    Kinematics2to3WithLabAngles.cpp
    ElectronKinematics.cpp
    LeptonKinematics.cpp
    LightConeKin2to2.cpp
    LightConeKin2to3.cpp
  LINK ToolsCore)
