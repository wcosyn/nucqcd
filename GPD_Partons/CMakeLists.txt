# libGPD: deuteron vector and transversity GPD convolutions, built on PARTONS.
# Only included with -DPARTONS_ENABLED=yes (the PARTONS libraries are found in the top level).
physcode_add_library(GPD
  SOURCES
    src/Deut_Conv_GPD_T.cpp
    src/Deut_Conv_GPD_V.cpp
    src/Deut_GPD_T_set.cpp
    src/Deut_GPD_V_set.cpp
    src/DiDVCS.cpp
    src/GPD_T_Nucl_grid.cpp
    src/GPD_V_Nucl_grid.cpp
    src/TransGPD_set.cpp
    src/TwoVector_Deut.cpp
    src/TwoVector_Nucl.cpp
  LINK partons_libs MePhys Integrator Glauber cuba)
# the headers are in include/, not next to this file
target_include_directories(GPD PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
if(STATIC_BUILD)
  target_include_directories(GPDST PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
endif()
