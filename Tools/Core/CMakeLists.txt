# libToolsCore: constants, errors, string and math utilities, the grid cache, run bookkeeping
# (RunLog), the Dirac/Lorentz algebra (FourVector, Matrix, GammaStructure, TLorentzQuaternion,
# TSpinor), Moshinsky brackets (recmosh), and the ROOT-free vector and Lorentz classes (Vector3,
# LorentzVector, Rotation, LorentzRotation; roadmap step 9, phase 4).

# constants.hpp (physical constants and nucqcd::share_dir = <repo>/share) is generated in the source
# directory, so that code compiled outside CMake can include it too; it is git-ignored.
set(SHARE_DIR "${PROJECT_SOURCE_DIR}/share")
configure_file(constants.hpp.in "${CMAKE_CURRENT_SOURCE_DIR}/constants.hpp")

# Build information for nucqcd::RunLog: the configuration (configure time) and the git state of the
# source tree (refreshed at every build by cmake/BuildInfo.cmake; the file only changes with the state)
set(NUCQCD_GENERATED_DIR "${PROJECT_BINARY_DIR}/generated")
set(NUCQCD_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
if(NOT NUCQCD_BUILD_TYPE)
  set(NUCQCD_BUILD_TYPE "none")
endif()
string(TOUPPER "${NUCQCD_BUILD_TYPE}" build_type_upper)
set(NUCQCD_FLAGS_STRING "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${build_type_upper}} -Wall -Wextra")
if(NUCQCD_BUILD_TYPE STREQUAL "Release")
  list(JOIN NUCQCD_RELEASE_FLAGS " " release_flags)
  string(APPEND NUCQCD_FLAGS_STRING " ${release_flags}")
endif()
string(STRIP "${NUCQCD_FLAGS_STRING}" NUCQCD_FLAGS_STRING)
string(REGEX REPLACE " +" " " NUCQCD_FLAGS_STRING "${NUCQCD_FLAGS_STRING}")
set(NUCQCD_LIBRARIES_STRING "GSL ${GSL_VERSION}")
if(ROOT_FOUND)
  string(PREPEND NUCQCD_LIBRARIES_STRING "ROOT ${ROOT_VERSION}, ")  # (the libraries do not use it)
endif()
if(LHAPDF_FOUND)
  string(APPEND NUCQCD_LIBRARIES_STRING ", LHAPDF ${LHAPDF_VERSION}")
endif()
configure_file("${PROJECT_SOURCE_DIR}/cmake/nucqcd_build_config.inc.in" "${NUCQCD_GENERATED_DIR}/nucqcd_build_config.inc" @ONLY)
add_custom_target(nucqcd_git_info
  COMMAND ${CMAKE_COMMAND} -DSRC=${PROJECT_SOURCE_DIR} -DOUT=${NUCQCD_GENERATED_DIR}/nucqcd_git_info.inc
          -P ${PROJECT_SOURCE_DIR}/cmake/BuildInfo.cmake
  BYPRODUCTS ${NUCQCD_GENERATED_DIR}/nucqcd_git_info.inc
  COMMENT "git state of the source tree for nucqcd::RunLog")
set_property(SOURCE RunLog.cpp APPEND PROPERTY INCLUDE_DIRECTORIES "${NUCQCD_GENERATED_DIR}")

# solve_3x3 has to round as ROOT's TDecompLU did (TYukawaPWF D-wave coefficients): no fused multiply-adds
set_source_files_properties(MathUtil.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)
# The vector and Lorentz classes reproduce ROOT's TVector3 & co. bitwise: the functions that ROOT
# compiles into libPhysics (built without fused multiply-adds) are in these files, so no FMAs here either
set_source_files_properties(Vector3.cpp LorentzVector.cpp Rotation.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)

nucqcd_add_library(ToolsCore
  SOURCES
    Utilfunctions.cpp
    RunLog.cpp
    MathUtil.cpp
    GridCache.cpp
    Matrix.cpp
    GammaStructure.cpp
    FourVector.cpp
    TSpinor.cpp
    TLorentzQuaternion.cpp
    recmosh.cpp
    Vector3.cpp
    LorentzVector.cpp
    Rotation.cpp
  LINK GSL::gsl GSL::gslcblas)
add_dependencies(ToolsCore nucqcd_git_info)
if(STATIC_BUILD)
  add_dependencies(ToolsCoreST nucqcd_git_info)
endif()
