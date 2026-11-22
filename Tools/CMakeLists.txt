# libMePhys: general physics and math (kinematics, spinors, deuteron wave functions,
# nucleon/nuclear structure functions and PDFs, partly in Fortran).

# constants.hpp (physical constants and HOMEDIR = <repo>/share) is generated in the source
# directory, so that code compiled outside CMake can include it too; it is git-ignored.
set(SHARE_DIR "${PROJECT_SOURCE_DIR}/share")
configure_file(constants.hpp.in "${CMAKE_CURRENT_SOURCE_DIR}/constants.hpp")

set(mpi_sources)
if(MPI_ENABLED)
  set(mpi_sources TMPI.cpp)
endif()

# solve_3x3 has to round as ROOT's TDecompLU did (TYukawaPWF D-wave coefficients): no fused multiply-adds
set_source_files_properties(MathUtil.cpp PROPERTIES COMPILE_OPTIONS -ffp-contract=off)

nucqcd_add_library(MePhys
  SOURCES
    NuclStructure.cpp
    NucleonStructure.cpp
    DeuteronStructure.cpp
    Utilfunctions.cpp
    MathUtil.cpp
    GridCache.cpp
    LightConeKin2to2.cpp
    LightConeKin2to3.cpp
    He3wf.cpp
    TCstPWF.cpp
    TDeuteron.cpp
    TGrossPWF.cpp
    TInterpolatingWavefunction.cpp
    TWavefunctionImplementation.cpp
    TYukawaPWF.cpp
    Matrix.cpp
    GammaStructure.cpp
    FourVector.cpp
    TSpinor.cpp
    TLorentzQuaternion.cpp
    TKinematics.cpp
    TKinematics2to2.cpp
    TElectronKinematics.cpp
    recmosh.cpp
    TKinematics2to3.cpp
    TKinematics2to3WithLabAngles.cpp
    TLeptonKinematics.cpp
    mstwpdf.cpp
    OldDeuteron.cpp
    ${mpi_sources}
    # Fortran
    alekhin.f
    F1F209_test.f
    cteq.f
    spline.f
    clebsch.f
    getwavemod_tm.f
    grsv2000pdf_g1.f
    R199x.f
    hmrs-b.f
  LINK root_libs GSL::gsl GSL::gslcblas)
if(MPI_ENABLED)
  target_link_libraries(MePhys PUBLIC MPI::MPI_CXX)
endif()
