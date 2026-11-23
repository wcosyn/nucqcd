# Test suite (Catch2 v3 + CTest). See tests/README.md.
#
#   cd build && cmake .. && make -j8 && ctest -LE slow        # default run
#   ctest -L slow                                            # Glauber-grid tests (minutes)
#   NUCQCD_UPDATE_REFERENCE=1 ctest -L "regression|slow"   # re-baseline golden values
#
# Labels: unit (analytic/exact checks), regression (golden values captured from master),
# slow (regression tests that compute Glauber grids). Test executables are not installed.

include(FetchContent)
FetchContent_Declare(Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.7.1
  GIT_SHALLOW    TRUE)
FetchContent_MakeAvailable(Catch2)
list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
include(Catch)

set(NUCQCD_TEST_SHARE "${CMAKE_CURRENT_BINARY_DIR}/share")   # sandbox copy of share/
set(NUCQCD_TEST_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/work")  # cwd of the tests (Fortran writes here)
set(NUCQCD_TEST_REFDIR "${CMAKE_CURRENT_SOURCE_DIR}/reference")
set(NUCQCD_TEST_OUTDIR "${CMAKE_CURRENT_BINARY_DIR}/bin")
file(MAKE_DIRECTORY "${NUCQCD_TEST_WORKDIR}")

# ---- support library: main(), reference store, sandbox location --------------------------------
add_library(nucqcd_test_support STATIC support/test_support.cpp support/test_main.cpp)
target_include_directories(nucqcd_test_support PUBLIC support)
target_link_libraries(nucqcd_test_support PUBLIC Catch2::Catch2)
target_compile_definitions(nucqcd_test_support PRIVATE
  NUCQCD_TEST_SHARE_DEFAULT="${NUCQCD_TEST_SHARE}"
  NUCQCD_TEST_REFDIR_DEFAULT="${NUCQCD_TEST_REFDIR}")

# ---- fixture: (re)create the sandbox share dir before any test that needs it --------------------
# Read-only data are symlinked; input/ is copied (the code writes into it);
# grids/ is emptied, so Glauber grids are always recomputed and the user's share/grids/ is
# never read or written.
add_test(NAME share_sandbox
  COMMAND ${CMAKE_COMMAND} -DSRC=${PROJECT_SOURCE_DIR}/share -DDEST=${NUCQCD_TEST_SHARE}
          -P ${CMAKE_CURRENT_SOURCE_DIR}/support/make_test_share.cmake)
set_tests_properties(share_sandbox PROPERTIES FIXTURES_SETUP share_sandbox LABELS "unit;regression;slow")

# nucqcd_add_test(<name> LABEL <unit|regression|slow> LINK <libs...> [TIMEOUT s])
#   Executable <name> from <name>.cpp in unit/ or regression/; each Catch2 TEST_CASE becomes
#   its own CTest test (own process: isolates Fortran COMMON blocks and other static state).
function(nucqcd_add_test name)
  cmake_parse_arguments(ARG "" "LABEL;TIMEOUT;DIR" "LINK;SOURCES" ${ARGN})
  if(NOT ARG_TIMEOUT)
    set(ARG_TIMEOUT 300)
  endif()
  add_executable(${name} ${ARG_DIR}/${name}.cpp ${ARG_SOURCES})
  target_link_libraries(${name} PRIVATE nucqcd_test_support ${ARG_LINK})
  # link against the libraries in build/lib, not the installed copies in <repo>/lib
  set_target_properties(${name} PROPERTIES
    BUILD_WITH_INSTALL_RPATH OFF
    RUNTIME_OUTPUT_DIRECTORY "${NUCQCD_TEST_OUTDIR}")
  catch_discover_tests(${name}
    TEST_PREFIX "${name}:"
    WORKING_DIRECTORY "${NUCQCD_TEST_WORKDIR}"
    DISCOVERY_MODE PRE_TEST
    PROPERTIES LABELS ${ARG_LABEL} TIMEOUT ${ARG_TIMEOUT} FIXTURES_REQUIRED share_sandbox)
endfunction()

# ---- unit tests: exact / analytic results ------------------------------------------------------
nucqcd_add_test(unit_numint          DIR unit LABEL unit LINK Integrator cuba)
nucqcd_add_test(unit_cuba            DIR unit LABEL unit LINK Integrator cuba)
nucqcd_add_test(unit_lorentz_gamma   DIR unit LABEL unit LINK MePhys cuba)
nucqcd_add_test(unit_deuteron_wf     DIR unit LABEL unit LINK MePhys cuba)
nucqcd_add_test(unit_kinematics      DIR unit LABEL unit LINK MePhys cuba)
nucqcd_add_test(unit_mathutil        DIR unit LABEL unit LINK MePhys cuba)
nucqcd_add_test(unit_lhapdf          DIR unit LABEL unit LINK MePhys cuba)
nucqcd_add_test(unit_he3wf           DIR unit LABEL unit LINK DIS cuba)
# Pair lives with its only user in progs/attic/twonucl (unfinished A(e,e'NN) code)
set(NUCQCD_TWONUCL_DIR ${PROJECT_SOURCE_DIR}/progs/attic/twonucl)
nucqcd_add_test(unit_pair            DIR unit LABEL unit LINK MePhys cuba
                SOURCES ${NUCQCD_TWONUCL_DIR}/pair.cpp ${NUCQCD_TWONUCL_DIR}/newcoef.cpp)
target_include_directories(unit_pair PRIVATE ${NUCQCD_TWONUCL_DIR})
nucqcd_add_test(unit_constants       DIR unit LABEL unit LINK MePhys cuba)
nucqcd_add_test(unit_meanfield       DIR unit LABEL unit LINK Glauber cuba)
nucqcd_add_test(unit_ownership       DIR unit LABEL unit LINK DIS cuba)
nucqcd_add_test(unit_errors          DIR unit LABEL unit LINK Glauber cuba)
nucqcd_add_test(unit_grid_bins       DIR unit LABEL unit LINK Glauber cuba)
nucqcd_add_test(unit_moshinsky       DIR unit LABEL unit LINK MePhys cuba)
nucqcd_add_test(unit_fortran_string  DIR unit LABEL unit LINK MePhys cuba)
nucqcd_add_test(unit_grid_cache      DIR unit LABEL unit LINK MePhys cuba)

# ---- regression tests: golden values from master (tests/reference/*.ref) ---------------------
nucqcd_add_test(reg_nucleon_structure DIR regression LABEL regression LINK MePhys cuba)
nucqcd_add_test(reg_deuteron_wf       DIR regression LABEL regression LINK MePhys cuba)
nucqcd_add_test(reg_meanfield         DIR regression LABEL regression LINK Glauber cuba)
nucqcd_add_test(reg_poldeut           DIR regression LABEL regression LINK DIS cuba)
nucqcd_add_test(reg_inclusive_dis     DIR regression LABEL regression LINK DIS cuba)
nucqcd_add_test(reg_deuteron_tagged   DIR regression LABEL regression LINK DIS cuba)
nucqcd_add_test(reg_eep_planewave     DIR regression LABEL regression LINK Knockout cuba)
nucqcd_add_test(reg_weakqe            DIR regression LABEL regression LINK Knockout cuba)
nucqcd_add_test(reg_glauber_points    DIR regression LABEL regression LINK Glauber cuba)
nucqcd_add_test(reg_glauber_grid      DIR regression LABEL slow LINK Knockout cuba TIMEOUT 1800)
