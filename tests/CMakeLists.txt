# Test suite (Catch2 v3 + CTest). See tests/README.md.
#
#   cd build && cmake .. && make -j8 && ctest -LE slow        # default run
#   ctest -L slow                                            # Glauber-grid tests (minutes)
#   PHYSCODE_UPDATE_REFERENCE=1 ctest -L "regression|slow"   # re-baseline golden values
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

set(PHYSCODE_TEST_SHARE "${CMAKE_CURRENT_BINARY_DIR}/share")   # sandbox copy of share/
set(PHYSCODE_TEST_WORKDIR "${CMAKE_CURRENT_BINARY_DIR}/work")  # cwd of the tests (Fortran writes here)
set(PHYSCODE_TEST_REFDIR "${CMAKE_CURRENT_SOURCE_DIR}/reference")
set(PHYSCODE_TEST_OUTDIR "${CMAKE_CURRENT_BINARY_DIR}/bin")
file(MAKE_DIRECTORY "${PHYSCODE_TEST_WORKDIR}")

# ---- support library: main(), reference store, sandbox location --------------------------------
add_library(physcode_test_support STATIC support/test_support.cpp support/test_main.cpp)
target_include_directories(physcode_test_support PUBLIC support)
target_link_libraries(physcode_test_support PUBLIC Catch2::Catch2)
target_compile_definitions(physcode_test_support PRIVATE
  PHYSCODE_TEST_SHARE_DEFAULT="${PHYSCODE_TEST_SHARE}"
  PHYSCODE_TEST_REFDIR_DEFAULT="${PHYSCODE_TEST_REFDIR}")

# ---- fixture: (re)create the sandbox share dir before any test that needs it --------------------
# Read-only data are symlinked; input/ is copied (the code writes into it);
# grids/ is emptied, so Glauber grids are always recomputed and the user's share/grids/ is
# never read or written.
add_test(NAME share_sandbox
  COMMAND ${CMAKE_COMMAND} -DSRC=${PROJECT_SOURCE_DIR}/share -DDEST=${PHYSCODE_TEST_SHARE}
          -P ${CMAKE_CURRENT_SOURCE_DIR}/support/make_test_share.cmake)
set_tests_properties(share_sandbox PROPERTIES FIXTURES_SETUP share_sandbox LABELS "unit;regression;slow")

# physcode_add_test(<name> LABEL <unit|regression|slow> LINK <libs...> [TIMEOUT s])
#   Executable <name> from <name>.cpp in unit/ or regression/; each Catch2 TEST_CASE becomes
#   its own CTest test (own process: isolates Fortran COMMON blocks and other static state).
function(physcode_add_test name)
  cmake_parse_arguments(ARG "" "LABEL;TIMEOUT;DIR" "LINK" ${ARGN})
  if(NOT ARG_TIMEOUT)
    set(ARG_TIMEOUT 300)
  endif()
  add_executable(${name} ${ARG_DIR}/${name}.cpp)
  target_link_libraries(${name} PRIVATE physcode_test_support ${ARG_LINK})
  # link against the libraries in build/lib, not the installed copies in <repo>/lib
  set_target_properties(${name} PROPERTIES
    BUILD_WITH_INSTALL_RPATH OFF
    RUNTIME_OUTPUT_DIRECTORY "${PHYSCODE_TEST_OUTDIR}")
  catch_discover_tests(${name}
    TEST_PREFIX "${name}:"
    WORKING_DIRECTORY "${PHYSCODE_TEST_WORKDIR}"
    DISCOVERY_MODE PRE_TEST
    PROPERTIES LABELS ${ARG_LABEL} TIMEOUT ${ARG_TIMEOUT} FIXTURES_REQUIRED share_sandbox)
endfunction()

# ---- unit tests: exact / analytic results ------------------------------------------------------
physcode_add_test(unit_numint          DIR unit LABEL unit LINK Integrator cuba)
physcode_add_test(unit_cuba            DIR unit LABEL unit LINK Integrator cuba)
physcode_add_test(unit_lorentz_gamma   DIR unit LABEL unit LINK MePhys cuba)
physcode_add_test(unit_deuteron_wf     DIR unit LABEL unit LINK MePhys cuba)
physcode_add_test(unit_kinematics      DIR unit LABEL unit LINK MePhys cuba)
physcode_add_test(unit_constants       DIR unit LABEL unit LINK MePhys cuba)
physcode_add_test(unit_meanfield       DIR unit LABEL unit LINK Glauber cuba)
physcode_add_test(unit_moshinsky       DIR unit LABEL unit LINK MePhys cuba)

# ---- regression tests: golden values from master (tests/reference/*.ref) ---------------------
physcode_add_test(reg_nucleon_structure DIR regression LABEL regression LINK MePhys cuba)
physcode_add_test(reg_deuteron_wf       DIR regression LABEL regression LINK MePhys cuba)
physcode_add_test(reg_meanfield         DIR regression LABEL regression LINK Glauber cuba)
physcode_add_test(reg_poldeut           DIR regression LABEL regression LINK DIS cuba)
physcode_add_test(reg_inclusive_dis     DIR regression LABEL regression LINK DIS cuba)
physcode_add_test(reg_deuteron_tagged   DIR regression LABEL regression LINK DIS cuba)
physcode_add_test(reg_eep_planewave     DIR regression LABEL regression LINK Knockout cuba)
physcode_add_test(reg_weakqe            DIR regression LABEL regression LINK Knockout cuba)
physcode_add_test(reg_glauber_points    DIR regression LABEL regression LINK Glauber cuba)
physcode_add_test(reg_glauber_grid      DIR regression LABEL slow LINK Knockout cuba TIMEOUT 1800)
