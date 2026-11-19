# physics-code: libraries (Integrator, MePhys, Glauber, Knockout, DIS, Evolution, GPD) and programs.
#
#   cd build
#   cmake -DCMAKE_BUILD_TYPE=Release [-DBUILD_PROGS=yes] [-DPARTONS_ENABLED=yes] ..
#   make install -j4        # libraries -> <repo>/lib, programs -> <repo>/bin
#   make doc                # Doxygen documentation
#
# Options: BUILD_PROGS, PARTONS_ENABLED, MPI_ENABLED, STATIC_BUILD (all OFF by default).
# ROOT (with MathMore and Minuit2) is found through its CMake config: source thisroot.sh,
# or ROOT in /opt/root is used.
cmake_minimum_required(VERSION 3.16)
project(physics-code LANGUAGES C CXX Fortran)

option(BUILD_PROGS "Build the programs in progs/ with 'make' and install them" OFF)
option(PARTONS_ENABLED "Build GPD_Partons and progs/GPD (needs PARTONS, see scripts/install_partons.sh)" OFF)
option(MPI_ENABLED "Enable MPI (TMPI class and MPI programs)" OFF)
option(STATIC_BUILD "Also build static versions of the libraries" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # compile_commands.json for VS Code

# ---- output and install layout -------------------------------------------------------------
# Build products go to build/lib and build/bin; 'make install' puts them into the source tree
# (<repo>/lib, <repo>/bin) unless another CMAKE_INSTALL_PREFIX is given. Binaries get an RPATH
# to <prefix>/lib and to the directories of the external libraries they link (ROOT, PARTONS...).
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
# (/usr/local, CMake's default, is also replaced: build directories configured with the old
# CMakeLists have it in their cache, where it is no longer flagged as the default)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR CMAKE_INSTALL_PREFIX STREQUAL "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" CACHE PATH "Install prefix" FORCE)
endif()
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)

# ---- compiler flags (all languages: C, C++, Fortran) ------------------------------------------
set(PHYSCODE_RELEASE_FLAGS)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(PHYSCODE_RELEASE_FLAGS -O3 -march=native -mtune=native -pipe)
  # -fopenmp only for C/C++ (which contain no OpenMP directives, as before). The Fortran in
  # DIS/evolution does contain !$omp directives but was always built without OpenMP; enabling
  # it there would also need OpenMP at link time (find_package(OpenMP), OpenMP::OpenMP_Fortran).
  add_compile_options("$<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:C,CXX>>:-fopenmp>")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")     # classic icpc
  set(PHYSCODE_RELEASE_FLAGS -O3 -shared-intel -fp-model source -unroll -pthread -ipo -xavx)
  if(MPI_ENABLED)
    list(APPEND PHYSCODE_RELEASE_FLAGS -parallel)
    add_link_options(-parallel)
  endif()
  find_program(XIAR xiar)                          # archiver needed for -ipo static libraries
  if(XIAR)
    set(CMAKE_AR "${XIAR}")
  endif()
  mark_as_advanced(XIAR)
endif()
add_compile_options(-Wall -Wextra "$<$<CONFIG:Release>:${PHYSCODE_RELEASE_FLAGS}>" "$<$<CONFIG:Debug>:-g>")

# ---- external dependencies -------------------------------------------------------------------
find_package(ROOT 6 CONFIG REQUIRED COMPONENTS MathMore Minuit Minuit2
             HINTS "$ENV{ROOTSYS}/cmake" /opt/root/cmake)
if(NOT ROOT_mathmore_FOUND)
  message(FATAL_ERROR "ROOT has to be built with MathMore (ROOT in ${ROOT_DIR})")
endif()
message(STATUS "ROOT ${ROOT_VERSION} in ${ROOT_DIR}")
# ROOT's libraries are linked by path (ROOT_LIBRARIES) rather than through the ROOT:: targets,
# because those require ROOT's own C++ standard (C++23 for ROOT 6.40) and would switch the
# whole code base to it.
add_library(root_libs INTERFACE)
target_include_directories(root_libs SYSTEM INTERFACE ${ROOT_INCLUDE_DIRS})
target_link_libraries(root_libs INTERFACE ${ROOT_LIBRARIES})

find_package(GSL REQUIRED)

if(MPI_ENABLED)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

if(PARTONS_ENABLED)
  # PARTONS and its dependencies are installed in ~/.local by scripts/install_partons.sh;
  # appended so that a user-supplied CMAKE_PREFIX_PATH takes precedence
  list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/.local")
  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/GPD_Partons/cmake/Modules")
  find_package(SFML 2.6 CONFIG REQUIRED COMPONENTS system)
  find_package(ElementaryUtils REQUIRED)
  find_package(NumA++ REQUIRED)
  find_package(PARTONS REQUIRED)
  add_library(partons_libs INTERFACE)
  target_include_directories(partons_libs SYSTEM INTERFACE
    ${ElementaryUtils_INCLUDE_DIR} ${NumA++_INCLUDE_DIR} ${PARTONS_INCLUDE_DIR})
  target_link_libraries(partons_libs INTERFACE
    ${PARTONS_LIBRARIES} ${NumA++_LIBRARIES} ${ElementaryUtils_LIBRARIES} sfml-system)
endif()

find_package(Doxygen)

# ---- the code ----------------------------------------------------------------------------------
include(cmake/PhysicsCode.cmake)

add_subdirectory(Integrator)    # also builds the bundled Cuba library (target 'cuba')
add_subdirectory(Tools)         # libMePhys
add_subdirectory(Glauber)
add_subdirectory(Knockout)
add_subdirectory(DIS)           # also DIS/evolution (libEvolution)
if(PARTONS_ENABLED)
  add_subdirectory(GPD_Partons) # libGPD
endif()
add_subdirectory(progs)

if(DOXYGEN_FOUND)
  configure_file(Doxyfile.in "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY)
  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} "${PROJECT_BINARY_DIR}/Doxyfile"
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif()
