# libIntegrator: numint (templated quadratures); the adaptive cubature of S. G. Johnson is downloaded and compiled in (below).
nucqcd_add_library(Integrator
  SOURCES
    numint/cube_romb.cpp
    numint/cube_tr.cpp
    numint/fubini.cpp
    numint/quad_cc.cpp
    numint/quad_romb.cpp
    numint/quad_tr.cpp
    numint/vfubini.cpp)

# ---- cubature (not bundled) -------------------------------------------------------------------
# hcubature of Steven G. Johnson's cubature 1.0.4 (https://github.com/stevengj/cubature, GPL-2-or-later), the engine of
# numint::cube_adaptive. It is downloaded at configure time (URL and SHA256 pinned; -DCUBATURE_TARBALL=<file> uses a local
# copy of v1.0.4.tar.gz, for machines without network), patched (cmake/patches/cubature-1.0.4-nucqcd.patch: hcubature_minEval, the
# minimum number of function evaluations that the 2010 adapt_integrate, which nucqcd bundled until 2026-11, had, and the 2010 rule to
# choose the dimension to split, which keeps all results unchanged) and its one source file is compiled into libIntegrator.
# -DCUBATURE_UPSTREAM_SPLIT=yes uses the rule of cubature 1.0.4 instead (results of integrands with steps change, see the ROADMAP).
# A cubature installed on the system (libcubature-dev) cannot be used: it has no minEval. The unpacked and patched source is redone
# when the patch changes. <cubature/cubature.h> is provided from <build>/Integrator/cubature-include.
set(CUBATURE_VERSION 1.0.4)
set(CUBATURE_SHA256 cd4899de0b047a9d220cfb751a8bdbb8fd0c97c1c894d07523b75168e6426f60)
set(CUBATURE_TARBALL "" CACHE FILEPATH "Local cubature v1.0.4.tar.gz to use instead of downloading it")
option(CUBATURE_UPSTREAM_SPLIT "Use the rule of cubature 1.0.4 to choose the dimension to split (changes results), not the 2010 one" OFF)
set(CUBATURE_PATCH "${PROJECT_SOURCE_DIR}/cmake/patches/cubature-${CUBATURE_VERSION}-nucqcd.patch")
set(CUBATURE_SRC_DIR "${CMAKE_CURRENT_BINARY_DIR}/cubature-src")
set(CUBATURE_STAMP "${CUBATURE_SRC_DIR}/.nucqcd-patch-sha256")
set(CUBATURE_INCLUDE_STAGE "${CMAKE_CURRENT_BINARY_DIR}/cubature-include")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CUBATURE_PATCH}")
file(SHA256 "${CUBATURE_PATCH}" CUBATURE_PATCH_SHA256)
set(CUBATURE_UP_TO_DATE FALSE)
if(EXISTS "${CUBATURE_SRC_DIR}/hcubature.c" AND EXISTS "${CUBATURE_STAMP}")
  file(READ "${CUBATURE_STAMP}" CUBATURE_STAMP_CONTENT)
  if(CUBATURE_STAMP_CONTENT STREQUAL "${CUBATURE_PATCH_SHA256}")
    set(CUBATURE_UP_TO_DATE TRUE)
  endif()
endif()
if(NOT CUBATURE_UP_TO_DATE)
  find_program(PATCH_EXECUTABLE patch REQUIRED)
  set(CUBATURE_ARCHIVE "${CMAKE_CURRENT_BINARY_DIR}/cubature-${CUBATURE_VERSION}.tar.gz")
  if(CUBATURE_TARBALL)
    file(SHA256 "${CUBATURE_TARBALL}" CUBATURE_TARBALL_SHA256)
    if(NOT CUBATURE_TARBALL_SHA256 STREQUAL CUBATURE_SHA256)
      message(FATAL_ERROR "cubature: ${CUBATURE_TARBALL} has SHA256 ${CUBATURE_TARBALL_SHA256}, expected ${CUBATURE_SHA256}")
    endif()
    configure_file("${CUBATURE_TARBALL}" "${CUBATURE_ARCHIVE}" COPYONLY)
  else()
    message(STATUS "cubature: getting ${CUBATURE_VERSION} from github.com/stevengj/cubature (download is skipped if it is there)")
    file(DOWNLOAD "https://github.com/stevengj/cubature/archive/refs/tags/v${CUBATURE_VERSION}.tar.gz" "${CUBATURE_ARCHIVE}"
         EXPECTED_HASH SHA256=${CUBATURE_SHA256} STATUS CUBATURE_DOWNLOAD_STATUS)
    list(GET CUBATURE_DOWNLOAD_STATUS 0 CUBATURE_DOWNLOAD_CODE)
    if(NOT CUBATURE_DOWNLOAD_CODE EQUAL 0)
      message(FATAL_ERROR "cubature: the download failed (${CUBATURE_DOWNLOAD_STATUS}); download v${CUBATURE_VERSION}.tar.gz by hand "
                          "and configure with -DCUBATURE_TARBALL=<file>")
    endif()
  endif()
  set(CUBATURE_UNPACK_DIR "${CMAKE_CURRENT_BINARY_DIR}/cubature-unpack")
  file(REMOVE_RECURSE "${CUBATURE_UNPACK_DIR}" "${CUBATURE_SRC_DIR}")
  file(MAKE_DIRECTORY "${CUBATURE_UNPACK_DIR}")
  execute_process(COMMAND "${CMAKE_COMMAND}" -E tar xzf "${CUBATURE_ARCHIVE}" WORKING_DIRECTORY "${CUBATURE_UNPACK_DIR}"
                  RESULT_VARIABLE CUBATURE_UNPACK_RESULT)
  if(NOT CUBATURE_UNPACK_RESULT EQUAL 0 OR NOT EXISTS "${CUBATURE_UNPACK_DIR}/cubature-${CUBATURE_VERSION}/hcubature.c")
    message(FATAL_ERROR "cubature: could not unpack ${CUBATURE_ARCHIVE}")
  endif()
  file(RENAME "${CUBATURE_UNPACK_DIR}/cubature-${CUBATURE_VERSION}" "${CUBATURE_SRC_DIR}")
  execute_process(COMMAND "${PATCH_EXECUTABLE}" -p1 -i "${CUBATURE_PATCH}"
                  WORKING_DIRECTORY "${CUBATURE_SRC_DIR}" RESULT_VARIABLE CUBATURE_PATCH_RESULT)
  if(NOT CUBATURE_PATCH_RESULT EQUAL 0)
    message(FATAL_ERROR "cubature: applying ${CUBATURE_PATCH} failed")
  endif()
  file(WRITE "${CUBATURE_STAMP}" "${CUBATURE_PATCH_SHA256}")
endif()
configure_file("${CUBATURE_SRC_DIR}/cubature.h" "${CUBATURE_INCLUDE_STAGE}/cubature/cubature.h" COPYONLY)
set(CUBATURE_COMPILE_FLAGS "-w")  # third-party code
if(CUBATURE_UPSTREAM_SPLIT)
  set(CUBATURE_COMPILE_FLAGS "${CUBATURE_COMPILE_FLAGS} -DNUCQCD_UPSTREAM_SPLIT")
  message(STATUS "cubature: using the split rule of cubature 1.0.4 (results differ from the 2010 rule)")
endif()
set_source_files_properties("${CUBATURE_SRC_DIR}/hcubature.c" PROPERTIES COMPILE_FLAGS "${CUBATURE_COMPILE_FLAGS}")
target_sources(Integrator PRIVATE "${CUBATURE_SRC_DIR}/hcubature.c")
target_include_directories(Integrator PUBLIC "${CUBATURE_INCLUDE_STAGE}")
if(STATIC_BUILD)
  target_sources(IntegratorST PRIVATE "${CUBATURE_SRC_DIR}/hcubature.c")
  target_include_directories(IntegratorST PUBLIC "${CUBATURE_INCLUDE_STAGE}")
endif()

# ---- Cuba (not bundled) --------------------------------------------------------------------
# Exposed as the INTERFACE target 'cuba' (library + the <Cuba/cuba.h> header); link 'cuba' to use it.
# 1. An installed Cuba is used if CMake finds one: -DCUBA_ROOT=<prefix>, $CUBA_ROOT, ~/.local or the
#    system (libcuba, cuba.h). scripts/install_cuba.sh installs the patched version into ~/.local.
#    Cuba 4.2.1 and 4.2.2 have a bug in Suave (see cmake/patches/cuba-4.2.1-nucqcd.patch), which only
#    shows with a flatness below 5; an installed Cuba that does not carry the patch (a system
#    package) is used as it is, and CMake says so. unit_cuba tests Suave and catches it.
# 2. Otherwise (or with -DCUBA_DOWNLOAD=yes) Cuba 4.2.1 is downloaded (URL and SHA256 pinned; -DCUBA_TARBALL=<file> uses a local
#    copy, for machines without network), patched and built in the build tree as an external project.
#    Cuba is old C code (typedef enum {false,true} bool): gcc >= 15 defaults to C23, where these are
#    keywords, hence -std=gnu17; it does not build in parallel (make -j1).
# <Cuba/cuba.h> is provided from <build>/Integrator/cuba-include.
set(CUBA_ROOT "" CACHE PATH "Prefix of an installed Cuba (lib/libcuba, include/cuba.h)")
option(CUBA_DOWNLOAD "Download and build Cuba even if an installed one is found" OFF)
set(CUBA_TARBALL "" CACHE FILEPATH "Local Cuba-4.2.1.tar.gz to use instead of downloading it")
set(CUBA_INCLUDE_STAGE "${CMAKE_CURRENT_BINARY_DIR}/cuba-include")
set(CUBA_SEARCH_HINTS ${CUBA_ROOT} $ENV{CUBA_ROOT} "$ENV{HOME}/.local")
if(NOT CUBA_DOWNLOAD)
  find_library(CUBA_FOUND_LIBRARY NAMES cuba HINTS ${CUBA_SEARCH_HINTS} PATH_SUFFIXES lib lib64)
  find_path(CUBA_FOUND_INCLUDE NAMES cuba.h HINTS ${CUBA_SEARCH_HINTS} PATH_SUFFIXES include)
endif()
add_library(cuba INTERFACE)
target_include_directories(cuba INTERFACE "${CUBA_INCLUDE_STAGE}")
if(NOT CUBA_DOWNLOAD AND CUBA_FOUND_LIBRARY AND CUBA_FOUND_INCLUDE)
  message(STATUS "Cuba: using the installed ${CUBA_FOUND_LIBRARY} (${CUBA_FOUND_INCLUDE}/cuba.h)")
  if(NOT CUBA_FOUND_LIBRARY MATCHES "^$ENV{HOME}/.local/")
    message(STATUS "Cuba: a Cuba outside ~/.local may lack the Suave patch (cmake/patches); "
                   "run scripts/install_cuba.sh, or configure with -DCUBA_DOWNLOAD=yes to get the patched one")
  endif()
  configure_file("${CUBA_FOUND_INCLUDE}/cuba.h" "${CUBA_INCLUDE_STAGE}/Cuba/cuba.h" COPYONLY)
  target_link_libraries(cuba INTERFACE "${CUBA_FOUND_LIBRARY}")
else()
  find_program(PATCH_EXECUTABLE patch REQUIRED)
  set(CUBA_VERSION 4.2.1)
  set(CUBA_SHA256 3789c458193440e73a603f6649acbe8f1903d61203905e1b63c00ef0db8580a8)
  if(CUBA_TARBALL)
    set(CUBA_URL "${CUBA_TARBALL}")
  else()
    set(CUBA_URL "https://feynarts.de/cuba/Cuba-${CUBA_VERSION}.tar.gz")
  endif()
  message(STATUS "Cuba: no installed Cuba found; ${CUBA_VERSION} will be downloaded from ${CUBA_URL} and built")
  set(CUBA_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/Cuba-install")
  set(CUBA_LIBRARY "${CUBA_INSTALL_DIR}/lib/libcuba.a")
  include(ExternalProject)
  ExternalProject_Add(cuba_external
    URL "${CUBA_URL}"
    URL_HASH SHA256=${CUBA_SHA256}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/Cuba-src"
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/Cuba-build"
    INSTALL_DIR "${CUBA_INSTALL_DIR}"
    PATCH_COMMAND "${PATCH_EXECUTABLE}" -p1 -i "${PROJECT_SOURCE_DIR}/cmake/patches/cuba-4.2.1-nucqcd.patch"
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env "CFLAGS=-O3 -std=gnu17"
                      "<SOURCE_DIR>/configure" --prefix=<INSTALL_DIR>
    BUILD_COMMAND make -j1 lib
    INSTALL_COMMAND make -j1 install
            COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/cuba.h" "${CUBA_INCLUDE_STAGE}/Cuba/cuba.h"
    BUILD_BYPRODUCTS "${CUBA_LIBRARY}" "${CUBA_INCLUDE_STAGE}/Cuba/cuba.h")
  target_link_libraries(cuba INTERFACE "${CUBA_LIBRARY}")
  add_dependencies(cuba cuba_external)
  install(FILES "${CUBA_LIBRARY}" DESTINATION lib)
endif()

# numint.hpp includes numint2Cuba.hpp and so <Cuba/cuba.h>: everything that uses numint needs the
# Cuba header (and libcuba.a, whose objects are only linked in where Cuba is actually called).
target_link_libraries(Integrator PUBLIC cuba)
if(STATIC_BUILD)
  target_link_libraries(IntegratorST PUBLIC cuba)
endif()
