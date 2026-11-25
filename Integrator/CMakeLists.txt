# libIntegrator: numint (templated quadratures) and adaptive cubature.
nucqcd_add_library(Integrator
  SOURCES
    numint/cube_romb.cpp
    numint/cube_tr.cpp
    numint/fubini.cpp
    numint/quad_cc.cpp
    numint/quad_romb.cpp
    numint/quad_tr.cpp
    numint/vfubini.cpp
    adaptive/cubature.c)

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
