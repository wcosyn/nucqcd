# Minuit2 (standalone, no ROOT) for the fit programs: provides the INTERFACE target minuit2_libs (headers as
# <Minuit2/FCNBase.h>, libMinuit2.a, libMinuit2Math.a). Programs that link minuit2_libs are skipped when it
# does not exist (MINUIT2_ENABLED=no). Standalone Minuit2 is the ROOT source of Minuit2 with its own CMake
# project (https://github.com/GooFit/Minuit2, same license as ROOT, LGPL-2.1).
# 1. An installed one is used if found: -DMINUIT2_ROOT=<prefix>, $MINUIT2_ROOT, ~/.local (cmake --install of
#    that project), or the system.
# 2. Otherwise (or with -DMINUIT2_DOWNLOAD=yes) the pinned release is downloaded (SHA256) and built in the
#    build tree (a few seconds; -DMINUIT2_TARBALL=<file> for a local copy, without network).
option(MINUIT2_ENABLED "Use Minuit2 for the fit programs (installed, or downloaded and built)" ON)
option(MINUIT2_DOWNLOAD "Download and build Minuit2 even if an installed one is found" OFF)
set(MINUIT2_ROOT "" CACHE PATH "Prefix of an installed standalone Minuit2 (lib/libMinuit2.a, include/Minuit2/Minuit2/*.h)")
set(MINUIT2_TARBALL "" CACHE FILEPATH "Local Minuit2-v6-40-02.tar.gz to use instead of downloading it")

if(MINUIT2_ENABLED)
  find_package(Threads REQUIRED)
  set(MINUIT2_HINTS ${MINUIT2_ROOT} $ENV{MINUIT2_ROOT} "$ENV{HOME}/.local")
  if(NOT MINUIT2_DOWNLOAD)
    # Only the standalone layout is accepted (<prefix>/include/Minuit2/Minuit2/FCNBase.h): ROOT's own
    # libMinuit2 needs ROOT's libraries, so the search must not look at ROOT's directories (they are on the
    # CMake prefix path after thisroot.sh). First the explicit prefixes, then the system directories.
    foreach(mode HINTS_ONLY SYSTEM)
      if(NOT (MINUIT2_LIBRARY AND MINUIT2_MATH_LIBRARY AND MINUIT2_INCLUDE_DIR))
        unset(MINUIT2_LIBRARY CACHE)
        unset(MINUIT2_MATH_LIBRARY CACHE)
        unset(MINUIT2_INCLUDE_DIR CACHE)
        if(mode STREQUAL "HINTS_ONLY")
          set(search_options HINTS ${MINUIT2_HINTS} NO_DEFAULT_PATH)
        else()
          set(search_options NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
        endif()
        find_library(MINUIT2_LIBRARY NAMES Minuit2 ${search_options} PATH_SUFFIXES lib lib64)
        find_library(MINUIT2_MATH_LIBRARY NAMES Minuit2Math ${search_options} PATH_SUFFIXES lib lib64)
        find_path(MINUIT2_INCLUDE_DIR NAMES Minuit2/FCNBase.h ${search_options} PATH_SUFFIXES include/Minuit2)
      endif()
    endforeach()
  endif()
  add_library(minuit2_libs INTERFACE)
  if(NOT MINUIT2_DOWNLOAD AND MINUIT2_LIBRARY AND MINUIT2_MATH_LIBRARY AND MINUIT2_INCLUDE_DIR)
    message(STATUS "Minuit2: using the installed ${MINUIT2_LIBRARY} (${MINUIT2_INCLUDE_DIR})")
    target_include_directories(minuit2_libs SYSTEM INTERFACE "${MINUIT2_INCLUDE_DIR}")
    target_link_libraries(minuit2_libs INTERFACE "${MINUIT2_LIBRARY}" "${MINUIT2_MATH_LIBRARY}" Threads::Threads)
  else()
    set(MINUIT2_VERSION v6-40-02)  # = ROOT 6.40.02
    set(MINUIT2_SHA256 72049214d0c6cdc5389ab7d8acf474638b68eaa1e66dd01389b5e7d70657fef8)
    if(MINUIT2_TARBALL)
      set(MINUIT2_URL "${MINUIT2_TARBALL}")
    else()
      set(MINUIT2_URL "https://github.com/GooFit/Minuit2/archive/refs/tags/${MINUIT2_VERSION}.tar.gz")
    endif()
    message(STATUS "Minuit2: no installed Minuit2 found; ${MINUIT2_VERSION} will be downloaded from ${MINUIT2_URL} and built")
    set(MINUIT2_INSTALL_DIR "${CMAKE_BINARY_DIR}/Minuit2-install")
    include(ExternalProject)
    ExternalProject_Add(minuit2_external
      URL "${MINUIT2_URL}"
      URL_HASH SHA256=${MINUIT2_SHA256}
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
      SOURCE_DIR "${CMAKE_BINARY_DIR}/Minuit2-src"
      BINARY_DIR "${CMAKE_BINARY_DIR}/Minuit2-build"
      INSTALL_DIR "${MINUIT2_INSTALL_DIR}"
      CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_CXX_STANDARD=17
                 -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target Minuit2 Minuit2Math --parallel
      BUILD_BYPRODUCTS "${MINUIT2_INSTALL_DIR}/lib/libMinuit2.a" "${MINUIT2_INSTALL_DIR}/lib/libMinuit2Math.a")
    # (the headers are installed at build time; the directory need not exist when configuring)
    target_include_directories(minuit2_libs SYSTEM INTERFACE "${MINUIT2_INSTALL_DIR}/include/Minuit2")
    target_link_libraries(minuit2_libs INTERFACE "${MINUIT2_INSTALL_DIR}/lib/libMinuit2.a"
                          "${MINUIT2_INSTALL_DIR}/lib/libMinuit2Math.a" Threads::Threads)
    add_dependencies(minuit2_libs minuit2_external)
  endif()
else()
  message(STATUS "Minuit2 disabled: the fit programs are not built")
endif()
