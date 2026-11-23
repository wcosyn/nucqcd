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

# ---- Cuba (bundled, autotools) -------------------------------------------------------------
# Built as an external project into the build tree and exposed as the target 'cuba'
# (static libcuba.a + the <Cuba/cuba.h> header); link 'cuba' to use it. libcuba.a is also
# installed into <prefix>/lib.
# Cuba is old C code (typedef enum {false,true} bool): gcc >= 15 defaults to C23, where these
# are keywords, hence -std=gnu17; it does not build in parallel (make -j1).
include(ExternalProject)
get_filename_component(CUBA_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Cuba" REALPATH)  # symlink
set(CUBA_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/Cuba-install")
set(CUBA_LIBRARY "${CUBA_INSTALL_DIR}/lib/libcuba.a")
ExternalProject_Add(cuba_external
  SOURCE_DIR "${CUBA_SOURCE_DIR}"
  BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/Cuba-build"
  INSTALL_DIR "${CUBA_INSTALL_DIR}"
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env "CFLAGS=-O3 -std=gnu17"
                    "${CUBA_SOURCE_DIR}/configure" --prefix=<INSTALL_DIR>
  BUILD_COMMAND make -j1 lib
  INSTALL_COMMAND make -j1 install
  BUILD_BYPRODUCTS "${CUBA_LIBRARY}")
add_library(cuba INTERFACE)
target_include_directories(cuba INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(cuba INTERFACE "${CUBA_LIBRARY}")
add_dependencies(cuba cuba_external)
install(FILES "${CUBA_LIBRARY}" DESTINATION lib)

# numint.hpp includes numint2Cuba.hpp and so <Cuba/cuba.h>: everything that uses numint needs the
# Cuba header (and libcuba.a, whose objects are only linked in where Cuba is actually called).
target_link_libraries(Integrator PUBLIC cuba)
if(STATIC_BUILD)
  target_link_libraries(IntegratorST PUBLIC cuba)
endif()
