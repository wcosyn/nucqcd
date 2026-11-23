# Helper functions shared by all CMakeLists.txt of nucqcd.
#
#   nucqcd_add_library(<name> SOURCES <src>... [LINK <lib>...])
#       Shared library lib<name>.so (version 1.0.0), installed into <prefix>/lib. The public
#       include directory is the repository root, so headers are included by their path from
#       there (#include <Tools/Core/FourVector.hpp>, <Glauber/FastParticle.hpp>, ...). LINK libraries are public dependencies, and
#       must resolve every symbol the library uses (linked with --no-undefined on Linux).
#       With STATIC_BUILD also a static lib<name>.a (target <name>ST).
#
#   nucqcd_add_program(<name> [SOURCES <src>...] [LINK <lib>...] [DESTINATION <dir>])
#       Executable <name> (default source <name>.cpp). With BUILD_PROGS it is built by
#       'make' and installed into <prefix>/<DESTINATION> (default: bin); otherwise it is
#       only built on request ('make <name>').
#
#   nucqcd_add_programs(LINK <lib>... PROGRAMS <name>... [DESTINATION <dir>])
#       Shortcut for several single-source programs with the same libraries.

include_guard(GLOBAL)

function(nucqcd_add_library name)
  cmake_parse_arguments(ARG "" "" "SOURCES;LINK" ${ARGN})
  add_library(${name} SHARED ${ARG_SOURCES})
  target_include_directories(${name} PUBLIC "${PROJECT_SOURCE_DIR}")
  target_link_libraries(${name} PUBLIC ${ARG_LINK})
  # Every symbol must be resolved by the library itself or by the libraries it links: this keeps the
  # LINK lists (and so the layering of the libraries) complete. macOS does this by default.
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_options(${name} PRIVATE "LINKER:--no-undefined")
  endif()
  # In the build tree a library must find its dependencies in build/lib, not the installed copies
  # in <prefix>/lib (otherwise the tests load stale libraries until 'make install'); 'make install'
  # replaces this RPATH with CMAKE_INSTALL_RPATH.
  set_target_properties(${name} PROPERTIES VERSION 1.0.0 SOVERSION 1 BUILD_WITH_INSTALL_RPATH OFF)
  install(TARGETS ${name} LIBRARY DESTINATION lib)

  if(STATIC_BUILD)
    get_target_property(all_sources ${name} SOURCES)
    add_library(${name}ST STATIC ${all_sources})
    target_include_directories(${name}ST PUBLIC "${PROJECT_SOURCE_DIR}")
    target_link_libraries(${name}ST PUBLIC ${ARG_LINK})
    set_target_properties(${name}ST PROPERTIES OUTPUT_NAME ${name})
    install(TARGETS ${name}ST ARCHIVE DESTINATION lib)
  endif()
endfunction()

function(nucqcd_add_program name)
  cmake_parse_arguments(ARG "" "DESTINATION" "SOURCES;LINK" ${ARGN})
  if(NOT ARG_SOURCES)
    set(ARG_SOURCES ${name}.cpp)
  endif()
  if(NOT ARG_DESTINATION)
    set(ARG_DESTINATION bin)
  endif()
  add_executable(${name} ${ARG_SOURCES})
  target_link_libraries(${name} PRIVATE ${ARG_LINK})
  if(BUILD_PROGS)
    install(TARGETS ${name} RUNTIME DESTINATION ${ARG_DESTINATION})
  else()
    set_target_properties(${name} PROPERTIES EXCLUDE_FROM_ALL ON)
  endif()
endfunction()

function(nucqcd_add_programs)
  cmake_parse_arguments(ARG "" "DESTINATION" "LINK;PROGRAMS" ${ARGN})
  foreach(prog IN LISTS ARG_PROGRAMS)
    nucqcd_add_program(${prog} LINK ${ARG_LINK} DESTINATION "${ARG_DESTINATION}")
  endforeach()
endfunction()
