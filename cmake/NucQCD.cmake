# Helper functions shared by all CMakeLists.txt of nucqcd.
#
#   nucqcd_add_library(<name> SOURCES <src>... [LINK <lib>...] [LINKDEFS <XLinkDef.h>...])
#       Shared library lib<name>.so (version 1.0.0), its directory as public include
#       directory, installed into <prefix>/lib. LINK libraries are public dependencies.
#       LINKDEFS: one ROOT dictionary per LinkDef file (see nucqcd_root_dictionary).
#       With STATIC_BUILD also a static lib<name>.a (target <name>ST).
#
#   nucqcd_add_program(<name> [SOURCES <src>...] [LINK <lib>...] [DESTINATION <dir>])
#       Executable <name> (default source <name>.cpp). With BUILD_PROGS it is built by
#       'make' and installed into <prefix>/<DESTINATION> (default: bin); otherwise it is
#       only built on request ('make <name>').
#
#   nucqcd_add_programs(LINK <lib>... PROGRAMS <name>... [DESTINATION <dir>])
#       Shortcut for several single-source programs with the same libraries.
#
#   nucqcd_root_dictionary(<target> <XLinkDef.h>...)
#       Runs rootcling for each X: X.h or X.hpp + XLinkDef.h -> XDict.cxx (+ XDict_rdict.pcm,
#       copied next to the library in the build tree and installed next to it), and adds
#       XDict.cxx to <target>.

include_guard(GLOBAL)

function(nucqcd_root_dictionary target)
  set(incflags "-I$<JOIN:$<TARGET_PROPERTY:${target},INCLUDE_DIRECTORIES>,;-I>")
  foreach(linkdef IN LISTS ARGN)
    string(REGEX REPLACE "LinkDef\\.h$" "" cls "${linkdef}")
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${cls}.h")
      set(header "${cls}.h")
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${cls}.hpp")
      set(header "${cls}.hpp")
    else()
      message(FATAL_ERROR "No header ${cls}.h(pp) for ${linkdef}")
    endif()
    set(dict "${CMAKE_CURRENT_BINARY_DIR}/${cls}Dict.cxx")
    set(pcm "${CMAKE_CURRENT_BINARY_DIR}/${cls}Dict_rdict.pcm")
    # ROOT looks for the .pcm next to the library, so also copy it to the build tree's lib/
    # (for tests and programs run from the build directory); install() covers <repo>/lib.
    set(libpcm "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${cls}Dict_rdict.pcm")
    add_custom_command(
      OUTPUT "${dict}" "${pcm}" "${libpcm}"
      COMMAND ${ROOT_rootcling_CMD} -f "${dict}" "${incflags}" "${header}" "${linkdef}"
      COMMAND ${CMAKE_COMMAND} -E copy "${pcm}" "${libpcm}"
      DEPENDS "${header}" "${linkdef}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      COMMENT "Generating ROOT dictionary ${cls}Dict"
      COMMAND_EXPAND_LISTS VERBATIM)
    target_sources(${target} PRIVATE "${dict}")
    install(FILES "${pcm}" DESTINATION lib)
  endforeach()
endfunction()

function(nucqcd_add_library name)
  cmake_parse_arguments(ARG "" "" "SOURCES;LINK;LINKDEFS" ${ARGN})
  add_library(${name} SHARED ${ARG_SOURCES})
  target_include_directories(${name} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(${name} PUBLIC ${ARG_LINK})
  # In the build tree a library must find its dependencies in build/lib, not the installed copies
  # in <prefix>/lib (otherwise the tests load stale libraries until 'make install'); 'make install'
  # replaces this RPATH with CMAKE_INSTALL_RPATH.
  set_target_properties(${name} PROPERTIES VERSION 1.0.0 SOVERSION 1 BUILD_WITH_INSTALL_RPATH OFF)
  if(ARG_LINKDEFS)
    nucqcd_root_dictionary(${name} ${ARG_LINKDEFS})
  endif()
  install(TARGETS ${name} LIBRARY DESTINATION lib)

  if(STATIC_BUILD)
    get_target_property(all_sources ${name} SOURCES)   # includes the dictionaries
    add_library(${name}ST STATIC ${all_sources})
    target_include_directories(${name}ST PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
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
