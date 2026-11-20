# Creates the sandbox share directory used by the tests (run as: cmake -DSRC=.. -DDEST=.. -P).
#   - read-only data directories: symlinks into <repo>/share
#   - input/: copied once (FsiCorrelator writes CorrGrid files into it), so the repository copy
#     is never modified
#   - grids/: emptied on every run, so FSI grids are computed by the code under test instead
#     of being read from a cache (and share/grids/ of the user is never touched)
if(NOT SRC OR NOT DEST)
  message(FATAL_ERROR "usage: cmake -DSRC=<repo>/share -DDEST=<sandbox> -P make_test_share.cmake")
endif()

file(MAKE_DIRECTORY "${DEST}")
file(GLOB entries RELATIVE "${SRC}" "${SRC}/*")
foreach(entry IN LISTS entries)
  if(entry STREQUAL "grids" OR NOT IS_DIRECTORY "${SRC}/${entry}")
    continue()
  endif()
  if(entry STREQUAL "input")
    if(NOT EXISTS "${DEST}/${entry}/.copied")
      file(REMOVE_RECURSE "${DEST}/${entry}")
      file(COPY "${SRC}/${entry}" DESTINATION "${DEST}")
      file(TOUCH "${DEST}/${entry}/.copied")
    endif()
  elseif(NOT EXISTS "${DEST}/${entry}")
    file(CREATE_LINK "${SRC}/${entry}" "${DEST}/${entry}" SYMBOLIC)
  endif()
endforeach()

file(REMOVE_RECURSE "${DEST}/grids")
file(MAKE_DIRECTORY "${DEST}/grids")
message(STATUS "test share sandbox ready: ${DEST}")
