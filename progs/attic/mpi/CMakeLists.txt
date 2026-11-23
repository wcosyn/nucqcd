# TMPI (from strangecalc: one object per program initializes MPI, distributes independent
# calculations over the ranks and gathers the results) and the two programs that used it: rhod
# (deuteron rho cross sections on a grid in nu and t or z) and rhottest (rho transparency test).
# Neither was used for actual calculations. Moved here from Tools and progs/RhoTransparency in
# 2026-11 (roadmap step 6, "Retire MPI"); the single-process programs in progs/RhoTransparency
# (rhoset, rhosingle, rhodint) do not need MPI. Built only when MPI is found.
find_package(MPI QUIET COMPONENTS CXX)
if(MPI_CXX_FOUND)
  add_library(attic_tmpi STATIC TMPI.cpp)
  target_link_libraries(attic_tmpi PUBLIC MPI::MPI_CXX)
  if(NOT BUILD_PROGS)
    set_target_properties(attic_tmpi PROPERTIES EXCLUDE_FROM_ALL ON)
  endif()
  nucqcd_add_programs(LINK Knockout Glauber MePhys attic_tmpi PROGRAMS rhod rhottest)
endif()
