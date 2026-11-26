# Old, one-off and personal programs. They are kept compiling (with -DBUILD_PROGS=yes, like the
# rest of progs/), so they don't rot silently, but are not maintained or documented.
add_subdirectory(twonucl)          # old work in progress on a 2N knockout model with correlated pairs
add_subdirectory(LDRD_test)        # LDRD project test, tagging with He3, never properly used
add_subdirectory(integrator)       # test of the integration routines (volume of part of a sphere)
add_subdirectory(gammatest)        # gamma-matrix identities and trace checks
add_subdirectory(DQEinc)           # inclusive deuteron QE with FSI, work in progress
add_subdirectory(cuba_test)        # tests of the Cuba library
add_subdirectory(camille)          # Camille's test programs
add_subdirectory(rhodeuteron)      # rho production on the deuteron (RhoDeuteron, rhodint): an abandoned calculation
add_subdirectory(mpi)              # TMPI (MPI helper from strangecalc) and the MPI rho programs rhod, rhottest
# progs/attic/threeres (InclusiveCrossRes, threeres, increstest: resonance model for inclusive DIS on the deuteron) is not
# built: it does not compile any more, see the CMakeLists.txt there.
