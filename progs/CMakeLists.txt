# Programs. They are only built by 'make' (and installed into <repo>/bin) with -DBUILD_PROGS=yes;
# otherwise build one with 'make <name>'. See cmake/PhysicsCode.cmake for physcode_add_program(s).
add_subdirectory(Obstest)          # A(e,e'p) observables
add_subdirectory(twonucl)          # old work in progress on a 2N knockout model with correlated pairs
# add_subdirectory(threeres)       # old, never properly used: resonance model for inclusive DIS
add_subdirectory(onenucl)          # A(e,e'p) cross sections and transparencies
add_subdirectory(SCX_test)         # checks of charge-exchange densities (Camille; roundoff problems in gsl)
add_subdirectory(RhoTransparency)  # rho (and pion) transparencies
add_subdirectory(poldeut)          # light-front deuteron spectator tagging
add_subdirectory(pionprog)         # writes out pion Glauber parameters
add_subdirectory(pionestimate)     # uncertainty from rho -> 2 pion decay in rho transparencies
add_subdirectory(neutrinoQE)       # neutrino QE double/single differential cross sections (RMF, RFG; nucleon, nucleus)
add_subdirectory(momdistr)         # momentum distributions with relativistic mean-field wave functions
add_subdirectory(LDRD_test)        # LDRD project test, tagging with He3, never properly used
add_subdirectory(integrator)       # test of the integration routines (volume of part of a sphere)
add_subdirectory(incldis)          # FSI in inclusive deuteron DIS (unpolarized and Azz)
add_subdirectory(gammatest)        # gamma-matrix identities and trace checks
add_subdirectory(DQEinc)           # inclusive deuteron QE with FSI, work in progress
add_subdirectory(cuba_test)        # tests of the Cuba library
add_subdirectory(camille)          # Camille's test programs
add_subdirectory(2bodymom/src)     # Camille's programs for the 2-nucleon knockout calculations
add_subdirectory(deepsdis)         # spectator tagging (DEEPS/BONuS) and deuteron tensor structure
if(PARTONS_ENABLED)
  add_subdirectory(GPD)            # transversity and vector GPDs (spin 1), diDVCS
endif()
add_subdirectory(F2_largex_fit)    # x > 1 nuclear F2 fit, cross-checks with Adam Freese
