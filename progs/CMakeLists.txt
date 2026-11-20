# Programs. They are only built by 'make' (and installed into <repo>/bin) with -DBUILD_PROGS=yes;
# otherwise build one with 'make <name>'. See cmake/NucQCD.cmake for nucqcd_add_program(s).
add_subdirectory(Obstest)          # A(e,e'p) observables
# add_subdirectory(threeres)       # old, never properly used: resonance model for inclusive DIS
add_subdirectory(onenucl)          # A(e,e'p) cross sections and transparencies
add_subdirectory(SCX_test)         # checks of charge-exchange densities (Camille; roundoff problems in gsl)
add_subdirectory(RhoTransparency)  # rho (and pion) transparencies
add_subdirectory(poldeut)          # light-front deuteron spectator tagging
add_subdirectory(pionprog)         # writes out pion Glauber parameters
add_subdirectory(pionestimate)     # uncertainty from rho -> 2 pion decay in rho transparencies
add_subdirectory(neutrinoQE)       # neutrino QE double/single differential cross sections (RMF, RFG; nucleon, nucleus)
add_subdirectory(momdistr)         # momentum distributions with relativistic mean-field wave functions
add_subdirectory(incldis)          # FSI in inclusive deuteron DIS (unpolarized and Azz)
add_subdirectory(2bodymom/src)     # Camille's programs for the 2-nucleon knockout calculations
add_subdirectory(deepsdis)         # spectator tagging (DEEPS/BONuS) and deuteron tensor structure
if(PARTONS_ENABLED)
  add_subdirectory(GPD)            # transversity and vector GPDs (spin 1), diDVCS
endif()
add_subdirectory(F2_largex_fit)    # x > 1 nuclear F2 fit, cross-checks with Adam Freese
add_subdirectory(attic)            # old, one-off and personal programs, still built (see attic/CMakeLists.txt)
