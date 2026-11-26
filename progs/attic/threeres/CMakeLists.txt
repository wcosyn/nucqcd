# NOT BUILT (no add_subdirectory in progs/attic/CMakeLists.txt): the resonance model for inclusive DIS on the deuteron
# (W. Cosyn, W. Melnitchouk and M. Sargsian, Phys. Rev. C 89, 014612 (2014)). InclusiveCrossRes used only the nested
# Romberg routine, was taken out of the DIS library at some point and no longer compiles (DeuteronStructure has a
# different constructor now); the programs are old and "never properly used". Archived in 2026-11 (it was already
# unbuilt); its history is in git. To revive it, port it to the current DeuteronStructure and to numint.
nucqcd_add_programs(LINK DIS MePhys PROGRAMS increstest threeres)
