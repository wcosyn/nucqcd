# fits with Minuit2 (cmake/Minuit2.cmake)
nucqcd_add_programs(LINK DIS Evolution MePhys minuit2_libs PROGRAMS F2_fit F2_bootstrap)
nucqcd_add_program(EMC_fit LINK MePhys minuit2_libs)
