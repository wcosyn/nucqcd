# fits with ROOT's TMinuit: only built with ROOT
nucqcd_add_programs(LINK DIS Evolution MePhys root_libs PROGRAMS F2_fit F2_bootstrap)
nucqcd_add_program(EMC_fit LINK MePhys root_libs)
