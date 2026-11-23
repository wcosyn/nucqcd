# GPD programs (need PARTONS). PARTONS programs only start with partons.properties next to the
# executable, so they are installed into <repo>/GPD_Partons/bin, where
# scripts/install_partons.sh creates that file.
nucqcd_add_programs(LINK GPD DESTINATION GPD_Partons/bin PROGRAMS
  two_vector_N
  two_vector_electro_N
  two_vector_D
  Nucl_GPD
  diDVCS
  diDVCS_Q2
  diDVCS_gen
  diDVCS_gen_grid)
nucqcd_add_program(conv_V_test LINK GPD MePhys Glauber DESTINATION GPD_Partons/bin)
nucqcd_add_program(diDVCS_gen_gridonly LINK ToolsCore DESTINATION GPD_Partons/bin)  # without PARTONS: reads the grids only
nucqcd_add_program(conv_T_test LINK GPD MePhys)    # transversity GPD test, installed into bin/
