# GPD programs (need PARTONS). PARTONS programs only start with partons.properties next to the
# executable, so they are installed into <repo>/GPD_Partons/bin, where
# scripts/install_partons.sh creates that file.
physcode_add_programs(LINK GPD DESTINATION GPD_Partons/bin PROGRAMS
  two_vector_N
  two_vector_electro_N
  two_vector_D
  Nucl_GPD
  diDVCS
  diDVCS_Q2
  diDVCS_gen
  diDVCS_gen_grid)
physcode_add_program(conv_V_test LINK GPD MePhys Glauber DESTINATION GPD_Partons/bin)
physcode_add_program(diDVCS_gen_gridonly LINK root_libs DESTINATION GPD_Partons/bin)  # ROOT only
physcode_add_program(conv_T_test LINK GPD MePhys)    # transversity GPD test, installed into bin/
