physcode_add_programs(LINK Knockout Glauber MePhys PROGRAMS
  rhotest              # tests rho grids (Glauber and momentum distributions)
  rhosingle            # single-node rho transparencies, never used for actual calculations
  rhoset               # single-node rho transparencies with the Glauber sigma set: used for the actual calculations
  rhodint              # single-node deuteron rho cross sections, never used for actual calculations
  pionT_calc12)        # pion transparencies

if(MPI_ENABLED)
  physcode_add_programs(LINK Knockout Glauber MePhys PROGRAMS
    rhottest           # tests rho transparency calculations, never used for proper calculations
    rhod)              # deuteron rho cross section, never used for actual calculations
endif()
