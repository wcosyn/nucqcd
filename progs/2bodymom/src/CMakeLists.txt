# Camille's programs for the 2-nucleon knockout calculations.

# sources shared by all programs below: compiled once
add_library(twobody_common OBJECT
  sampler_HALLA.cpp sampler_PRL99_072501.cpp sampler_HALLB.cpp parser.cpp event.cpp
  GlauberGridThick_SCX.cpp GlauberGridThick_SEL.cpp ClassGridThickNucleus_SCX.cpp)
target_include_directories(twobody_common PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(twobody_common PUBLIC Knockout Glauber MePhys cuba)
if(NOT BUILD_PROGS)
  set_target_properties(twobody_common PROPERTIES EXCLUDE_FROM_ALL ON)
endif()

physcode_add_program(sampler LINK twobody_common)
physcode_add_program(sampler_test LINK twobody_common)
physcode_add_program(2bodymom_new SOURCES 2bodymom.cpp LINK twobody_common)
physcode_add_program(SCX_test_new SOURCES SCX_test.cpp LINK twobody_common)
physcode_add_program(class_SCX LINK twobody_common)

# wrapper library used by MeanFieldNucleus.py in this directory (loaded from <repo>/lib)
add_library(MeanFieldNucleusWrapper SHARED MeanFieldNucleusWrapper.cpp)
target_link_libraries(MeanFieldNucleusWrapper PRIVATE Glauber)
if(BUILD_PROGS)
  install(TARGETS MeanFieldNucleusWrapper LIBRARY DESTINATION lib)
else()
  set_target_properties(MeanFieldNucleusWrapper PROPERTIES EXCLUDE_FROM_ALL ON)
endif()
