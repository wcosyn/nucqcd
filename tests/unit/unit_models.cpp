// Model names as enumerations (Tools/WaveFunctions/WaveFunctionModels.hpp,
// Tools/StructureFunctions/StructureFunctionModel.hpp): parsing in one place, early errors that list the
// valid names, and the overloads that take the enumeration.
#include <test_support.hpp>

#include <catch2/matchers/catch_matchers_string.hpp>

#include <Tools/StructureFunctions/NucleonStructure.hpp>
#include <Tools/WaveFunctions/OldDeuteron.hpp>
#include <Tools/WaveFunctions/TDeuteron.hpp>
#include <Tools/WaveFunctions/WaveFunctionModels.hpp>

#include <memory>
#include <stdexcept>
#include <string>
using namespace nucqcd;

using Catch::Matchers::ContainsSubstring;

TEST_CASE("DeuteronModel: names round trip, case insensitive, unknown names list the valid ones", "[models]") {
  CHECK(deuteron_models().size() == 11);
  for (DeuteronModel m : deuteron_models()) {
    CHECK(parse_deuteron_model(deuteron_model_name(m)) == m);
    std::string upper = deuteron_model_name(m);
    for (char &c : upper) c = char(std::toupper(c));
    CHECK(parse_deuteron_model(upper) == m);
  }
  CHECK(parse_deuteron_model("Paris") == DeuteronModel::Paris);
  CHECK_THROWS_WITH(parse_deuteron_model("Nijmegen"), ContainsSubstring("unknown wave function \"Nijmegen\"") &&
                                                         ContainsSubstring("known: cdbonn, cdbonn87, paris, av18, av18b, nijm3, "
                                                                           "grossiib, grossl4, grosswjc1, grosswjc2, test"));
}

TEST_CASE("TDeuteron factory: the name and the enumeration give the same wave function", "[models]") {
  for (DeuteronModel m : deuteron_models()) {
    std::unique_ptr<TDeuteron::Wavefunction> by_name(TDeuteron::Wavefunction::CreateWavefunction(deuteron_model_name(m)));
    std::unique_ptr<TDeuteron::Wavefunction> by_model(TDeuteron::Wavefunction::CreateWavefunction(m));
    REQUIRE(by_name);
    REQUIRE(by_model);
    CHECK(std::string(by_name->GetName()) == by_model->GetName());
    CHECK(by_name->GetUp(0.2) == by_model->GetUp(0.2));
  }
}

TEST_CASE("OldDeuteron: unknown names throw (until 2026-11 they gave an empty parametrization in a release build)", "[models]") {
  CHECK_NOTHROW(OldDeuteron("Paris"));
  CHECK_NOTHROW(OldDeuteron(DeuteronModel::AV18));
  CHECK_THROWS_AS(OldDeuteron("nonsense"), std::invalid_argument);
  CHECK_THROWS_WITH(OldDeuteron(DeuteronModel::GrossIIb), ContainsSubstring("grossiib") && ContainsSubstring("available: paris"));
}

TEST_CASE("He3Model: names", "[models]") {
  CHECK(he3_models().size() == 5);
  for (He3Model m : he3_models()) CHECK(parse_he3_model(he3_model_name(m)) == m);
  CHECK(parse_he3_model("AV18-UIX") == He3Model::AV18_UIX);
  CHECK_THROWS_AS(parse_he3_model("av18"), std::invalid_argument);  // case sensitive, as before
  CHECK_THROWS_WITH(parse_he3_model("nonsense"), ContainsSubstring("invalid name for the he3 wave function: nonsense") &&
                                                    ContainsSubstring("known: AV18, AV18-UIX"));
}

TEST_CASE("NucleonStructure: a misspelled plain name fails at construction", "[models]") {
  // until 2026-11 this only threw when a structure function was first asked for
  CHECK_THROWS_WITH(NucleonStructure("SLACK"), ContainsSubstring("invalid name for structure function \"SLACK\"") &&
                                                  ContainsSubstring("known: SLAC, Alekhin, CB, CTEQ, HMRS, MSTW"));
  CHECK_THROWS_AS(NucleonStructure("cb"), std::invalid_argument);  // case sensitive, as before
  CHECK_THROWS_AS(NucleonStructure("CB+SLACK"), std::invalid_argument);
  CHECK_NOTHROW(NucleonStructure("CB+SLAC@9"));
  for (StructureFunctionModel m : {StructureFunctionModel::SLAC, StructureFunctionModel::Alekhin, StructureFunctionModel::CB,
                                   StructureFunctionModel::CTEQ, StructureFunctionModel::HMRS})
    CHECK(parse_structure_function_model(structure_function_model_name(m)) == m);
  NucleonStructure slac(StructureFunctionModel::SLAC), by_name("SLAC");
  CHECK(slac.getName() == "SLAC");
  CHECK(slac.getF2_xQ(true, 0.3, 4.e6) == by_name.getF2_xQ(true, 0.3, 4.e6));
  CHECK(slac.isF2Fit(1.e6));
  CHECK_FALSE(NucleonStructure(StructureFunctionModel::CTEQ).isF2Fit(1.e6));
  NucleonStructure copy(slac);
  CHECK(copy.getF2_xQ(true, 0.3, 4.e6) == slac.getF2_xQ(true, 0.3, 4.e6));  // the model is copied
}
