// Errors in the libraries are exceptions (Tools/Core/ErrorMessage.hpp); until 2026-11 they printed a
// message and called exit(1), which ended the program and could not be tested. One or two cases of
// each kind: std::invalid_argument, std::out_of_range, std::logic_error, std::runtime_error.
#include <test_support.hpp>

#include <catch2/matchers/catch_matchers_string.hpp>

#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/IntegratorFlag.hpp>
#include <Glauber/GlauberDecayGridThick.hpp>
#include <Glauber/OneGlauberGrid.hpp>
#include <Glauber/FastParticle.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Tools/Core/Matrix.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Tools/StructureFunctions/NucleonStructure.hpp>
#include <Tools/Kinematics/Kinematics2to3.hpp>
#include <Tools/Core/constants.hpp>

#include <stdexcept>
#include <string>
using namespace nucqcd;

using Catch::Matchers::ContainsSubstring;

TEST_CASE("ErrorMessage builds the message with stream syntax", "[errors]") {
  const std::string msg = nucqcd::ErrorMessage() << "r out of range: " << 12.5 << " fm";
  CHECK(msg == "r out of range: 12.5 fm");
  CHECK_THROWS_WITH(throw std::out_of_range(nucqcd::ErrorMessage() << "x = " << 3), "x = 3");
}

TEST_CASE("invalid options throw std::invalid_argument", "[errors]") {
  const std::string dir = nucqcd_test::share_dir();
  CHECK_THROWS_AS(MeanFieldNucleus(42, dir), std::invalid_argument);
  double F1 = 0., F2 = 0.;
  CHECK_THROWS_AS(NucleonStructure("nonsense").getF_xQ(F1, F2, true, 0.3, 2.e06), std::invalid_argument);
  CHECK_THROWS_AS(FastParticle(77, 0, 1000., 0., 0., 1., 0., 1., 1., dir), std::invalid_argument);
  CHECK_THROWS_AS(Kinematics2to3(1, Kinematics2to3::kK, "qsquared:wlab:pk", 1.e06, 600., 300., 0.5, 0.5, 0.5),
                  std::invalid_argument);
}

TEST_CASE("values outside a grid or parametrization throw std::out_of_range", "[errors]") {
  const std::string dir = nucqcd_test::share_dir();
  const MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  CHECK_NOTHROW(carbon.getWave_F(0, 0.5 * carbon.getRange()));
  CHECK_THROWS_AS(carbon.getWave_F(0, carbon.getRange() + 1.), std::out_of_range);
  CHECK_THROWS_AS(carbon.getProtonDensity(carbon.getRange() + 1.), std::out_of_range);
  // Christy-Bosted is valid up to W^2 = 25 GeV^2
  double F1 = 0., F2 = 0.;
  const NucleonStructure cb("CB");
  CHECK_NOTHROW(cb.getF_xQ(F1, F2, true, 0.5, 2.e06));
  CHECK_THROWS_WITH(cb.getF_xQ(F1, F2, true, 0.01, 2.e06), ContainsSubstring("Christy & Bosted"));
  CHECK_THROWS_AS(cb.getF_xQ(F1, F2, true, 0.01, 2.e06), std::out_of_range);
  Matrix<2, 2> m;
  CHECK_THROWS_AS(m(2, 0), std::out_of_range);
}

TEST_CASE("calls in the wrong order throw std::logic_error", "[errors]") {
  const std::string dir = nucqcd_test::share_dir();
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
  GlauberGridThick grid(8, 4, 2, &carbon, 1.E-03, IntegrationMethod::Adaptive, dir);
  CHECK_THROWS_AS(grid.getFsiGridFull_interp3(1., 0.5, 0.), std::logic_error);  // not filled yet
}

TEST_CASE("missing input throws std::runtime_error", "[errors]") {
  CHECK_THROWS_AS(MeanFieldNucleus(MeanFieldNucleus::C, "/nonexistent/share/dir"), std::runtime_error);
}

TEST_CASE("the integrator option is 1 or 2; 0 (the old rombergerN) was removed", "[errors][integrator]") {
  // the grids take an IntegrationMethod since 2026-11, so an invalid integer can no longer reach a constructor: it is
  // rejected when it is converted (before 2026-11 a constructor argument 0 silently gave grids without FSI)
  CHECK(integration_method(1) == IntegrationMethod::Romberg);
  CHECK(integration_method(2) == IntegrationMethod::Adaptive);
  CHECK(numint::method_flag(IntegrationMethod::Romberg) == 1);
  CHECK(numint::method_flag(IntegrationMethod::Adaptive) == 2);
  CHECK_THROWS_WITH(integration_method(0), ContainsSubstring("rombergerN"));
  CHECK_THROWS_AS(integration_method(3), std::invalid_argument);
  CHECK(integration_method("romberg") == IntegrationMethod::Romberg);
  CHECK(integration_method("Adaptive") == IntegrationMethod::Adaptive);
  CHECK(integration_method("2") == IntegrationMethod::Adaptive);
  CHECK_THROWS_AS(integration_method("simpson"), std::invalid_argument);
  CHECK_THROWS_AS(integration_method(""), std::invalid_argument);
  // the Knockout classes still take an integer and validate it with checked_integrator_flag
  CHECK(checked_integrator_flag(1) == 1);
  CHECK(checked_integrator_flag(2) == 2);
  CHECK_THROWS_AS(checked_integrator_flag(0), std::invalid_argument);
}
