// MePhys: Tools/Core/StringUtil.hpp, the std::string helpers that replace TString (roadmap step 9, phase 3).
#include <test_support.hpp>

#include <StringUtil.hpp>
#include <TDeuteron.h>

#include <memory>
#include <string>
#include <vector>

TEST_CASE("StringUtil: lower case, case-insensitive comparison and search, split without empty tokens", "[string]") {
  CHECK(nucqcd::to_lower("QSquared:WLab") == "qsquared:wlab");
  CHECK(nucqcd::iequals("CDBonn", "cdbonn"));
  CHECK_FALSE(nucqcd::iequals("cdbonn87", "cdbonn"));
  CHECK(nucqcd::icontains("qsquared:wlab:STOT", "stot"));
  CHECK_FALSE(nucqcd::icontains("qsquared:wlab", "stot"));
  // as TString::Tokenize: empty tokens are left out
  CHECK(nucqcd::split("a::b:", ':') == std::vector<std::string>{"a", "b"});
  CHECK(nucqcd::split(":a", ':') == std::vector<std::string>{"a"});
  CHECK(nucqcd::split("", ':').empty());
  CHECK(nucqcd::split("abc", ':') == std::vector<std::string>{"abc"});
}

TEST_CASE("TDeuteron wave-function names are case insensitive", "[string]") {
  std::unique_ptr<TDeuteron::Wavefunction> a(TDeuteron::Wavefunction::CreateWavefunction("AV18"));
  std::unique_ptr<TDeuteron::Wavefunction> b(TDeuteron::Wavefunction::CreateWavefunction("av18"));
  CHECK(a->GetUp(100.) == b->GetUp(100.));
  CHECK(std::string(a->GetName()) == "AV18");
}
