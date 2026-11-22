// Ownership of resources (docs/ROADMAP.md, step 5): classes that own memory either copy it correctly
// (std::vector / rule of zero) or cannot be copied at all. Until 2026-11 several classes owned raw
// pointers without copy operations, so a copy would have deleted the memory twice, and
// MeanFieldNucleus[Thick] had a copy constructor but no assignment.
#include <test_support.hpp>

#include <ClassGridThick_SCX.hpp>
#include <DQEinclusive.hpp>
#include <Deut_Conv_EMTLF.hpp>
#include <FsiCorrelator.hpp>
#include <GlauberDecayGridThick.hpp>
#include <GlauberGrid.hpp>
#include <GlauberGridThick.hpp>
#include <InclusiveCross.hpp>
#include <InclusiveCrossRes.hpp>
#include <MeanFieldNucleusThick.hpp>
#include <Model.hpp>
#include <OneGlauberGrid.hpp>
#include <PionTCross.hpp>
#include <ROMEAGrid.hpp>
#include <Resonance.hpp>
#include <RhoTCross.hpp>

#include <complex>
#include <memory>
#include <type_traits>

// grids own raw arrays: not copyable (and neither is Model, which holds grids)
static_assert(!std::is_copy_constructible_v<GlauberGrid> && !std::is_copy_assignable_v<GlauberGrid>);
static_assert(!std::is_copy_constructible_v<GlauberGridThick> && !std::is_copy_assignable_v<GlauberGridThick>);
static_assert(!std::is_copy_constructible_v<GlauberDecayGridThick>);
static_assert(!std::is_copy_constructible_v<OneGlauberGrid>);
static_assert(!std::is_copy_constructible_v<ROMEAGrid>);
static_assert(!std::is_copy_constructible_v<Model>);
// their grids point to a member nucleus: neither copyable nor movable
static_assert(!std::is_copy_constructible_v<PionTCross> && !std::is_move_constructible_v<PionTCross>);
static_assert(!std::is_copy_constructible_v<RhoTCross> && !std::is_move_constructible_v<RhoTCross>);
// a deuteron wave function held by std::unique_ptr: movable, not copyable
static_assert(!std::is_copy_constructible_v<InclusiveCross> && std::is_move_constructible_v<InclusiveCross>);
static_assert(!std::is_copy_constructible_v<InclusiveCrossRes> && std::is_move_constructible_v<InclusiveCrossRes>);
static_assert(!std::is_copy_constructible_v<DQEinclusive> && std::is_move_constructible_v<DQEinclusive>);
static_assert(!std::is_copy_constructible_v<Deut_Conv_EMTLF> && std::is_move_constructible_v<Deut_Conv_EMTLF>);
// std::vector members: copyable and assignable (rule of zero)
static_assert(std::is_copy_constructible_v<MeanFieldNucleusThick> && std::is_copy_assignable_v<MeanFieldNucleusThick>);
static_assert(std::is_copy_constructible_v<FsiCorrelator> && std::is_copy_assignable_v<FsiCorrelator>);
static_assert(std::is_copy_constructible_v<ClassGridThick_SCX>);

TEST_CASE("MeanFieldNucleusThick: copies and assignments are independent and equal", "[ownership]") {
  const std::string dir = nucqcd_test::share_dir();
  auto copy = std::make_unique<MeanFieldNucleusThick>(MeanFieldNucleus::C, dir);
  MeanFieldNucleusThick assigned(MeanFieldNucleus::O, dir);
  {
    const MeanFieldNucleusThick carbon(MeanFieldNucleus::C, dir);
    *copy = MeanFieldNucleusThick(carbon);
    assigned = carbon;  // did not compile before 2026-11 (no operator=), now a deep copy
    REQUIRE(assigned.getA() == 12);
  }  // the original is gone: the copies must own their data
  const MeanFieldNucleusThick fresh(MeanFieldNucleus::C, dir);
  for (const MeanFieldNucleusThick *nucl : {copy.get(), &assigned}) {
    CHECK(nucl->getTotalLevels() == fresh.getTotalLevels());
    for (double r : {0.5, 1.7, 3.2}) {
      CHECK(nucl->getProtonDensity(r) == fresh.getProtonDensity(r));
      CHECK(nucl->getNeutronDensity(r) == fresh.getNeutronDensity(r));
      for (int shell = 0; shell < fresh.getTotalLevels(); shell++) {
        CHECK(nucl->getWave_F(shell, r) == fresh.getWave_F(shell, r));
        CHECK(nucl->getWave_G(shell, r) == fresh.getWave_G(shell, r));
      }
    }
    std::complex<double> wave[4], wave_fresh[4];
    nucl->getWaveFunction(wave, 1, 1, 0.9, 0.3, 0.7);
    fresh.getWaveFunction(wave_fresh, 1, 1, 0.9, 0.3, 0.7);
    for (int i = 0; i < 4; i++) CHECK(wave[i] == wave_fresh[i]);
  }
}
