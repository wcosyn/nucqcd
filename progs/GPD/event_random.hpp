// Random numbers for the diDVCS event generators (diDVCS_gen, diDVCS_gen_grid, diDVCS_gen_gridonly).
// Replaces ROOT's TRandom3 and TF1::GetRandom (roadmap step 9, phase 6) with <random>:
//  - one std::mt19937_64 engine per program; seeded, as ROOT's SetSeed(0) was, differently for every run
//    (from std::random_device), unless the environment variable NUCQCD_SEED is set to a number, which
//    makes a run reproducible. The seed is written to stderr ("# random seed: ..."), so that stdout,
//    which holds the events, is unchanged;
//  - uniform(a, b) and breit_wigner(mean, gamma) with ROOT's formulas (TRandom::Uniform,
//    TRandom::BreitWigner), using a uniform number in (0, 1] as TRandom3::Rndm;
//  - sin2_angle() and cos2_angle(): angles in [0, 2pi) with densities proportional to sin^2 and cos^2,
//    sampled exactly by accept-reject (TF1::GetRandom used an interpolated inverse of a binned CDF).
// The event samples therefore differ from the ROOT ones event by event, but not in distribution.
#ifndef NUCQCD_PROGS_GPD_EVENT_RANDOM_HPP
#define NUCQCD_PROGS_GPD_EVENT_RANDOM_HPP

#include <Tools/Core/constants.hpp>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

namespace event_random {

inline std::mt19937_64 &engine() {
  static std::mt19937_64 gen = [] {
    std::uint64_t seed;
    if (const char *s = std::getenv("NUCQCD_SEED"))
      seed = std::stoull(s);
    else
      seed = (std::uint64_t(std::random_device{}()) << 32) ^ std::random_device{}();
    std::cerr << "# random seed: " << seed << " (set NUCQCD_SEED to repeat the run)" << std::endl;
    return std::mt19937_64(seed);
  }();
  return gen;
}

/// uniform in (0, 1], as TRandom3::Rndm
inline double rndm() { return 1.0 - std::uniform_real_distribution<double>(0.0, 1.0)(engine()); }

/// uniform in (a, b], as TRandom::Uniform(a, b)
inline double uniform(double a, double b) { return a + (b - a) * rndm(); }

/// Breit-Wigner (Cauchy) with peak mean and full width gamma, as TRandom::BreitWigner
inline double breit_wigner(double mean, double gamma) {
  double rval = 2 * rndm() - 1;
  return mean + 0.5 * gamma * std::tan(rval * nucqcd::pi / 2);
}

/// angle in [0, 2pi) with density proportional to sin^2
inline double sin2_angle() {
  for (;;) {
    double x = 2 * nucqcd::pi * (1.0 - rndm());
    double s = std::sin(x);
    if (rndm() <= s * s) return x;
  }
}

/// angle in [0, 2pi) with density proportional to cos^2
inline double cos2_angle() {
  for (;;) {
    double x = 2 * nucqcd::pi * (1.0 - rndm());
    double c = std::cos(x);
    if (rndm() <= c * c) return x;
  }
}

}  // namespace event_random

#endif
