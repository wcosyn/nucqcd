// Cuba integrators through the numint wrappers (Integrator/numint/numint2Cuba.hpp).
// The main() of the tests sets CUBACORES=0 (serial), so runs with a fixed seed are reproducible.
#include <test_support.hpp>

#include <numint/numint.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
const double pi = std::acos(-1.);
// smooth 4D Gaussian on [-1,1]^4: exact (sqrt(pi) erf(1))^4
void gauss4(const numint::array<double, 4> &x, void *, double &ret) {
  ret = std::exp(-(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]));
}
const double gauss4_exact = std::pow(std::sqrt(pi) * std::erf(1.), 4);
// indicator of the 4D unit ball on [-1,1]^4: exact pi^2/2
void ball4(const numint::array<double, 4> &x, void *, double &ret) {
  ret = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] <= 1.) ? 1. : 0.;
}
const numint::array<double, 4> lower = {{-1., -1., -1., -1.}};
const numint::array<double, 4> upper = {{1., 1., 1., 1.}};

double run_vegas(numint::mdfunction<double, 4> &f, double &err) {
  int neval = 0, fail = 0;
  double result = 0., prob = 0.;
  numint::vegas(f, lower, upper, 1, 1e-3, 1e-12, 0, 9875661, 10000, 400000, 1000, 10000, 20000, 0, nullptr,
                neval, fail, result, err, prob);
  return result;
}
}  // namespace

TEST_CASE("Cuhre on a smooth 4D Gaussian", "[cuba]") {
  numint::mdfunction<double, 4> f;
  f.func = &gauss4;
  int nregions = 0, neval = 0, fail = 0;
  double result = 0., err = 0., prob = 0.;
  numint::cuhre(f, lower, upper, 1, 1e-8, 1e-12, 0, 1000, 1000000, 9, nullptr, nregions, neval, fail, result, err,
                prob);
  CHECK(fail == 0);
  CHECK_THAT(result, WithinRel(gauss4_exact, 1e-7));
}

TEST_CASE("Vegas on a smooth 4D Gaussian, reproducible with fixed seed", "[cuba]") {
  numint::mdfunction<double, 4> f;
  f.func = &gauss4;
  double err1 = 0., err2 = 0.;
  const double r1 = run_vegas(f, err1);
  const double r2 = run_vegas(f, err2);
  CHECK_THAT(r1, WithinAbs(gauss4_exact, 5. * err1));
  CHECK_THAT(r1, WithinRel(gauss4_exact, 5e-3));
  CHECK(r1 == r2);  // same seed, serial: bit-identical
}

TEST_CASE("Suave on the 4D unit ball (discontinuous integrand)", "[cuba]") {
  numint::mdfunction<double, 4> f;
  f.func = &ball4;
  int nregions = 0, neval = 0, fail = 0;
  double result = 0., err = 0., prob = 0.;
  numint::suave(f, lower, upper, 1, 1e-3, 1e-12, 0, 6667, 10000, 1000000, 5000, 2, 25., nullptr, nregions, neval,
                fail, result, err, prob);
  CHECK_THAT(result, WithinRel(pi * pi / 2., 1e-2));
}

TEST_CASE("Divonne on the 4D unit ball (discontinuous integrand)", "[cuba]") {
  numint::mdfunction<double, 4> f;
  f.func = &ball4;
  int nregions = 0, neval = 0, fail = 0;
  double result = 0., err = 0., prob = 0.;
  numint::divonne(f, lower, upper, 1, 1e-3, 1e-12, 0, 1235488, 10000, 1000000, 47, 1, 1, 5, 0., 10., 0.25, 0, 4,
                  nullptr, 0, nullptr, nullptr, nregions, neval, fail, result, err, prob);
  CHECK_THAT(result, WithinRel(pi * pi / 2., 1e-2));
}

// Hidden test (not run by ctest; run by hand with: unit_cuba "[known-bug]"): Suave with a small
// flatness parameter (0.25, as used in progs/cuba_test/*.cpp; 0.2 in progs/2bodymom/src/SCX_test.cpp)
// overflows its fluctuation estimate (minfluct = inf in src/suave/Integrate.c), computes a negative
// sample count and writes out of bounds: "buffer overflow detected" / segfault. Cuba's demo uses 25.
TEST_CASE("Known bug: Suave crashes for small flatness", "[.][cuba][known-bug]") {
  numint::mdfunction<double, 4> f;
  f.func = &gauss4;
  int nregions = 0, neval = 0, fail = 0;
  double result = 0., err = 0., prob = 0.;
  numint::suave(f, lower, upper, 1, 1e-3, 1e-12, 0, 6667, 10000, 1000000, 5000, 2, 0.25, nullptr, nregions, neval,
                fail, result, err, prob);
  CHECK_THAT(result, WithinRel(gauss4_exact, 1e-2));
}
