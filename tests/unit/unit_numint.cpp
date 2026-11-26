// Integrator library: numint quadratures and adaptive cubature against exact integrals.
#include <test_support.hpp>

#include <Integrator/adaptive/cubature.h>
#include <Integrator/numint/cube.hpp>
#include <Integrator/numint/numint.hpp>

#include <cmath>
#include <sstream>
#include <vector>

using Catch::Matchers::WithinRel;

namespace {
const double pi = std::acos(-1.);

void sin_1d(double x, void *, double &ret) { ret = std::sin(x); }
void exp_1d(double x, void *, double &ret) { ret = std::exp(x); }
void spherical_volume(const numint::array<double, 3> &x, void *, double &ret) {
  ret = x[0] * x[0] * std::sin(x[1]);  // r^2 sin(theta) on [0,1]x[0,pi]x[0,2pi]
}
// integrand of progs/attic/integrator/inttest.cpp: r^2 cos^2(phi) + sqrt(phi) on [0,1]x[-1,1]x[0,2pi]
void cubature_integrand(unsigned, const double *x, void *, unsigned, double *ret) {
  ret[0] = x[0] * x[0] * std::cos(x[2]) * std::cos(x[2]) + std::sqrt(x[2]);
}
void cubature_integrand_md(const numint::array<double, 3> &x, void *, numint::vector_d &ret) {
  ret = numint::vector_d(1, x[0] * x[0] * std::cos(x[2]) * std::cos(x[2]) + std::sqrt(x[2]));
}
const double cubature_exact = 2. * pi / 3. + 4. / 3. * std::pow(2. * pi, 1.5);
}  // namespace

TEST_CASE("quad_romb integrates sin over [0,pi]", "[numint]") {
  numint::function<double> f;
  f.func = &sin_1d;
  double result = 0.;
  unsigned count = 0;
  numint::quad_romb(f, 0., pi, 1e-12, 1e-10, result, count);
  CHECK_THAT(result, WithinRel(2., 1e-9));
  CHECK(count > 0);
}

TEST_CASE("quad_cc integrates exp over [0,1]", "[numint]") {
  numint::function<double> f;
  f.func = &exp_1d;
  double result = 0.;
  unsigned count = 0;
  numint::quad_cc(f, 0., 1., 1e-12, 1e-10, result, count);
  CHECK_THAT(result, WithinRel(std::exp(1.) - 1., 1e-9));
}

TEST_CASE("cube_romb gives the volume of the unit sphere", "[numint]") {
  numint::mdfunction<double, 3> f;
  f.func = &spherical_volume;
  numint::array<double, 3> lower = {{0., 0., 0.}};
  numint::array<double, 3> upper = {{1., pi, 2. * pi}};
  double result = 0.;
  unsigned count = 0;
  numint::cube_romb(f, lower, upper, 1e-8, 1e-7, result, count);
  CHECK_THAT(result, WithinRel(4. * pi / 3., 1e-6));
}

TEST_CASE("adapt_integrate (adaptive cubature) on the inttest integrand", "[cubature]") {
  double xmin[3] = {0., -1., 0.}, xmax[3] = {1., 1., 2. * pi};
  double res[1], err[1];
  const int error = adapt_integrate(1, cubature_integrand, nullptr, 3, xmin, xmax, 100, 0, 0., 1e-8, res, err);
  CHECK(error == 0);
  CHECK_THAT(res[0], WithinRel(cubature_exact, 1e-7));
}

TEST_CASE("numint::cube_adaptive wrapper around the adaptive cubature", "[cubature]") {
  numint::mdfunction<numint::vector_d, 3> f;
  f.func = &cubature_integrand_md;
  numint::array<double, 3> lower = {{0., -1., 0.}};
  numint::array<double, 3> upper = {{1., 1., 2. * pi}};
  std::vector<double> result(1, 0.);
  unsigned count = 0;
  numint::cube_adaptive(f, lower, upper, 1e-12, 1e-8, 100, 2000000, result, count, 0);
  CHECK_THAT(result[0], WithinRel(cubature_exact, 1e-7));
}

TEST_CASE("numint::cube_adaptive reports integrals that stop at maxEval", "[cubature][adaptive_stats]") {
  numint::mdfunction<numint::vector_d, 3> f;
  f.func = &cubature_integrand_md;
  numint::array<double, 3> lower = {{0., -1., 0.}};
  numint::array<double, 3> upper = {{1., 1., 2. * pi}};
  auto &stats = numint::adaptive_stats();
  stats.reset();
  stats.warn = false;  // no stderr noise in the test log
  std::vector<double> result(1, 0.);
  unsigned count = 0;

  // converges: counted, not capped
  numint::cube_adaptive(f, lower, upper, 1e-12, 1e-6, 100, 2000000, result, count, 0);
  CHECK(stats.calls == 1);
  CHECK(stats.capped == 0);

  // 1000 evaluations cannot give 1e-10: capped, with a large estimated error; return value stays 0
  int code = numint::cube_adaptive(f, lower, upper, 1e-14, 1e-10, 100, 1000, result, count, 0);
  CHECK(code == 0);
  CHECK(stats.calls == 2);
  CHECK(stats.capped == 1);
  CHECK(stats.worst_rel_error > 1e-10);
  CHECK(stats.max_eval_seen == 1000);

  std::ostringstream out;
  CHECK(numint::adaptive_report(out));
  CHECK(out.str().find("# cube_adaptive: 2 integrals, 1 stopped at maxEval") == 0);
  stats.reset();
  stats.warn = true;
  std::ostringstream none;
  CHECK_FALSE(numint::adaptive_report(none));
}

TEST_CASE("numint::cube calls cube_romb and cube_adaptive with the arguments of its options", "[numint][cubature]") {
  // double and vector results, as used by the Glauber grids: the dispatch must give bit-identical results and counts
  numint::mdfunction<double, 3> f;
  f.func = &spherical_volume;
  numint::array<double, 3> lower = {{0., 0., 0.}};
  numint::array<double, 3> upper = {{1., pi, 2. * pi}};
  {
    double direct = 0., viacube = 0.;
    unsigned count_direct = 0, count_cube = 0;
    numint::cube_romb(f, lower, upper, 1e-8, 1e-7, direct, count_direct);
    numint::cube({numint::IntegrationMethod::Romberg, 1e-8, 1e-7, 0, 0}, f, lower, upper, viacube, count_cube);
    CHECK(viacube == direct);
    CHECK(count_cube == count_direct);
    CHECK_THAT(viacube, WithinRel(4. * pi / 3., 1e-6));
  }
  {
    double direct = 0., viacube = 0.;
    unsigned count_direct = 0, count_cube = 0;
    numint::cube_adaptive(f, lower, upper, 1e-8, 1e-7, 100, 200000, direct, count_direct, 0);
    numint::cube({numint::IntegrationMethod::Adaptive, 1e-8, 1e-7, 100, 200000}, f, lower, upper, viacube, count_cube);
    CHECK(viacube == direct);
    CHECK(count_cube == count_direct);
  }
  numint::mdfunction<numint::vector_d, 3> g;
  g.func = &cubature_integrand_md;
  numint::array<double, 3> lower2 = {{0., -1., 0.}};
  numint::array<double, 3> upper2 = {{1., 1., 2. * pi}};
  for (const auto method : {numint::IntegrationMethod::Romberg, numint::IntegrationMethod::Adaptive}) {
    std::vector<double> direct(1, 0.), viacube(1, 0.);
    unsigned count_direct = 0, count_cube = 0;
    if (method == numint::IntegrationMethod::Romberg)
      numint::cube_romb(g, lower2, upper2, 1e-12, 1e-4, direct, count_direct, 0);
    else
      numint::cube_adaptive(g, lower2, upper2, 1e-12, 1e-4, 100, 2000000, direct, count_direct, 0);
    numint::cube({method, 1e-12, 1e-4, 100, 2000000}, g, lower2, upper2, viacube, count_cube, 0);
    CHECK(viacube[0] == direct[0]);
    CHECK(count_cube == count_direct);
  }
}
