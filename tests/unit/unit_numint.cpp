// Integrator library: numint quadratures and adaptive cubature against exact integrals.
#include <test_support.hpp>

#include <adaptive/cubature.h>
#include <numint/numint.hpp>

#include <cmath>
#include <vector>

using Catch::Matchers::WithinRel;

namespace {
const double pi = std::acos(-1.);

void sin_1d(double x, void *, double &ret) { ret = std::sin(x); }
void exp_1d(double x, void *, double &ret) { ret = std::exp(x); }
void spherical_volume(const numint::array<double, 3> &x, void *, double &ret) {
  ret = x[0] * x[0] * std::sin(x[1]);  // r^2 sin(theta) on [0,1]x[0,pi]x[0,2pi]
}
// integrand of progs/integrator/inttest.cpp: r^2 cos^2(phi) + sqrt(phi) on [0,1]x[-1,1]x[0,2pi]
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
