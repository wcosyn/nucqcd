// MePhys: Tools/MathUtil.hpp, the replacements of ROOT's MathMore polynomial roots and of TDecompLU
// (roadmap step 9, phase 2). The reference values (hex floats) were computed with ROOT 6.40:
// ROOT::Math::Polynomial(a,b,c).FindRealRoots() and TDecompLU (tolerance 1e-16) followed by Solve.
// Before the switch both replacements were also compared with ROOT on millions of random inputs (no
// difference), and the deuteron wave functions and TKinematics2to3(WithLabAngles) were bitwise the same
// (solve_3x3 gives the TYukawaPWF D-wave coefficients; their values are pinned by reg_deuteron_wf).
#include <test_support.hpp>

#include <MathUtil.hpp>

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("quadratic_real_roots: simple cases", "[mathutil]") {
  CHECK(nucqcd::quadratic_real_roots(1., -3., 2.) == std::vector<double>{1., 2.});
  CHECK(nucqcd::quadratic_real_roots(-1., 3., -2.) == std::vector<double>{1., 2.});  // ascending order
  CHECK(nucqcd::quadratic_real_roots(1., 2., 1.) == std::vector<double>{-1., -1.});  // double root twice
  CHECK(nucqcd::quadratic_real_roots(1., 0., 1.).empty());                          // complex roots
  CHECK(nucqcd::quadratic_real_roots(0., 2., -1.) == std::vector<double>{0.5});      // linear
  CHECK_THROWS_AS(nucqcd::quadratic_real_roots(0., 0., 1.), std::invalid_argument);
}

TEST_CASE("quadratic_real_roots: bitwise equal to ROOT::Math::Polynomial::FindRealRoots", "[mathutil]") {
  struct Case { std::array<double, 3> abc; std::vector<double> roots; };
  const Case cases[] = {
      // b = 0: gsl_poly_solve_quadratic would give a different last bit
      {{0x1.3eb55ef1fddecp+2, 0x0p+0, -0x1.f0299d883ba34p+1}, {-0x1.c3b88a6d9aad3p-1, 0x1.c3b88a6d9aad3p-1}},
      {{0x1.3333333333333p-2, -0x1.c666666666666p+2, 0x1.7333333333333p+1}, {0x1.a9bb7e65dc3ebp-2, 0x1.7403bcb11339bp+4}},
      {{-0x1.4p+1, 0x1.4p+0, 0x1.bp+2}, {-0x1.697de111f5b1ep+0, 0x1.e97de111f5b1ep+0}},
      {{0x1.0624dd2f1a9fcp-10, 0x1.4p+2, -0x1p+1}, {-0x1.38866644d9d6fp+12, 0x1.99913675babd7p-2}},
      {{0x1.b333333333333p+0, -0x1.2666666666666p+1, 0x0p+0}, {0x0p+0, 0x1.5a5a5a5a5a5a5p+0}},
  };
  for (const auto &c : cases) CHECK(nucqcd::quadratic_real_roots(c.abc[0], c.abc[1], c.abc[2]) == c.roots);
}

TEST_CASE("solve_3x3: bitwise equal to ROOT's TDecompLU", "[mathutil]") {
  // matrices of the TYukawaPWF boundary conditions {1/m_i^2, 1, m_i^2} for three masses m_i
  struct Case { std::array<double, 9> a; std::array<double, 3> b, x; };
  const Case cases[] = {
      {{0x1.b609b0ca8719bp-8, 0x1.03c1695dbe564p-7, 0x1.390306b4eb91bp-7, 0x1p+0, 0x1p+0, 0x1p+0, 0x1.2b39ac075a0a4p+7,
        0x1.f898fa0dcce3dp+6, 0x1.a2be9c0ce5b33p+6},
       {-0x1.930be0ded288dp-7, 0x1.d2f1a9fbe76c9p-2, -0x1.f8f5c28f5c28fp+2},
       {-0x1.38e9ea45295f3p+5, 0x1.3d203911772afp+6, -0x1.3db0a489cd27cp+5}},
      {{0x1.3c0ca4587e6b7p+0, 0x1.3c0ca4587e6b7p-2, 0x1.18eecaf953edbp-3, 0x1p+0, 0x1p+0, 0x1p+0, 0x1.9eb851eb851ecp-1,
        0x1.9eb851eb851ecp+1, 0x1.d28f5c28f5c2ap+2},
       {0x1p+0, -0x1p+1, 0x1.8p+1},
       {0x1.39f093a8847d5p+1, -0x1.087ac763a4e09p+3, 0x1.e7fa89e60f04fp+1}},
      {{0x1.c71c71c71c71cp-4, 0x1p+2, 0x1.c71c71c71c71cp-2, 0x1p+0, 0x1p+0, 0x1p+0, 0x1.2p+3, 0x1p-2, 0x1.2p+1},
       {0x0p+0, 0x0p+0, 0x1p+0},
       {0x1.3813813813814p-3, 0x1.d41d41d41d41bp-7, -0x1.5555555555555p-3}},
  };
  for (const auto &c : cases) {
    const auto x = nucqcd::solve_3x3(c.a, c.b);
    CHECK(x == c.x);
    for (int i = 0; i < 3; i++)  // and it is a solution
      CHECK_THAT(c.a[3 * i] * x[0] + c.a[3 * i + 1] * x[1] + c.a[3 * i + 2] * x[2] - c.b[i], WithinAbs(0., 1e-12));
  }
  CHECK_THROWS_AS(nucqcd::solve_3x3({1., 2., 3., 2., 4., 6., 0., 1., 1.}, {1., 1., 1.}), std::runtime_error);
}
