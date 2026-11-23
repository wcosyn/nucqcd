// Tools/Core/recmosh: Moshinsky brackets <(n l, N Lambda) L | (n1 l1, n2 l2) L> (equal masses), computed
// recursively. They form an orthogonal transformation within each oscillator shell E and total L.
#include <test_support.hpp>

#include <recmosh.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <vector>

using Catch::Matchers::WithinAbs;

TEST_CASE("Moshinsky brackets: closed-form values", "[moshinsky]") {
  RecMosh m0000(0, 0, 0, 0, "");
  CHECK_THAT(m0000.getCoefficient(0, 0, 0, 0, 0), WithinAbs(1., 1e-14));

  // two 0p particles (values and phase convention as in the former share/recmosh/recmosh0101.dat)
  RecMosh m0101(0, 1, 0, 1, "");
  const double s = 1. / std::sqrt(2.);
  CHECK_THAT(m0101.getCoefficient(0, 0, 1, 0, 0), WithinAbs(s, 1e-14));   // L = 0
  CHECK_THAT(m0101.getCoefficient(1, 0, 0, 0, 0), WithinAbs(-s, 1e-14));
  CHECK_THAT(m0101.getCoefficient(0, 1, 0, 1, 1), WithinAbs(1., 1e-14));  // L = 1
  CHECK_THAT(m0101.getCoefficient(0, 0, 0, 2, 2), WithinAbs(s, 1e-14));   // L = 2
  CHECK_THAT(m0101.getCoefficient(0, 2, 0, 0, 2), WithinAbs(-s, 1e-14));
  // energy is conserved: brackets between different shells vanish
  CHECK(m0101.getCoefficient(0, 0, 0, 0, 0) == 0.);
}

// Regression for the former share/recmosh/ cache, whose integer key 100000*n + 10000*l + ...
// collided once l or Lambda reached 10 (shells E >= 10): checking E up to 12 covers that.
TEST_CASE("Moshinsky brackets: orthonormal within each shell", "[moshinsky]") {
  const int Emax = 12;
  double worst = 0.;
  for (int E = 0; E <= Emax; ++E) {
    for (int L = 0; L <= E; ++L) {
      std::vector<std::array<int, 4>> states;  // (n1 l1 n2 l2) with 2n1+l1+2n2+l2 = E
      for (int n1 = 0; 2 * n1 <= E; ++n1)
        for (int l1 = 0; 2 * n1 + l1 <= E; ++l1)
          for (int n2 = 0; 2 * n1 + l1 + 2 * n2 <= E; ++n2) {
            int l2 = E - 2 * n1 - l1 - 2 * n2;
            if (l2 >= std::abs(l1 - L) && l2 <= l1 + L) states.push_back({n1, l1, n2, l2});
          }
      std::vector<RecMosh> mosh;
      for (auto& s : states) mosh.emplace_back(s[0], s[1], s[2], s[3], "");
      for (size_t a = 0; a < states.size(); ++a)
        for (size_t b = a; b < states.size(); ++b) {
          double sum = 0.;
          for (int n = 0; 2 * n <= E; ++n)
            for (int l = 0; 2 * n + l <= E; ++l)
              for (int N = 0; 2 * n + l + 2 * N <= E; ++N) {
                int Lambda = E - 2 * n - l - 2 * N;
                if (Lambda < std::abs(l - L) || Lambda > l + L) continue;
                sum += mosh[a].getCoefficient(n, l, N, Lambda, L) * mosh[b].getCoefficient(n, l, N, Lambda, L);
              }
          worst = std::max(worst, std::fabs(sum - (a == b ? 1. : 0.)));
        }
    }
  }
  CHECK(worst < 1e-12);
}
