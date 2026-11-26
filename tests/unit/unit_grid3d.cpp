// Tools/Core/Grid3D.hpp: the owning 3D array that replaced the `T ***` arrays of the Glauber grids.
#include <test_support.hpp>

#include <Tools/Core/Grid3D.hpp>
#include <Tools/Core/Utilfunctions.hpp>

#include <complex>
#include <vector>
using namespace nucqcd;

using Catch::Matchers::WithinAbs;

TEST_CASE("Grid3D is zero-initialized and indexed like a C array", "[grid3d]") {
  Grid3D<std::complex<double> > grid(3, 4, 5);
  CHECK(grid.size0() == 3);
  CHECK(grid.size1() == 4);
  CHECK(grid.size2() == 5);
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 4; j++)
      for (std::size_t k = 0; k < 5; k++) CHECK(grid[i][j][k] == std::complex<double>(0., 0.));
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 4; j++)
      for (std::size_t k = 0; k < 5; k++) grid[i][j][k] = std::complex<double>(100. * i + 10. * j + k, -1.);
  // the last index runs fastest, so the iterators give the order of the loops i, j, k (the order of the cache files)
  std::size_t n = 0;
  for (std::size_t i = 0; i < 3; i++)
    for (std::size_t j = 0; j < 4; j++)
      for (std::size_t k = 0; k < 5; k++, n++) CHECK(*(grid.begin() + n) == std::complex<double>(100. * i + 10. * j + k, -1.));
  CHECK(grid.end() - grid.begin() == 60);
  // grid[i][j] is a pointer to the row of fixed i and j
  CHECK(&grid[1][2][0] + 3 == &grid[1][2][3]);
  const Grid3D<std::complex<double> > &constgrid = grid;
  CHECK(constgrid[2][3][4] == std::complex<double>(234., -1.));
}

TEST_CASE("Grid3D can be copied and moved, and the copy is independent", "[grid3d]") {
  Grid3D<double> grid(2, 2, 2);
  grid[1][1][1] = 7.;
  Grid3D<double> copy = grid;
  copy[1][1][1] = 9.;
  CHECK(grid[1][1][1] == 7.);
  CHECK(copy[1][1][1] == 9.);
  Grid3D<double> moved = std::move(copy);
  CHECK(moved[1][1][1] == 9.);
  // the arrays of the Glauber grids: [variant][nucleon] of Grid3D
  std::vector<std::vector<Grid3D<double> > > sets(2, std::vector<Grid3D<double> >(2, Grid3D<double>(2, 3, 4)));
  sets[1][0][1][2][3] = 5.;
  CHECK(sets[1][0][1][2][3] == 5.);
  CHECK(sets[0][0][1][2][3] == 0.);  // not shared with the others
  CHECK(sets[1][1][1][2][3] == 0.);
}

TEST_CASE("Interp3d on a Grid3D is the trilinear interpolation", "[grid3d]") {
  // a function that is linear in each variable is reproduced exactly: f = 1 + 2 i + 3 j - k + 0.5 i j k (complex: i f)
  Grid3D<std::complex<double> > grid(4, 4, 4);
  auto f = [](double i, double j, double k) { return 1. + 2. * i + 3. * j - k + 0.5 * i * j * k; };
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      for (int k = 0; k < 4; k++) grid[i][j][k] = std::complex<double>(f(i, j, k), 2. * f(i, j, k));
  const double s = 0.25, t = 0.5, u = 0.75;
  const std::complex<double> value = Interp3d(grid, s, t, u, 1. - s, 1. - t, 1. - u, 1, 2, 0);
  CHECK_THAT(value.real(), WithinAbs(f(1. + s, 2. + t, 0. + u), 1e-12));
  CHECK_THAT(value.imag(), WithinAbs(2. * f(1. + s, 2. + t, 0. + u), 1e-12));
}

TEST_CASE("LUsolve for std::vector and for pointers solve the same linear system", "[grid3d][lusolve]") {
  // A x = b with x = (1, -2, 3)
  const double A[3][3] = {{4., -1., 2.}, {-1., 5., 1.}, {2., 1., 6.}};
  std::vector<std::vector<double> > a(3, std::vector<double>(3));
  std::vector<double> b(3, 0.);
  double *rows[3];
  double raw[3][3];
  double braw[3] = {0., 0., 0.};
  const double x[3] = {1., -2., 3.};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a[i][j] = raw[i][j] = A[i][j];
      b[i] += A[i][j] * x[j];
      braw[i] += A[i][j] * x[j];
    }
    rows[i] = raw[i];
  }
  LUsolve(a, b);
  LUsolve(rows, braw, 3);
  for (int i = 0; i < 3; i++) {
    CHECK_THAT(b[i], WithinAbs(x[i], 1e-12));
    CHECK(b[i] == braw[i]);  // the vector form calls the pointer form: the same numbers
  }
}
