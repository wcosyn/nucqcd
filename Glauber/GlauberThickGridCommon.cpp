#include <Glauber/GlauberThickGridCommon.hpp>

#include <Tools/Core/ErrorMessage.hpp>
#include <Integrator/numint/cube.hpp>
#include <Tools/Core/constants.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

namespace nucqcd {

GlauberThickGridCommon::GlauberThickGridCommon(int nvariants, int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnuclthick,
                                               double prec, IntegrationMethod integrator, std::string dir)
  : AbstractFsiCTGrid(r_grid, cth_grid, phi_grid, pnuclthick, prec, integrator, dir),
    variants(nvariants), pnucleusthick(pnuclthick), fsicorrelator(pnuclthick, 90, 90, dir), error(0.) {
  number_of_grids = 2 * nvariants;
  // [variant][nucleon: 0 neutron, 1 proton] of an array of (r, cos theta, phi) points, set to zero
  const Grid3D<complex<double> > points(getRgrid() + 1, getCthgrid() + 1, getPhigrid() + 1);
  fsi_grid.assign(variants, std::vector<Grid3D<complex<double> > >(2, points));
  fsi_ct_grid.assign(variants, std::vector<Grid3D<complex<double> > >(2, points));
}

GlauberThickGridCommon::~GlauberThickGridCommon() {
}

int GlauberThickGridCommon::getGridNumber(GridVariant variant) const {
  const int bits = static_cast<int>(variant);
  const int srcdecay = bits & 3;  // bit 0 SRC, bit 1 decay
  if (srcdecay >= variants) return -1;
  return ((bits & 4) ? variants : 0) + srcdecay;  // bit 2: CT, the CT grids follow the others
}

void GlauberThickGridCommon::setFilenames(string dir) {
  AbstractFsiCTGrid::setFilenames(dir + "Thick.");
}

complex<double> GlauberThickGridCommon::interpolateVariant(const std::vector<std::vector<Grid3D<std::complex<double> > > > &array, int variant) {
  if (!filledallgrid) {
    throw std::logic_error(nucqcd::ErrorMessage() << "You have to fill the grids first!");
  }
  return pow(getInterp(array[variant][0]), getPnucleus()->getN() - getNeutronKnockout())
    * pow(getInterp(array[variant][1]), getPnucleus()->getZ() - getProtonKnockout());
}

complex<double> GlauberThickGridCommon::getFsiGridN_interp(int grid) {
  if (grid < 0 || grid >= 2 * variants) {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "grid number not supported");
  }
  return grid < variants ? interpolateVariant(fsi_grid, grid) : interpolateVariant(fsi_ct_grid, grid - variants);
}

void GlauberThickGridCommon::print_grid(int grid) {
  printKnockout();
  for (int i = 0; i <= getRgrid(); i++) {
    double r = float(i) * getPnucleus()->getRange() / getRgrid();
    for (int j = 0; j <= getCthgrid(); j++) {
      double costheta = -2. * j / getCthgrid() + 1.;
      for (int k = 0; k <= getPhigrid(); k++) {
        double phi = (getAllinplane() ? 1. : 2.) * PI * k / getPhigrid();
        if (std::isnan(phi)) phi = 0.;
        setinterp(r, costheta, phi);
        complex<double> value = getFsiGridN_interp(grid);
        cout << r << " " << costheta << " " << phi * nucqcd::rad_to_deg << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }
}

bool GlauberThickGridCommon::closeToOne(bool ctOnly, int i, int j, int k) const {
  const MeanFieldNucleus *nucleus = getPnucleus();
  for (int v = 0; v < variants; v++) {
    if (!(abs(pow(fsi_ct_grid[v][0][i][j][k], nucleus->getN()) * pow(fsi_ct_grid[v][1][i][j][k], nucleus->getZ())) > 0.999999)) return false;
    if (!ctOnly && !(abs(pow(fsi_grid[v][0][i][j][k], nucleus->getN()) * pow(fsi_grid[v][1][i][j][k], nucleus->getZ())) > 0.999999))
      return false;
  }
  return true;
}

void GlauberThickGridCommon::constructAllGrids() { constructGrids(false); }
void GlauberThickGridCommon::constructCtGrid() { constructGrids(true); }

void GlauberThickGridCommon::constructGrids(bool ctOnly) {
  // The loops over the grid points. A point where all factors are within 1e-6 of one (the "threshold") makes
  // the points at larger r with the same angles one too, which saves computing time; the point at r = 0
  // and the points along the z axis (cos theta = +-1) are copied to all angles.
  std::vector<std::vector<int> > treshold(getCthgrid() + 1, std::vector<int>(getPhigrid() + 1, 0));

  for (int i = 0; i <= getRgrid(); i++) {
    r_hit = float(i) * getPnucleus()->getRange() / getRgrid();
    for (int j = 0; j <= getCthgrid(); j++) {
      costheta_hit = -2. * j / getCthgrid() + 1.;
      sintheta_hit = sqrt(1. - costheta_hit * costheta_hit);
      for (int k = 0; k <= getPhigrid(); k++) {
        phi_hit = (getAllinplane() ? 1. : 2.) * PI * k / getPhigrid();
        if (std::isnan(phi_hit)) phi_hit = 0.;
        sincos(phi_hit, &sinphi_hit, &cosphi_hit);
        if (treshold[j][k] == 0) {  // if treshold is true, we're close enough to one and we don't have to compute everything again and again
          for (size_t it = 0; it < getParticles().size(); it++) {
            getParticles()[it].setHitcoord(r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit);
          }
          if (ctOnly) calcGlauberphasesCt(i, j, k);
          else calcGlauberphasesBoth(i, j, k);

          if (closeToOne(ctOnly, i, j, k)) {
            treshold[j][k] = 1;
            if (j == 0 || j == getCthgrid()) {
              for (int ll = 1; ll <= getPhigrid(); ll++) treshold[j][ll] = 1;
            }
          }
          // r=0 symmetry shortcut
          // (j and k are the variables of the loops above: the loops end after this, as they always did)
          if (i == 0) {
            for (j = 0; j <= getCthgrid(); j++) {
              for (k = 0; k <= getPhigrid(); k++) {
                for (int src = 0; src < variants; src++) {
                  for (int proton = 0; proton < 2; proton++) {
                    if (!ctOnly) fsi_grid[src][proton][i][j][k] = fsi_grid[src][proton][0][0][0];
                    fsi_ct_grid[src][proton][i][j][k] = fsi_ct_grid[src][proton][0][0][0];
                  }
                }
              }
            }
          }
          // theta=0 or Pi symmetry shortcut
          else if (j == 0 || j == getCthgrid()) {
            for (k = 1; k <= getPhigrid(); k++) {
              for (int src = 0; src < variants; src++) {
                for (int proton = 0; proton < 2; proton++) {
                  if (!ctOnly) fsi_grid[src][proton][i][j][k] = fsi_grid[src][proton][i][j][0];
                  fsi_ct_grid[src][proton][i][j][k] = fsi_ct_grid[src][proton][i][j][0];
                }
              }
            }
          }
        }
        // treshold reached, fill with 1s
        else {
          for (int s = 0; s < variants; s++) {
            for (int mm = 0; mm < 2; mm++) {
              if (!ctOnly) fsi_grid[s][mm][i][j][k] = 1.;
              fsi_ct_grid[s][mm][i][j][k] = 1.;
            }
          }
        }
      }
    }
  }
}

namespace {
// reading and writing a 5D array in the order of the cache files: variant, nucleon, r, cos theta, phi
template <class F>
void for_all_points(std::vector<std::vector<Grid3D<std::complex<double> > > > &grid, F f) {
  for (auto &variant : grid)
    for (auto &nucleon : variant)
      for (auto &value : nucleon) f(value);  // r, cos theta, phi: phi runs fastest
}
}  // namespace

void GlauberThickGridCommon::readinFsiGrid(ifstream &infile) {
  if (infile.peek() == std::ifstream::traits_type::eof()) {
    cerr << "File is empty " << getFsi_Filename() << endl;
    assert(1 == 0);
  }
  for_all_points(fsi_grid,
                 [&infile](complex<double> &v) { infile.read(reinterpret_cast<char *>(&v), sizeof(complex<double>)); });
}

void GlauberThickGridCommon::readinFsiCtGrid(ifstream &infile) {
  if (infile.peek() == std::ifstream::traits_type::eof()) {
    cerr << "File is empty " << getFsi_Ct_Filename() << endl;
    assert(1 == 0);
  }
  for_all_points(fsi_ct_grid,
                 [&infile](complex<double> &v) { infile.read(reinterpret_cast<char *>(&v), sizeof(complex<double>)); });
}

void GlauberThickGridCommon::writeoutFsiGrid(ofstream &outfile) {
  for_all_points(fsi_grid,
                 [&outfile](complex<double> &v) { outfile.write(reinterpret_cast<char *>(&v), sizeof(complex<double>)); });
}

void GlauberThickGridCommon::writeoutFsiCtGrid(ofstream &outfile) {
  for_all_points(fsi_ct_grid,
                 [&outfile](complex<double> &v) { outfile.write(reinterpret_cast<char *>(&v), sizeof(complex<double>)); });
}


namespace {
/*! integrator functor of the tilted frame */
struct TiltedFtor {
  static void exec(const numint::array<double,3> &x, void *param, numint::vector_d &ret) {
    TiltedFtor &p = *static_cast<TiltedFtor *>(param);
    p.f(ret, x[0], x[1], x[2], *p.grid, p.proton, p.frame);
  }
  GlauberThickGridCommon *grid;
  int proton;
  GlauberThickGridCommon::TiltedFrame frame;
  GlauberThickGridCommon::TiltedIntegrand f;
};
}  // namespace

std::complex<double> GlauberThickGridCommon::tiltedPhase(int proton, const std::vector<double> &ret, int index, double src) {
  const std::complex<double> front = getParticles()[0].getScatterfront(proton);
  if (index % 2) return (ret.back() - front * ret[index]) * src;
  return 1. - front * ret[index];
}

void GlauberThickGridCommon::srcNormIntegrand(numint::vector_d &results, double b, double z, double phi, GlauberThickGridCommon &grid,
                                              int proton, const TiltedFrame &frame) {
  results = numint::vector_d(1, 0.);
  TiltedPoint p;
  if (tiltedPoint(grid, b, z, phi, proton, frame, p)) results[0] = p.dens * p.src;
}

bool GlauberThickGridCommon::tiltedPoint(GlauberThickGridCommon &grid, double b, double z, double phi, int proton,
                                         const TiltedFrame &frame, TiltedPoint &point) {
  const double r = std::sqrt(b * b + z * z);
  if (r > grid.getPnucleus()->getRange()) return false;
  FastParticle &part = grid.getParticles()[0];
  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);
  const double hbn = part.getHitbnorm();
  point.dens = grid.getPnucleusthick()->getDensity(r, proton) / (r * r) * b;
  point.fsi = point.dens *
    std::exp(-(b * b + hbn * hbn - 2. * b * hbn * cosphi) / (2. * part.getBetasq(proton)));
  // lab coordinates of the point: b cos(phi) u + b sin(phi) w + z e, and its distance to the hit point
  const double px = b * cosphi * frame.u[0] + b * sinphi * frame.w[0] + z * part.getEx();
  const double py = b * cosphi * frame.u[1] + b * sinphi * frame.w[1] + z * part.getEy();
  const double pz = b * cosphi * frame.u[2] + b * sinphi * frame.w[2] + z * part.getEz();
  const double dx = px - grid.getR_hit() * grid.getSintheta_hit() * grid.getCosphi_hit();
  const double dy = py - grid.getR_hit() * grid.getSintheta_hit() * grid.getSinphi_hit();
  const double dz = pz - grid.getR_hit() * grid.getCostheta_hit();
  point.src = grid.getFsiCorrelator().correlation(std::sqrt(dx * dx + dy * dy + dz * dz))
              * grid.getFsiCorrelator().getCorrGrid_interp(r, pz / r, proton);
  return true;
}

void GlauberThickGridCommon::integrateTilted(int proton, TiltedIntegrand f, int nresults, int maxeval, std::vector<double> &ret) {
  const double range = getPnucleusthick()->getRange();
  FastParticle &part = getParticles()[0];
  ret.assign(nresults, 0.);
  // the frame
  TiltedFtor F;
  F.grid = this;
  F.proton = proton;
  F.f = f;
  const double e[3] = {part.getEx(), part.getEy(), part.getEz()};
  const double *hb = part.getHitbvec();
  const double norm = part.getHitbnorm();
  if (norm > 1.E-10) {
    for (int i = 0; i < 3; i++) F.frame.u[i] = hb[i] / norm;
  }
  else {  // any direction perpendicular to e
    double a[3] = {1., 0., 0.};
    if (std::fabs(e[0]) > 0.9) {a[0] = 0.; a[1] = 1.;}
    const double dot = a[0] * e[0] + a[1] * e[1] + a[2] * e[2];
    double n2 = 0.;
    for (int i = 0; i < 3; i++) {F.frame.u[i] = a[i] - dot * e[i]; n2 += F.frame.u[i] * F.frame.u[i];}
    for (int i = 0; i < 3; i++) F.frame.u[i] /= std::sqrt(n2);
  }
  F.frame.w[0] = e[1] * F.frame.u[2] - e[2] * F.frame.u[1];
  F.frame.w[1] = e[2] * F.frame.u[0] - e[0] * F.frame.u[2];
  F.frame.w[2] = e[0] * F.frame.u[1] - e[1] * F.frame.u[0];
  // the particle scatters on the nucleons it still has to pass (z > hitz), or has passed (an incoming one)
  numint::array < double, 3 > lower = { {1.E-06, part.getIncoming() ? -range : part.getHitz(), 0.} };
  numint::array < double, 3 > upper = { {range, part.getIncoming() ? part.getHitz() : range, 2. * PI} };
  unsigned count = 0;
  numint::mdfunction < numint::vector_d, 3 > mdf;
  mdf.func = &TiltedFtor::exec;
  mdf.param = &F;
  if (lower[1] < upper[1])
    numint::cube({integrator, 1.E-12, prec, 2000, maxeval}, mdf, lower, upper, ret, count, 0);
  evalCount = count;
  // the SRC normalization: the same integration over all z of dens * src
  F.f = srcNormIntegrand;
  lower[1] = -range;
  upper[1] = range;
  std::vector<double> srcnorm(1, 0.);
  // N is a narrow peak at the hit point (the correlation function): at the precision of the other integrals the
  // integrator stops too early (errors of 3e-3 at prec 1e-3), at prec/100 the error is ~1e-4. The cost is small.
  numint::cube({integrator, 1.E-12, 0.01 * prec, 2000, maxeval}, mdf, lower, upper, srcnorm, count, 0);
  evalCount += count;
  ret.push_back(srcnorm[0]);
}

}  // namespace nucqcd
