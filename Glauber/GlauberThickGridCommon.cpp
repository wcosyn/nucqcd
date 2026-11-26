#include <Glauber/GlauberThickGridCommon.hpp>

#include <Tools/Core/ErrorMessage.hpp>
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
    variants(nvariants), pnucleusthick(pnuclthick), fsicorrelator(pnuclthick, 90, 90, dir),
    fsi_grid(NULL), fsi_ct_grid(NULL), treshold(NULL), error(0.) {
  number_of_grids = 2 * nvariants;
  auto allocate = [this]() {
    std::complex<double> *****grid = new complex<double> ****[variants];
    for (int i = 0; i < variants; i++) {
      grid[i] = new complex<double> ***[2];
      for (int j = 0; j < 2; j++) {
        grid[i][j] = new complex<double> **[getRgrid() + 1];
        for (int k = 0; k < (getRgrid() + 1); k++) {
          grid[i][j][k] = new complex<double> *[getCthgrid() + 1];
          for (int l = 0; l < (getCthgrid() + 1); l++) {
            grid[i][j][k][l] = new complex<double>[getPhigrid() + 1];
          }
        }
      }
    }
    return grid;
  };
  fsi_grid = allocate();
  fsi_ct_grid = allocate();
}

GlauberThickGridCommon::~GlauberThickGridCommon() {
  auto release = [this](std::complex<double> *****grid) {
    for (int i = 0; i < variants; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < (getRgrid() + 1); k++) {
          for (int l = 0; l < (getCthgrid() + 1); l++) delete[] grid[i][j][k][l];
          delete[] grid[i][j][k];
        }
        delete[] grid[i][j];
      }
      delete[] grid[i];
    }
    delete[] grid;
  };
  release(fsi_grid);
  release(fsi_ct_grid);
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

complex<double> GlauberThickGridCommon::interpolateVariant(std::complex<double> *****array, int variant) {
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
  if (treshold == NULL) {
    treshold = new int *[getCthgrid() + 1];
    for (int i = 0; i <= getCthgrid(); i++) {
      treshold[i] = new int[getPhigrid() + 1];
      for (int j = 0; j <= getPhigrid(); j++) treshold[i][j] = 0;
    }
  }
  else {
    throw std::logic_error(nucqcd::ErrorMessage() << "Treshold should only be initialized once!");
  }

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
  // mem cleanup
  for (int i = 0; i <= getCthgrid(); i++) delete[] treshold[i];
  delete[] treshold;
  treshold = NULL;
}

namespace {
// reading and writing a 5D array in the order of the cache files: variant, nucleon, r, cos theta, phi
template <class F>
void for_all_points(std::complex<double> *****grid, int variants, int nr, int nth, int nph, F f) {
  for (int i = 0; i < variants; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < nr; k++)
        for (int l = 0; l < nth; l++)
          for (int mm = 0; mm < nph; mm++) f(grid[i][j][k][l][mm]);
}
}  // namespace

void GlauberThickGridCommon::readinFsiGrid(ifstream &infile) {
  if (infile.peek() == std::ifstream::traits_type::eof()) {
    cerr << "File is empty " << getFsi_Filename() << endl;
    assert(1 == 0);
  }
  for_all_points(fsi_grid, variants, getRgrid() + 1, getCthgrid() + 1, getPhigrid() + 1,
                 [&infile](complex<double> &v) { infile.read(reinterpret_cast<char *>(&v), sizeof(complex<double>)); });
}

void GlauberThickGridCommon::readinFsiCtGrid(ifstream &infile) {
  if (infile.peek() == std::ifstream::traits_type::eof()) {
    cerr << "File is empty " << getFsi_Ct_Filename() << endl;
    assert(1 == 0);
  }
  for_all_points(fsi_ct_grid, variants, getRgrid() + 1, getCthgrid() + 1, getPhigrid() + 1,
                 [&infile](complex<double> &v) { infile.read(reinterpret_cast<char *>(&v), sizeof(complex<double>)); });
}

void GlauberThickGridCommon::writeoutFsiGrid(ofstream &outfile) {
  for_all_points(fsi_grid, variants, getRgrid() + 1, getCthgrid() + 1, getPhigrid() + 1,
                 [&outfile](complex<double> &v) { outfile.write(reinterpret_cast<char *>(&v), sizeof(complex<double>)); });
}

void GlauberThickGridCommon::writeoutFsiCtGrid(ofstream &outfile) {
  for_all_points(fsi_ct_grid, variants, getRgrid() + 1, getCthgrid() + 1, getPhigrid() + 1,
                 [&outfile](complex<double> &v) { outfile.write(reinterpret_cast<char *>(&v), sizeof(complex<double>)); });
}

}  // namespace nucqcd
