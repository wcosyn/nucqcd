#include <cstdlib>
#include <Integrator/numint/cube.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <iostream>
#include <fstream>
#include <complex>
#include <cassert>

using namespace std;

#include <Glauber/GlauberGrid.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/AuxFunction.hpp>

namespace nucqcd {


//constructor, calls abstractfsigrid's constructor
GlauberGrid::GlauberGrid(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleus * pnucl,
                         double prec, IntegrationMethod integrator, string dir):
AbstractFsiCTGrid(r_grid, cth_grid, phi_grid, pnucl, prec, integrator, dir) {

  //the grids, index [shelllevel][mindex (positive only due to symmetry)][rgrid][thgrid][phgrid]; the last shell holds the product of all
  const int shells = getPnucleus()->getTotalLevels() + 1;
  const Grid3D<complex<double> > points(getRgrid() + 1, getCthgrid() + 1, getPhigrid() + 1);
  fsi_grid.resize(shells);
  fsi_ct_grid.resize(shells);
  for (int i = 0; i < shells; i++) {
    const int mvalues = (i == getPnucleus()->getTotalLevels() ? 1 : ((getPnucleus()->getJ_array()[i] + 1) / 2));
    fsi_grid[i].assign(mvalues, points);
    fsi_ct_grid[i].assign(mvalues, points);
  }
}

//destructor
GlauberGrid::~GlauberGrid() {
}

//interpolation of grid number grid (0: RMSGA, 1: RMSGA+CT) after r,cth,phi have been set
complex < double >GlauberGrid::getFsiGridN_interp(int grid) {
  if (grid != 0 && grid != 1) {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "grid number not supported!" );
  }
  if (!filledallgrid) {
    throw std::logic_error(nucqcd::ErrorMessage() << "You have to fill the grids first!" );
  }
  auto &array = (grid == 0 ? fsi_grid : fsi_ct_grid);
  complex < double >result = getInterp(array[getPnucleus()->getTotalLevels()][0]);
  for (size_t it = 0; it < getKnockoutLevels().size(); it++) {
    result /= getInterp(array[getKnockoutLevels()[it]][(abs(getKnockoutM()[it]) - 1) / 2]);
  }
  return result;
}

int GlauberGrid::getGridNumber(GridVariant variant) const {
  if (variant == GridVariant::Fsi) return 0;
  if (variant == GridVariant::Ct) return 1;
  return -1;
}

void GlauberGrid::printFsi_grid() {
  cout << "Printing Glauberarray for " << getFsi_Filename() << endl;
  print_grid(0);
}
void GlauberGrid::printFsi_ct_grid() {
  cout << "Printing CT Glauberarray for " << getFsi_Ct_Filename() << endl;
  print_grid(1);
}

void GlauberGrid::print_grid(int gridindex) {
  printKnockout();
  for (int i = 0; i <= getRgrid(); i++) {
    double r = float (i) * getPnucleus()->getRange() / getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for (int j = 0; j <= getCthgrid(); j++) {
      double costheta = -2. * j / getCthgrid() + 1.;
      for (int k = 0; k <= getPhigrid(); k++) {
        double phi = (getAllinplane()? 1. : 2.) * PI * k / getPhigrid();
        if (std::isnan(phi))
          phi = 0.;
        setinterp(r,costheta,phi);
        complex < double >value = getFsiGridN_interp(gridindex);
        cout << r << " " << costheta << " " << phi * nucqcd::rad_to_deg << " " << real(value) << " " << imag(value) << endl;
      }
    }
  }
}

//set filenames, don't forget to add ending "/" to dir if necessary!!!
void GlauberGrid::setFilenames(string dir) {
  AbstractFsiCTGrid::setFilenames(dir + "Gl.");
}

//readin fsi grid
void GlauberGrid::readinFsiGrid(ifstream & infile) {
  if(infile.peek() == std::ifstream::traits_type::eof()){
    cerr << "File is empty " << getFsi_Filename() << endl;
    assert(1==0);
  }
  
  if (!fsi_grid.empty()) {
    for (int i = 0; i < (getPnucleus()->getTotalLevels() + 1); i++) {
      for (int j = 0; j < (i == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[i] + 1) / 2)); j++) {
        for (int k = 0; k < (getRgrid() + 1); k++) {
          for (int l = 0; l < (getCthgrid() + 1); l++) {
            for (int mm = 0; mm < (getPhigrid() + 1); mm++) {
              infile.read(reinterpret_cast < char *>(&fsi_grid[i][j][k][l][mm]), sizeof(complex < double >));
            }
          }
        }
      }
    }
  }
  else {
    throw std::logic_error(nucqcd::ErrorMessage() << "fsi_grid not initialized!!!!" );
  }
}

//readin fsi+ct grid
void GlauberGrid::readinFsiCtGrid(ifstream & infile) {
  if(infile.peek() == std::ifstream::traits_type::eof()){
    cerr << "File is empty " << getFsi_Ct_Filename() << endl;
    assert(1==0);
  }
  if (!fsi_ct_grid.empty()) {
    for (int i = 0; i < (getPnucleus()->getTotalLevels() + 1); i++) {
      for (int j = 0; j < (i == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[i] + 1) / 2)); j++) {
        for (int k = 0; k < (getRgrid() + 1); k++) {
          for (int l = 0; l < (getCthgrid() + 1); l++) {
            for (int mm = 0; mm < (getPhigrid() + 1); mm++) {
              infile.read(reinterpret_cast < char *>(&fsi_ct_grid[i][j][k][l][mm]), sizeof(complex < double >));
            }
          }
        }
      }
    }
  }
  else {
    throw std::logic_error(nucqcd::ErrorMessage() << "fsi_ct_grid not initialized!!!!" );
  }
}

//write fsi grid to file
void GlauberGrid::writeoutFsiGrid(ofstream & outfile) {
  for (int i = 0; i < (getPnucleus()->getTotalLevels() + 1); i++) {
    for (int j = 0; j < (i == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[i] + 1) / 2)); j++) {
      for (int k = 0; k < (getRgrid() + 1); k++) {
        for (int l = 0; l < (getCthgrid() + 1); l++) {
          for (int mm = 0; mm < (getPhigrid() + 1); mm++) {
            outfile.write(reinterpret_cast < char *>(&fsi_grid[i][j][k][l][mm]), sizeof(complex < double >));
          }
        }
      }
    }
  }
}

//write fsi+ct grid to file
void GlauberGrid::writeoutFsiCtGrid(ofstream & outfile) {
  for (int i = 0; i < (getPnucleus()->getTotalLevels() + 1); i++) {
    for (int j = 0; j < (i == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[i] + 1) / 2)); j++) {
      for (int k = 0; k < (getRgrid() + 1); k++) {
        for (int l = 0; l < (getCthgrid() + 1); l++) {
          for (int mm = 0; mm < (getPhigrid() + 1); mm++) {
            outfile.write(reinterpret_cast < char *>(&fsi_ct_grid[i][j][k][l][mm]), sizeof(complex < double >));
          }
        }
      }
    }
  }
}

// 

//   }
// }
// 
// 
// 
//   }
// }

}  // namespace nucqcd

