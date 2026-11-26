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
AbstractFsiCTGrid(r_grid, cth_grid, phi_grid, pnucl, prec, integrator, dir), fsi_grid(NULL), fsi_ct_grid(NULL), treshold(NULL) {

  //mem reservation for the grids, index [shelllevel][mindex(postive only due to symmetry][rgrid][thgrid][phgrid]
  //fsi grid
  if (fsi_grid == NULL) {
    //cout << "Mem reservation for fsi_grid" << endl;
    fsi_grid = new complex < double >****[getPnucleus()->getTotalLevels() + 1];
    for (int i = 0; i < (getPnucleus()->getTotalLevels() + 1); i++) {
      if (i == getPnucleus()->getTotalLevels())
        fsi_grid[i] = new complex < double >***[1];
      else
        fsi_grid[i] = new complex < double >***[(getPnucleus()->getJ_array()[i] + 1) / 2];
      for (int j = 0; j < (i == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[i] + 1) / 2)); j++) {
        fsi_grid[i][j] = new complex < double >**[getRgrid() + 1];
        for (int k = 0; k < (getRgrid() + 1); k++) {
          fsi_grid[i][j][k] = new complex < double >*[getCthgrid() + 1];
          for (int l = 0; l < (getCthgrid() + 1); l++) {
            fsi_grid[i][j][k][l] = new complex < double >[getPhigrid() + 1];
          }
        }
      }
    }
  }
  else {
    throw std::logic_error(nucqcd::ErrorMessage() << "fsi_grid should not be initialized more than once!!!!" );
  }
  //ct grid
  if (fsi_ct_grid == NULL) {
    //cout << "Mem reservation for fsi_ct_grid" << endl;
    fsi_ct_grid = new complex < double >****[getPnucleus()->getTotalLevels() + 1];
    for (int i = 0; i < (getPnucleus()->getTotalLevels() + 1); i++) {
      if (i == getPnucleus()->getTotalLevels())
        fsi_ct_grid[i] = new complex < double >***[1];
      else
        fsi_ct_grid[i] = new complex < double >***[(getPnucleus()->getJ_array()[i] + 1) / 2];
      for (int j = 0; j < (i == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[i] + 1) / 2)); j++) {
        fsi_ct_grid[i][j] = new complex < double >**[getRgrid() + 1];
        for (int k = 0; k < (getRgrid() + 1); k++) {
          fsi_ct_grid[i][j][k] = new complex < double >*[getCthgrid() + 1];
          for (int l = 0; l < (getCthgrid() + 1); l++) {
            fsi_ct_grid[i][j][k][l] = new complex < double >[getPhigrid() + 1];
          }
        }
      }
    }
  }
  else {
    throw std::logic_error(nucqcd::ErrorMessage() << "fsi_ct_grid should not be initialized more than once!!!!" );
  }
}

//destructor
GlauberGrid::~GlauberGrid() {
  if (fsi_grid != NULL) {
    for (int i = 0; i < (getPnucleus()->getTotalLevels() + 1); i++) {
      for (int j = 0; j < (i == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[i] + 1) / 2)); j++) {
        for (int k = 0; k < (getRgrid() + 1); k++) {
          for (int l = 0; l < (getCthgrid() + 1); l++) {
            delete[]fsi_grid[i][j][k][l];
            delete[]fsi_ct_grid[i][j][k][l];
          }
          delete[]fsi_grid[i][j][k];
          delete[]fsi_ct_grid[i][j][k];
        }
        delete[]fsi_grid[i][j];
        delete[]fsi_ct_grid[i][j];
      }
      delete[]fsi_grid[i];
      delete[]fsi_ct_grid[i];
    }
    delete[]fsi_grid;
    delete[]fsi_ct_grid;
  }
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

//calc both fsi and fsi+ct grid
void GlauberGrid::constructAllGrids() {
  //mem reservation for treshold array
  if (treshold == NULL) {
    //cout << "Mem reservation for treshold array" << endl;
    treshold = new int *[getCthgrid() + 1];
    for (int i = 0; i <= getCthgrid(); i++) {
      treshold[i] = new int[getPhigrid() + 1];
      for (int j = 0; j <= getPhigrid(); j++)
        treshold[i][j] = 0;
    }
  }
  else {
    throw std::logic_error(nucqcd::ErrorMessage() << "Treshold should only be initialized once!" );
  }
  
  //fill the arrays!
  for (int i = 0; i <= getRgrid(); i++) {
    r_hit = float (i) * getPnucleus()->getRange() / getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for (int j = 0; j <= getCthgrid(); j++) {
      costheta_hit = -2. * j / getCthgrid() + 1.;
      //sincos(theta_hit,&sintheta_hit,&costheta_hit);
      sintheta_hit = sqrt(1. - costheta_hit * costheta_hit);
      for (int k = 0; k <= getPhigrid(); k++) {
        phi_hit = (getAllinplane()? 1. : 2.) * PI * k / getPhigrid();
        if (std::isnan(phi_hit))
          phi_hit = 0.;
        sincos(phi_hit, &sinphi_hit, &cosphi_hit);
        if (treshold[j][k] == 0) {  //if treshold is true, we're close enough to one and we don't have to compute everything again and again
          for (size_t it = 0; it < getParticles().size(); it++) {
            getParticles()[it].setHitcoord(r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit);
          }
          calcGlauberphasesBoth(i, j, k);

          fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = 1.;
          fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = 1.;

          for (int s = 0; s < getPnucleus()->getPLevels() - 1; s++) {
            for (int mm = 0; mm < ((getPnucleus()->getJ_array()[s] + 1) / 2); mm++) {
              fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_grid[s][mm][i][j][k], 2.);
              fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[s][mm][i][j][k], 2.);
            }
          }
          //special care for last shell if it's not full
          //finalmproton was determined earlier!
          for (int mm = 0; mm < getPnucleus()->getFinalMProton(); mm++) {
            fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_grid[getPnucleus()->getPLevels() - 1][mm][i][j][k], 2.);
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[getPnucleus()->getPLevels() - 1][mm][i][j][k], 2.);
          }
          //odd # of protons in last shell
          if (getPnucleus()->getOnlyOneProton()) {
            fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= fsi_grid[getPnucleus()->getPLevels() - 1][getPnucleus()->getFinalMProton()][i][j][k];
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= fsi_ct_grid[getPnucleus()->getPLevels() - 1][getPnucleus()->getFinalMProton()][i][j][k];
          }
          //pair # of protons in last shell
          else {
            fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_grid[getPnucleus()->getPLevels() - 1][getPnucleus()->getFinalMProton()][i][j][k], 2.);
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[getPnucleus()->getPLevels() - 1][getPnucleus()->getFinalMProton()][i][j][k], 2.);
          }
          //multiply all neutron phases except last shell
          for (int s = getPnucleus()->getPLevels(); s < (getPnucleus()->getTotalLevels() - 1); s++) {
            for (int mm = 0; mm < ((getPnucleus()->getJ_array()[s] + 1) / 2); mm++) {
              fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_grid[s][mm][i][j][k], 2.);
              fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[s][mm][i][j][k], 2.);
            }
          }
          //special care for last shell if it's not full
          //getPnucleus()->getFinalMNeutron was determined earlier!
          for (int mm = 0; mm < getPnucleus()->getFinalMNeutron(); mm++) {
            fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_grid[getPnucleus()->getTotalLevels() - 1][mm][i][j][k], 2.);
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[getPnucleus()->getTotalLevels() - 1][mm][i][j][k], 2.);
          }
          //odd # of neutrons in last shell
          if (getPnucleus()->getOnlyOneNeutron()) {
            fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= fsi_grid[getPnucleus()->getTotalLevels() - 1][getPnucleus()->getFinalMNeutron()][i][j][k];
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= fsi_ct_grid[getPnucleus()->getTotalLevels() - 1][getPnucleus()->getFinalMNeutron()][i][j][k];
          }
          //pair # of neutrons in last shell
          else {
            fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_grid[getPnucleus()->getTotalLevels() - 1][getPnucleus()->getFinalMNeutron()][i][j][k], 2.);
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[getPnucleus()->getTotalLevels() - 1][getPnucleus()->getFinalMNeutron()][i][j][k], 2.);
          }
          //cout << r_hit << " " << theta_hit*nucqcd::rad_to_deg << " " << phi_hit*nucqcd::rad_to_deg << " " << fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] << endl;

          if (abs(fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k]) > 0.999999 && abs(fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k]) > 0.999999) {
            treshold[j][k] = 1;
            if (j == 0 || j == getCthgrid()) {
              for (int ll = 1; ll <= getPhigrid(); ll++)
                treshold[j][ll] = 1;
            }
          }
          //r=0 symmetry shortcut
          if (i == 0) {
            for (j = 0; j <= getCthgrid(); j++) {
              for (k = 0; k <= getPhigrid(); k++) {
                for (int level = 0; level < getPnucleus()->getTotalLevels(); level++) {
                  for (int mm = 0; mm < ((getPnucleus()->getJ_array()[level] + 1) / 2); mm++) {
                    fsi_grid[level][mm][i][j][k] = fsi_grid[level][mm][0][0][0];
                    fsi_ct_grid[level][mm][i][j][k] = fsi_ct_grid[level][mm][0][0][0];

                  }
                }
                fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = fsi_grid[getPnucleus()->getTotalLevels()][0][0][0][0];
                fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = fsi_ct_grid[getPnucleus()->getTotalLevels()][0][0][0][0];
                //cout << r_hit << " " << theta_hit*nucqcd::rad_to_deg << " " << phi_hit*nucqcd::rad_to_deg << " " << fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] << endl;
              }
            }
          }
          //theta=0 or Pi symmetry shortcut
          else if (j == 0 || j == getCthgrid()) {
            for (k = 1; k <= getPhigrid(); k++) {
              for (int level = 0; level < getPnucleus()->getTotalLevels(); level++) {
                for (int mm = 0; mm < ((getPnucleus()->getJ_array()[level] + 1) / 2); mm++) {
                  fsi_grid[level][mm][i][j][k] = fsi_grid[level][mm][i][j][0];
                  fsi_ct_grid[level][mm][i][j][k] = fsi_ct_grid[level][mm][i][j][0];
                }
              }
              fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][0];
              fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][0];
              //cout << r_hit << " " << theta_hit*nucqcd::rad_to_deg << " " << phi_hit*nucqcd::rad_to_deg << " " << fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] << endl;
            }
          }
        }
        //treshold reached, fill with 1s
        else {
          for (int s = 0; s < (getPnucleus()->getTotalLevels() + 1); s++) {
            for (int mm = 0; mm < (s == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[s] + 1) / 2)); mm++) {
              fsi_grid[s][mm][i][j][k] = 1.;
              fsi_ct_grid[s][mm][i][j][k] = 1.;
            }
          }
          //cout << r_hit << " " << theta_hit*nucqcd::rad_to_deg << " " << phi_hit*nucqcd::rad_to_deg << " " << fsi_grid[getPnucleus()->getTotalLevels()][0][i][j][k] << endl;
        }
      }
    }
  }
  //mem cleanup
  for (int i = 0; i <= getCthgrid(); i++)
    delete[]treshold[i];
  delete[]treshold;
  treshold = NULL;
}

//calc only ct grid
void GlauberGrid::constructCtGrid() {
  //mem reservation for treshold
  if (treshold == NULL) {
    treshold = new int *[getCthgrid() + 1];
    for (int i = 0; i <= getCthgrid(); i++) {
      treshold[i] = new int[getPhigrid() + 1];
      for (int j = 0; j <= getPhigrid(); j++)
        treshold[i][j] = 0;
    }
  }
  else {
    throw std::logic_error(nucqcd::ErrorMessage() << "Treshold should only be initialized once!" );
  }
  //fill ct array!
  for (int i = 0; i <= getRgrid(); i++) {
    r_hit = float (i) * getPnucleus()->getRange() / getRgrid();
    //if(r_hit==0.) r_hit=0.001;
    for (int j = 0; j <= getCthgrid(); j++) {
      costheta_hit = -2. * j / getCthgrid() + 1.;
      //sincos(theta_hit,&sintheta_hit,&costheta_hit);
      sintheta_hit = sqrt(1. - costheta_hit * costheta_hit);
      for (int k = 0; k <= getPhigrid(); k++) {
        phi_hit = (getAllinplane()? 1. : 2.) * PI * k / getPhigrid();
        if (std::isnan(phi_hit))
          phi_hit = 0.;
        sincos(phi_hit, &sinphi_hit, &cosphi_hit);
        if (treshold[j][k] == 0) {  //if treshold is true, we're close enough to one and we don't have to compute everything again and again
          //set z coord along all fast particles for coord of hard interaction
          for (size_t it = 0; it < getParticles().size(); it++) {
            getParticles()[it].setHitcoord(r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit);
          }
          calcGlauberphasesCt(i, j, k);
          fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = 1.;
          for (int s = 0; s < getPnucleus()->getPLevels() - 1; s++) {
            for (int mm = 0; mm < ((getPnucleus()->getJ_array()[s] + 1) / 2); mm++) {
              fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[s][mm][i][j][k], 2.);
            }
          }
          //special care for last shell if it's not full
          //finalmproton was determined earlier!
          for (int mm = 0; mm < getPnucleus()->getFinalMProton(); mm++) {
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[getPnucleus()->getPLevels() - 1][mm][i][j][k], 2.);
          }
          //odd # of protons in last shell
          if (getPnucleus()->getOnlyOneProton()) {
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= fsi_ct_grid[getPnucleus()->getPLevels() - 1][getPnucleus()->getFinalMProton()][i][j][k];
          }
          //pair # of protons in last shell
          else {
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[getPnucleus()->getPLevels() - 1][getPnucleus()->getFinalMProton()][i][j][k], 2.);
          }
          //multiply all neutron phases except last shell
          for (int s = getPnucleus()->getPLevels(); s < (getPnucleus()->getTotalLevels() - 1); s++) {
            for (int mm = 0; mm < ((getPnucleus()->getJ_array()[s] + 1) / 2); mm++) {
              fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[s][mm][i][j][k], 2.);
            }
          }
          //special care for last shell if it's not full
          //getPnucleus()->getFinalMNeutron was determined earlier!
          for (int mm = 0; mm < getPnucleus()->getFinalMNeutron(); mm++) {
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[getPnucleus()->getTotalLevels() - 1][mm][i][j][k], 2.);
          }
          //odd # of neutrons in last shell
          if (getPnucleus()->getOnlyOneNeutron()) {
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= fsi_ct_grid[getPnucleus()->getTotalLevels() - 1][getPnucleus()->getFinalMNeutron()][i][j][k];
          }
          //pair # of neutrons in last shell
          else {
            fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] *= pow(fsi_ct_grid[getPnucleus()->getTotalLevels() - 1][getPnucleus()->getFinalMNeutron()][i][j][k], 2.);
          }
          if (abs(fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k]) > 0.999999) {
            treshold[j][k] = 1;
            if (j == 0 || j == getCthgrid()) {
              for (int ll = 1; ll <= getPhigrid(); ll++)
                treshold[j][ll] = 1;
            }
          }
          //r=0 symmetry shortcut
          if (i == 0) {
            for (j = 0; j <= getCthgrid(); j++) {
              for (k = 0; k <= getPhigrid(); k++) {
                for (int level = 0; level < getPnucleus()->getTotalLevels(); level++) {
                  for (int mm = 0; mm < ((getPnucleus()->getJ_array()[level] + 1) / 2); mm++) {
                    fsi_ct_grid[level][mm][i][j][k] = fsi_ct_grid[level][mm][0][0][0];
                  }
                }
                fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = fsi_ct_grid[getPnucleus()->getTotalLevels()][0][0][0][0];
              }
            }
          }
          //theta=0 or Pi symmetry shortcut
          else if (j == 0 || j == getCthgrid()) {
            for (k = 1; k <= getPhigrid(); k++) {
              for (int level = 0; level < getPnucleus()->getTotalLevels(); level++) {
                for (int mm = 0; mm < ((getPnucleus()->getJ_array()[level] + 1) / 2); mm++) {
                  fsi_grid[level][mm][i][j][k] = fsi_grid[level][mm][i][j][0];
                  fsi_ct_grid[level][mm][i][j][k] = fsi_ct_grid[level][mm][i][j][0];
                }
              }
              fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][k] = fsi_ct_grid[getPnucleus()->getTotalLevels()][0][i][j][0];
            }
          }

        }
        //treshold reached, fill with 1s
        else {
          for (int s = 0; s < (getPnucleus()->getTotalLevels() + 1); s++) {
            for (int mm = 0; mm < (s == getPnucleus()->getTotalLevels()? 1 : ((getPnucleus()->getJ_array()[s] + 1) / 2)); mm++) {
              fsi_ct_grid[s][mm][i][j][k] = 1.;
            }
          }
        }
      }
    }
  }
  for (int i = 0; i <= getCthgrid(); i++)
    delete[]treshold[i];
  delete[]treshold;
  treshold = NULL;
}

//calc glauberphases for one gridpoint
void GlauberGrid::calcGlauberphasesBoth(const int i, const int j, const int k) {
  for (int level = 0; level < getPnucleus()->getTotalLevels(); level++) {
    for (int mm = 0; mm < ((getPnucleus()->getJ_array()[level] + 1) / 2); mm++) {
      int res = 90;
      unsigned count = 0;
      double deeserror = 0.;
      if (getParticles().size() == 1 && getParticles()[0].getCosTheta() == 1.) {
        throw std::invalid_argument(nucqcd::ErrorMessage() << "You should be using the OneGlauberGrid class boy!" );
      }
      else {
        numint::array < double, 3 > lower = { {0., -1., 0.} };
        numint::array < double, 3 > upper = { {getPnucleus()->getRange(), 1., 2. * PI} };

        GlauberGrid::Ftor_one F;
        F.grid = this;
        F.level = level;
        F.mm = mm;
        numint::mdfunction < numint::vector_z, 3 > mdf;
        mdf.func = &Ftor_one::exec;
        mdf.param = &F;
        vector<complex<double> > ret(2, 0.);
        F.f = klaas_one_bound;
        res = numint::cube({integrator, 1.E-12, prec, 2000, 2000000}, mdf, lower, upper, ret, count, 0);
        fsi_grid[level][mm][i][j][k] = ret[0];
        fsi_ct_grid[level][mm][i][j][k] = ret[1];

      }
//       cout << i << " " << j << " " << k << " " << level << " " << mm << fsi_grid[level][mm][i][j][k] << " " << fsi_ct_grid[level][mm][i][j][k] <<
//       " " << res << " " << count << endl;
    }
  }
}

//calc glauberphases for one gridpoint,only CT grid
void GlauberGrid::calcGlauberphasesCt(const int i, const int j, const int k) {
  for (int level = 0; level < getPnucleus()->getTotalLevels(); level++) {
    for (int mm = 0; mm < ((getPnucleus()->getJ_array()[level] + 1) / 2); mm++) {
      int res = 90;
      unsigned count = 0;
      double deeserror = 0.;
      if (getParticles().size() == 1 && getParticles()[0].getCosTheta() == 1.) {
        throw std::invalid_argument(nucqcd::ErrorMessage() << "You should be using the OneGlauberGrid class boy!" );
      }
      else {
        numint::array < double, 3 > lower = { {0., -1., 0.} };
        numint::array < double, 3 > upper = { {getPnucleus()->getRange(), 1., 2. * PI} };

        GlauberGrid::Ftor_one F;
        F.grid = this;
        F.level = level;
        F.mm = mm;
        numint::mdfunction < numint::vector_z, 3 > mdf;
        mdf.func = &Ftor_one::exec;
        mdf.param = &F;
        vector<complex<double> > ret(1, 0.);
        F.f = klaas_one_bound_ct;
        res = numint::cube({integrator, 1.E-12, prec, 2000, 2000000}, mdf, lower, upper, ret, count, 0);
        fsi_ct_grid[level][mm][i][j][k] = ret[2];

      }
//       cout << i << " " << j << " " << k << " " << level << " " << mm << fsi_grid[level][mm][i][j][k] << " " << fsi_ct_grid[level][mm][i][j][k] <<
//       " " << res << " " << count << endl;
    }
  }
}

//readin fsi grid
void GlauberGrid::readinFsiGrid(ifstream & infile) {
  if(infile.peek() == std::ifstream::traits_type::eof()){
    cerr << "File is empty " << getFsi_Filename() << endl;
    assert(1==0);
  }
  
  if (fsi_grid != NULL) {
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
  if (fsi_ct_grid != NULL) {
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

void GlauberGrid::klaas_one_bound(numint::vector_z & results, double r, double costheta, double phi, GlauberGrid & grid, int level, int m) {
  double sintheta = sqrt(1. - costheta * costheta);
  results = numint::vector_z(2, 1.);
  if (r > grid.getPnucleus()->getRange())
    return;
  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);
  
  for (size_t it = 0; it < grid.getParticles().size(); it++) {
    double zmom = grid.getParticles()[it].calcZ(r, costheta, sintheta, cosphi, sinphi);
    bool check = (zmom > grid.getParticles()[it].getHitz());
    if (grid.getParticles()[it].getIncoming())
      check = !check;
    if (check) {
      complex < double >temp = grid.getParticles()[it].getScatterfront(level, grid.getPnucleus())
        * exp(-grid.getParticles()[it].getBdist(r, costheta, sintheta, cosphi, sinphi, zmom)
              / (2. * grid.getParticles()[it].getBetasq(level, grid.getPnucleus())));
      results[0] *= 1. - temp;
      results[1] *= 1. - temp * grid.getParticles()[it].getCTsigma(zmom);
    }
  }
  double wf = (power(grid.getPnucleus()->getWave_F(level, r), 2) * grid.getPnucleus()->getYminkappacos(level, m, costheta)
               + power(grid.getPnucleus()->getWave_G(level, r), 2) * grid.getPnucleus()->getYkappacos(level, m, costheta));
  for (int i = 0; i < 2; ++i)
    results[i] *= wf;
  return;
}

void GlauberGrid::klaas_one_bound_ct(numint::vector_z & results, double r, double costheta, double phi, GlauberGrid & grid, int level, int m) {
  double sintheta = sqrt(1. - costheta * costheta);
  results = numint::vector_z(1, 1.);
  if (r > grid.getPnucleus()->getRange())
    return;
  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);
  
  for (size_t it = 0; it < grid.getParticles().size(); it++) {
    double zmom = grid.getParticles()[it].calcZ(r, costheta, sintheta, cosphi, sinphi);
    bool check = (zmom > grid.getParticles()[it].getHitz());
    if (grid.getParticles()[it].getIncoming())
      check = !check;
    if (check) {
      results[0] *= 1. - grid.getParticles()[it].getScatterfront(level, grid.getPnucleus())
        * exp(-grid.getParticles()[it].getBdist(r, costheta, sintheta, cosphi, sinphi, zmom)
              / (2. * grid.getParticles()[it].getBetasq(level, grid.getPnucleus())))
        * grid.getParticles()[it].getCTsigma(zmom);
    }
  }
  results[0] *= (power(grid.getPnucleus()->getWave_F(level, r), 2) * grid.getPnucleus()->getYminkappacos(level, m, costheta)
                 + power(grid.getPnucleus()->getWave_G(level, r), 2) * grid.getPnucleus()->getYkappacos(level, m, costheta));
  return;
}

}  // namespace nucqcd

