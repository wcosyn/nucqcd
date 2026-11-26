#include <cstdlib>
#include <Tools/Core/ErrorMessage.hpp>
#include <iostream>
#include <fstream>
#include <complex>
#include <Integrator/adaptive/cubature.h>
#include <cassert>

using namespace std;

#include <Glauber/GlauberGridThick.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/AuxFunction.hpp>

namespace nucqcd {


GlauberGridThick::GlauberGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick * pnucl,
                                   double prec, int integrator, string dir):AbstractFsiGrid(r_grid, cth_grid, phi_grid, pnucl, prec, integrator, dir),
AbstractFsiCTGrid(r_grid, cth_grid, phi_grid, pnucl, prec, integrator, dir),
AbstractFsiCTGridThick(r_grid, cth_grid, phi_grid, pnucl, prec, integrator, dir),
GlauberThickGridCommon(2, r_grid, cth_grid, phi_grid, pnucl, prec, integrator, dir) {
}

GlauberGridThick::~GlauberGridThick() {
}



//interpolation of the grid after r,th,phi have been set
complex < double >GlauberGridThick::getFsiGridFull_interp() {
  return interpolateVariant(fsi_grid, 0);
}

//interpolation of the grid after r,th,phi have been set
complex < double >GlauberGridThick::getFsiCtGridFull_interp() {
  return interpolateVariant(fsi_ct_grid, 0);
}

//interpolation of the grid after r,th,phi have been set
complex < double >GlauberGridThick::getFsiSrcGridFull_interp() {
  return interpolateVariant(fsi_grid, 1);
}

//interpolation of the grid after r,th,phi have been set
complex < double >GlauberGridThick::getFsiSrcCtGridFull_interp() {
  return interpolateVariant(fsi_ct_grid, 1);
}



void GlauberGridThick::printFsi_grid() {
  cout << "Printing Glauberarray for " << getFsi_Filename() << endl;
  print_grid(0);
}
void GlauberGridThick::printFsi_ct_grid() {
  print_grid(2);
}

void GlauberGridThick::printFsi_src_grid() {
  print_grid(1);
}
void GlauberGridThick::printFsi_src_ct_grid() {
  print_grid(3);
}

//set filenames, don't forget to add ending "/" to dir if necessary!!!
void GlauberGridThick::setFilenames(string dir) {
  AbstractFsiCTGridThick::setFilenames(dir + "Gl.");
}

void GlauberGridThick::getFsiphaseAll(vector < complex < double > >&phases, double rhit, double costhetahit, double phihit) {
  phases = vector < complex < double > >(4, 0.);
  r_hit = rhit;
  costheta_hit = costhetahit;
  sintheta_hit = sqrt(1. - costheta_hit * costheta_hit);
  phi_hit = phihit;
  sincos(phi_hit, &sinphi_hit, &cosphi_hit);
  vector < complex < double > >protonphases(4, 0.);
  vector < complex < double > >neutronphases(4, 0.);
  for (size_t it = 0; it < getParticles().size(); it++) {
    getParticles()[it].setHitcoord(r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit);
  }
  calcGlauberphasesBoth(protonphases, neutronphases);
  for (int i = 0; i < 4; i++)
    phases[i] = pow(protonphases[i], getPnucleusthick()->getZ() - getProtonKnockout()) * pow(neutronphases[i], getPnucleusthick()->getN() - getNeutronKnockout());
}

void GlauberGridThick::getFsiphaseCt(vector < complex < double > >&phases, double rhit, double costhetahit, double phihit) {
  phases = vector < complex < double > >(2, 0.);
  r_hit = rhit;
  costheta_hit = costhetahit;
  sintheta_hit = sqrt(1. - costheta_hit * costheta_hit);
  phi_hit = phihit;
  sincos(phi_hit, &sinphi_hit, &cosphi_hit);
  vector < complex < double > >protonphases(2, 0.);
  vector < complex < double > >neutronphases(2, 0.);
  for (size_t it = 0; it < getParticles().size(); it++) {
    getParticles()[it].setHitcoord(r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit);
  }
  calcGlauberphasesCt(protonphases, neutronphases);
  for (int i = 0; i < 2; i++)
    phases[i] = pow(protonphases[i], getPnucleusthick()->getZ() - getProtonKnockout()) * pow(neutronphases[i], getPnucleusthick()->getN() - getNeutronKnockout());
}



//calc glauberphases for one gridpoint
void GlauberGridThick::calcGlauberphasesBoth(const int i, const int j, const int k) {
  for (int proton = 0; proton < 2; proton++) {
    int res = 90;
    unsigned count = 0;
    double deeserror = 0.;
    double src = getFsiCorrelator().getCorrGrid_interp(r_hit, costheta_hit, proton);
    if (integrator == 1 || integrator == 2) {
      if (getParticles().size() == 1 && getParticles()[0].getCosTheta() == 1.) {
        numint::array < double, 3 > lower = { {1.E-06, getParticles()[0].getHitz(), 0.} };
        numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), getPnucleusthick()->getRange(), 2. * PI} };
        
        GlauberGridThick::Ftor_bound F;
        F.grid = this;
        F.proton = proton;
        numint::mdfunction < numint::vector_d, 3 > mdf;
        mdf.func = &Ftor_bound::exec;
        mdf.param = &F;
        vector < double >ret(4, 0.);
        F.f = klaas_int_bound;
        if (integrator == 1)
          res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
        else
          res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E03,2E06, ret, count, 0);
        fsi_grid[0][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[0];
        fsi_grid[1][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[1] * src;
        fsi_ct_grid[0][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[2];
        fsi_ct_grid[1][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[3] * src;
      }
      else {
        numint::array < double, 3 > lower = { {0., -1., 0.} };
        numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), 1., 2. * PI} };
        
        GlauberGridThick::Ftor_all F;
        F.grid = this;
        F.proton = proton;
        numint::mdfunction < numint::vector_z, 3 > mdf;
        mdf.func = &Ftor_all::exec;
        mdf.param = &F;
        vector < complex < double > >ret(4, 0.);
        F.f = klaas_int_all;
        if (integrator == 1)
          res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
        else
          res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E03,2E06, ret, count, 0);
        fsi_grid[0][proton][i][j][k] = ret[0];
        fsi_grid[1][proton][i][j][k] = ret[1] * src;
        fsi_ct_grid[0][proton][i][j][k] = ret[2];
        fsi_ct_grid[1][proton][i][j][k] = ret[3] * src;
      }
    }
    else {
      throw std::invalid_argument(nucqcd::ErrorMessage() << "integrator type not implemented" );
    }
//     cout << i << " " << j << " " << k << " " << proton << " "<< fsi_grid[0][proton][i][j][k] << " " <<
//     fsi_grid[1][proton][i][j][k] << " " << fsi_ct_grid[0][proton][i][j][k] << " " <<
//     fsi_ct_grid[1][proton][i][j][k] << 
//  " " << res << " " << count << " " << deeserror << endl;
//     
  }
  //delete[] results;
}

//calc glauberphases for one gridpoint,only CT grid
void GlauberGridThick::calcGlauberphasesCt(const int i, const int j, const int k) {
  for (int proton = 0; proton < 2; proton++) {
    int res = 90;
    unsigned count = 0;
    double deeserror = 0.;
    double src = getFsiCorrelator().getCorrGrid_interp(r_hit, costheta_hit, proton);
    if (integrator == 1 || integrator == 2) {
      if (getParticles().size() == 1 && getParticles()[0].getCosTheta() == 1.) {
        numint::array < double, 3 > lower = { {1.E-06, getParticles()[0].getHitz(), 0.} };
        numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), getPnucleusthick()->getRange(), 2. * PI} };
        
        GlauberGridThick::Ftor_bound F;
        F.grid = this;
        F.proton = proton;
        numint::mdfunction < numint::vector_d, 3 > mdf;
        mdf.func = &Ftor_bound::exec;
        mdf.param = &F;
        vector < double >ret(2, 0.);
        F.f = GlauberGridThick::klaas_int_bound_ct;
        if (integrator == 1)
          res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
        else
          res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E03,5E06, ret, count, 0);
        fsi_ct_grid[0][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[0];
        fsi_ct_grid[1][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[1] * src;
      }
      else {
        numint::array < double, 3 > lower = { {0., -1., 0.} };
        numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), 1., 2. * PI} };
        
        GlauberGridThick::Ftor_all F;
        F.grid = this;
        F.proton = proton;
        numint::mdfunction < numint::vector_z, 3 > mdf;
        mdf.func = &Ftor_all::exec;
        mdf.param = &F;
        vector < complex < double > >ret(2, 0.);
        F.f = GlauberGridThick::klaas_int_all_ct;
        if (integrator == 1)
          res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
        else
          res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E03,5E06, ret, count, 0);
        fsi_ct_grid[0][proton][i][j][k] = ret[0];
        fsi_ct_grid[1][proton][i][j][k] = ret[1] * src;
      }
    }
    else {
      throw std::invalid_argument(nucqcd::ErrorMessage() << "integrator type not implemented" );
    }
//     cout << i << " " << j << " " << k << " " << proton << " "<< fsi_grid[0][proton][i][j][k] << " " <<
//     fsi_grid[1][proton][i][j][k] << " " << fsi_ct_grid[0][proton][i][j][k] << " " <<
//     fsi_ct_grid[1][proton][i][j][k] << 
//  " " << res << " " << count << " " << deeserror << endl;
  }
  //delete[] results;
}

//calc glauberphases for one gridpoint
void GlauberGridThick::calcGlauberphasesBoth(std::vector < std::complex < double > >&protonphases, std::vector < std::complex < double > >&neutronphases) {
  protonphases = vector < complex < double > >(4, 0.);
  neutronphases = vector < complex < double > >(4, 0.);
  for (int proton = 0; proton < 2; proton++) {
    int res = 90;
    unsigned count = 0;
    double deeserror = 0.;
    double src = getFsiCorrelator().getCorrGrid_interp(r_hit, costheta_hit, proton);
    if (integrator == 1 || integrator == 2) {
      if (getParticles().size() == 1 && getParticles()[0].getCosTheta() == 1.) {
        numint::array < double, 3 > lower = { {1.E-06, getParticles()[0].getHitz(), 0.} };
        numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), getPnucleusthick()->getRange(), 2. * PI} };
        
        GlauberGridThick::Ftor_bound F;
        F.grid = this;
        F.proton = proton;
        numint::mdfunction < numint::vector_d, 3 > mdf;
        mdf.func = &Ftor_bound::exec;
        mdf.param = &F;
        vector < double >ret(4, 0.);
        F.f = klaas_int_bound;
        if (integrator == 1)
          res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
        else
          res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E03,2E06, ret, count, 0);
        (proton ? protonphases[0] : neutronphases[0]) = 1. - getParticles()[0].getScatterfront(proton) * ret[0];
        (proton ? protonphases[1] : neutronphases[1]) = 1. - getParticles()[0].getScatterfront(proton) * ret[1] * src;
        (proton ? protonphases[2] : neutronphases[2]) = 1. - getParticles()[0].getScatterfront(proton) * ret[2];
        (proton ? protonphases[3] : neutronphases[3]) = 1. - getParticles()[0].getScatterfront(proton) * ret[3] * src;
      }
      else {
        numint::array < double, 3 > lower = { {0., -1., 0.} };
        numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), 1., 2. * PI} };
        
        GlauberGridThick::Ftor_all F;
        F.grid = this;
        F.proton = proton;
        numint::mdfunction < numint::vector_z, 3 > mdf;
        mdf.func = &Ftor_all::exec;
        mdf.param = &F;
        vector < complex < double > >ret(4, 0.);
        F.f = klaas_int_all;
        if (integrator == 1)
          res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
        else
          res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E03,2E06, ret, count, 0);
        (proton ? protonphases[0] : neutronphases[0]) = ret[0];
        (proton ? protonphases[1] : neutronphases[1]) = ret[1] * src;
        (proton ? protonphases[2] : neutronphases[2]) = ret[2];
        (proton ? protonphases[3] : neutronphases[3]) = ret[3] * src;
      }
    }
    else {
      throw std::invalid_argument(nucqcd::ErrorMessage() << "integrator type not implemented" );
    }
//     cout << i << " " << j << " " << k << " " << proton << " "<< fsi_grid[0][proton][i][j][k] << " " <<
//     fsi_grid[1][proton][i][j][k] << " " << fsi_ct_grid[0][proton][i][j][k] << " " <<
//     fsi_ct_grid[1][proton][i][j][k] << 
//  " " << res << " " << count << " " << deeserror << endl;
  }
  //delete[] results;
}

//calc glauberphases for one gridpoint,only CT grid
void GlauberGridThick::calcGlauberphasesCt(std::vector < std::complex < double > >&protonphases, std::vector < std::complex < double > >&neutronphases) {
  protonphases = vector < complex < double > >(2, 0.);
  neutronphases = vector < complex < double > >(2, 0.);
  for (int proton = 0; proton < 2; proton++) {
    int res = 90;
    unsigned count = 0;
    double deeserror = 0.;
    double src = getFsiCorrelator().getCorrGrid_interp(r_hit, costheta_hit, proton);
    if (integrator == 1 || integrator == 2) {
      if (getParticles().size() == 1 && getParticles()[0].getCosTheta() == 1.) {
        numint::array < double, 3 > lower = { {1.E-06, getParticles()[0].getHitz(), 0.} };
        numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), getPnucleusthick()->getRange(), 2. * PI} };
        
        GlauberGridThick::Ftor_bound F;
        F.grid = this;
        F.proton = proton;
        numint::mdfunction < numint::vector_d, 3 > mdf;
        mdf.func = &Ftor_bound::exec;
        mdf.param = &F;
        vector < double >ret(2, 0.);
        F.f = GlauberGridThick::klaas_int_bound_ct;
        if (integrator == 1)
          res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
        else
          res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E03,5E06, ret, count, 0);
        (proton ? protonphases[0] : neutronphases[0]) = 1. - getParticles()[0].getScatterfront(proton) * ret[0];
        (proton ? protonphases[1] : neutronphases[1]) = 1. - getParticles()[0].getScatterfront(proton) * ret[1] * src;
      }
      else {
        numint::array < double, 3 > lower = { {0., -1., 0.} };
        numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), 1., 2. * PI} };
        
        GlauberGridThick::Ftor_all F;
        F.grid = this;
        F.proton = proton;
        numint::mdfunction < numint::vector_z, 3 > mdf;
        mdf.func = &Ftor_all::exec;
        mdf.param = &F;
        vector < complex < double > >ret(2, 0.);
        F.f = GlauberGridThick::klaas_int_all_ct;
        if (integrator == 1)
          res = numint::cube_romb(mdf, lower, upper, 1.E-12, prec, ret, count, 0);
        else
          res = numint::cube_adaptive(mdf, lower, upper, 1.E-12, prec, 2E03,5E06, ret, count, 0);
        (proton ? protonphases[0] : neutronphases[0]) = ret[0];
        (proton ? protonphases[1] : neutronphases[1]) = ret[1] * src;
      }
    }
    else {
      throw std::invalid_argument(nucqcd::ErrorMessage() << "integrator type not implemented" );
    }
//     cout << i << " " << j << " " << k << " " << proton << " "<< fsi_grid[0][proton][i][j][k] << " " <<
//     fsi_grid[1][proton][i][j][k] << " " << fsi_ct_grid[0][proton][i][j][k] << " " <<
//     fsi_ct_grid[1][proton][i][j][k] << 
//  " " << res << " " << count << " " << deeserror << endl;
  }
  //delete[] results;
}





// 
// 
//   }
//   double gr = getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi, r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit));
//   results[1] = results[0] * gr;
//   results[3] = results[2] * gr;
// }
// 
// 
// 
//   }
//   double gr = getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi, r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit));
//   results[1] = results[0] * gr;
// }

void GlauberGridThick::klaas_int_all(numint::vector_z & results, double r, double costheta, double phi, GlauberGridThick & grid, int proton) {
  double sintheta = sqrt(1. - costheta * costheta);
  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);
  results = numint::vector_z(4, 0.);
  results[0] = 1.;
  results[2] = 1.;
  
  for (size_t it = 0; it < grid.getParticles().size(); it++) {
    double zmom = grid.getParticles()[it].calcZ(r, costheta, sintheta, cosphi, sinphi);
    bool check = (zmom > grid.getParticles()[it].getHitz());
    if (grid.getParticles()[it].getIncoming())
      check = !check;
    if (check) {
      complex < double >temp = grid.getParticles()[it].getScatterfront(proton)
        * exp(-grid.getParticles()[it].getBdist(r, costheta, sintheta, cosphi, sinphi, zmom)
              / (2. * grid.getParticles()[it].getBetasq(proton)));
      results[0] *= 1. - temp;
      results[2] *= 1. - temp * grid.getParticles()[it].getCTsigma(zmom);
    }
  }
  double src = grid.getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi,
                                                         grid.getR_hit(), grid.getCostheta_hit(), grid.getSintheta_hit(),
                                                         grid.getCosphi_hit(), grid.getSinphi_hit())) * grid.getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
  results[1] = results[0] * src;
  results[3] = results[2] * src;
  double dens = grid.getPnucleusthick()->getDensity(r, proton);
  for (int i = 0; i < 4; i++) {
    results[i] *= dens;
  }
  return;
}

void GlauberGridThick::klaas_int_bound(numint::vector_d & results, double b, double z, double phi, GlauberGridThick & grid, int proton) {
  //double sintheta=sqrt(1.-costheta*costheta);
  results = numint::vector_d(4, 0.);
  double r = sqrt(b * b + z * z);
  if (r > grid.getPnucleus()->getRange())
    return;
  double costheta = z / r;
  double sintheta = b / r;
  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);
  
  results[0] = grid.getPnucleusthick()->getDensity(r, proton) / (r * r) * b *
    exp(-(b * b + grid.getParticles()[0].getHitbnorm() * grid.getParticles()[0].getHitbnorm() -
          2. * b * grid.getParticles()[0].getHitbnorm() * cosphi) / (2. * grid.getParticles()[0].getBetasq(proton)));
  double src = grid.getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi,
                                                         grid.getR_hit(), grid.getCostheta_hit(), grid.getSintheta_hit(),
                                                         grid.getCosphi_hit(), grid.getSinphi_hit())) * grid.getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
  results[2] = results[0] * grid.getParticles()[0].getCTsigma(z);
  results[1] = results[0] * src;
  results[3] = results[2] * src;
  return;
}

void GlauberGridThick::klaas_int_all_ct(numint::vector_z & results, double r, double costheta, double phi, GlauberGridThick & grid, int proton) {
  double sintheta = sqrt(1. - costheta * costheta);
  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);
  results = numint::vector_z(2, 0.);
  results[0] = 1.;
  
  for (size_t it = 0; it < grid.getParticles().size(); it++) {
    double zmom = grid.getParticles()[it].calcZ(r, costheta, sintheta, cosphi, sinphi);
    bool check = (zmom > grid.getParticles()[it].getHitz());
    if (grid.getParticles()[it].getIncoming())
      check = !check;
    if (check) {
      complex < double >temp = grid.getParticles()[it].getScatterfront(proton)
        * exp(-grid.getParticles()[it].getBdist(r, costheta, sintheta, cosphi, sinphi, zmom)
              / (2. * grid.getParticles()[it].getBetasq(proton)));
      results[0] *= 1. - temp * grid.getParticles()[it].getCTsigma(zmom);
    }
  }
  double src = grid.getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi,
                                                         grid.getR_hit(), grid.getCostheta_hit(), grid.getSintheta_hit(),
                                                         grid.getCosphi_hit(), grid.getSinphi_hit())) * grid.getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
  results[1] = results[0] * src;
  double dens = grid.getPnucleusthick()->getDensity(r, proton);
  for (int i = 0; i < 2; i++) {
    results[i] *= dens;
  }
  return;
}

void GlauberGridThick::klaas_int_bound_ct(numint::vector_d & results, double b, double z, double phi, GlauberGridThick & grid, int proton) {
  //double sintheta=sqrt(1.-costheta*costheta);
  results = numint::vector_d(2, 0.);
  double r = sqrt(b * b + z * z);
  if (r > grid.getPnucleus()->getRange())
    return;
  double costheta = z / r;
  double sintheta = b / r;
  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);
  
  results[0] = grid.getPnucleusthick()->getDensity(r, proton) / (r * r) * b *
    exp(-(b * b + grid.getParticles()[0].getHitbnorm() * grid.getParticles()[0].getHitbnorm() -
          2. * b * grid.getParticles()[0].getHitbnorm() * cosphi) / (2. * grid.getParticles()[0].getBetasq(proton)))
    * grid.getParticles()[0].getCTsigma(z);
  double src = grid.getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi,
                                                         grid.getR_hit(), grid.getCostheta_hit(), grid.getSintheta_hit(),
                                                         grid.getCosphi_hit(), grid.getSinphi_hit())) * grid.getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
  results[1] = results[0] * src;
  return;
}

//   for (int i = 0; i < 4; i++)
//     results[i] *= b * exp(-(b * b + getParticles()[0].getHitbnorm() * getParticles()[0].getHitbnorm())
//                           / (2. * getParticles()[0].getBetasq(proton)));
// }
// 
//   double costheta = z / r;
//   double sintheta = b / r;
//   double dens = getPnucleusthick()->getDensity(r, proton) / (r * r);
//   double src = getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
//   double intresults[2] = { 0., 0. };
//   //cout << r << " " << acos(costheta)*nucqcd::rad_to_deg << " " << getFsiCorrelator()->getRindex() << endl;
//   results[0] = intresults[0];
//   results[1] = intresults[1] * src;
//   results[2] = results[0] * getParticles()[0].getCTsigma(z);
//   results[3] = results[1] * getParticles()[0].getCTsigma(z);
//   for (int i = 0; i < 4; i++)
//     results[i] *= dens;
// }
// 
// 
//   for (int i = 0; i < 2; i++)
//     results[i] *= b * exp(-(b * b + getParticles()[0].getHitbnorm() * getParticles()[0].getHitbnorm())
//                           / (2. * getParticles()[0].getBetasq(proton)));
// }
// 
//   double costheta = z / r;
//   double sintheta = b / r;
//   double dens = getPnucleusthick()->getDensity(r, proton);
//   double src = getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
//   double intresults;
//   //cout << r << " " << acos(costheta)*nucqcd::rad_to_deg << " " << getFsiCorrelator()->getRindex() << endl;
//   results[0] = intresults * getParticles()[0].getCTsigma(z);
//   results[1] = intresults * getParticles()[0].getCTsigma(z) * src;
//   for (int i = 0; i < 2; i++)
//     results[i] *= dens / (r * r);
// }
// 

}  // namespace nucqcd

