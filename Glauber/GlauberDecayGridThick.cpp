#include <Integrator/numint/cube.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <iostream>
#include <fstream>
#include <complex>
#include <cassert>

using namespace std;

#include <Glauber/GlauberDecayGridThick.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/AuxFunction.hpp>

namespace nucqcd {


GlauberDecayGridThick::GlauberDecayGridThick(const int r_grid, const int cth_grid, const int phi_grid, MeanFieldNucleusThick * pnucl,
                                             double prec, IntegrationMethod integrator, string dir):
GlauberThickGridCommon(4, r_grid, cth_grid, phi_grid, pnucl, prec, integrator, dir) {
}

GlauberDecayGridThick::~GlauberDecayGridThick() {
}



void GlauberDecayGridThick::printFsi_grid() {
  cout << "Printing Glauberarray for " << getFsi_Filename() << endl;
  print_grid(0);
}

void GlauberDecayGridThick::printFsi_src_grid() {
  cout << "Printing SRC Glauberarray for " << getFsi_Filename() << endl;
  print_grid(1);
}

void GlauberDecayGridThick::printFsi_decay_grid() {
  cout << "Printing Decay Glauberarray for " << getFsi_Filename() << endl;
  print_grid(2);
}

void GlauberDecayGridThick::printFsi_src_decay_grid() {
  cout << "Printing SRC Decay Glauberarray for " << getFsi_Filename() << endl;
  print_grid(3);
}

void GlauberDecayGridThick::printFsi_ct_grid() {
  cout << "Printing CT Glauberarray for " << getFsi_Ct_Filename() << endl;
  print_grid(4);
}

void GlauberDecayGridThick::printFsi_src_ct_grid() {
  cout << "Printing SRC CT Glauberarray for " << getFsi_Ct_Filename() << endl;
  print_grid(5);
}

void GlauberDecayGridThick::printFsi_ct_decay_grid() {
  cout << "Printing CT Decay Glauberarray for " << getFsi_Ct_Filename() << endl;
  print_grid(6);
}

void GlauberDecayGridThick::printFsi_src_ct_decay_grid() {
  cout << "Printing SRC CT Decay Glauberarray for " << getFsi_Ct_Filename() << endl;
  print_grid(7);
}

//set filenames, don't forget to add ending "/" to dir if necessary!!!
void GlauberDecayGridThick::setFilenames(string dir) {
  GlauberThickGridCommon::setFilenames(dir + "Gl.");
  // the decay dilation of the (first) ejected particle is part of the file names
  const string decay = ".Decay" + to_string(lround(getParticles()[0].getDecay_dil()));
  fsi_ct_filename.insert(fsi_ct_filename.size() - 4, decay);
  fsi_filename.insert(fsi_filename.size() - 4, decay);
}



//calc glauberphases for one gridpoint
void GlauberDecayGridThick::calcGlauberphasesBoth(const int i, const int j, const int k){
  complex<double>* results= new complex<double>[8];
  for (int proton = 0; proton < 2; proton++) {
    int res = 90;
    unsigned count = 0;
    double deeserror = 0.;
    double src = getFsiCorrelator().getCorrGrid_interp(r_hit, costheta_hit, proton);
    if (getParticles().size() == 1 && getParticles()[0].getCosTheta() == 1.) {
      numint::array < double, 3 > lower = { {1.E-06, getParticles()[0].getHitz(), 0.} };
      numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), getPnucleusthick()->getRange(), 2. * PI} };

      GlauberDecayGridThick::Ftor_bound F;
      F.grid = this;
      F.proton = proton;
      numint::mdfunction < numint::vector_d, 3 > mdf;
      mdf.func = &Ftor_bound::exec;
      mdf.param = &F;
      vector < double >ret(8, 0.);
      F.f = klaas_int_bound;
      res = numint::cube({integrator, 1.E-12, prec, 2000, 2000000}, mdf, lower, upper, ret, count, 0);
      fsi_grid[0][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[0];
      fsi_grid[1][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[1] * src;
      fsi_grid[2][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[2];
      fsi_grid[3][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[3] * src;
      fsi_ct_grid[0][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[4];
      fsi_ct_grid[1][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[5] * src;
      fsi_ct_grid[2][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[6];
      fsi_ct_grid[3][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[7] * src;
    }
    else if (getParticles().size() == 1 && useTiltedFrame) {
      vector < double >ret;
      integrateTilted(proton, klaas_int_tilted, 8, 2000000, ret);
      for (int v = 0; v < 4; v++) {
        fsi_grid[v][proton][i][j][k] = tiltedPhase(proton, ret, v, src);
        fsi_ct_grid[v][proton][i][j][k] = tiltedPhase(proton, ret, 4 + v, src);
      }
    }
    else {
      numint::array < double, 3 > lower = { {0., -1., 0.} };
      numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), 1., 2. * PI} };
      GlauberDecayGridThick::Ftor_all F;
      F.grid = this;
      F.proton = proton;
      numint::mdfunction < numint::vector_z, 3 > mdf;
      mdf.func = &Ftor_all::exec;
      mdf.param = &F;
      vector < complex < double > >ret(8, 0.);
      F.f = klaas_int_all;
      res = numint::cube({integrator, 1.E-12, prec, 2000, 2000000}, mdf, lower, upper, ret, count, 0);
      fsi_grid[0][proton][i][j][k] = ret[0];
      fsi_grid[1][proton][i][j][k] = ret[1] * src;
      fsi_grid[2][proton][i][j][k] = ret[2];
      fsi_grid[3][proton][i][j][k] = ret[3] * src;
      fsi_ct_grid[0][proton][i][j][k] = ret[4];
      fsi_ct_grid[1][proton][i][j][k] = ret[5] * src;
      fsi_ct_grid[2][proton][i][j][k] = ret[6];
      fsi_ct_grid[3][proton][i][j][k] = ret[7] * src;
    }
//     cout << i << " " << j << " " << k << " " << proton << " "<< fsi_grid[0][proton][i][j][k] << " " <<
//     fsi_grid[1][proton][i][j][k] << " " << fsi_ct_grid[0][proton][i][j][k] << " " <<
//     fsi_ct_grid[1][proton][i][j][k] << 
//   " " << res << " " << count << " " << deeserror << endl;
  }
}

//calc glauberphases for one gridpoint,only CT grid
void GlauberDecayGridThick::calcGlauberphasesCt(const int i, const int j, const int k) {
  for (int proton = 0; proton < 2; proton++) {
    int res = 90;
    unsigned count = 0;
    double deeserror = 0.;
    double src = getFsiCorrelator().getCorrGrid_interp(r_hit, costheta_hit, proton);
    if (getParticles().size() == 1 && getParticles()[0].getCosTheta() == 1.) {
      numint::array < double, 3 > lower = { {1.E-06, getParticles()[0].getHitz(), 0.} };
      numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), getPnucleusthick()->getRange(), 2. * PI} };

      GlauberDecayGridThick::Ftor_bound F;
      F.grid = this;
      F.proton = proton;
      numint::mdfunction < numint::vector_d, 3 > mdf;
      mdf.func = &Ftor_bound::exec;
      mdf.param = &F;
      vector < double >ret(4, 0.);
      F.f = GlauberDecayGridThick::klaas_int_bound_ct;
      res = numint::cube({integrator, 1.E-12, prec, 2000, 5000000}, mdf, lower, upper, ret, count, 0);
      fsi_ct_grid[0][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[0];
      fsi_ct_grid[1][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[1] * src;
      fsi_ct_grid[2][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[2];
      fsi_ct_grid[3][proton][i][j][k] = 1. - getParticles()[0].getScatterfront(proton) * ret[3] * src;
    }
    else if (getParticles().size() == 1 && useTiltedFrame) {
      vector < double >ret;
      integrateTilted(proton, klaas_int_tilted_ct, 4, 5000000, ret);
      for (int v = 0; v < 4; v++)
        fsi_ct_grid[v][proton][i][j][k] = tiltedPhase(proton, ret, v, src);
    }
    else {
      numint::array < double, 3 > lower = { {0., -1., 0.} };
      numint::array < double, 3 > upper = { {getPnucleusthick()->getRange(), 1., 2. * PI} };

      GlauberDecayGridThick::Ftor_all F;
      F.grid = this;
      F.proton = proton;
      numint::mdfunction < numint::vector_z, 3 > mdf;
      mdf.func = &Ftor_all::exec;
      mdf.param = &F;
      vector < complex < double > >ret(4, 0.);
      F.f = GlauberDecayGridThick::klaas_int_all_ct;
      res = numint::cube({integrator, 1.E-12, prec, 2000, 5000000}, mdf, lower, upper, ret, count, 0);
      fsi_ct_grid[0][proton][i][j][k] = ret[0];
      fsi_ct_grid[1][proton][i][j][k] = ret[1] * src;
      fsi_ct_grid[2][proton][i][j][k] = ret[2];
      fsi_ct_grid[3][proton][i][j][k] = ret[3] * src;
    }
//     cout << i << " " << j << " " << k << " " << proton << " "<< fsi_grid[0][proton][i][j][k] << " " <<
//     fsi_grid[1][proton][i][j][k] << " " << fsi_ct_grid[0][proton][i][j][k] << " " <<
//     fsi_ct_grid[1][proton][i][j][k] << 
//   " " << res << " " << count << " " << deeserror << endl;
  }
}





// 
// 
//   }
//   double gr = getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi, r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit));
//   results[1] = results[0] * gr;
//   results[3] = results[2] * gr;
//   results[5] = results[4] * gr;
//   results[7] = results[6] * gr;
// }
// 
// 
// 
//   }
//   double gr = getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi, r_hit, costheta_hit, sintheta_hit, cosphi_hit, sinphi_hit));
//   results[1] = results[0] * gr;
//   results[3] = results[2] * gr;
// }

void GlauberDecayGridThick::klaas_int_tilted(numint::vector_d & results, double b, double z, double phi, GlauberThickGridCommon & grid,
                                             int proton, const TiltedFrame &frame) {
  results = numint::vector_d(8, 0.);
  TiltedPoint p;
  if (!tiltedPoint(grid, b, z, phi, proton, frame, p)) return;
  FastParticle &part = grid.getParticles()[0];
  results[0] = p.fsi;
  results[2] = p.fsi * part.getDecay_sigma(z, proton);
  results[4] = p.fsi * part.getCTsigma(z);
  results[6] = results[2] * part.getCTsigma(z);
  for (int v = 1; v < 8; v += 2) results[v] = results[v - 1] * p.src;
}

void GlauberDecayGridThick::klaas_int_tilted_ct(numint::vector_d & results, double b, double z, double phi, GlauberThickGridCommon & grid,
                                                int proton, const TiltedFrame &frame) {
  results = numint::vector_d(4, 0.);
  TiltedPoint p;
  if (!tiltedPoint(grid, b, z, phi, proton, frame, p)) return;
  FastParticle &part = grid.getParticles()[0];
  results[0] = p.fsi * part.getCTsigma(z);
  results[2] = results[0] * part.getDecay_sigma(z, proton);
  results[1] = results[0] * p.src;
  results[3] = results[2] * p.src;
}

void GlauberDecayGridThick::klaas_int_all(numint::vector_z & results, double r, double costheta, double phi, GlauberDecayGridThick & grid, int proton) {
  double sintheta = sqrt(1. - costheta * costheta);
  double cosphi, sinphi;
  sincos(phi, &sinphi, &cosphi);
  results = numint::vector_z(8, 0.);
  results[0] = 1.;
  results[2] = 1.;
  results[4] = 1.;
  results[6] = 1.;

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
      results[2] *= 1. - temp * grid.getParticles()[it].getDecay_sigma(zmom, proton);
      results[4] *= 1. - temp * grid.getParticles()[it].getCTsigma(zmom);
      results[6] *= 1. - temp * grid.getParticles()[it].getCTsigma(zmom) * grid.getParticles()[it].getDecay_sigma(zmom, proton);
    }
  }
  double src = grid.getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi,
                                                         grid.getR_hit(), grid.getCostheta_hit(), grid.getSintheta_hit(),
                                                         grid.getCosphi_hit(), grid.getSinphi_hit())) * grid.getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
  results[1] = results[0] * src;
  results[3] = results[2] * src;
  results[5] = results[4] * src;
  results[7] = results[6] * src;
  double dens = grid.getPnucleusthick()->getDensity(r, proton);
  for (int i = 0; i < 8; i++) {
    results[i] *= dens;
  }
  return;

}
void GlauberDecayGridThick::klaas_int_bound(numint::vector_d & results, double b, double z, double phi, GlauberDecayGridThick & grid, int proton) {
  //double sintheta=sqrt(1.-costheta*costheta);
  results = numint::vector_d(8, 0.);
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
  results[2] = results[0] * grid.getParticles()[0].getDecay_sigma(z, proton);
  results[4] = results[0] * grid.getParticles()[0].getCTsigma(z);
  results[6] = results[0] * grid.getParticles()[0].getDecay_sigma(z, proton) * grid.getParticles()[0].getCTsigma(z);
  results[1] = results[0] * src;
  results[3] = results[2] * src;
  results[5] = results[4] * src;
  results[7] = results[6] * src;
  return;
}

void GlauberDecayGridThick::klaas_int_all_ct(numint::vector_z & results, double r, double costheta, double phi, GlauberDecayGridThick & grid, int proton) {
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
      results[0] *= 1. - temp * grid.getParticles()[it].getCTsigma(zmom);
      results[2] *= 1. - temp * grid.getParticles()[it].getCTsigma(zmom) * grid.getParticles()[it].getDecay_sigma(zmom, proton);
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

void GlauberDecayGridThick::klaas_int_bound_ct(numint::vector_d & results, double b, double z, double phi, GlauberDecayGridThick & grid, int proton) {
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
          2. * b * grid.getParticles()[0].getHitbnorm() * cosphi) / (2. * grid.getParticles()[0].getBetasq(proton)))
    * grid.getParticles()[0].getCTsigma(z);
  results[2] = results[0] * grid.getParticles()[0].getDecay_sigma(z, proton);
  double src = grid.getFsiCorrelator().correlation(normr(r, costheta, sintheta, cosphi, sinphi,
                                                         grid.getR_hit(), grid.getCostheta_hit(), grid.getSintheta_hit(),
                                                         grid.getCosphi_hit(), grid.getSinphi_hit())) * grid.getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
  results[1] = results[0] * src;
  results[3] = results[2] * src;
  return;
}

//   for (int i = 0; i < 8; i++)
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
//   double ct = getParticles()[0].getCTsigma(z);
//   double decay = getParticles()[0].getDecay_sigma(z, proton);
//   results[0] = intresults[0];
//   results[1] = intresults[1] * src;
//   results[2] = results[0] * decay;
//   results[3] = results[1] * decay;
//   results[4] = results[0] * ct;
//   results[5] = results[1] * ct;
//   results[6] = results[2] * ct;
//   results[7] = results[3] * ct;
//   for (int i = 0; i < 8; i++)
//     results[i] *= dens;
// }
// 
// 
//   for (int i = 0; i < 4; i++)
//     results[i] *= b * exp(-(b * b + getParticles()[0].getHitbnorm() * getParticles()[0].getHitbnorm())
//                           / (2. * getParticles()[0].getBetasq(proton)));
// }
// 
//   double costheta = z / r;
//   double sintheta = b / r;
//   double dens = getPnucleusthick()->getDensity(r, proton);
//   double src = getFsiCorrelator().getCorrGrid_interp(r, costheta, proton);
//   double intresults = 0.;
//   //cout << r << " " << acos(costheta)*nucqcd::rad_to_deg << " " << getFsiCorrelator()->getRindex() << endl;
//   double ct = getParticles()[0].getCTsigma(z);
//   double decay = getParticles()[0].getDecay_sigma(z, proton);
//   results[0] = intresults * ct;
//   results[1] = results[0] * src;
//   results[2] = results[0] * decay;
//   results[3] = results[1] * decay;
//   for (int i = 0; i < 4; i++)
//     results[i] *= dens / (r * r);
// }
// 

}  // namespace nucqcd

