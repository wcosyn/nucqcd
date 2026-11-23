#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include <Glauber/AbstractFsiCTGridThick.hpp>
#include <Tools/Core/constants.hpp>
#include <Tools/Core/Utilfunctions.hpp>

AbstractFsiCTGridThick::AbstractFsiCTGridThick(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleusThick *pnucl, 
					      double prec, int integrator, string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator, homedir),
AbstractFsiGridThick(r_grid,cth_grid,phi_grid,pnucl,prec, integrator,homedir), 
AbstractFsiCTGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir){
  number_of_grids=4;
}  

AbstractFsiCTGridThick::~AbstractFsiCTGridThick(){
  
}

//interpolation functions
complex<double> AbstractFsiCTGridThick::getFsiSrcCtGridFull_interpvec(const nucqcd::Vector3 &rvec){
  return getFsiSrcCtGridFull_interp3(rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiCTGridThick::getFsiSrcCtGridFull_interp3(const double r, const double costheta, const double phi){
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}

complex<double> AbstractFsiCTGridThick::getFsiSrcCtGridFull_interp2(const double costheta, const double phi){
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}


complex<double> AbstractFsiCTGridThick::getFsiSrcCtGridFull_interp1(const double phi){
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
  return getFsiSrcCtGridFull_interp();
}


void AbstractFsiCTGridThick::setFilenames(string homedir){
  AbstractFsiCTGrid::setFilenames(homedir+"Thick.");
}
  