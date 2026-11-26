#include <fstream>
#include <Tools/Core/ErrorMessage.hpp>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <Glauber/FsiCorrelator.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/GridCache.hpp>
#include <Tools/Core/AuxFunction.hpp>

namespace nucqcd {

using namespace std;


//constructor
FsiCorrelator::FsiCorrelator(MeanFieldNucleusThick *inputnucleus, const int rgrid, const int cthgrid, const string dir):
corr_rgrid(rgrid),corr_cthgrid(cthgrid),corr_phigrid(6),pnucleus(inputnucleus)
{
  //cout << "Initializing Object for correlated glauber" << endl;
  //initialization
  invrstep=corr_rgrid/pnucleus->getRange();
  invcthstep=corr_cthgrid/2.;
  setCorrFilename(dir);
  corrgridfull.assign(getCorr_rgrid(),std::vector<double>(getCorr_cthgrid(),0.));
  corrgridproton.assign(getCorr_rgrid(),std::vector<double>(getCorr_cthgrid(),0.));
  corrgridneutron.assign(getCorr_rgrid(),std::vector<double>(getCorr_cthgrid(),0.));
  ifstream infile(corrfilename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  bool read=false;
  if(infile.is_open()){
    if(infile.peek() == std::ifstream::traits_type::eof()){
      cerr << "File is empty " << corrfilename << endl;
      assert(1==0);
    }
    //cout << "Reading in correlation grids from memory: " << corrfilename << endl;
    for(int i=0;i<corr_rgrid;i++){
      for(int j=0;j<corr_cthgrid;j++){
	infile.read(reinterpret_cast<char *>(&corrgridfull[i][j]),sizeof(double));
	infile.read(reinterpret_cast<char *>(&corrgridproton[i][j]),sizeof(double));
	infile.read(reinterpret_cast<char *>(&corrgridneutron[i][j]),sizeof(double));
      }
    }
    read=!infile.fail(); // a truncated file sets failbit
    infile.close();
    if(!read) cerr << "Correlation grid file is incomplete, recomputing it: " << corrfilename << endl;
    //printCorrGridAll();
  }
  //else construct it now!
  if(!read){
    cout << "Constructing grids" << endl;
    constructCorrgrid(pnucleus->getTotalDensity(),corrgridfull);
    constructCorrgrid(pnucleus->getProtonDensity(),corrgridproton);
    constructCorrgrid(pnucleus->getNeutronDensity(),corrgridneutron);
    //printCorrGridAll();
    bool written=nucqcd::write_file_atomically(corrfilename,[this](ofstream &outfile){
      for(int i=0;i<corr_rgrid;i++){
	for(int j=0;j<corr_cthgrid;j++){
	  outfile.write(reinterpret_cast<char *>(&corrgridfull[i][j]),sizeof(double));
	  outfile.write(reinterpret_cast<char *>(&corrgridproton[i][j]),sizeof(double));
	  outfile.write(reinterpret_cast<char *>(&corrgridneutron[i][j]),sizeof(double));
	}
      }
    });
    if(written) cout << "Wrote grids to file :" << corrfilename << endl;
  }

}


//get corrfilename
const string FsiCorrelator::getCorrFilename() const{
  return corrfilename;
}

//get corr_rgrid
int FsiCorrelator::getCorr_rgrid() const{
  return corr_rgrid;
}

//get corr_thgrid
int FsiCorrelator::getCorr_cthgrid() const{
  return corr_cthgrid;
}

//get rstep
double FsiCorrelator::getInvRstep() const{
  return invrstep;
}

//get thstep
double FsiCorrelator::getInvCthstep() const{
  return invcthstep;
}

//get s_interp
double FsiCorrelator::getS_interp() const{
  return s_interp;
}

//get t_interp
double FsiCorrelator::getT_interp() const{
  return t_interp;
}

//get comp_s_interp
double FsiCorrelator::getComp_s_interp() const{
  return comp_s_interp;
}

//get comp_t_interp
double FsiCorrelator::getComp_t_interp() const{
  return comp_t_interp;
}

//get rindex
int FsiCorrelator::getRindex() const{
  return rindex;
}

//get thindex
int FsiCorrelator::getCthindex() const{
  return cthindex;
}

//set up inteprolation help variables for r
void FsiCorrelator::setRinterp(const double r){
  if(r>pnucleus->getRange()){
    throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range in getWave_Rpart: " << r );
  }
  rindex = int(floor(r*getInvRstep()));
  if(rindex==getCorr_rgrid()) rindex-=2;
  if(rindex==(getCorr_rgrid()-1)) rindex-=1;
  s_interp = (r*getInvRstep() - (rindex));
  comp_s_interp=1.-s_interp;
}

//set up interpolation help variables for theta
void FsiCorrelator::setCthinterp(double costheta){
  //if(costheta<0.) costheta=-costheta;
  //if(costheta>1.) costheta=1.;
  cthindex = int(floor((1.-costheta)*getInvCthstep()));
  if(cthindex==getCorr_cthgrid()) cthindex-=2;
  if(cthindex==(getCorr_cthgrid()-1)) cthindex-=1;
  t_interp = ((1.-costheta)*getInvCthstep() - (cthindex));  
  comp_t_interp=1.-t_interp;
}

//set r, theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridFull_interp(const double r, const double costheta){
  setRinterp(r); 
  setCthinterp(costheta);
  return getCorrGridFull_interp();
}

//set r, theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridProton_interp(const double r, const double costheta){
  setRinterp(r); 
  setCthinterp(costheta);
  return getCorrGridProton_interp();
}


//set r, theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridNeutron_interp(const double r, const double costheta){
  setRinterp(r); 
  setCthinterp(costheta);
  return getCorrGridNeutron_interp();
}

//set theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridFull_interp(const double costheta){
  setCthinterp(costheta);
  return getCorrGridFull_interp();
}

//set theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridProton_interp(const double costheta){
  setCthinterp(costheta);
  return getCorrGridProton_interp();
}


//set theta and then get interpolation value of the grid
double FsiCorrelator::getCorrGridNeutron_interp(const double costheta){
  setCthinterp(costheta);
  return getCorrGridNeutron_interp();
}

double FsiCorrelator::getCorrGrid_interp(const double r, const double costheta, bool proton){
 return proton? getCorrGridProton_interp(r,costheta) : getCorrGridNeutron_interp(r,costheta);
}

double FsiCorrelator::getCorrGrid_interp(const double costheta, bool proton){
 return proton? getCorrGridProton_interp(costheta) : getCorrGridNeutron_interp(costheta);
}

// get interpolation value of the grid
double FsiCorrelator::getCorrGridFull_interp() const{
  return comp_t_interp*(comp_s_interp*corrgridfull[rindex][cthindex]
	    + s_interp*corrgridfull[rindex+1][cthindex])
	    + t_interp*(comp_s_interp*corrgridfull[rindex][cthindex+1]
	    + s_interp*corrgridfull[rindex+1][cthindex+1]);  
}


// get interpolation value of the grid
double FsiCorrelator::getCorrGridProton_interp() const{
  return comp_t_interp*(comp_s_interp*corrgridproton[rindex][cthindex]
	    + s_interp*corrgridproton[rindex+1][cthindex])
	    + t_interp*(comp_s_interp*corrgridproton[rindex][cthindex+1]
	    + s_interp*corrgridproton[rindex+1][cthindex+1]);  
}


// get interpolation value of the grid
double FsiCorrelator::getCorrGridNeutron_interp() const{
//  cout << rindex << " " << thindex << " " << s_interp << " " << t_interp << endl;
  return comp_t_interp*(comp_s_interp*corrgridneutron[rindex][cthindex]
	    + s_interp*corrgridneutron[rindex+1][cthindex])
	    + t_interp*(comp_s_interp*corrgridneutron[rindex][cthindex+1]
	    + s_interp*corrgridneutron[rindex+1][cthindex+1]);  
}
  
double FsiCorrelator::getCorrGrid_interp(bool proton) const{
 return proton? getCorrGridProton_interp() : getCorrGridNeutron_interp();
}
  
  
  
  
//set filename where the class object is stored or written to
void FsiCorrelator::setCorrFilename(const string dir){
  corrfilename=dir+"/input/CorrGrid."+pnucleus->getNucleusName()+".tot."+to_string(corr_rgrid)+".c"+to_string(corr_cthgrid)+".dat";  
}

void FsiCorrelator::printCorrGridFull() const{
  for(int i=0;i<corr_rgrid;i++){
    for(int j=0;j<corr_cthgrid;j++){
      cout << pnucleus->getRange()*i/corr_rgrid << " " << 1.-j/invcthstep << " " << corrgridfull[i][j] << endl;
    }
  }
  cout << endl << endl;
}
  
void FsiCorrelator::printCorrGridProton() const{
  for(int i=0;i<corr_rgrid;i++){
    for(int j=0;j<corr_cthgrid;j++){
      cout << pnucleus->getRange()*i/corr_rgrid << " " << 1.-j/invcthstep << " " << corrgridproton[i][j] << endl;
    }
  }
  cout << endl << endl;
}
  
void FsiCorrelator::printCorrGridNeutron() const{
  for(int i=0;i<corr_rgrid;i++){
    for(int j=0;j<corr_cthgrid;j++){
      cout << pnucleus->getRange()*i/corr_rgrid << " " << 1.-j/invcthstep << " " << corrgridneutron[i][j] << endl;
    }
  }
  cout << endl << endl;
}
  
void FsiCorrelator::printCorrGridAll() const{
  for(int i=0;i<corr_rgrid;i++){
    for(int j=0;j<corr_cthgrid;j++){
      cout << pnucleus->getRange()*i/corr_rgrid << " " << 1.-j/invcthstep << " " << corrgridfull[i][j]
      << " " << corrgridproton[i][j] << " " << corrgridneutron[i][j] << endl;
    }
  }
  cout << endl << endl;
}
  




}  // namespace nucqcd



