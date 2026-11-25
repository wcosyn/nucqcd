#include <iostream>
#include <Tools/Core/ErrorMessage.hpp>
#include <cstdlib>


#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Tools/Core/Utilfunctions.hpp>

namespace nucqcd {


using namespace std;

//constructor
MeanFieldNucleusThick::MeanFieldNucleusThick(const int nucleus, const string & dir)
 : MeanFieldNucleus(nucleus,dir) {
   //cout << "Initializing nucleus with thickness..." << endl;
   constructDensity();
   
}


/**
 * Construct the mean densities
 * note that these are normed to 1!
 * so \f$ \int d^{2} \Omega \int \textrm{d} r^{2} \rho(r) = 1 \f$
 */
void MeanFieldNucleusThick::constructDensity(){
  if(!protondensity.empty()||!neutrondensity.empty()||!totaldensity.empty()){
    cerr<< "Warning: densities already initialized, skipping..." << endl;
    return;
  }
  
  double norm=2.*PI;
  protondensity.assign(getWF_r_lines(),0.);
  neutrondensity.assign(getWF_r_lines(),0.);
  totaldensity.assign(getWF_r_lines(),0.);
  
  for(int j=0;j<getWF_r_lines();j++){
    protondensity[j]=0.;
    neutrondensity[j]=0.;
  }
  int counter=0; //keep track of #nucleons
  
  //protondensity
  for(int i=0;i<getPLevels();i++){
    if((counter+2*abs(getKappas()[i]))<=getZ()){
      for(int j=0;j<getWF_r_lines();j++){
	protondensity[j]+=(getF()[i][j]*getF()[i][j]+getG()[i][j]*getG()[i][j])*abs(getKappas()[i]); //2j+1=2*abs(kappa), cancelled a 2 with the 4pi from the denom
      }
    }
    else{
      for(int j=0;j<getWF_r_lines();j++){
	protondensity[j]+=(getF()[i][j]*getF()[i][j]+getG()[i][j]*getG()[i][j])*(getZ()-counter)/2.; //modified for not-full-shells 
      }
    }
    counter+=2*abs(getKappas()[i]);
  }
  for(int j=0;j<getWF_r_lines();j++){
    totaldensity[j]=protondensity[j];
    protondensity[j]/=norm*getZ(); /**< this norms protondensity to one! **/
  }
  if(getZ()==0){
    for(int j=0;j<getWF_r_lines();j++){
      protondensity[j]=0.;
    }
  }

  counter =0;
  //neutron density needed
  for(int i=getPLevels();i<getTotalLevels();i++){
    if((counter+2*abs(getKappas()[i]))<=getN()){
      for(int j=0;j<getWF_r_lines();j++){
	neutrondensity[j]+=(getF()[i][j]*getF()[i][j]+getG()[i][j]*getG()[i][j])*abs(getKappas()[i]); //2j+1=2*abs(kappa)
      }
    }
    else{
      for(int j=0;j<getWF_r_lines();j++){
	neutrondensity[j]+=(getF()[i][j]*getF()[i][j]+getG()[i][j]*getG()[i][j])*(getN()-counter)/2.; //not-full shell	  
      }
    }
    counter+=2*abs(getKappas()[i]);
  }
  for(int j=0;j<getWF_r_lines();j++){
    totaldensity[j]+=neutrondensity[j];
    neutrondensity[j]/=norm*getN(); /**< this norms neutron density to one! **/
    totaldensity[j]/=norm*getA(); /**< this norms total density to one! **/
  }    
  if(getN()==0){
    for(int j=0;j<getWF_r_lines();j++){
      neutrondensity[j]=0.;
    }
  }

  

}


double MeanFieldNucleusThick::getProtonDensity(const double r) const{
  if(r>getRange()){
    throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range in getWave_Rpart: " << r );
  }
  return interpolate(protondensity.data(),r,getWF_r_step(),getWF_r_lines(),1); 
}

double MeanFieldNucleusThick::getNeutronDensity(const double r) const{
  if(r>getRange()){
    throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range in getWave_Rpart: " << r );
  }
  return interpolate(neutrondensity.data(),r,getWF_r_step(),getWF_r_lines(),1); 
}

double MeanFieldNucleusThick::getTotalDensity(const double r) const{
  if(r>getRange()){
    throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range in getWave_Rpart: " << r );
  }
  return interpolate(totaldensity.data(),r,getWF_r_step(),getWF_r_lines(),1); 
}

double MeanFieldNucleusThick::getDensity(const double r, bool proton) const{

  return proton? getProtonDensity(r):getNeutronDensity(r);
}


double MeanFieldNucleusThick::getDensity(const double r, const double *density) const{ //get density for r
  if(r>getRange()){
    throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range in getWave_Rpart: " << r );
  }
  return interpolate(density,r,getWF_r_step(),getWF_r_lines(),1); 
}

}  // namespace nucqcd



