/*
 * FastParticle_ScatteringData.cpp
 * FastParticle: the scattering parameters of the ejected particle from the data
 * (nucleons: setGlauberParameters; pions: setPionGlauberData, interpPionGlauberData).
 * Part of the implementation of FastParticle, split from FastParticle.cpp in 2026-11 (the functions were moved verbatim).
 */

#include <iostream>
#include <Tools/Core/ErrorMessage.hpp>
#include <cstdlib>
#include <fstream>

using namespace std; 

#include <Glauber/FastParticle.hpp>
#include <Tools/Core/constants.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <cassert>
#include <vector>
#include <cmath>

namespace nucqcd {

namespace {
// A table of measured points: three columns of the same length (the logarithm of the momentum, a value and its sigma); table[j] is a
// pointer to column j, so that table[j][k] and the pointer arguments of fitdata work as they did for the raw arrays (which leaked)
struct DataTable {
  std::vector<double> column[3];
  void resize(int points) { for (auto &c : column) c.assign(points, 0.); }
  double *operator[](int j) { return column[j].data(); }
};
}  // namespace


using namespace std;

/////////////////////////////////////////////
//calculate parameters for profile function//
/////////////////////////////////////////////
void FastParticle::setGlauberParameters(double mom, double &sigmapp, double &beta2pp, double &epspp, 
					double &sigmapn, double &beta2pn, double &epspn){

  double p_gev = mom/1000.;
  double sigmappel, sigmapnel;

  // sigmapp

  if (p_gev > 1.4) sigmapp = 48 + 0.522*pow(log(p_gev),2) - 4.51*log(p_gev);
  else {
    if (p_gev > 0.9) sigmapp = 23 + 47.1 * (p_gev - 0.9);
    else {
      if (p_gev > 0.5 ) sigmapp = 25 - 5 * (p_gev - 0.5);
      else  sigmapp = 8.667/p_gev/p_gev - 23.114/p_gev + 36.559;
    }
  }	
  
  sigmapp /= 10;

  // sigmapp elastic

  if (p_gev > 1.7) 
    {
      sigmappel = 11.9 + 26.9 * pow(p_gev,-1.21) + 0.169 * pow(log(p_gev),2) - 
	1.85 * log(p_gev);
    }
  else{
    if (p_gev > 0.9) sigmappel = 23 + 2.65 * (p_gev - 0.9);
    else {
      if (p_gev > 0.5) sigmappel = 25 - 5 * (p_gev - 0.5);
      else  sigmappel = 8.667/p_gev/p_gev - 23.114/p_gev + 36.559;
    }
  }
  
  sigmappel /= 10;

  // sigmapn

  if (p_gev > 2) sigmapn = 40;
  else{
    if (p_gev > 1) sigmapn = 33.7 + 6.3 * (p_gev - 1);
    else{
      if (p_gev > 0.4) sigmapn = 19.7/p_gev/p_gev - 34.083/p_gev + 48.083;
      else  sigmapn = 26.426/p_gev/p_gev - 48.595/p_gev + 42.325;
    }
  }
  sigmapn /= 10;

  // sigmapnel
 
  if (p_gev > 1) sigmapnel = 31 - 10.85 * log(p_gev);
  else{
    if (p_gev > 0.4) sigmapnel = 19.7/p_gev/p_gev - 34.083/p_gev + 48.083;
    else  sigmapnel = 26.426/p_gev/p_gev - 48.595/p_gev + 42.325;
  }
  sigmapnel /= 10;

  // epspp

  if (p_gev > 0.32) epspp = -0.18/p_gev/p_gev + 1.45/p_gev - 0.77;
  else epspp = 0.8 + 6.49 * (p_gev - 0.135);
  
  // epspn

  if (p_gev > 3) epspn = 6.61/p_gev/p_gev - 3.51/p_gev - 0.0655;
  else{
    if (p_gev > 0.50) epspn = 0.27/p_gev/p_gev + 0.39/p_gev - 0.66;
    else  epspn = 0.25 + 1.04 * (log(p_gev) + 1.61);
  }
  // betapp

  beta2pp = (pow(sigmapp,2)*(1+pow(epspp,2)))/(16*PI*sigmappel);
  beta2pn = (pow(sigmapn,2)*(1+pow(epspn,2)))/(16*PI*sigmapnel);
}


void FastParticle::setPionGlauberData(double mom, double &sigmatotp, double &beta2p, double &epsp, 
		     double &sigmatotn, double &beta2n, double &epsn, const string dir){
  double logp=log10(mom/1000.);
  double sigmaelp, sigmaeln;
  double betap2, betan2;
  //read in all the data
  const string files[6]={dir+"/input/pimp_elastic.dat",dir+"/input/pimp_elastic.reim",dir+"/input/pimp_total.dat",
			       dir+"/input/pipp_elastic.dat",dir+"/input/pipp_elastic.reim",dir+"/input/pipp_total.dat"};
  std::vector<int> lines(6);
  std::vector<DataTable> data(6);  //dim1 = 6 files, dim2 = 3 (momentum, y-value, sigma), dim3 = grid
  double dummy,dummy2; //used to read in redundant info etc
  for(int i=0;i<6;i++){
    ifstream file(files[i].c_str(),ios::in);
    
    if(file.is_open()){

      //ignore first 10 lines, no usefull information
      for(int j=0;j<10;j++) file.ignore(500,'\n');

      //read in number of data points
      file >> lines[i];
      data[i].resize(lines[i]);
      
      //read in data points
      for(int k=0;k<lines[i];k++){
	file >> dummy;
	file >> dummy;
	data[i][0][k]=log10(dummy); //momentum
	file >> dummy; file >> dummy;
	file >> dummy2;
	if(i==1||i==4)	data[i][1][k]=dummy2; //value(total,elastic or epsilon)
	else data[i][1][k]=log10(dummy2);
	file >> dummy;
	data[i][2][k]=1.;//dummy/log(10); //sigma for value
	file.ignore(500,'\n');
      }
      file.close();
    }
    
    else{
      throw std::runtime_error(nucqcd::ErrorMessage() << "File " << files[i] << " could not be opened" );
    }
  }
  

  //new data for epsilon! extra input
  const string files2[2]={dir+"/input/epspim.inp",dir+"/input/epspip.inp"};
  
  std::vector<DataTable> epsdata(2);
  int lines2[2]={164,134};
  for(int i=0;i<2;i++){
    
    ifstream file(files2[i].c_str(),ios::in);    
    if(file.is_open()){
      epsdata[i].resize(lines2[i]);
      for(int k=0;k<lines2[i];k++){      
	file >> dummy;
	epsdata[i][0][k] = log10(dummy/1000.);
	file >> epsdata[i][1][k];
	epsdata[i][2][k] = 1.;
      }
      file.close();
    }
    
    else{
      throw std::runtime_error(nucqcd::ErrorMessage() << "File " << files2[i] << " could not be opened" );
    }
  }


  //new data for beta, extra input
  const string files3[2]={dir+"/input/betapim.inp",dir+"/input/betapip.inp"};

  std::vector<DataTable> betadata(2);
  int lines3[2]={130,96};
  for(int i=0;i<2;i++){
    
    ifstream file(files3[i].c_str(),ios::in);    
    if(file.is_open()){
      betadata[i].resize(lines3[i]);
      for(int k=0;k<lines3[i];k++){      
	file >> dummy;
	betadata[i][0][k] = log10(dummy/1000.);
	file >> betadata[i][1][k];
	betadata[i][2][k] = 1.;
	file.ignore(500,'\n');
      }
      file.close();
    }
    
    else{
      throw std::runtime_error(nucqcd::ErrorMessage() << "File " << files3[i] << " could not be opened" );
    }
  }

  double *a;
  double chi2;

  //pi- scattering
    
  //sigmatot in fm^2
  if(logp<=-0.19){
    a=fitdata(data[2][0],data[2][1],data[2][2],0,195,9,chi2);
    sigmatotn=0.;
    for(int i=0;i<9;i++) sigmatotn+=a[i]*pow(logp,i);
    sigmatotn=pow(10.,sigmatotn)/10;
    delete [] a;
  }
  else if(logp<=0.03){
    a=fitdata(data[2][0],data[2][1],data[2][2],183,330,8,chi2);
    sigmatotn=0.;
    for(int i=0;i<8;i++) sigmatotn+=a[i]*pow(logp,i);
    sigmatotn=pow(10.,sigmatotn)/10;
      delete [] a;
    }
  else if(logp<=0.47){
    a=fitdata(data[2][0],data[2][1],data[2][2],315,450,10,chi2);
    sigmatotn=0.;
    for(int i=0;i<10;i++) sigmatotn+=a[i]*pow(logp,i);
    sigmatotn=pow(10.,sigmatotn)/10;
    delete [] a;
  }
  else {
    a=fitdata(data[2][0],data[2][1],data[2][2],440,lines[2]-1,5,chi2);
    sigmatotn=0.;
    for(int i=0;i<5;i++) sigmatotn+=a[i]*pow(logp,i);
    sigmatotn=pow(10.,sigmatotn)/10;
    delete [] a;
    }

   //sigmael
  if(logp<=-0.19){
    a=fitdata(data[0][0],data[0][1],data[0][2],0,55,9,chi2);
    sigmaeln=0.;
    for(int i=0;i<9;i++) sigmaeln+=a[i]*pow(logp,i);
    sigmaeln=pow(10.,sigmaeln)/10;
    delete [] a;
  }
  else if(logp<=0.1){
    a=fitdata(data[0][0],data[0][1],data[0][2],45,120,10,chi2);
    sigmaeln=0.;
    for(int i=0;i<10;i++) sigmaeln+=a[i]*pow(logp,i);
    sigmaeln=pow(10.,sigmaeln)/10;
    delete [] a;
  }
  else {
    a=fitdata(data[0][0],data[0][1],data[0][2],110,lines[0]-1,4,chi2);
    sigmaeln=0.;
    for(int i=0;i<4;i++) sigmaeln+=a[i]*pow(logp,i);
    sigmaeln=pow(10.,sigmaeln)/10;
    delete [] a;
  }
    


  //eps
  /*a=fitdata(data[1][0],data[1][1],data[1][2],0,lines[1]-1,4,chi2);
  epsn=0.;
  for(int i=0;i<4;i++) epsn+=a[i]*pow(logp,i);
  delete [] a;*/

  //eps2, fit to new data
  if(logp<=0.04){
    std::vector<double> tempdata(41);
    for(int i=0;i<41;i++) tempdata[i]=pow(10.,epsdata[0][0][i])*1000.;
    a=fitdata(tempdata.data(),epsdata[0][1],epsdata[0][2],0,40,20,chi2);
    epsn=0.;
    for(int i=0;i<20;i++) epsn+=a[i]*pow(mom,i);
    delete [] a;
  }

  else if(logp<=0.57){
    a=fitdata(epsdata[0][0],epsdata[0][1],epsdata[0][2],35,130,12,chi2);
    epsn=0.;
    for(int i=0;i<12;i++) epsn+=a[i]*pow(logp,i);
    delete [] a;
  }
  else{
    a=fitdata(epsdata[0][0],epsdata[0][1],epsdata[0][2],120,lines2[0]-1,4,chi2);
    epsn=0.;
    for(int i=0;i<4;i++) epsn+=a[i]*pow(logp,i);
    delete [] a;
  }

  //beta, 2 possibilites:
  // 1. perform fit
  // 2. calc from other parameters
  if(logp<=0.09){
    a=fitdata(betadata[0][0],betadata[0][1],betadata[0][2],0,47,11,chi2);
    beta2n=0.;
    for(int i=0;i<11;i++) beta2n+=a[i]*pow(logp,i);
    beta2n = beta2n*pow(nucqcd::hbarc/1000.,2.);
    delete [] a;
  }
  else if(logp<=0.6){
    a=fitdata(betadata[0][0],betadata[0][1],betadata[0][2],40,105,11,chi2);
    beta2n=0.;
    for(int i=0;i<11;i++) beta2n+=a[i]*pow(logp,i);
    beta2n = beta2n*pow(nucqcd::hbarc/1000.,2.);
    delete [] a;
  }
  else{
    a=fitdata(betadata[0][0],betadata[0][1],betadata[0][2],102,lines3[0]-1,5,chi2);
    beta2n=0.;
    for(int i=0;i<5;i++) beta2n+=a[i]*pow(logp,i);
    beta2n = beta2n*pow(nucqcd::hbarc/1000.,2.);
    delete [] a;
  }    

  betan2=sigmatotn*sigmatotn*(1+epsn*epsn)/(16*PI*sigmaeln);

  //epsn=sqrt((16*PI*sigmaeln*betan*betan)/sigmatotn/sigmatotn-1);
  
  //pi+ scattering
  //sigmatot in fm^2
  if(logp<=-0.22){
    a=fitdata(data[5][0],data[5][1],data[5][2],0,135,10,chi2);
    sigmatotp=0.;
    for(int i=0;i<10;i++) sigmatotp+=a[i]*pow(logp,i);
    sigmatotp=pow(10.,sigmatotp)/10;
    delete [] a;
  }
  else if(logp<=0.32){
    a=fitdata(data[5][0],data[5][1],data[5][2],120,283,11,chi2);
    sigmatotp=0.;
    for(int i=0;i<11;i++) sigmatotp+=a[i]*pow(logp,i);
    sigmatotp=pow(10.,sigmatotp)/10;
    delete [] a;
  }
  else {
    a=fitdata(data[5][0],data[5][1],data[5][2],277,lines[5]-1,10,chi2);
    sigmatotp=0.;
    for(int i=0;i<10;i++) sigmatotp+=a[i]*pow(logp,i);
    sigmatotp=pow(10.,sigmatotp)/10;
    delete [] a;
  }

  //sigmael
  if(logp<=-0.18){
    /*a=fitdata(data[3][0],data[3][1],data[3][2],0,40,8,chi2);
      sigmael=0.;
      for(int i=0;i<8;i++) sigmael+=a[i]*pow(logp,i);
      sigmael=pow(10.,sigmael);
      delete [] a;*/
    sigmaelp=sigmatotp; //data coincide, so no need to do second fit
    }
  else if(logp<=0.26){
    a=fitdata(data[3][0],data[3][1],data[3][2],30,85,10,chi2);
    sigmaelp=0.;
    for(int i=0;i<10;i++) sigmaelp+=a[i]*pow(logp,i);
    sigmaelp=pow(10.,sigmaelp)/10;
    delete [] a;
  }
  else {
    a=fitdata(data[3][0],data[3][1],data[3][2],77,lines[3]-1,7,chi2);
    sigmaelp=0.;
    for(int i=0;i<7;i++) sigmaelp+=a[i]*pow(logp,i);
    sigmaelp=pow(10.,sigmaelp)/10;
    delete [] a;
  }

  
  //eps
  /*a=fitdata(data[4][0],data[4][1],data[4][2],0,lines[4]-1,2,chi2);
  epsp=0.;
  for(int i=0;i<2;i++) epsp+=a[i]*pow(logp,i);
  delete [] a;*/

  //eps2, fit to new data
  if(logp<=0.15){
    a=fitdata(epsdata[1][0],epsdata[1][1],epsdata[1][2],0,57,10,chi2);
    epsp=0.;
    for(int i=0;i<10;i++) epsp+=a[i]*pow(logp,i);
    delete [] a;
  }

  else if(logp<=0.50){
    a=fitdata(epsdata[1][0],epsdata[1][1],epsdata[1][2],50,123,10,chi2);
    epsp=0.;
    for(int i=0;i<10;i++) epsp+=a[i]*pow(logp,i);
    delete [] a;
  }
  else{
    a=fitdata(epsdata[1][0],epsdata[1][1],epsdata[1][2],105,lines2[1]-1,3,chi2);
    epsp=0.;
    for(int i=0;i<3;i++) epsp+=a[i]*pow(logp,i);
    delete [] a;
  }

  
  //beta
  //again 2 possibilities
  if(logp<=0.255){
    a=fitdata(betadata[1][0],betadata[1][1],betadata[1][2],0,52,8,chi2);
    beta2p=0.;
    for(int i=0;i<8;i++) beta2p+=a[i]*pow(logp,i);
    beta2p=beta2p*pow(nucqcd::hbarc/1000.,2.);
    delete [] a;
  }
  else{
    a=fitdata(betadata[1][0],betadata[1][1],betadata[1][2],46,lines3[1]-1,8,chi2);
    beta2p=0.;
    for(int i=0;i<8;i++) beta2p+=a[i]*pow(logp,i);
    beta2p=beta2p*pow(nucqcd::hbarc/1000.,2.);
    delete [] a;
  }

  betap2=sigmatotp*sigmatotp*(1+epsp*epsp)/(16*PI*sigmaelp);
  //epsp=sqrt((16*PI*sigmaelp*betap*betap)/sigmatotp/sigmatotp-1);
      
}

void FastParticle::interpPionGlauberData(int particletype, double mom, double &sigmap, double &beta2p, double &epsilonp, 
		     double &sigman, double &beta2n, double &epsilonn){
  switch(particletype){
  case(2):
    //cout <<"pi+" << endl;
    sigmap = interpolate(sigmap_array,log10(mom),0.01,151,250.);
    beta2p = interpolate(beta2p_array,log10(mom),0.01,151,250.);
    epsilonp = interpolate(epsp_array,log10(mom),0.01,151,250.);
    sigman = interpolate(sigman_array,log10(mom),0.01,151,250.);
    beta2n = interpolate(beta2n_array,log10(mom),0.01,151,250.);
    epsilonn = interpolate(epsn_array,log10(mom),0.01,151,250.);
    break;
  case(3):
    //cout <<"pi-" << endl;
    // isospin: pi- p = pi+ n and pi- n = pi+ p, so the proton and neutron tables are swapped
    sigmap = interpolate(sigman_array,log10(mom),0.01,151,250.);
    beta2p = interpolate(beta2n_array,log10(mom),0.01,151,250.);
    epsilonp = interpolate(epsn_array,log10(mom),0.01,151,250.);
    sigman = interpolate(sigmap_array,log10(mom),0.01,151,250.);
    beta2n = interpolate(beta2p_array,log10(mom),0.01,151,250.);
    epsilonn = interpolate(epsp_array,log10(mom),0.01,151,250.);
    break;
  case(4):
    //cout <<"rho0" << endl;
    sigman=sigmap = (interpolate(sigman_array,log10(mom),0.01,151,250.)+interpolate(sigmap_array,log10(mom),0.01,151,250.))*0.5;
    beta2n=beta2p = (interpolate(beta2n_array,log10(mom),0.01,151,250.)+interpolate(beta2p_array,log10(mom),0.01,151,250.))*0.5;
    epsilonn = epsilonp = (interpolate(epsn_array,log10(mom),0.01,151,250.)+interpolate(epsp_array,log10(mom),0.01,151,250.))*0.5;
    break;
  default:
    throw std::invalid_argument(nucqcd::ErrorMessage() << "Particle type is not yet supported!!: " << particletype );
  }
  return;
}


}  // namespace nucqcd

  
