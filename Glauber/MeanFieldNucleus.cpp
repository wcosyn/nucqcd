#include <iostream>
#include <Tools/Core/ErrorMessage.hpp>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <gsl/gsl_sf_legendre.h>

#include <Glauber/MeanFieldNucleus.hpp>
// #include "MeanFieldNucleusWrapper.hpp"
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/WaveFunctions/Deuteron.hpp>

namespace nucqcd {


using namespace std;

const std::map<std::string,MeanFieldNucleus::Type> MeanFieldNucleus::TypeNames = MeanFieldNucleus::initTypeNames();
//constructor
MeanFieldNucleus::MeanFieldNucleus(const int nucleus, const string & dir){
  
  setNucleusName(nucleus);
  setInputfile(nucleus, dir);  //set inputfile location
  //read in operation
  ifstream file(inputfile.c_str(),ios::in);
  if(file.is_open()){
    //cout << "Reading in nucleus properties from " << inputfile << endl;
    file >> A >> Z;
    N=A-Z;
    file.ignore(500,'\n');
   
    
    file >> massA;
    file.ignore(500,'\n');
    file >> massA_min_proton;
    file.ignore(500,'\n');
    file >> massA_min_neutron;
    file.ignore(500,'\n');
    file >> massA_min_pp;
    file.ignore(500,'\n');
    file >> massA_min_pn;
    file.ignore(500,'\n');
    file >> massA_min_nn;
    file.ignore(500,'\n');
       
    
    file >> plevels;
    file >> nlevels;
    totallevels=plevels+nlevels;
    file.ignore(500,'\n');
    excitation.assign(totallevels,0);
    for(int i=0;i<plevels;i++) file >> excitation[i]; //read in proton excitation energies
    file.ignore(500,'\n');
    for(int i=plevels;i<totallevels;i++) file >> excitation[i]; //neutron excitation energies
    for(int j=0;j<3;j++) file.ignore(500,'\n');
    n_array.assign(totallevels,0);
    for(int i=0;i<plevels;i++) file >> n_array[i]; //read in kappa for proton levels
    file.ignore(500,'\n');
    for(int i=plevels;i<totallevels;i++) file >> n_array[i]; //read in kappa for neutron levels
    file.ignore(500,'\n');
    kappas.assign(totallevels,0);
    for(int i=0;i<plevels;i++) file >> kappas[i]; //read in kappa for proton levels
    file.ignore(500,'\n');
    for(int i=plevels;i<totallevels;i++) file >> kappas[i]; //read in kappa for neutron levels
    file.close();
    l_array.assign(totallevels,0);
    lbar_array.assign(totallevels,0);
    j_array.assign(totallevels,0);
    for(int i=0;i<totallevels;i++){
      l_array[i]=kappas[i]>0? kappas[i]:-kappas[i]-1;
      lbar_array[i]=kappas[i]>0? kappas[i]-1:-kappas[i];
      j_array[i]=2*abs(kappas[i])-1;      
    }
    int counter=0;
    for(int i=0;i<getPLevels();i++)    counter+=2*abs(getKappas()[i]);  //# of protons in final shell
    onlyoneproton = ((counter-Z)%2)==0? 0:1; //redundant protons in final shell pair number?
    finalmproton = abs(kappas[getPLevels()-1])-1-(counter-Z-onlyoneproton)/2; //max m value we have to include in the summation over m to have the correct number of protons
    counter=0;  //same as above now for neutrons
    for (int i=getPLevels();i<getTotalLevels();i++) counter+=2*abs(kappas[i]);
    onlyoneneutron = ((counter-N)%2)==0? 0:1;
    finalmneutron = abs(kappas[getTotalLevels()-1])-1-(counter-N-onlyoneneutron)/2;

  }
  else{
    throw std::runtime_error(nucqcd::ErrorMessage() << "nucleus file could not be opened: " << inputfile );
  }
  
  //read in F&G grids
  readRgrids(nucleus, dir);
  
  //read in Ykappa grids
  constructThetaArray();
}





//set inputfile location
void MeanFieldNucleus::setInputfile(const int nucleus, const string &dir){
  inputfile=dir+"/nuclei/";
  switch(nucleus){
  case 0: //helium
    inputfile+="he4.inp";
    break;
  case 1: //carbon
    inputfile+="c12.inp";
    break;
  case 2: //oxygen
    inputfile+="o16.inp";
    break;
  case 3: //iron
    inputfile+="fe56.inp";
    break;
  case 4: //lead
    inputfile+="pb208.inp";
    break;
  case 5: //aluminium
    inputfile+="al27.inp";
    break;
  case 6: //cupper
    inputfile+="cu63.inp";
    break;
  case 7: //gold
    inputfile+="au197.inp";
    break;
  case 8: //fe54
    inputfile+="fe54.inp";
    break;
  case 9: //gold
    inputfile+="ca40.inp";
    break;
  case 10: //gold
    inputfile+="ca48.inp";
    break;
  default:
    throw std::invalid_argument(nucqcd::ErrorMessage() << "invalid nucleus input" );
  }
}

//set location of F&G grids file
void MeanFieldNucleus::setRgridfile(const int nucleus, const string &dir){
  rgridfile=dir+"/input/";
  switch(nucleus){
  case 0: //helium
    rgridfile+="he4.radial2.bound";
    break;
  case 1: //carbon
    rgridfile+="c12.radial2.bound";
    break;
  case 2: //oxygen
    rgridfile+="o16.radial2.bound";
    break;
  case 3: //iron
    rgridfile+="fe56.radial.bound";
    break;
  case 4: //lead
    rgridfile+="pb208.radial.bound";
    break;
  case 5: //lead
    rgridfile+="al27.radial2.bound";
    break;
  case 6: //lead
    rgridfile+="cu63.radial2.bound";
    break;
  case 7: //lead
    rgridfile+="au197.radial.bound";
    break;
  case 8: //fe54
    rgridfile+="fe56.radial.bound";
    break;
  case 9: //ca40
    rgridfile+="ca40.radial2.bound";
    break;
  case 10: //ca48
    rgridfile+="ca48.radial2.bound";
    break;
  default:
    throw std::runtime_error(nucqcd::ErrorMessage() << "invalid nucleus input for F&G grid file" );
  }  
}

void MeanFieldNucleus::setNucleusName(const int nucleus){
    switch(nucleus){
  case 0: //He
    nucleusname= "He";
    break;
  case 1: //C
    nucleusname= "C";
    break;
   case 2: //O
    nucleusname= "O";
    break;
  case 3: //Fe
    nucleusname= "Fe";
    break;
  case 4: //Pb
    nucleusname= "Pb";
    break;
  case 5: //Al
    nucleusname= "Al";
    break;
  case 6: //Cu
    nucleusname= "Cu";
    break;
  case 7: //Au
    nucleusname= "Au";
    break;
  case 8: //fe54
    nucleusname= "Fe54";
    break;
  case 9: //ca40
    nucleusname= "Ca40";
    break;
  case 10: //ca48
    nucleusname= "Ca48";
    break;
  default:
    throw std::invalid_argument(nucqcd::ErrorMessage() << "Implement Other Nucleus!" );
  }
}


void MeanFieldNucleus::readRgrids(const int nucleus, const string &dir){
  if(!F.empty()||!G.empty()){
     cerr << "Warning: F & G already initialized!!" << endl << "skipping initialization..." << endl;
     return;
  }

  setRgridfile(nucleus, dir);  //set inputfile location
  //read in operation
  ifstream file(rgridfile.c_str(),ios::in);
  
  if(file.is_open()){
    //cout << "Reading in F&G grids..." << endl;
    F.assign(totallevels,std::vector<double>());
    G.assign(totallevels,std::vector<double>());
    file >> range;
    file >> wf_r_step;
    file >> wf_r_lines;
    wf_r_lines--;
    int dummy;
    for(int i=0;i<totallevels;i++){
      F[i].assign(wf_r_lines,0.); G[i].assign(wf_r_lines,0.);
      file >> dummy;
      for(int j=0;j<wf_r_lines;j++){
	file >> dummy;
	file >> G[i][j];
	file >> F[i][j];
      }
    }
  }
  else{
    throw std::runtime_error(nucqcd::ErrorMessage() << "rgrids file could not be opened: " << rgridfile );
  }
}
/////////////////////////////////////////////////////////
//Spherical Harmonics without the phi part exp(i*m*phi)//
/////////////////////////////////////////////////////////
double MeanFieldNucleus::Spher_Harm(const int l, const int m, const double x){
  
  // This subroutine calculates the spherical harmonics minus the factor
  // exp(i*m*phi). Cos(theta) is x.

  if (abs(m) > l) return 0.;

  if (m < 0) return Spher_Harm(l,-m,x) * pow(-1.,m);

  else
    {
      return gsl_sf_legendre_sphPlm(l,m,x);
    }
}


///////////////////////////////////////////////////////////////////
//construct arrays with spherical harmonics  
//used in calculation of the glauberphase
///////////////////////////////////////////////////////////////////
void MeanFieldNucleus::constructThetaArray(){
  
  //first index: shell
  //second index: m quantum number.  Only positive m needed, so abs(kappa) in total.  0 is 1/2, 1 is 3/2 and so on
  //third index: grid in costheta

  //memory reservation
  //calculation of the grid
  if(!Yminkappa.empty()||!Ykappa.empty()){
    cerr << "Ykappa grids have already been initialized, skipping..." << endl;
    return;
  }
  //cout << "Initializing Ykappa grids..." << endl;
  Yminkappa.assign(totallevels,std::vector<std::vector<double>>());
  Ykappa.assign(totallevels,std::vector<std::vector<double>>());

  for(int i=0;i<totallevels;i++){
    Yminkappa[i].assign(abs(kappas[i]),std::vector<double>(GRIDP,0.));
    Ykappa[i].assign(abs(kappas[i]),std::vector<double>(GRIDP,0.));

    for(int j=0;j<abs(kappas[i]);j++){
      int m=2*j+1;  //2*m_j

      for(int k=0;k<GRIDP;k++){
	Yminkappa[i][j][k]=0.;
	Ykappa[i][j][k]=0.;
	for(int ms=-1;ms<=1;ms+=2){
	  Yminkappa[i][j][k]+=pow(Deuteron::Wavefunction::ClebschGordan(2*l_array[i],m-ms,1,ms, j_array[i],m)*Spher_Harm(l_array[i], (m-ms)/2, float(k)/float((GRIDP-1)/2)-1.),2);
	  Ykappa[i][j][k]+=pow(Deuteron::Wavefunction::ClebschGordan(2*lbar_array[i],m-ms,1,ms, j_array[i],m)*Spher_Harm(lbar_array[i], (m-ms)/2, float(k)/float((GRIDP-1)/2)-1.),2);
	}
      }
    }
  }
    
}

void MeanFieldNucleus::getWaveFunction(complex<double> *wave, const int shellindex, const int m, 
							       const double r, const double costheta, const double phi) const{
  // j, ms and m are inputted as 2*j, 2*ms and 2*m respectively
  double Gr=getWave_G(shellindex,r);
  double Fr=getWave_F(shellindex,r);
  complex<double> spinup=getWave_Phipart(m,1,phi);
  complex<double> spindown=getWave_Phipart(m,-1,phi);
  wave[0]=nucqcd::i_unit*Gr*spinup*Deuteron::Wavefunction::ClebschGordan(2*l_array[shellindex],m-1,1,1, j_array[shellindex],m)*Spher_Harm(l_array[shellindex], (m-1)/2, costheta);
  wave[1]=nucqcd::i_unit*Gr*spindown*Deuteron::Wavefunction::ClebschGordan(2*l_array[shellindex],m+1,1,-1, j_array[shellindex],m)*Spher_Harm(l_array[shellindex], (m+1)/2, costheta);
  wave[2]=-Fr*spinup*Deuteron::Wavefunction::ClebschGordan(2*lbar_array[shellindex],m-1,1,1, j_array[shellindex],m)*Spher_Harm(lbar_array[shellindex], (m-1)/2, costheta);
  wave[3]=-Fr*spindown*Deuteron::Wavefunction::ClebschGordan(2*lbar_array[shellindex],m+1,1,-1, j_array[shellindex],m)*Spher_Harm(lbar_array[shellindex], (m+1)/2, costheta);
}


double MeanFieldNucleus::getWave_F(const int shellindex, const double r) const{
  if(r>range){
    throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range in getWave_F: " << r << " " << r-range );
  }
  return interpolate(getF()[shellindex].data(),r,wf_r_step,wf_r_lines,1);  
}


double MeanFieldNucleus::getWave_G(const int shellindex, const double r) const{
  if(r>range){
    throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range in getWave_G: " << r << " " << r-range );
  }
  return interpolate(getG()[shellindex].data(),r,wf_r_step,wf_r_lines,1);
}
  

//interpolation of the Ykappa grid
double MeanFieldNucleus::getYkappacos(const int shellindex, const int m, const double costheta) const{
 return  interpolate(getYkappa()[shellindex][m].data(),costheta,double(2./(GRIDP-1)),GRIDP,-(GRIDP-1)/2);
}

//interpolation of the Yminkappa grid
double MeanFieldNucleus::getYminkappacos(const int shellindex, const int m, const double costheta) const{
 return  interpolate(getYminkappa()[shellindex][m].data(),costheta,double(2./(GRIDP-1)),GRIDP,-(GRIDP-1)/2);
}
//get exp(I*phi/2*m_l)
complex<double> MeanFieldNucleus::getWave_Phipart(const int m, const int spin, const double phi) const{
  // j, ms and m are inputted as 2*j, 2*ms and 2*m respectively
  return exp(nucqcd::i_unit*0.5*phi*double(m-spin));  
}

}  // namespace nucqcd


