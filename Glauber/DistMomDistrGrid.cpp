#include <Tools/Core/IntegratorFlag.hpp>
#include <Integrator/numint/cube.hpp>
#include <cstdlib>
#include <Tools/Core/ErrorMessage.hpp>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include <cassert>

using namespace std;

#include <Glauber/DistMomDistrGrid.hpp>
#include <Glauber/AbstractFsiCTGrid.hpp>
#include <Tools/Core/constants.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/GridCache.hpp>
#include <Tools/Core/AuxFunction.hpp>
#include <Tools/Core/FourVector.hpp>
#include <Tools/Core/Matrix.hpp>
#include <Tools/Core/Spinor.hpp>
#include <Tools/Core/Rotation.hpp>
#include <Glauber/MeanFieldSpinor.hpp>

namespace nucqcd {


DistMomDistrGrid::DistMomDistrGrid():
pgrid(0),
cthgrid(0),
phigrid(0),
pfsigrid(NULL),
dir(""){
//   cerr << "This DistMomDistrGrid default constructor shouldn't be called!" << endl;
//   exit(1);
}

DistMomDistrGrid::DistMomDistrGrid(const int shell, const double p_max, const int p_grid, const int cth_grid, const int phi_grid,
				   AbstractFsiCTGrid *pfsi_grid, const double precision, const IntegrationMethod integr, 
				   const int max_Eval, const double theta_rot, const std::string homedir):
shellindex(shell),
pmax(p_max),
filledgrid(0),
filledctgrid(0),
filledallgrid(0),
dir(nucqcd::grid_cache_dir(homedir)),
pgrid(p_grid),
cthgrid(cth_grid),
phigrid(phi_grid),
pfsigrid(pfsi_grid),
prec(precision),
integrator(integr),
maxEval(max_Eval),
thetarot(theta_rot),
number_of_Grids(pfsi_grid->getNumber_of_grids()){
  
  invpstep=pgrid/p_max;
  invcthstep=cthgrid/2.;
  invphistep=0.5*phigrid/PI;
  if(pfsigrid!=NULL){
    mass = getShellindex()<getPfsigrid()->getPnucleus()->getPLevels()? nucqcd::mass_p:nucqcd::mass_n;
    const Grid3D<double> points(getPgrid()+1,getCthgrid()+1,getPhigrid()+1);
    rhogrid.assign(number_of_Grids/2,points);
    rhoctgrid.assign(number_of_Grids/2,points);
    rhopwgrid.assign(getPgrid()+1,0.);
    constructpwGrid();
  }
   //setFilenames(dir);
  //fillGrids();
  
  
}  
  
//set up inteprolation help variables for r
void DistMomDistrGrid::setPinterp(const double p){
  if(p>getPmax()){
    throw std::out_of_range(nucqcd::ErrorMessage() << "p out of range: " << p );
  }
  pindex = int(floor(p*getInvPstep()));
  if(pindex==getPgrid()) pindex-=1;
  s_interp = (p*getInvPstep() - (pindex));
  comp_s_interp=1.-s_interp;
}

//set up interpolation help variables for theta
void DistMomDistrGrid::setCthinterp(const double costheta){
  cthindex = int(floor((-costheta+1.)*getInvCthstep()));
  if(cthindex==getCthgrid()) cthindex-=1;
  t_interp = ((-costheta+1.)*getInvCthstep() - (cthindex));  
  comp_t_interp=1.-t_interp;
}

//set up interpolation help variables for phi
void DistMomDistrGrid::setPhiinterp(const double phi){
  phiindex = int(floor(phi*getInvPhistep()));
  if(phiindex==getPhigrid()) phiindex-=1;
  u_interp = (phi*getInvPhistep() - (phiindex));  
  comp_u_interp=1.-u_interp;
}

//interpolation functions
double DistMomDistrGrid::getRhoGridFull_interpvec(int gridindex, const nucqcd::Vector3 &pvec){
  return getRhoGridFull_interp3(gridindex, pvec.Mag(),pvec.CosTheta(),pvec.Phi());
}

double DistMomDistrGrid::getRhoGridFull_interp3(int gridindex, const double p, const double costheta, const double phi){
  setPinterp(p);
  setCthinterp(costheta);
  setPhiinterp(phi);
  return getRhoGridFull_interp(gridindex);
}

double DistMomDistrGrid::getRhoGridFull_interp2(int gridindex, const double costheta, const double phi){
  setCthinterp(costheta);
  setPhiinterp(phi);
  return getRhoGridFull_interp(gridindex);
}

double DistMomDistrGrid::getRhoGridFull_interp1(int gridindex, const double phi){
  setPhiinterp(phi);
  return getRhoGridFull_interp(gridindex);
}

double DistMomDistrGrid::getRhoGridFull_interp(int gridindex){
  return Interp3d((gridindex < (number_of_Grids/2))? 
		    rhogrid[gridindex] : rhoctgrid[gridindex-number_of_Grids/2],
		  getS_interp(), getT_interp(), getU_interp(), 
		  getComp_s_interp(), getComp_t_interp(), getComp_u_interp(), 
		  getPindex(), getCthindex(), getPhiindex());  
}



double DistMomDistrGrid::getRhopwGridFull_interp(double p){
  if(p>getPmax()){
    throw std::out_of_range(nucqcd::ErrorMessage() << "p out of range in : getRhopwGridFull" << p );
  }
  return interpolate(rhopwgrid.data(),p,1./getInvPstep(),getPgrid()+1,0); 
  
}



//set filenames
void DistMomDistrGrid::setFilenames(string homedir){
  rho_filename=getPfsigrid()->getFsi_Filename();
  rhoct_filename=getPfsigrid()->getFsi_Ct_Filename();
  size_t found;
  found=rho_filename.rfind("/");
  if (found!=string::npos){
    string addendum = "Rhogrid.l"+to_string(getShellindex())+".throt"+to_string(int(round(thetarot*10.)))
		    +".p"+to_string(getPgrid())
		    +".cth"+to_string(getCthgrid())+".phi"+to_string(getPhigrid())+".pmax"+
		    to_string(getPmax())+".prec"+to_string(getPrec()*1.E05)+".integr"+to_string(numint::method_flag(integrator))+".";
    rho_filename.insert(found+1,addendum);
    rhoct_filename.insert(found+1,addendum);
  }
}
  

void DistMomDistrGrid::printRhopw_grid(){
  cout << "Printing plane-wave momentum array for " << getRho_Filename() << endl;
//   double tot = 0.;
  for(int i=0; i<=getPgrid(); i++){
    double p = float(i)/getInvPstep();
    cout << p << " " << getRhopwGridFull_interp(p) << endl;    
//     tot +=p*p*getRhopwGridFull_interp(p);
  }
//   tot*=4.*PI/getInvPstep()*pow(nucqcd::inv_hbarc,3.);
//   cout << "tot " << tot << endl;
}
  
void DistMomDistrGrid::printRho_grid(int gridindex){
  cout << "Printing Distorted momentum array for " << getRho_Filename() << endl << "Gridindex " << gridindex << endl;
  for(int i=0; i<=getPgrid(); i++){
    double p = float(i)/getInvPstep();
    for(int j=0;j<=getCthgrid();j++){
      double costheta = -2.*j/getCthgrid()+1.;
      for(int k=0;k<=getPhigrid();k++){
	double phi = k/getInvPhistep();	
	cout << p << " " << costheta << " " << phi*nucqcd::rad_to_deg << " " << getRhoGridFull_interp3(gridindex,p,costheta,phi) << endl;
      }
    }
  }

  
}
  
//readin fsi grid
void DistMomDistrGrid::readinRhoGrid(ifstream &infile){
  if(infile.peek() == std::ifstream::traits_type::eof()){
    cerr << "File is empty " << rho_filename << endl;
    assert(1==0);
  }
  
  for(int i=0;i<number_of_Grids/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  infile.read(reinterpret_cast<char *>(&rhogrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}
  
void DistMomDistrGrid::readinRhoCtGrid(ifstream &infile){
  if(infile.peek() == std::ifstream::traits_type::eof()){
    cerr << "File is empty " << rhoct_filename << endl;
    assert(1==0);
  }
  
  for(int i=0;i<number_of_Grids/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  infile.read(reinterpret_cast<char *>(&rhoctgrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}
  
  
//write out rho grid
void DistMomDistrGrid::writeoutRhoGrid(ofstream &outfile){
  for(int i=0;i<number_of_Grids/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  outfile.write(reinterpret_cast<char *>(&rhogrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}


//write out rho grid
void DistMomDistrGrid::writeoutRhoCtGrid(ofstream &outfile){
  for(int i=0;i<number_of_Grids/2;i++){
    for(int j=0;j<(getPgrid()+1);j++){
      for(int k=0;k<(getCthgrid()+1);k++){
        for(int l=0;l<(getPhigrid()+1);l++) 
	  outfile.write(reinterpret_cast<char *>(&rhoctgrid[i][j][k][l]),sizeof(double));	    	 
      }
    }
  }
}


void DistMomDistrGrid::fillGrids(nucqcd::Rotation & rot){
  setFilenames(dir);
  ifstream infile(rho_filename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  bool read=false;
  if(infile.is_open()){
    //cout << "Reading in FSI grid from memory: " << rho_filename << endl;
    readinRhoGrid(infile);
    read=!infile.fail(); // a truncated file sets failbit
    infile.close();
    if(read) filledgrid=1;
    else cerr << "Grid file is incomplete, recomputing it: " << rho_filename << endl;
  }
  if(!read){
    cout << "Constructing all grids" << endl;
    constructAllGrids(rot);
    filledgrid=filledallgrid=filledctgrid=1;
    if(nucqcd::write_file_atomically(rho_filename,[this](ofstream &outfile){ writeoutRhoGrid(outfile); }))
      cout << "Wrote FSI grid: " << rho_filename << endl;
    if(nucqcd::write_file_atomically(rhoct_filename,[this](ofstream &outfile){ writeoutRhoCtGrid(outfile); }))
      cout << "Wrote FSI+CT grid: " << rhoct_filename << endl;
    return;
  } 
  readinOrConstructCtGrid(rot);
}

void DistMomDistrGrid::readinOrConstructCtGrid(nucqcd::Rotation & rot){
  ifstream infile2(rhoct_filename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  if(infile2.is_open()){
    //cout << "Reading in FSI+CT grid from memory: " << rhoct_filename << endl;
    readinRhoCtGrid(infile2);
    bool read=!infile2.fail(); // a truncated file sets failbit
    infile2.close();
    if(read){
      filledctgrid=1;
      return;
    }
    cerr << "Grid file is incomplete, recomputing it: " << rhoct_filename << endl;
  }
  cout << "Constructing FSI+CT grid" << endl;
  constructCtGrid(rot);
  filledctgrid=1;
  if(nucqcd::write_file_atomically(rhoct_filename,[this](ofstream &outfile){ writeoutRhoCtGrid(outfile); }))
    cout << "Wrote FSI+CT grid: " << rhoct_filename << endl;
}


void DistMomDistrGrid::updateGrids(AbstractFsiCTGrid *pfsi_grid, int shell, nucqcd::Rotation & rot){
  nucqcd::Vector3 axis;
  rot.AngleAxis(thetarot,axis);
  thetarot*=-SIGN(axis[1]);
  if(!filledgrid){
  fillGrids(rot);
//      cout << "Grids still empty " << endl;
  }
  else{
    string old_rho_filename = rho_filename;
    string old_rhoct_filename = rhoct_filename;
    pfsigrid = pfsi_grid;
    shellindex=shell;
    setFilenames(dir);
  
    if(rho_filename.compare(old_rho_filename)){
//        cout << "fsi grid not equal" << rho_filename << endl << old_rho_filename << endl;
      fillGrids(rot);      
    }
//      else cout << "fsi grid equal to the earlier one, doing nothing" << endl << rho_filename << endl << old_rho_filename << endl;
    if(old_rhoct_filename.compare(rhoct_filename)){
//      cout << "fsict grid not equal " << endl << rhoct_filename << endl << old_rhoct_filename << endl;
    
      readinOrConstructCtGrid(rot);
    }
//   else cout << "fsict grid equal to the earlier one, doing nothing" << endl << rhoct_filename << endl << old_rhoct_filename << endl;

  }
}
 
 
/*! integrandum function (clean ones)*/
}  // namespace nucqcd





