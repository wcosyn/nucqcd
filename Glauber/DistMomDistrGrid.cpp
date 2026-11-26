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


//calc both fsi and fsi+ct grid
void DistMomDistrGrid::constructAllGrids(nucqcd::Rotation & rot){
  //fill the arrays!
  for(int i=0; i<=getPgrid(); i++){
    p_hit = float(i)/getInvPstep();
    for(int j=0;j<=getCthgrid();j++){
      costheta_hit = -2.*j/getCthgrid()+1.;
      sintheta_hit=sqrt(1.-costheta_hit*costheta_hit);
      for(int k=0;k<=getPhigrid();k++){
	phi_hit = 2.*PI*k/getPhigrid();
	sincos(phi_hit,&sinphi_hit,&cosphi_hit);
	pvec_hit= nucqcd::Vector3(p_hit*sintheta_hit*cosphi_hit,p_hit*sintheta_hit*sinphi_hit,p_hit*costheta_hit);
	pvec_hit=pvec_hit.Transform(rot);	
	for(int l=0;l<number_of_Grids/2;l++) rhogrid[l][i][j][k]=rhoctgrid[l][i][j][k]=0.;
	//symmetry for (-m,-ms)(m,ms), so only half needed, multiply everything by 2 in the end.
	double absprec=1.E-12;
	for(int m=1;m<=getPfsigrid()->getPnucleus()->getJ_array()[shellindex];m+=2){
	  getPfsigrid()->clearKnockout();
	  getPfsigrid()->addKnockout(getShellindex(),m);
	  for(int ms=0;ms<=1;ms++){
	    FourVector<double> pf(sqrt(p_hit*p_hit+mass*mass),pvec_hit.X(),pvec_hit.Y(),pvec_hit.Z());
	    if(ms) Upm_bar=Spinor::Bar(Spinor(pf,mass,Spinor::Polarization(0.,0.,Spinor::Polarization::kDown),Spinor::kUnity));
	    else Upm_bar=Spinor::Bar(Spinor(pf,mass,Spinor::Polarization(0.,0.,Spinor::Polarization::kUp),Spinor::kUnity));
	    int res=90;
	    unsigned count=0;
	    double deeserror=0.;
	      numint::array<double,3> lower = {{0.,-1.,0.}};
	      numint::array<double,3> upper = {{getPfsigrid()->getPnucleus()->getRange(),1.,2.*PI}};
	      DistMomDistrGrid::Ftor_one F;
	      F.grid = this;
	      F.m = m;
	      numint::mdfunction<numint::vector_z,3> mdf;
	      mdf.func = &DistMomDistrGrid::Ftor_one::exec;
	      mdf.param = &F;
	      vector<complex<double> > ret(number_of_Grids,0.);
	      F.f=DistMomDistrGrid::klaas_distint;
	      res = numint::cube({integrator, absprec, prec, 2000, maxEval}, mdf, lower, upper, ret, count, 0);
	      if(abs(ret[i])*prec*0.1>absprec) absprec=abs(ret[i])*prec*0.1;
	      for(int l=0;l<number_of_Grids/2;l++){
          rhogrid[l][i][j][k]+=norm(ret[l]);
          rhoctgrid[l][i][j][k]+=norm(ret[l+number_of_Grids/2]);
	      }
//   	    cout << i << " " << j << " " << k << " " << m << " " << ms << " " << rhogrid[0][i][j][k] << " " << rhogrid[1][i][j][k]<< " " << rhogrid[2][i][j][k]<< " " << rhoctgrid[0][i][j][k]<< " " << res << " " << count << endl;
	    
	  }
	}
	for(int l=0;l<number_of_Grids/2;l++){
	  //extra factor of 2 due to symmetry
	  rhogrid[l][i][j][k]/=pow(2.*PI,3.)/2.;
	  rhoctgrid[l][i][j][k]/=pow(2.*PI,3.)/2.;
	}
	//r=0 symmetry shortcut
// 	cout << p_hit << " " << costheta_hit << " " << phi_hit*nucqcd::rad_to_deg << " " << rhogrid[0][i][j][k] << " " << rhogrid[1][i][j][k]<< " " << rhogrid[2][i][j][k]<< " " << rhoctgrid[0][i][j][k]<< endl;
	if(i==0){
	  for(j=0;j<=getCthgrid();j++){
	    for(k=0;k<=getPhigrid();k++){
	      for(int l=0;l<number_of_Grids/2;l++){
		rhogrid[l][i][j][k]=rhogrid[l][0][0][0];
		rhoctgrid[l][i][j][k]=rhoctgrid[l][0][0][0];
	      }
	    }
	  }
	}
	//theta=0 or Pi symmetry shortcut
	else if(j==0||j==getCthgrid()){
	  for(k=1;k<=getPhigrid();k++){
	    for(int l=0;l<number_of_Grids/2;l++){
	      rhogrid[l][i][j][k]=rhogrid[l][i][j][0];
	      rhoctgrid[l][i][j][k]=rhoctgrid[l][i][j][0];
	    }
	  }
	}
		    

	
	
      }
    }
  }
}

//calc both fsi and fsi+ct grid
void DistMomDistrGrid::constructpwGrid(){
//   cout << "Constructing pw momentum distribution" << endl;
  //fill the arrays!
  for(int i=0; i<=getPgrid(); i++){
    p_hit = float(i)/getInvPstep();
    // radial integral of r G(r) j_l(p r): a 1D adaptive cubature to 1e-7 relative
    Ftor_pw F;
    F.grid = this;
    numint::mdfunction<numint::vector_d,1> mdf;
    mdf.func = &Ftor_pw::exec;
    mdf.param = &F;
    numint::array<double,1> lower = {{0.}};
    numint::array<double,1> upper = {{getPfsigrid()->getPnucleus()->getRange()}};
    numint::vector_d ret(1,0.);
    unsigned count=0;
    numint::cube_adaptive(mdf,lower,upper,1.E-14,1.E-07,1E02,1E05,ret,count,0);
    double result = ret[0];
    rhopwgrid[i] = result*result*2.*getMass()/(sqrt(getMass()*getMass()+p_hit*p_hit)+getMass())*
		    (getPfsigrid()->getPnucleus()->getJ_array()[getShellindex()]+1)/(4.*PI)*(2./PI);  //(2j+1)/(4pi)*(2/pi)
  }
}

//calc both fsi and fsi+ct grid
void DistMomDistrGrid::constructCtGrid(nucqcd::Rotation & rot){
  //fill the arrays!
  for(int i=0; i<=getPgrid(); i++){
    p_hit = float(i)/getInvPstep();
    for(int j=0;j<=getCthgrid();j++){
      costheta_hit = -2.*j/getCthgrid()+1.;
      sintheta_hit=sqrt(1.-costheta_hit*costheta_hit);
      for(int k=0;k<=getPhigrid();k++){
	phi_hit = 2.*PI*k/getPhigrid();
	sincos(phi_hit,&sinphi_hit,&cosphi_hit);
	pvec_hit= nucqcd::Vector3(p_hit*sintheta_hit*cosphi_hit,p_hit*sintheta_hit*sinphi_hit,p_hit*costheta_hit);	
	pvec_hit=pvec_hit.Transform(rot);
	for(int l=0;l<number_of_Grids/2;l++) rhoctgrid[l][i][j][k]=0.;
	double absprec=1.E-12;
	//symmetry in m_j,m_s final proton needs only half the values
	for(int m=1;m<=getPfsigrid()->getPnucleus()->getJ_array()[shellindex];m+=2){
	  getPfsigrid()->clearKnockout();
	  getPfsigrid()->addKnockout(getShellindex(),m);
	  for(int ms=0;ms<=1;ms++){
	    FourVector<double> pf(sqrt(p_hit*p_hit+mass*mass),pvec_hit.X(),pvec_hit.Y(),pvec_hit.Z());
	    if(ms) Upm_bar=Spinor::Bar(Spinor(pf,mass,Spinor::Polarization(0.,0.,Spinor::Polarization::kDown),Spinor::kUnity));
	    else Upm_bar=Spinor::Bar(Spinor(pf,mass,Spinor::Polarization(0.,0.,Spinor::Polarization::kUp),Spinor::kUnity));
	    
    	    int res=90;
	    unsigned count=0;
	    double deeserror=0.;
	      numint::array<double,3> lower = {{0.,-1.,0.}};
	      numint::array<double,3> upper = {{getPfsigrid()->getPnucleus()->getRange(),1.,2.*PI}};
	      DistMomDistrGrid::Ftor_one F;
	      F.grid = this;
	      F.m = m;
	      numint::mdfunction<numint::vector_z,3> mdf;
	      mdf.func = &Ftor_one::exec;
	      mdf.param = &F;
	      vector<complex<double> > ret(number_of_Grids/2,0.);
	      F.f=DistMomDistrGrid::klaas_distint_ct;
	      res = numint::cube({integrator, absprec, prec, 2000, maxEval}, mdf, lower, upper, ret, count, 0);
	      if(abs(ret[i])*prec*0.1>absprec) absprec=abs(ret[i])*prec*0.1;
	      for(int l=0;l<number_of_Grids/2;l++){
		rhoctgrid[l][i][j][k]+=norm(ret[l]);
	      }
// 	    cout << i << " " << j << " " << k << " " << m << " " << ms << " " << rhoctgrid[0][i][j][k] << " " << rhoctgrid[1][i][j][k]<< " " << rhoctgrid[2][i][j][k]<< " " << res << " " << count << endl;
	  }
	}
	for(int l=0;l<number_of_Grids/2;l++){
	  rhoctgrid[l][i][j][k]/=pow(2.*PI,3.)/2.;  //factor 2 because of symmetry
	}
	//r=0 symmetry shortcut
	if(i==0){
	  for(j=0;j<=getCthgrid();j++){
	    for(k=0;k<=getPhigrid();k++){
	      for(int l=0;l<number_of_Grids/2;l++){
		rhoctgrid[l][i][j][k]=rhoctgrid[l][0][0][0];
	      }
	    }
	  }
	}
	//theta=0 or Pi symmetry shortcut
	else if(j==0||j==getCthgrid()){
	  for(k=1;k<=getPhigrid();k++){
	    for(int l=0;l<number_of_Grids/2;l++){
	      rhoctgrid[l][i][j][k]=rhoctgrid[l][i][j][0];
	    }
	  }
	}
		    

	
	
      }
    }
  }
}

// void DistMomDistrGrid::intRhoR(const double r, complex<double> *results, va_list ap){
//   
//   int m = va_arg(ap,int);
//   double *pthetaestimate = va_arg(ap,double*);
//   double *pphiestimate = va_arg(ap,double*);
//   
//   
//   getPfsigrid()->setRinterp(r);
//   rombergerN(this,&DistMomDistrGrid::intRhoCosTheta,-1.,1.,number_of_Grids,
// 	     results,getPrec(),3,6,pthetaestimate, m,  r, pphiestimate);
//   for(int i=0;i<number_of_Grids;i++) results[i]*=r;
//   return;
// }
// 
// void DistMomDistrGrid::intRhoRCT(const double r, complex<double> *results, va_list ap){
//   
//   int m = va_arg(ap,int);
//   double *pthetaestimate = va_arg(ap,double*);
//   double *pphiestimate = va_arg(ap,double*);
// 
//   getPfsigrid()->setRinterp(r);
//   rombergerN(this,&DistMomDistrGrid::intRhoCosThetaCT,-1.,1.,number_of_Grids/2,
// 	     results,getPrec(),3,6,pthetaestimate, m,  r, pphiestimate);
//   for(int i=0;i<number_of_Grids/2;i++) results[i]*=r;
//   return;
// }
// 
double DistMomDistrGrid::intRhoRpw(const double r) const{
  return r*getPfsigrid()->getPnucleus()->getWave_G(getShellindex(),r)
	*gsl_sf_bessel_jl(getPfsigrid()->getPnucleus()->getL_array()[getShellindex()],p_hit*r*nucqcd::inv_hbarc);
}
// void DistMomDistrGrid::intRhoCosTheta(const double costheta, complex<double> *results, va_list ap){
//   
//   int m = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double *pphiestimate = va_arg(ap,double*);
//   double sintheta=sqrt(1.-costheta*costheta);
//   getPfsigrid()->setCthinterp(costheta);
//   rombergerN(this,&DistMomDistrGrid::intRhoPhi,0.,2.*PI,number_of_Grids,
// 	     results,getPrec(),3,6,pphiestimate, m, r, costheta, sintheta);
//   return;
// }
// 
// void DistMomDistrGrid::intRhoCosThetaCT(const double costheta, complex<double> *results, va_list ap){
//   
//   int m = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double *pphiestimate = va_arg(ap,double*);
//   double sintheta=sqrt(1.-costheta*costheta);
//   getPfsigrid()->setCthinterp(costheta);
//   rombergerN(this,&DistMomDistrGrid::intRhoPhiCT,0.,2.*PI,number_of_Grids/2,
// 	     results,getPrec(),3,6,pphiestimate, m, r, costheta, sintheta);
//   return;
// }
// 
// 
// void DistMomDistrGrid::intRhoPhi(const double phi, complex<double> *results, va_list ap){
// 
//   int m = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double costheta = va_arg(ap,double);
//   double sintheta = va_arg(ap,double);
//   double cosphi,sinphi;
//   sincos(phi,&sinphi,&cosphi);
//   
//   MeanFieldSpinor wave(*getPfsigrid()->getPnucleus(),getShellindex(),m,r,costheta,phi);
//   complex<double> exp_pr=exp(-nucqcd::inv_hbarc*pvec_hit*nucqcd::Vector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*nucqcd::i_unit)
// 			  *(Upm_bar*wave);  
// 
//   for(int i=0;i<number_of_Grids;i++) results[i]=exp_pr*getPfsigrid()->getFsiGridN_interp1(i,phi);
//   return;
// }
// 
// void DistMomDistrGrid::intRhoPhiCT(const double phi, complex<double> *results, va_list ap){
// 
//   int m = va_arg(ap,int);
//   double r = va_arg(ap,double);
//   double costheta = va_arg(ap,double);
//   double sintheta = va_arg(ap,double);
//   double cosphi,sinphi;
//   sincos(phi,&sinphi,&cosphi);
//   MeanFieldSpinor wave(*getPfsigrid()->getPnucleus(),getShellindex(),m,r,costheta,phi);
//   
//   
//   complex<double> exp_pr=exp(-nucqcd::inv_hbarc*pvec_hit*nucqcd::Vector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta)*nucqcd::i_unit)
// 			  *(Upm_bar*wave);  
// 
//   for(int i=0;i<number_of_Grids/2;i++) results[i]=exp_pr*getPfsigrid()->getFsiGridN_interp1(i+number_of_Grids/2,phi);
//   return;
// }


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
void DistMomDistrGrid::klaas_distint(numint::vector_z &results, double r, double costheta, double phi, DistMomDistrGrid & grid, 
			  int m){
  results=numint::vector_z(grid.number_of_Grids,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  MeanFieldSpinor wave(*grid.getPfsigrid()->getPnucleus(),grid.getShellindex(),m,r,costheta,phi);
  complex<double> exp_pr=exp(-nucqcd::inv_hbarc*((grid.getPvec_hit())*nucqcd::Vector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta))*nucqcd::i_unit)
			  *((grid.getUpm_bar())*wave);  
  for(int i=0;i<grid.number_of_Grids;i++) results[i]=exp_pr*grid.getPfsigrid()->getFsiGridN_interp3(i,r,costheta,phi);
  for(int i=0;i<grid.number_of_Grids;i++) results[i]*=r;
  return;
}
/*! integrandum function (clean ones), only CT*/
void DistMomDistrGrid::klaas_distint_ct(numint::vector_z &results, double r, double costheta, double phi, DistMomDistrGrid & grid, 
			      int m){

  results=numint::vector_z(grid.number_of_Grids/2,0.);
  double sintheta=sqrt(1.-costheta*costheta);
  double cosphi,sinphi;
  sincos(phi,&sinphi,&cosphi);
  MeanFieldSpinor wave(*grid.getPfsigrid()->getPnucleus(),grid.getShellindex(),m,r,costheta,phi);
  complex<double> exp_pr=exp(-nucqcd::inv_hbarc*((grid.getPvec_hit())*nucqcd::Vector3(r*sintheta*cosphi,r*sintheta*sinphi,r*costheta))*nucqcd::i_unit)
			  *((grid.getUpm_bar())*wave);  

  for(int i=0;i<grid.number_of_Grids/2;i++) results[i]=exp_pr*
  grid.getPfsigrid()->getFsiGridN_interp3(i+grid.number_of_Grids/2,r,costheta,phi);
  for(int i=0;i<grid.number_of_Grids/2;i++) results[i]*=r;
  return;

}

}  // namespace nucqcd





