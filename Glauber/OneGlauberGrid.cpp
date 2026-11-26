#include <complex>
#include <Integrator/numint/cube.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <iostream>
#include <cstdlib>

using namespace std;

#include <Glauber/OneGlauberGrid.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/AuxFunction.hpp>

namespace nucqcd {


//constructor, calls glaubergrid's constructor
OneGlauberGrid::OneGlauberGrid(const int r_grid, const int cth_grid, MeanFieldNucleus *pnucl, 
			       double prec, IntegrationMethod integrator, string dir):
AbstractFsiGrid(r_grid,cth_grid,0,pnucl,prec,integrator,dir),
GlauberGrid(r_grid,cth_grid,0,pnucl,prec,integrator,dir){
}

OneGlauberGrid::~OneGlauberGrid(){
}

//2d interpolation of array
complex<double> OneGlauberGrid::getInterp(complex<double> ***grid){
  return (
    complex<double>(getComp_t_interp()*(getComp_s_interp()*real(grid[getRindex()][getCthindex()][0])
	+ getS_interp()*real(grid[getRindex()+1][getCthindex()][0]))
	+ getT_interp()*(getComp_s_interp()*real(grid[getRindex()][getCthindex()+1][0])
	+ getS_interp()*real(grid[getRindex()+1][getCthindex()+1][0])),
	getComp_t_interp()*(getComp_s_interp()*imag(grid[getRindex()][getCthindex()][0])
	+ getS_interp()*imag(grid[getRindex()+1][getCthindex()][0])) 
	+ getT_interp()*(getComp_s_interp()*imag(grid[getRindex()][getCthindex()+1][0])
	+ getS_interp()*imag(grid[getRindex()+1][getCthindex()+1][0]))
 		  ));
}

void OneGlauberGrid::addParticle(FastParticle& newparticle){
 if(getParticles().size()!=0){
   throw std::invalid_argument(nucqcd::ErrorMessage() << "You cannot add more than one particle to the Glauber phase for one particle!!!" );
 }
 else{
   if(newparticle.getPhi()!=0.||newparticle.getTheta()!=0.){
     throw std::invalid_argument(nucqcd::ErrorMessage() << "phi and theta angles of the particle should be 0 in this special case!!!" );
   }
   else AbstractFsiGrid::addParticle(newparticle);
 }
}


//set filenames, don't forget to add ending "/" to dir if necessary!!!
void OneGlauberGrid::setFilenames(string dir){
 AbstractFsiCTGrid::setFilenames(dir+"OneGl."); 
}




//calc glauberphases for one gridpoint
void OneGlauberGrid::calcGlauberphasesBoth(const int i, const int j, const int k){
  double results[2]={0.,0.};
  for(int level=0;level<getPnucleus()->getTotalLevels();level++){
    for(int mm=0; mm<((getPnucleus()->getJ_array()[level]+1)/2);mm++){
      int res=90;
      unsigned count=0;
      double deeserror=0.;
	numint::array<double,3> lower = {{1.E-06,getParticles()[0].getHitz(),0.}};
	numint::array<double,3> upper = {{getPnucleus()->getRange(),getPnucleus()->getRange(),0.5*PI}};
	OneGlauberGrid::Ftor_one F;
	F.grid = this;
	F.level = level;
	F.mm = mm;
	numint::mdfunction<numint::vector_d,3> mdf;
	mdf.func = &Ftor_one::exec;
	mdf.param = &F;
	vector<double> ret(2,0.);
	F.f=OneGlauberGrid::klaas_one_bound;
	res = numint::cube({integrator, 1.E-12, prec, 2000, 2000000}, mdf, lower, upper, ret, count, 0);
	results[0]=ret[0]; results[1]=ret[1];
	
      
      fsi_grid[level][mm][i][j][k]=1.-getParticles()[0].getScatterfront(level,getPnucleus())*results[0];
      fsi_ct_grid[level][mm][i][j][k]=1.-getParticles()[0].getScatterfront(level,getPnucleus())*results[1];
//       cout << i << " " << j << " " << k << " " << level << " " << mm << fsi_grid[level][mm][i][j][k] << " " << fsi_ct_grid[level][mm][i][j][k] << 
//       " " << res << " " << count << endl;
    }
  }   
  
}

//calc glauberphases for one gridpoint,only CT grid
void OneGlauberGrid::calcGlauberphasesCt(const int i, const int j, const int k){
  double results=0.;
  for(int level=0;level<getPnucleus()->getTotalLevels();level++){
    for(int mm=0; mm<((getPnucleus()->getJ_array()[level]+1)/2);mm++){
      int res=90;
      unsigned count=0;
      double deeserror=0.;
	numint::array<double,3> lower = {{1.E-06,getParticles()[0].getHitz(),0.}};
	numint::array<double,3> upper = {{getPnucleus()->getRange(),getPnucleus()->getRange(),0.5*PI}};
	OneGlauberGrid::Ftor_one F;
	F.grid = this;
	F.level = level;
	F.mm = mm;
	numint::mdfunction<numint::vector_d,3> mdf;
	mdf.func = &Ftor_one::exec;
	mdf.param = &F;
	vector<double> ret(2,0.);
	F.f=OneGlauberGrid::klaas_one_bound_ct;
	res = numint::cube({integrator, 1.E-12, prec, 2000, 2000000}, mdf, lower, upper, ret, count, 0);
	results=ret[0];
	
      
      fsi_ct_grid[level][mm][i][j][k]=1.-getParticles()[0].getScatterfront(level,getPnucleus())*results;
//       cout << i << " " << j << " " << k << " " << level << " " << mm << fsi_grid[level][mm][i][j][k] << " " << fsi_ct_grid[level][mm][i][j][k] << 
//       " " << res << " " << count << endl;
    }
  }   
  
}








void OneGlauberGrid::klaas_one_bound(numint::vector_d & results, double b, double z, double phi, OneGlauberGrid & grid, 
		   int level, int m){
  //double sintheta=sqrt(1.-costheta*costheta);
  results=numint::vector_d(2,0.);
  double r=sqrt(b*b+z*z);
  if(r>grid.getPnucleus()->getRange()) return;
  double costheta=z/r;
  double sintheta=b/r;
  double cosphi, sinphi;
  sincos(phi,&sinphi,&cosphi);
  double zmom=grid.getParticles()[0].calcZ(r,costheta,sintheta,cosphi,sinphi);
  
  results[0]=(power(grid.getPnucleus()->getWave_F(level,r)/r,2)*grid.getPnucleus()->getYminkappacos(level,m,costheta)
	      +power(grid.getPnucleus()->getWave_G(level,r)/r,2)*grid.getPnucleus()->getYkappacos(level,m,costheta))
	      *4.*b*exp(-(power(b-grid.getParticles()[0].getHitbnorm(),2)
	      +b*grid.getParticles()[0].getHitbnorm()*4.*sinphi*sinphi)
	      /(2.*grid.getParticles()[0].getBetasq(level,grid.getPnucleus())));
  results[1]=results[0]*grid.getParticles()[0].getCTsigma(zmom);
  return;
  
}

void OneGlauberGrid::klaas_one_bound_ct(numint::vector_d & results, double b, double z, double phi, OneGlauberGrid & grid, 
		   int level, int m){
  //double sintheta=sqrt(1.-costheta*costheta);
  results=numint::vector_d(1,0.);
  double r=sqrt(b*b+z*z);
  if(r>grid.getPnucleus()->getRange()) return;
  double costheta=z/r;
  double sintheta=b/r;
  double cosphi, sinphi;
  sincos(phi,&sinphi,&cosphi);
  double zmom=grid.getParticles()[0].calcZ(r,costheta,sintheta,cosphi,sinphi);
  
  results[0]=(power(grid.getPnucleus()->getWave_F(level,r)/r,2)*grid.getPnucleus()->getYminkappacos(level,m,costheta)
	      +power(grid.getPnucleus()->getWave_G(level,r)/r,2)*grid.getPnucleus()->getYkappacos(level,m,costheta))
	      *4.*b*exp(-(power(b-grid.getParticles()[0].getHitbnorm(),2)
	      +b*grid.getParticles()[0].getHitbnorm()*4.*sinphi*sinphi)
	      /(2.*grid.getParticles()[0].getBetasq(level,grid.getPnucleus())))*grid.getParticles()[0].getCTsigma(zmom);
  return;
  
}

}  // namespace nucqcd




