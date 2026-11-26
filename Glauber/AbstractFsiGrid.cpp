#include <Tools/Core/IntegratorFlag.hpp>
#include <cstdlib>
#include <Tools/Core/ErrorMessage.hpp>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include <Glauber/AbstractFsiGrid.hpp>
#include <Tools/Core/constants.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/GridCache.hpp>

namespace nucqcd {


AbstractFsiGrid::AbstractFsiGrid(int r_grid, int cth_grid, int phi_grid, 
				 MeanFieldNucleus *pnucl, double precision, IntegrationMethod integr, std::string homedir):
filledgrid(0),filledallgrid(0),dir(nucqcd::grid_cache_dir(homedir)),rgrid(r_grid),cthgrid(cth_grid),
phigrid(phi_grid),allinplane(0),totalprotonout(0), totalneutronout(0),
protonknockout(0), neutronknockout(0),pnucleus(pnucl),number_of_grids(1),prec(precision),integrator(integr){
  
  invrstep=rgrid/pnucleus->getRange();
  invcthstep=cthgrid/2.;
  invphistep=0.5*phigrid/PI;
  
  
}  
  
AbstractFsiGrid::~AbstractFsiGrid(){
  //cout << "Deleting FSI object" << endl;
  
}
  
//add particle to the particles vector
void AbstractFsiGrid::addParticle(FastParticle& particle)
{
  // the grid is computed for the centre of the cache bin of the particle, see FastParticle::gridBinCentre
  FastParticle newparticle=particle.gridBinCentre();
  if (particles.empty()){
    particles.push_back(newparticle);
    if(newparticle.getParticletype()==0){
      if(newparticle.getIncoming()) totalprotonout--;
      else totalprotonout++;
    }
    if(newparticle.getParticletype()==1){
      if(newparticle.getIncoming()) totalneutronout--;
      else totalneutronout++;
    }
  } else {
    if(newparticle.getIncoming()){
      if(particles[0].getIncoming()){
        cerr << "Not more than one incoming particle possible in the particle list!!!" << endl <<
        "Replacing old incoming particle with this one" << endl;
        if(particles[0].getParticletype()==0) totalprotonout++;
        if(particles[0].getParticletype()==1) totalneutronout++;
        particles[0]=newparticle;
        if(particles[0].getParticletype()==0) totalprotonout--;
        if(particles[0].getParticletype()==1) totalneutronout--;
      }
      else particles.insert(particles.begin(),newparticle);
      if(particles[0].getParticletype()==0) totalprotonout--;
      if(particles[0].getParticletype()==1) totalneutronout--;
    }
    else particles.push_back(newparticle);
    if(newparticle.getParticletype()==0) totalprotonout--;
    if(newparticle.getParticletype()==1) totalneutronout--;
  }
}

//print the particles in the particles vector
void AbstractFsiGrid::clearParticles(){
 particles.clear();
 totalprotonout=0;
 totalneutronout=0;
 return;
}


//print the particles in the particles vector
void AbstractFsiGrid::printParticles() const{
  //cout << endl << endl << "Printing contents of particle vector in FSI grid" << endl;
  for(size_t i=0;i<particles.size();i++){
    cout <<"Particle " << i+1 << endl;
    particles[i].printParticle();
  }
}

//add a nucleon that gets knocked out
void AbstractFsiGrid::addKnockout(const int level, const int m){
  if(level<0||level>pnucleus->getTotalLevels()){
    throw std::invalid_argument(nucqcd::ErrorMessage() << "invalid knockout level: " << level );
  }
  if(abs(m)>pnucleus->getJ_array()[level]){
    throw std::invalid_argument(nucqcd::ErrorMessage() << "invalid m level: " <<pnucleus->getJ_array()[level] << " " << m );
  }
  for(size_t i=0;i<knockoutlevels.size();i++){
    if((level==knockoutlevels[i])&&(m==knockoutm[i])) {
      throw std::logic_error(nucqcd::ErrorMessage() << "knockoutlevel is already present!!!" << "\n" << i+1 << " level: " << knockoutlevels[i] << ", m: " << knockoutm[i] );
    }
  }
  if(level<pnucleus->getPLevels()) protonknockout++;
  else neutronknockout++;
  knockoutlevels.push_back(level);
  knockoutm.push_back(m);
  return;
}

//empty knockout vectors
void AbstractFsiGrid::clearKnockout(){
 knockoutlevels.clear();
 knockoutm.clear();
 protonknockout=0;
 neutronknockout=0;
 return;
}

//print knockout vectors
void AbstractFsiGrid::printKnockout() const{
  cout << endl << endl << "Printing contents of knockout vector in FSI grid" << endl;
  for(size_t i=0;i<knockoutlevels.size();i++){
    cout <<"knockout " << i+1 << endl << "level: " << knockoutlevels[i] << ", m: " << knockoutm[i] << endl;
    
  }
}


//set up inteprolation help variables for r
void AbstractFsiGrid::setRinterp(const double r){
  if(r>pnucleus->getRange()){
    throw std::out_of_range(nucqcd::ErrorMessage() << "r out of range in : AbstractFsiGrid::setRinterp" << r );
  }
  rindex = int(floor(r*getInvRstep()));
  if(rindex==getRgrid()) rindex-=1; //take care at high end
  s_interp = (r*getInvRstep() - (rindex));
  comp_s_interp=1.-s_interp;
}

//set up interpolation help variables for theta
void AbstractFsiGrid::setCthinterp(const double costheta){
  cthindex = int(floor((-costheta+1.)*getInvCthstep()));
  if(cthindex==getCthgrid()) cthindex-=1;
  t_interp = ((-costheta+1.)*getInvCthstep() - (cthindex));  
  comp_t_interp=1.-t_interp;
}

//set up interpolation help variables for phi
void AbstractFsiGrid::setPhiinterp(const double phi){
  phiindex = int(floor(phi*getInvPhistep()));
  if(phiindex==getPhigrid()) phiindex-=1;
  u_interp = (phi*getInvPhistep() - (phiindex));  
  comp_u_interp=1.-u_interp;
}

void AbstractFsiGrid::setinterp(const double r, const double costheta, const double phi)
{
  setRinterp(r);
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
}

void AbstractFsiGrid::setinterp(const double costheta, const double phi)
{
  setCthinterp(costheta);
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
}

void AbstractFsiGrid::setinterp(const double phi)
{
  if(getAllinplane()&&phi>PI) setPhiinterp(2.*PI-phi);
  else setPhiinterp(phi);
}

//interpolation functions
int AbstractFsiGrid::variantNumber(GridVariant variant) const{
  const int number=getGridNumber(variant);
  if(number<0){
    static const char *names[]={"FSI","FSI+SRC","FSI+decay","FSI+SRC+decay","FSI+CT","FSI+SRC+CT","FSI+CT+decay","FSI+SRC+CT+decay"};
    throw std::logic_error(nucqcd::ErrorMessage() << "this grid class has no " << names[static_cast<int>(variant)] << " grid");
  }
  return number;
}

complex<double> AbstractFsiGrid::getFsiGridN_interpvec(const int grid, const nucqcd::Vector3 &rvec){
  return getFsiGridN_interp3(grid, rvec.Mag(),rvec.CosTheta(),rvec.Phi());
}

complex<double> AbstractFsiGrid::getFsiGridN_interp3(const int grid, const double r, const double costheta, const double phi){
  setinterp(r,costheta,phi);
  return getFsiGridN_interp(grid);
}

complex<double> AbstractFsiGrid::getFsiGridN_interp2(const int grid, const double costheta, const double phi){
  setinterp(costheta,phi);
  return getFsiGridN_interp(grid);
}

complex<double> AbstractFsiGrid::getFsiGridN_interp1(const int grid, const double phi){
  setinterp(phi);
  return getFsiGridN_interp(grid);
}


//3d interpolation of array
complex<double> AbstractFsiGrid::getInterp(complex<double> ***grid){
//   return getComp_u_interp()*(getComp_t_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()][getPhiindex()]
// 	+ getS_interp()*grid[getRindex()+1][getCthindex()][getPhiindex()]) 
// 	+ getT_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()+1][getPhiindex()] 
// 	+ getS_interp()*grid[getRindex()+1][getCthindex()+1][getPhiindex()]))
// 	+  getU_interp()*(getComp_t_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()][getPhiindex()+1]
// 	+ getS_interp()*grid[getRindex()+1][getCthindex()][getPhiindex()+1]) 
// 	+ getT_interp()*(getComp_s_interp()*grid[getRindex()][getCthindex()+1][getPhiindex()+1]
// 	+ getS_interp()*grid[getRindex()+1][getCthindex()+1][getPhiindex()+1]));

//way faster with doubles!
  return Interp3d(grid, getS_interp(), getT_interp(), getU_interp(), 
		  getComp_s_interp(), getComp_t_interp(), getComp_u_interp(), 
		  getRindex(), getCthindex(), getPhiindex());
}




/**set filenames
 *
 * The file name is the cache key. It contains the particle properties the grid depends on, on the
 * lattice of the bins of FastParticle::gridBinCentre (the particles were moved to their bin
 * centres in addParticle, so the grid in a file is the same whichever kinematics computed it):
 * the momentum in units of 10 MeV, theta and phi in units of 0.1 rad, and sigma in mb; for
 * particles with user-set scattering parameters also beta^2 and epsilon in units of 0.01. Then the
 * grid sizes, the integrator and its precision, and the version of the key. \n
 * Key versions: none (until 2026-11): p truncated to 10 MeV, grid computed at the first
 * kinematics in the bin, no momentum for user-set parameters; key2: bin centres, cos theta in
 * units of 0.1 (".th"); key3: theta in units of 0.1 rad (".tht").
 */
void AbstractFsiGrid::setFilenames(string cachedir){
  fsi_filename=cachedir+pnucleus->getNucleusName();
  //add all particle properties
  for(size_t i=0;i<particles.size();i++){
    const FastParticle &part=particles[i];
    if(part.getIncoming()) fsi_filename+=".inc";
    fsi_filename+=".Ptype"+to_string(part.getParticletype())+".p"+to_string(lround(part.getP()/10.))
		  +".tht"+to_string(lround(part.getTheta()*10))+".ph"+to_string(lround(part.getPhi()*10))
		  +".sigma"+to_string(lround(part.getSigmap()*10.));
    if(part.getUserset()) fsi_filename+=".beta"+to_string(lround(part.getBeta2p()*100.))
		  +".eps"+to_string(lround(part.getEpsilonp()*100.));
  }
  fsi_filename+=".r"+to_string(getRgrid())+".cth"+to_string(getCthgrid())
		+".phi"+to_string(getPhigrid())+".prec"+to_string(getPrec()*1.E05)+".intgr"+to_string(numint::method_flag(getIntegrator()))+".key3.dat";
//     cout << fsi_filename << endl;
}
  
  
  
//fills the grids, uses pure virtual functions!!!!
void AbstractFsiGrid::fillGrids(){
  invphistep=0.5*phigrid/PI;
  for(size_t i=0;i<getParticles().size();i++){
    allinplane=1;
    if((getParticles()[i].getPhi()!=0.)&&(getParticles()[i].getPhi()!=PI)){
      allinplane=0;
      break;
    }
  }
  if(allinplane){
    //cout <<"All particles lie in the phi=0 or phi=Pi plane, special case activated" << endl;
    //exploit symmetry in phi
    invphistep*=2.;
  }
  setFilenames(dir);
  ifstream infile(fsi_filename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  bool read=false;
  if(infile.is_open()){
    //cout << "Reading in FSI grid from memory: " << fsi_filename << endl;
    readinFsiGrid(infile); // ! pure virtual
    read=!infile.fail(); // a truncated file sets failbit
    infile.close();
    // only the FSI grid was read: the other grids (FSI+CT) in memory may belong to earlier particles,
    // so filledallgrid is cleared and AbstractFsiCTGrid::fillGrids reads or computes them
    if(read){ filledgrid=1; filledallgrid=0; }
    else cerr << "FSI grid file is incomplete, recomputing it: " << fsi_filename << endl;
  }
  if(!read){
//     cout << "Constructing all grids" << endl;
    constructAllGrids(); // ! pure virtual
    filledgrid=filledallgrid=1;
    nucqcd::write_file_atomically(fsi_filename,[this](ofstream &outfile){ writeoutFsiGrid(outfile); }); // ! pure virtual
  } 
}
  
//fills the grids, uses pure virtual functions!!!!
void AbstractFsiGrid::updateGrids(){
  if(!filledgrid){
    fillGrids();
//     cout << "Grids still empty " << endl;
  }
  else{
    string old_fsi_filename=fsi_filename;
    setFilenames(dir);
    if(fsi_filename.compare(old_fsi_filename)){
//       cout << "fsi grid not equal" << fsi_filename << endl << old_fsi_filename << endl;
      fillGrids();      
    }
//     else cout << "fsi grid equal to the earlier one, doing nothing" << endl << fsi_filename << endl << old_fsi_filename << endl;
  }
}

}  // namespace nucqcd

