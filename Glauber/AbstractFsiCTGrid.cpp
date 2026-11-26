#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include <Glauber/AbstractFsiCTGrid.hpp>
#include <Tools/Core/constants.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <Tools/Core/GridCache.hpp>

namespace nucqcd {


AbstractFsiCTGrid::AbstractFsiCTGrid(int r_grid, int cth_grid, int phi_grid, MeanFieldNucleus *pnucl, 
				     double prec, IntegrationMethod integrator, string homedir):
AbstractFsiGrid(r_grid,cth_grid,phi_grid,pnucl,prec,integrator,homedir),filledctgrid(0){
  number_of_grids=2;
}  
  
AbstractFsiCTGrid::~AbstractFsiCTGrid(){
  //cout << "Deleting FSI object" << endl;
  
}

//get corrfilename
const string AbstractFsiCTGrid::getFsi_Ct_Filename() const{
  return fsi_ct_filename;
}


//set filenames
void AbstractFsiCTGrid::setFilenames(string cachedir){
  AbstractFsiGrid::setFilenames(cachedir);
  //cout << "bla " << getFsi_Filename() << endl;
  fsi_ct_filename=getFsi_Filename();
  fsi_ct_filename.insert(fsi_ct_filename.size()-4,".CT");
  for(size_t i=0;i<getParticles().size();i++){
    // bin centres, see FastParticle::gridBinCentre
    fsi_ct_filename.insert(fsi_ct_filename.size()-4,("."+to_string(lround(getParticles()[i].getHardScale()*10))
                                                    +"."+to_string(lround(getParticles()[i].get_lc_mod()*100))
                                                    +"."+to_string(lround(getParticles()[i].get_nkt_mod()*100))));
  }
  //cout << fsi_ct_filename << endl;
}
  
//reads the fsi+ct grid from its cache file if that exists and is complete, otherwise computes and writes it
void AbstractFsiCTGrid::readinOrConstructCtGrid(){
  ifstream infile2(fsi_ct_filename.c_str(),ios::in|ios::binary);
  //check if object has been created sometime earlier and read it in
  if(infile2.is_open()){
    // cout << "Reading in FSI+CT grid from memory: " << fsi_ct_filename << endl;
    readinFsiCtGrid(infile2);
    bool read=!infile2.fail(); // a truncated file sets failbit
    infile2.close();
    if(read){
      filledctgrid=filledallgrid=1;
      return;
    }
    cerr << "FSI+CT grid file is incomplete, recomputing it: " << fsi_ct_filename << endl;
  }
  //we have to calc the grid
//   cout << "Constructing FSI+CT grid" << endl;
  constructCtGrid(); // ! pure virtual
  filledctgrid=filledallgrid=1;
  //write it out now too
  nucqcd::write_file_atomically(fsi_ct_filename,[this](ofstream &outfile){ writeoutFsiCtGrid(outfile); });
}

//fills the grids, uses pure virtual functions!!!!
void AbstractFsiCTGrid::fillGrids(){
  AbstractFsiGrid::fillGrids(); // fill regular FSI (RMSGA) grid, this function is implemented in AbstractFsiGrid but it uses pure virtual functions
  setFilenames(getDir());
  if(filledallgrid){
    //ct grid was filled while constructing RMSGA grid
    filledctgrid=1;
    //write to file
    nucqcd::write_file_atomically(fsi_ct_filename,[this](ofstream &outfile){ writeoutFsiCtGrid(outfile); });
  }
  //grid is still empty
  else readinOrConstructCtGrid();
}

//updates the grids, uses pure virtual functions!!!!
void AbstractFsiCTGrid::updateGrids(){
  AbstractFsiGrid::updateGrids(); // implemented virtual function in AbstractFsiGrid but uses pure virtual functions.
                                    // the class that will eventually call this will have inherited the virtual function
                                    // updateGrids(), which uses the [former] pure virtual functions that are now implemented
                                    // within this class and hence will be called, as virtual functions are called using the 
                                    // object rather then handle type
  string old_fsi_ct_filename = fsi_ct_filename;
  setFilenames(getDir());
  //check is something is different so we needto update the grids
  if(old_fsi_ct_filename.compare(fsi_ct_filename)){
//     cout << "fsict grid not equal " << endl << fsi_ct_filename << endl << old_fsi_ct_filename << endl;
    readinOrConstructCtGrid();
  }
//   else cout << "fsict grid equal to the earlier one, doing nothing" << endl << fsi_ct_filename << endl << old_fsi_ct_filename << endl;
  
}

}  // namespace nucqcd

