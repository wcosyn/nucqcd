/*
 * InterpolatingWavefunction.cpp
 *
 * Author: Pieter Vancraeyveld
 *         pieter.vancraeyveld@ugent.be
 *
 * Worked started February,16 2009
 *
 */

////////////////////////////////////////////////////////////////////////
//                                                                    //
// InterpolatingWavefunction                                         //
//                                                                    //
// InterpolatingWavefunction implements the abstract interface       //
// WavefunctionImplementation.                                       //
//                                                                    //
// This class calculates (non-)relativistic wavefunctions by          //
// interpolating (linear) between a list of tabulated values.         //
//                                                                    //
// WavefunctionImplementation fixes the following conventions        //
// concerning the units:                                              //
// * position in [fm]                                                 //
// * momentum in [MeV]                                                //
// * r-space wavefunctions in [fm^-1/2]                               //
// * p-space wavefunctions in [MeV^-3/2]                              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <Tools/WaveFunctions/InterpolatingWavefunction.hpp>
#include <iostream>
#include <Tools/WaveFunctions/Deuteron.hpp>

namespace nucqcd {


using std::cout; using std::cerr; using std::endl;
using std::map; using std::pair;


//_____________________________________________________________________
InterpolatingWavefunction::InterpolatingWavefunction()
: WavefunctionImplementation(),
  fUr(0), fWr(0), fVTr(0), fVSr(0), fUp(0), fWp(0), fVTp(0), fVSp(0), inclS(1), inclD(1)
{
  // Constructor
}

//_____________________________________________________________________
InterpolatingWavefunction::InterpolatingWavefunction(const InterpolatingWavefunction& rhs)
: WavefunctionImplementation(),
  fUr(0), fWr(0), fVTr(0), fVSr(0), fUp(0), fWp(0), fVTp(0), fVSp(0)
{
  // Copy constructor

  if(rhs.fUr) fUr = new map<double,double>(*rhs.fUr);
  if(rhs.fWr) fWr = new map<double,double>(*rhs.fWr);
  if(rhs.fVTr) fVTr = new map<double,double>(*rhs.fVTr);
  if(rhs.fVSr) fVSr = new map<double,double>(*rhs.fVSr);

  if(rhs.fUp) fUp = new map<double,double>(*rhs.fUp);
  if(rhs.fWp) fWp = new map<double,double>(*rhs.fWp);
  if(rhs.fVTp) fVTp = new map<double,double>(*rhs.fVTp);
  if(rhs.fVSp) fVSp = new map<double,double>(*rhs.fVSp);

  inclS=rhs.inclS;
  inclD=rhs.inclD;

}

//_____________________________________________________________________
InterpolatingWavefunction& InterpolatingWavefunction::operator=(const InterpolatingWavefunction& rhs)
{
  // Assignment

  if(this!=&rhs) { // avoid self-assignment

    if(rhs.fUr) {
      if(fUr) *fUr = *rhs.fUr;
      else fUr = new map<double,double>(*rhs.fUr);
    } else {
      delete fUr;
      fUr = 0;
    }

    if(rhs.fWr) {
      if(fWr) *fWr = *rhs.fWr;
      else fWr = new map<double,double>(*rhs.fWr);
    } else {
      delete fWr;
      fWr = 0;
    }

    if(rhs.fVTr) {
      if(fVTr) *fVTr = *rhs.fVTr;
      else fVTr = new map<double,double>(*rhs.fVTr);
    } else {
      delete fVTr;
      fVTr = 0;
    }

    if(rhs.fVSr) {
      if(fVSr) *fVSr = *rhs.fVSr;
      else fVSr = new map<double,double>(*rhs.fVSr);
    } else {
      delete fVSr;
      fVSr = 0;
    }

    if(rhs.fUp) {
      if(fUp) *fUp = *rhs.fUp;
      else fUp = new map<double,double>(*rhs.fUp);
    } else {
      delete fUp;
      fUp = 0;
    }

    if(rhs.fWp) {
      if(fWp) *fWp = *rhs.fWp;
      else fWp = new map<double,double>(*rhs.fWp);
    } else {
      delete fWp;
      fWp = 0;
    }

    if(rhs.fVTp) {
      if(fVTp) *fVTp = *rhs.fVTp;
      else fVTp = new map<double,double>(*rhs.fVTp);
    } else {
      delete fVTp;
      fVTp = 0;
    }

    if(rhs.fVSp) {
      if(fVSp) *fVSp = *rhs.fVSp;
      else fVSp = new map<double,double>(*rhs.fVSp);
    } else {
      delete fVSp;
      fVSp = 0;
    }
  }

  return *this;
}

//_____________________________________________________________________
InterpolatingWavefunction::~InterpolatingWavefunction()
{
  // Destructor

  delete fUr;
  delete fWr;
  delete fVTr;
  delete fVSr;
  
  delete fUp;
  delete fWp;
  delete fVTp;
  delete fVSp;
}

//_____________________________________________________________________
void InterpolatingWavefunction::RemoveUr()
{
  // Clear list of tabulated values
  if(fUr) fUr->clear();
}

//_____________________________________________________________________
void InterpolatingWavefunction::RemoveWr()
{
  // Clear list of tabulated values
  if(fWr) fWr->clear();
}

//_____________________________________________________________________
void InterpolatingWavefunction::RemoveVTr()
{
  // Clear list of tabulated values
  if(fVTr) fVTr->clear();
}

//_____________________________________________________________________
void InterpolatingWavefunction::RemoveVSr()
{
  // Clear list of tabulated values
  if(fVSr) fVSr->clear();
}

//_____________________________________________________________________
void InterpolatingWavefunction::RemoveUp()
{
  // Clear list of tabulated values
  if(fUp) fUp->clear();
}

//_____________________________________________________________________
void InterpolatingWavefunction::RemoveWp()
{
  // Clear list of tabulated values
  if(fWp) fWp->clear();
}

//_____________________________________________________________________
void InterpolatingWavefunction::RemoveVTp()
{
  // Clear list of tabulated values
  if(fVTp) fVTp->clear();
}

//_____________________________________________________________________
void InterpolatingWavefunction::RemoveVSp()
{
  // Clear list of tabulated values
  if(fVSp) fVSp->clear();
}

//_____________________________________________________________________
void InterpolatingWavefunction::SetUr(int n, double *r, double *list)
{
  // Add 'n' L=0 wavefunction values.
  // r in [fm]
  // list in [fm^-1/2]

  for(int i=0; i<n; ++i)
    AddUr(r[i],list[i]);  
}

//_____________________________________________________________________
void InterpolatingWavefunction::SetWr(int n, double *r, double *list)
{
  // Add 'n' L=2 wavefunction values.
  // r in [fm]
  // list in [fm^-1/2]

  for(int i=0; i<n; ++i)
    AddWr(r[i],list[i]);  
}

//_____________________________________________________________________
void InterpolatingWavefunction::SetVTr(int n, double *r, double *list)
{
  // Add 'n' L=1 triplet wavefunction values.
  // r in [fm]
  // list in [fm^-1/2]

  for(int i=0; i<n; ++i)
    AddVTr(r[i],list[i]);  
}

//_____________________________________________________________________
void InterpolatingWavefunction::SetVSr(int n, double *r, double *list)
{
  // Add 'n' L=1 singlet wavefunction values.
  // r in [fm]
  // list in [fm^-1/2]

  for(int i=0; i<n; ++i)
    AddVSr(r[i],list[i]);  
}

//_____________________________________________________________________
void InterpolatingWavefunction::SetUp(int n, double *p, double *list)
{
  // Add 'n' L=0 wavefunction values.
  // p in [MeV]
  // list in [MeV^-3/2]

  for(int i=0; i<n; ++i)
    AddUp(p[i],list[i]);  
}

//_____________________________________________________________________
void InterpolatingWavefunction::SetWp(int n, double *p, double *list)
{
  // Add 'n' L=2 wavefunction values.
  // p in [MeV]
  // list in [MeV^-3/2]

  for(int i=0; i<n; ++i)
    AddWp(p[i],list[i]);  
}

//_____________________________________________________________________
void InterpolatingWavefunction::SetVTp(int n, double *p, double *list)
{
  // Add 'n' L=1 triplet wavefunction values.
  // p in [MeV]
  // list in [MeV^-3/2]

  for(int i=0; i<n; ++i)
    AddVTp(p[i],list[i]);  
}

//_____________________________________________________________________
void InterpolatingWavefunction::SetVSp(int n, double *p, double *list)
{
  // Add 'n' L=1 singlet wavefunction values.
  // p in [MeV]
  // list in [MeV^-3/2]

  for(int i=0; i<n; ++i)
    AddVSp(p[i],list[i]);  
}

//_____________________________________________________________________
void InterpolatingWavefunction::AddUr(double r, double value)
{
  // Add L=0 wavefunction value.
  // r in [fm]
  // wavefunction in [fm^-1/2]

  if(!fUr) fUr = new map<double,double>;

  fUr->insert( pair<double,double>(r,value) );
}

//_____________________________________________________________________
void InterpolatingWavefunction::AddWr(double r, double value)
{
  // Add L=2 wavefunction value.
  // r in [fm]
  // wavefunction in [fm^-1/2]

  if(!fWr) fWr = new map<double,double>;

  fWr->insert( pair<double,double>(r,value) );
}

//_____________________________________________________________________
void InterpolatingWavefunction::AddVTr(double r, double value)
{
  // Add L=1 triplet wavefunction value.
  // r in [fm]
  // wavefunction in [fm^-1/2]

  if(!fVTr) fVTr = new map<double,double>;

  fVTr->insert( pair<double,double>(r,value) );
}

//_____________________________________________________________________
void InterpolatingWavefunction::AddVSr(double r, double value)
{
  // Add L=1 singlet wavefunction value.
  // r in [fm]
  // wavefunction in [fm^-1/2]

  if(!fVSr) fVSr = new map<double,double>;

  fVSr->insert( pair<double,double>(r,value) );
}

//_____________________________________________________________________
void InterpolatingWavefunction::AddUp(double p, double value)
{
  // Add L=0 wavefunction value.
  // p in [MeV]
  // wavefunction in [MeV^-3/2]

  if(!fUp) fUp = new map<double,double>;

  fUp->insert( pair<double,double>(p,value) );
}

//_____________________________________________________________________
void InterpolatingWavefunction::AddWp(double p, double value)
{
  // Add L=2 wavefunction value.
  // p in [MeV]
  // wavefunction in [MeV^-3/2]

  if(!fWp) fWp = new map<double,double>;

  fWp->insert( pair<double,double>(p,value) );
}

//_____________________________________________________________________
void InterpolatingWavefunction::AddVTp(double p, double value)
{
  // Add L=1 triplet wavefunction value.
  // p in [MeV]
  // wavefunction in [MeV^-3/2]

  if(!fVTp) fVTp = new map<double,double>;

  fVTp->insert( pair<double,double>(p,value) );
}

//_____________________________________________________________________
void InterpolatingWavefunction::AddVSp(double p, double value)
{
  // Add L=1 singlet wavefunction value.
  // p in [MeV]
  // wavefunction in [MeV^-3/2]

  if(!fVSp) fVSp = new map<double,double>;

  fVSp->insert( pair<double,double>(p,value) );
}

//_____________________________________________________________________
double InterpolatingWavefunction::GetUr(double r) const
{
  // L=0 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  return Interpolate(fUr,r);
}

//_____________________________________________________________________
double InterpolatingWavefunction::GetWr(double r) const
{
  // L=2 wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  return Interpolate(fWr,r);
}

//_____________________________________________________________________
double InterpolatingWavefunction::GetVTr(double r) const
{
  // L=1 triplet wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  return Interpolate(fVTr,r);
}

//_____________________________________________________________________
double InterpolatingWavefunction::GetVSr(double r) const
{
  // L=1 singlet wave in configuration space
  // r = relative position in [fm]
  // Return value in units [fm^-1/2]

  return Interpolate(fVSr,r);
}

//_____________________________________________________________________
double InterpolatingWavefunction::GetUp(double p) const
{
  // L=0 wave in momentum space
  // p = relative position in [MeV]
  // Return value in units [MeV^-3/2]

  return Interpolate(fUp,p);
}

//_____________________________________________________________________
double InterpolatingWavefunction::GetWp(double p) const
{
  // L=2 wave in momentum space
  // p = relative position in [MeV]
  // Return value in units [MeV^-3/2]

  return Interpolate(fWp,p);
}

//_____________________________________________________________________
double InterpolatingWavefunction::GetVTp(double p) const
{
  // L=1 triplet wave in momentum space
  // p = relative position in [MeV]
  // Return value in units [MeV^-3/2]

  return Interpolate(fVTp,p);
}

//_____________________________________________________________________
double InterpolatingWavefunction::GetVSp(double p) const
{
  // L=1 singlet wave in momentum space
  // p = relative position in [MeV]
  // Return value in units [MeV^-3/2]

  return Interpolate(fVSp,p);
}

double InterpolatingWavefunction::GetUpoff(const nucqcd::Vector3& p) const{
 
  cerr << "Off-shell wave function not implemented in Interpolating wave function class" << endl;
  return 0.;
  
}

double InterpolatingWavefunction::GetWpoff1(const nucqcd::Vector3& p) const{
 
  cerr << "Off-shell wave function not implemented in Interpolating wave function class" << endl;
  return 0.;
  
}

double InterpolatingWavefunction::GetWpoff2(double pperp2) const{
 
  cerr << "Off-shell wave function not implemented in Interpolating wave function class" << endl;
  return 0.;
  
}


//_____________________________________________________________________
double InterpolatingWavefunction::Interpolate(map<double,double>* list,
					       double q) const
{
  // Find the wavefunction value corresponding to 'q', by interpolating
  // between the known values found in the ordered 'list'
  if(std::isnan(q)) return 0.;
  if( !list || list->size()<2 ) {
    cerr << "WARNING in InterpolatingWavefunction::Interpolate(map*,double): "
	 << "nothing to interpolate!\n";
    return 0;
  }

  // Find point v2=(q2,wf2) such that q < q2
  map<double,double>::iterator v2 = list->upper_bound(q);
  
  // check whether it is within bounds
  if( v2==list->end() ) {
    --v2; // move iterator one back
    if( v2->first!=q ) {
      cerr << q << " WARNING in InterpolatingWavefunction::Interpolate(map*,double): "
	   << "requested point outside range of tabulated values.\n";
      return 0;
    }
  }
  
  if( v2==list->begin() ) {
    cerr << q << " WARNING in InterpolatingWavefunction::Interpolate(map*,double): "
	 << "requested point outside range of tabulated values.\n";
    return 0;
  }
  
  // Find point v1=(q1,wf1) such that q1 <= q
  map<double,double>::iterator v1 = v2;
  --v1;

  // Now interpolation is children's play
  return v1->second + (q-v1->first) 
    * (v2->second-v1->second) / (v2->first-v1->first);
}



std::complex<double> InterpolatingWavefunction::DeuteronPState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const nucqcd::Vector3& p) const
{
  // Calculate the state of a deuteron in momentum space. 
  // deuteronPol is the deuteron polarization (times two)
  // nucleon1pol fixes the polarization of the first nucleon (the one that interacts with the photon). (times two)
  //nucleon2pol fixes spectator nucleon polarization (times two)
  
  // The Dnp vertex matrix element is
  // BEGIN_LATEX
  // #frac{1}{2E_{1}} #bar{u}(p_{2},m_{2})u(#dot{p}_{1},m_{1})#Gamma_{#mu}#phi^{#mu}_{D}(m_{D})
  // = #sum_{L} u_{L}(|#vec{p}_{1}|) #sum_{m_{L},m_{S}} <Lm_{L}Sm_{S}|1m_{D}><#frac{1}{2}m_{1}#frac{1}{2} m_{2}|1m_{S}> Y_{Lm_{L}}(#hat{p}_{1})
  // END_LATEX
  // with BEGIN_LATEX #dot{p}_{1} = ( #sqrt{#vec{p}_{1}^{2}+M_{1}^{2}},#vec{p}_{1}) END_LATEX the on-shell 4-momentum of particle 1.
  // The labelling particle 1 and 2 is arbitrary in the non-relativistic case (You can do the math).

  
  std::complex<double> state=0.;
// L=2, mL=mD-mS, S=1, mS=m1+m2
  // for(int mS=-1; mS<=1; ++mS)
  if(inclD)
    state -=
      Deuteron::Wavefunction::ClebschGordan(4,deuteronPol-nucleon1Pol-nucleon2Pol,2,nucleon1Pol+nucleon2Pol,2,deuteronPol)
      *Deuteron::Wavefunction::ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,nucleon1Pol+nucleon2Pol)
      *Deuteron::Wavefunction::SphericalHarmonicCos(2,(deuteronPol-nucleon1Pol-nucleon2Pol)/2,p.CosTheta(),p.Phi());
  state*=GetWp(p.Mag());
  // L=0, mL=0, S=1, mS=mD
  if(inclS)
  state += 
    Deuteron::Wavefunction::ClebschGordan(0,0,2,deuteronPol,2,deuteronPol)
    *Deuteron::Wavefunction::ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,deuteronPol)
    *Deuteron::Wavefunction::SphericalHarmonicCos(0,0,p.CosTheta(),p.Phi())
    *GetUp(p.Mag());

  return state;
}


std::complex<double> InterpolatingWavefunction::DeuteronRState(int deuteronPol, int nucleon2Pol,
						       int nucleon1Pol, 
						       const nucqcd::Vector3& r) const
{
  std::complex<double> state=0.;
// L=2, mL=mD-mS, S=1, mS=m1+m2
  // for(int mS=-1; mS<=1; ++mS)
  if(inclD)
    state +=
      Deuteron::Wavefunction::ClebschGordan(4,deuteronPol-nucleon1Pol-nucleon2Pol,2,nucleon1Pol+nucleon2Pol,2,deuteronPol)
      *Deuteron::Wavefunction::ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,nucleon1Pol+nucleon2Pol)
      *Deuteron::Wavefunction::SphericalHarmonicCos(2,(deuteronPol-nucleon1Pol-nucleon2Pol)/2,r.CosTheta(),r.Phi());
  state*=GetWr(r.Mag());
  // L=0, mL=0, S=1, mS=mD
  if(inclS)
  state += 
    Deuteron::Wavefunction::ClebschGordan(0,0,2,deuteronPol,2,deuteronPol)
    *Deuteron::Wavefunction::ClebschGordan(1,nucleon1Pol,1,nucleon2Pol,2,deuteronPol)
    *Deuteron::Wavefunction::SphericalHarmonicCos(0,0,r.CosTheta(),r.Phi())
    *GetUr(r.Mag());

  return state;
}

}  // namespace nucqcd

