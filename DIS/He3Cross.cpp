#include <DIS/He3Cross.hpp>
#include <fstream>
#include <Tools/Core/FourVector.hpp>
#include <Tools/Core/Vector3.hpp>

namespace nucqcd {


using namespace std;

He3Cross::He3Cross(std::string wfname, std::string inputdir, std::string strucname, bool proton_): 
  proton(proton_),
  wavefunction(wfname,inputdir+"/he3"),
  strucfunc(proton_,strucname){
}

He3Cross::~He3Cross(){
}

He3Cross::He3Cross(const He3Cross& rhs){
  proton=rhs.proton;
  wavefunction=rhs.wavefunction;
  strucfunc=rhs.strucfunc;
}

He3Cross& He3Cross::operator=(const He3Cross& rhs){
  if(this!=&rhs) { // avoid self-assignment
    proton=rhs.proton;
    wavefunction=rhs.wavefunction;
    strucfunc=rhs.strucfunc;
  }
  return *this;  
}



double He3Cross::getMomentumDistribution(const nucqcd::Vector3 &sp1, const nucqcd::Vector3 &sp2){
  // sum over the spins of all three nucleons (ms3, the struck one, independent of the others; until
  // 2026-11 it was fixed to mA-ms1-ms2, which left out part of the D and P waves) and over the
  // isospins of the spectators
  double result=0.;
  const int mA=-1;  // the density summed over the nucleon spins is the same for mA = +1 and -1
  for(int ms1=-1;ms1<=1;ms1+=2){
    for(int ms2=-1;ms2<=1;ms2+=2){
      for(int ms3=-1;ms3<=1;ms3+=2){
	for(int mt1=(proton?-1:1);mt1<=1;mt1+=2){
	  int mt2=1-mt1-(proton?1:-1);
	  result+=norm(wavefunction.getWF(sp1,sp2,ms1,ms2,ms3,mt1,mt2,mA));
	}
      }
    }
  }
  // A = 3: the He3wf wave function is normalized over all isospin arrangements, the note's per struck
  // nucleon (see the class description); missing until 2026-11
  return 3.*result;
}

double He3Cross::getDensity(LightConeKin2to3 &kin){
  nucqcd::Vector3 sp1(kin.getSp1_mu()[1]-kin.getA_mu()[1]/3.,
		kin.getSp1_mu()[2]-kin.getA_mu()[2]/3.,
		kin.getSp1_mu()[3]-kin.getA_mu()[3]/3.);
  nucqcd::Vector3 sp2(kin.getSp2_mu()[1]-kin.getA_mu()[1]/3.,
		kin.getSp2_mu()[2]-kin.getA_mu()[2]/3.,
		kin.getSp2_mu()[3]-kin.getA_mu()[3]/3.);
  return getMomentumDistribution(sp1,sp2)*sqrt((nucqcd::mass_nucleon*nucqcd::mass_nucleon+sp1.Mag2())*(nucqcd::mass_nucleon*nucqcd::mass_nucleon+sp2.Mag2())); //He3 spin average not needed, symmetry used in m_A
}

double He3Cross::getCross(LightConeKin2to3 &kin){
  return 2.*pow(kin.getYA()*nucqcd::alpha_em/kin.getQ2(),2.)*getDensity(kin)*strucfunc.getStructureLC(kin)
      *nucqcd::hbarc*nucqcd::hbarc*1.E25;  //conversion from MeV-8 to nb/GeV-6
}

}  // namespace nucqcd
