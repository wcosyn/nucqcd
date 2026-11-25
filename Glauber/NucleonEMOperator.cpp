#include <Glauber/NucleonEMOperator.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>

namespace nucqcd {

using namespace std;

const FourVector<GammaStructure> NucleonEMOperator::gamma_mu=FourVector<GammaStructure>(GammaStructure(0.,0.,1.),
											GammaStructure(0.,0.,0.,1.),
				      GammaStructure(0.,0.,0.,0.,1.),GammaStructure(0.,0.,0.,0.,0.,1.));
const GammaStructure NucleonEMOperator::Id=GammaStructure(1.);


NucleonEMOperator::NucleonEMOperator(){
}

NucleonEMOperator::NucleonEMOperator(const double q2, const bool prot, const int para):
proton(prot),parametrization(para),Q2(q2){
  tau=Q2/(4.*(proton? nucqcd::mass_p*nucqcd::mass_p:nucqcd::mass_n*nucqcd::mass_n));
  GE_null=proton? 1.:1.91304;
  GM_null=proton? 2.79:-1.91304;
  setGM(); //GM first because GE can depend on GM for Q^2 above 6 GeV^2
  setGE();
  setF1();
  setF2();
  Q2index=int(floor(Q2*1.E-06/0.05));
  if(Q2index>=60) Q2index=59;
  Q2_interp = (Q2*1.E-06/0.05 - (Q2index));
}

NucleonEMOperator::~NucleonEMOperator(){
}

NucleonEMOperator::NucleonEMOperator(const NucleonEMOperator &rhs){
  Q2=rhs.Q2;
  proton=rhs.proton;
  parametrization=rhs.parametrization;
  tau=rhs.tau;
  GE_null=rhs.GE_null;
  GM_null=rhs.GM_null;
  GE=rhs.GE;
  GM=rhs.GM;
  F1=rhs.F1;
  F2=rhs.F2;
  Q2index=rhs.Q2index;
  Q2_interp=rhs.Q2_interp;
}

NucleonEMOperator & NucleonEMOperator::operator=(const NucleonEMOperator &rhs){
  if(this!=&rhs) { // avoid self-assignment
    Q2=rhs.Q2;
    proton=rhs.proton;
    parametrization=rhs.parametrization;
    tau=rhs.tau;
    GE_null=rhs.GE_null;
    GM_null=rhs.GM_null;
    GE=rhs.GE;
    GM=rhs.GM;
    F1=rhs.F1;
    F2=rhs.F2;
    Q2index=rhs.Q2index;
    Q2_interp=rhs.Q2_interp;
  }
  return *this;
  
}


FourVector<GammaStructure> NucleonEMOperator::getCC1(const FourVector<double> &pi, const FourVector<double> &pf, 
						     const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  if(medium) return getGM(r,medium,nucleus)*gamma_mu-getF2(r,medium,nucleus)/(2.*(proton?nucqcd::mass_p:nucqcd::mass_n))*(pi+pf)*Id;
  else return getGM()*gamma_mu-getF2()/(2.*(proton?nucqcd::mass_p:nucqcd::mass_n))*(pi+pf)*Id;
}
FourVector<GammaStructure> NucleonEMOperator::getCC2(const FourVector<double> &q, const double r, const int medium,
						     const MeanFieldNucleusThick &nucleus) const{
  if(medium) return gamma_mu*getF1(r,medium,nucleus)
    +getF2(r,medium,nucleus)/(4.*(proton?nucqcd::mass_p:nucqcd::mass_n))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
  else return gamma_mu*getF1()+getF2()/(4.*(proton?nucqcd::mass_p:nucqcd::mass_n))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
  
}
FourVector<GammaStructure> NucleonEMOperator::getCC3(const FourVector<double> &q, const FourVector<double> &pi,
						     const FourVector<double> &pf, const double r, const int medium,
						     const MeanFieldNucleusThick &nucleus) const{
  if(medium) return getF1(r,medium,nucleus)/(2.*(proton?nucqcd::mass_p:nucqcd::mass_n))*(pi+pf)*Id 
    + getGM(r,medium,nucleus)/(4.*(proton?nucqcd::mass_p:nucqcd::mass_n))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q)); 
  else return getF1()/(2.*(proton?nucqcd::mass_p:nucqcd::mass_n))*(pi+pf)*Id 
    + getGM()/(4.*(proton?nucqcd::mass_p:nucqcd::mass_n))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q)); 
}

FourVector<GammaStructure> NucleonEMOperator::getCC(const int current, const FourVector<double>& q, 
						    const FourVector<double>& pi, const FourVector<double>& pf, 
						    const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
switch(current){
    case(1):
      return getCC1(pi,pf,r,medium,nucleus);
      break;
    case(2):
      return getCC2(q,r,medium,nucleus);
      break;
    case(3):
      return getCC3(q,pi,pf,r,medium,nucleus);
      break;
    default:
      throw std::invalid_argument(nucqcd::ErrorMessage() << "Current operator not supported " << current );
  }
}


FourVector<GammaStructure> NucleonEMOperator::getCC1(const FourVector<double> &pi, const FourVector<double> &pf) const{
  return getGM()*gamma_mu-getF2()/(2.*(proton?nucqcd::mass_p:nucqcd::mass_n))*(pi+pf)*Id;
}
FourVector<GammaStructure> NucleonEMOperator::getCC2(const FourVector<double> &q) const{
  return gamma_mu*getF1()+getF2()/(4.*(proton?nucqcd::mass_p:nucqcd::mass_n))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q));
  
}
FourVector<GammaStructure> NucleonEMOperator::getCC3(const FourVector<double> &q, const FourVector<double> &pi,
						     const FourVector<double> &pf) const{
  return getF1()/(2.*(proton?nucqcd::mass_p:nucqcd::mass_n))*(pi+pf)*Id 
    + getGM()/(4.*(proton?nucqcd::mass_p:nucqcd::mass_n))*((gamma_mu*q)*gamma_mu-gamma_mu*(gamma_mu*q)); 
}

FourVector<GammaStructure> NucleonEMOperator::getCC(const int current, const FourVector<double>& q, 
						    const FourVector<double>& pi, const FourVector<double>& pf) const{
switch(current){
    case(1):
      return getCC1(pi,pf);
      break;
    case(2):
      return getCC2(q);
      break;
    case(3):
      return getCC3(q,pi,pf);
      break;
    default:
      throw std::invalid_argument(nucqcd::ErrorMessage() << "Current operator not supported " << current );
  }
}



void NucleonEMOperator::setGE(){
  if(parametrization==0){ //BBA
    if(!proton) GE=GE_null*pow(1+Q2/.71E06,-2)*tau*0.942/(1+4.61*tau); //neutron
    else{
      if(Q2<6.*1.E06) GE=GE_null/(1. + 3.253 * Q2*1.E-06 + 1.422 * pow(Q2*1.E-06,2) + 0.08582 * pow(Q2*1.E-06,3) 
	    + 0.3318 * pow(Q2*1.E-06,4) +  (-0.09371)*pow(Q2*1.E-06,5) + 
	    0.01076*pow(Q2*1.E-06,6));
      else{ //above 6GeV^2 scale like GM
	NucleonEMOperator form6(6.,proton,parametrization);
	GE=GM/form6.getGM()*form6.getGE();
      }
    }
    return;
  }
  else if(parametrization==1){ //SLAC dipole
    GE=proton? Get_Gdipole(Q2)*GE_null : GE_null*Get_Gdipole(Q2)*tau/(1+5.6*tau);
    return;
  }
  else if(parametrization==2){ //JJ Kelly
    GE = proton? GE_null*(1.-0.24*tau)/(1.+10.98*tau+12.82*tau*tau+21.97*pow(tau,3.)):
          1.7*tau/(1.+3.3*tau)*Get_Gdipole(Q2);                                  
    return;
  }
  else if(parametrization==3){ //BBBA 
    GE = proton? GE_null*(1.-0.0578*tau)/(1.+11.1*tau+13.6*tau*tau+33.0*pow(tau,3.)):
    (1.25*tau+1.3*tau*tau)/(1.-9.86*tau+305.*tau*tau-758.0*pow(tau,3.)+802.0*pow(tau,4.));
    return;
  }
}

void NucleonEMOperator::setGM(){
  if(parametrization==0){ //BBA
     GM=proton? GM_null/(1. + 3.104 *Q2*1.E-06 + 1.428 * pow(Q2*1.E-06,2) + 0.1112 * pow(Q2*1.E-06,3)
	       +(-0.006981)*pow(Q2*1.E-06,4) + 0.0003705*pow(Q2*1.E-06,5) + 
	       (-0.000007063)*pow(Q2*1.E-06,6)):
	GM_null/(1. + 3.043 *Q2*1.E-06 + 0.8548 * pow(Q2*1.E-06,2) +
		0.6806 * pow(Q2*1.E-06,3)
		+(-0.1287)*pow(Q2*1.E-06,4) + 0.008912*pow(Q2*1.E-06,5));
     return;
  }
  else if(parametrization==1){ //SLAC dipole 
    GM=proton? GM_null*Get_Gdipole(Q2) : GM_null*Get_Gdipole(Q2);
    return;
  }
  else if(parametrization==2){ //JJ Kelly
    GM = proton? GM_null*(1.+0.12*tau)/(1.+10.97*tau+18.86*tau*tau+6.55*pow(tau,3.)):
	  GM_null*(1.+2.33*tau)/(1.+14.72*tau+24.20*tau*tau+84.1*pow(tau,3.));
    return;
  }
  else if(parametrization==3){ //BBBA 
    GM = proton? GM_null*(1.+0.15*tau)/(1.+11.1*tau+19.6*tau*tau+7.54*pow(tau,3.)):
      GM_null*(1.+1.81*tau)/(1.+14.1*tau+20.7*tau*tau+68.7*pow(tau,3.));

    return;
  }
}


double NucleonEMOperator::Get_Gdipole(double Q2)  
{
  return pow(1+Q2/0.71E06,-2);
}


void NucleonEMOperator::setF1(){
  F1=(tau*GM+GE)/(1.+tau);
}
void NucleonEMOperator::setF2(){
  F2=(GM-GE)/(1.+tau);
}


double NucleonEMOperator::getF1(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  return (tau*getGM(r,medium,nucleus)+getGE(r,medium,nucleus))/(1.+tau);
}

double NucleonEMOperator::getF2(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  return (getGM(r,medium,nucleus)-getGE(r,medium,nucleus))/(1.+tau);
}

double NucleonEMOperator::getEmod(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  double dens=nucleus.getDensity(r,proton)/(r*r)*nucleus.getA();
  if(r==0.) dens=nucleus.getTotalDensity(1.E-02)/(1.E-04)*nucleus.getA();
  int densindex = int(floor(dens/0.016));
  if(densindex>=(medium==1?14:15)) densindex=(medium==1?14:15); //take care at high end
  double dens_interp = (dens/0.016 - (densindex));
  if(medium==1) return dens_interp*(Q2_interp*QMCGE[Q2index][densindex]+(1.-Q2_interp)*QMCGE[Q2index+1][densindex])
    +(1.-dens_interp)*(Q2_interp*QMCGE[Q2index][densindex+1]+(1.-Q2_interp)*QMCGE[Q2index+1][densindex+1]);
  else return dens_interp*(Q2_interp*CQSMGE[Q2index][densindex]+(1.-Q2_interp)*CQSMGE[Q2index+1][densindex])
    +(1.-dens_interp)*(Q2_interp*CQSMGE[Q2index][densindex+1]+(1.-Q2_interp)*CQSMGE[Q2index+1][densindex+1]);
  
}
double NucleonEMOperator::getMmod(const double r, const int medium, const MeanFieldNucleusThick &nucleus) const{
  double dens=nucleus.getDensity(r,proton)/(r*r)*nucleus.getA();
  if(r==0.) dens=nucleus.getTotalDensity(1.E-02)/(1.E-04)*nucleus.getA();
  int densindex = int(floor(dens/0.016));
  if(densindex>=(medium==1?14:15)) densindex=(medium==1?14:15); //take care at high end
  double dens_interp = (dens/0.016 - (densindex));
  
  if(medium==1) return dens_interp*(Q2_interp*QMCGM[Q2index][densindex]+(1.-Q2_interp)*QMCGM[Q2index+1][densindex])
    +(1.-dens_interp)*(Q2_interp*QMCGM[Q2index][densindex+1]+(1.-Q2_interp)*QMCGM[Q2index+1][densindex+1]);
  else return dens_interp*(Q2_interp*CQSMGM[Q2index][densindex]+(1.-Q2_interp)*CQSMGM[Q2index+1][densindex])
    +(1.-dens_interp)*(Q2_interp*CQSMGM[Q2index][densindex+1]+(1.-Q2_interp)*CQSMGM[Q2index+1][densindex+1]);
  
}





}  // namespace nucqcd
