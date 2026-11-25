#include <Tools/Kinematics/LightConeKin2to3.hpp>

#include <cassert>
#include <Tools/Core/constants.hpp>

namespace nucqcd {


using namespace std;



LightConeKin2to3::LightConeKin2to3(double massA, double Q2, double massSp1, double massSp2, 
		   nucqcd::Vector3 &vecA, nucqcd::Vector3 &vecq, nucqcd::Vector3 &vecSp1, nucqcd::Vector3 &vecSp2, nucqcd::Vector3 &vecBeam, nucqcd::Vector3 &vecSpcm):
		   LightConeKin2to2(massA,Q2,
				   sqrt(pow(sqrt(massSp1*massSp1+vecSp1.Mag2())+sqrt(massSp2*massSp2+vecSp1.Mag2()),2.)
					-(vecSpcm).Mag2()),
					vecA,vecq,vecSpcm,vecBeam){

  Sp1_mu=FourVector<double>(sqrt(massSp1*massSp1+vecSp1.Mag2()),vecSp1.X(),vecSp1.Y(),vecSp1.Z());
  Sp2_mu=FourVector<double>(sqrt(massSp2*massSp2+vecSp2.Mag2()),vecSp2.X(),vecSp2.Y(),vecSp2.Z());
		     
}

LightConeKin2to3::LightConeKin2to3(double massA, double Q2, double massSp1, double massSp2, 
		   double pA, double vecq, nucqcd::Vector3 &vecSp1, nucqcd::Vector3 &vecSp2, nucqcd::Vector3 &vecBeam, nucqcd::Vector3 &vecSpcm):
		    LightConeKin2to2(massA,Q2, 
				 sqrt(pow(sqrt(massSp1*massSp1+vecSp1.Mag2())+sqrt(massSp2*massSp2+vecSp1.Mag2()),2.)
					-(vecSpcm).Mag2()),     
			pA, vecq, vecSpcm, vecBeam){
  Sp1_mu=FourVector<double>(sqrt(massSp1*massSp1+vecSp1.Mag2()),vecSp1.X(),vecSp1.Y(),vecSp1.Z());
  Sp2_mu=FourVector<double>(sqrt(massSp2*massSp2+vecSp2.Mag2()),vecSp2.X(),vecSp2.Y(),vecSp2.Z());		     
}
		      
LightConeKin2to3::LightConeKin2to3(LightConeKin2to3 &rhs): LightConeKin2to2(rhs){
  Sp1_mu=rhs.Sp1_mu;
  Sp2_mu=rhs.Sp2_mu;
}
  
  
LightConeKin2to3& LightConeKin2to3::operator=(LightConeKin2to3 &rhs){
 if(this!=&rhs) { // avoid self-assignment
    Sp1_mu=rhs.Sp1_mu;
    Sp2_mu=rhs.Sp2_mu;
  }
  return *this;
  
}

}  // namespace nucqcd

