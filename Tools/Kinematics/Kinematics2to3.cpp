/*
 * Kinematics2to3.cpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 * 
 */

///////////////////////////////////////////////////////////////////////////
/* Begin_Html
   <center><h2>Kinematics2to3</h2></center>
   Kinematics2to3 allows to calculate kinematics for g D -> K Y N.
   <p>
   All quantities are in MeV or rad.
   <p>


   <h3>Reference frame</h3>
   This class solves the kinematics in the LAB frame, where the deuteron
   is at rest. The <b>axis convention</b>, however, can be chosen in a number of
   different ways. The user can specify this with the <b>EAxisConvention</b> type.
   <p>
   The z axis is always along the photon momentum and we define x ~ y x z.
   This leaves us the choice of y axis:
   <ul>
   <li><b>kN</b>:&nbsp  y ~ p<sub>g</sub> x p<sub>N</sub>
   <li><b>kK</b>:&nbsp  y ~ p<sub>g</sub> x p<sub>K</sub>
   <li><b>kY</b>:&nbsp  y ~ p<sub>g</sub> x p<sub>Y</sub>
   </ul>
   Internally, Kinematics2to3 uses the kN convention. Getters that depend
   on the type of axis convention, take an EAxisConvention as argument.  
   The constant <b>kS</b> is defined as the axis convention chosen by the user  
   in the constructor. This allows the user to omit the EAxisConvention  
   argument in getters, which defaults to kS.
   <p>
   Some getters return quantities in the centre-of-mass (CM) of a 2-particle
   subsystem. These getters don't use the same axis as the LAB frame getters.
   The y axis is still fixed with the axis convention, but the x and z axis are
   rotated such that the z axis lies along the total momentum of the 2 particles.


   <h3>Isospin channels</h3>
   Kinematics2to3 is designed to work with a fixed number of final states, confusingly labeled <b>isospin channels</b>.
   <ul>
   <li> (1)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>+</sup> + L<sup>0</sup> + n
   <li> (2)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>+</sup> + S<sup>0</sup> + n
   <li> (3)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>0</sup> + S<sup>+</sup> + n
   <li> (4)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>0</sup> + L<sup>0</sup> + p
   <li> (5)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>0</sup> + S<sup>0</sup> + p
   <li> (6)&nbsp&nbsp   ===> &nbsp&nbsp  g + d &nbsp  
   -> &nbsp K<sup>+</sup> + S<sup>-</sup> + p
   </ul>


   <a name="indepvar"></a><h3>Independent variables</h3>
   The user needs to specify 6 variables to unambiguously fix the kinematics.
   Obviously, he/she needs to made sure the variables of choice are independent. 
   Kinematics2to3 is capable of dealing with 2 situations.
   <a name="indepvarCaseA"></a><h3><h4>Case A: Mother with children</h4>
   We consider the final state as the sum of one particle (mother) and the subsystem
   of the two other particles (children). The mother should be the same particle
   as the one that fixes the x and y axis (axis convention).
   <ul>
   <li>The user provides info about the energy and the theta angle (either cos or angle) of the mother.
   <li>A 3rd variable should either specify
   <ul><li>the photon's energy
   <li>the invariant mass of the total system
   <li>the invariant mass of the children subsystem</ul>
   <li>Finally, the user specifies the solid angle (theta+phi) of one of the children in the CM of the subsystem.
   <li>The virtualiy Q<sup>2</sup> of the photon is always one of the independent variables.
   </ul> 
   Note: the user can provide the solid angle of one of the children in the LAB frame 
   as well. Such a set of independent variables, however, can have two solutions. This
   class is not capable of dealing with this and will print a warning before arbitrarily
   choosing one of both solutions. See Kinematics2to3WithLabAngles for a kinematics
   class that is designed to deal with this situation.
   <h4>Case B: Papa Mandelstam</h4>
   All independent variables are Mandelstam variables, with the exception of Q<sup>2</sup> of course.
   Again, we define particle 1 as the one that fixes the x and y axis.
   <ul><li>s<sub>tot</sub> should always be given
   <li>For the other variables the user has 7 options:
   <ul><li>s<sub>12</sub>, s<sub>23</sub>, t<sub>g1</sub>, t<sub>d3</sub>
   <li>s<sub>12</sub>, s<sub>23</sub>, t<sub>g1</sub>, t<sub>d2</sub>
   <li>s<sub>12</sub>, s<sub>23</sub>, t<sub>d2</sub>, t<sub>d3</sub>
   <li>s<sub>12</sub>, t<sub>d2</sub>, t<sub>g1</sub>, t<sub>d3</sub>
   <li>s<sub>12</sub>, s<sub>23</sub>, t<sub>g1</sub>, t<sub>g2</sub>
   <li>s<sub>12</sub>, s<sub>23</sub>, t<sub>g2</sub>, t<sub>d3</sub>
   <li>t<sub>g2</sub>, s<sub>23</sub>, t<sub>g1</sub>, t<sub>d3</sub>
   </ul></ul>
   <p>
   End_Html */
///////////////////////////////////////////////////////////////////////////

#include <Tools/Kinematics/Kinematics2to3.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/LorentzQuaternion.hpp>
#include <Tools/Core/constants.hpp>
#include <Tools/Core/MathUtil.hpp>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <Tools/Core/StringUtil.hpp>
#include <algorithm>
#include <optional>
#include <vector>
#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/Vector3.hpp>

namespace nucqcd {

// #include <numtoa.h>
using std::cout; using std::cerr;
using std::vector;
using std::strcmp;


//_____________________________________________________________________
const double Kinematics2to3::kUflow = 1.e-7;

//_____________________________________________________________________
Kinematics2to3::Kinematics2to3(int isospin, EAxisConvention ac, 
				 const char *format,
				 double var1, double var2, double var3,
				 double var4, double var5, double var6,
				 double Md, double Mn, double Mk, double My)
  : fKYsystem(0), fKNsystem(0), fYNsystem(0),
    fKvector(0), fYvector(0), fNvector(0), fGvector(0), fDvector(0),
    fIsospin(0), fMd(Md), fMk(Mk), fMy(My), fMn(Mn), fIsPhysical(false),
    fQsquared(0), fWlab(0), fKlab(0), fEk(0), fEy(0), fEn(0),
    fPk(0), fPy(0), fPn(0), fThk(0), fThy(0), fThn(0), fCosthk(0),
    fCosthy(0), fCosthn(0), fPhik(0), fPhiy(0), fPhin(0), fWtot(0),
    fWky(0), fWkn(0), fWyn(0), fStot(0), fSky(0), fSkn(0), fSyn(0), fTgk(0),
    fTgy(0), fTgn(0), fKYcosthkcm(0), fKYcosthycm(0), fKNcosthkcm(0),
    fKNcosthncm(0), fYNcosthycm(0), fYNcosthncm(0),
    fVar(0), fE1Var(0), fTh1Var(0), fE2Var(0), fTh2Var(0), fPhi2Var(0),
    fTdk(0), fTdy(0), fTdn(0), fKYphikcm(0), fKYphiycm(0), fKNphikcm(0),
    fKNphincm(0), fYNphiycm(0), fYNphincm(0), fKYthkcm(0), fKYthycm(0),
    fKNthkcm(0), fKNthncm(0), fYNthycm(0), fYNthncm(0),
    fEVar(0), fIsVar(0), fLowerLimit(0), fUpperLimit(0),
    fStepSize(0), fNumberOfSteps(0), fStep(0), fFormat(0),
    fAC(ac), fKinSolution(0), fNrOfPhysicals(-1)
{
  // Constructor
  //
  // See SetFormat(const char*) for info about what constitutes a valid
  // format string. 
  //
  // var1-6 will set the initial value of the independent variables in the 
  // same order as listed in the format string.
  //
  // Make sure the independent variables you choose in the format string are
  // compatible with your axis convention.

  if( ac==kS ) {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in Kinematics2to3::Kinematics2to3(..): "
	 << "Kinematics2to3::kS is NOT a valid axis convention.");
  }

  ChangeIsospinChannel(isospin);

  // allocate memory for fFormat, fVar and f*vector
  fFormat = new char[100];
  fKinSolution = new char[4];
  fVar    = new double*[6];
  fKvector = new nucqcd::LorentzVector;
  fYvector = new nucqcd::LorentzVector;
  fNvector = new nucqcd::LorentzVector;
  fGvector = new nucqcd::LorentzVector;
  fDvector = new nucqcd::LorentzVector;

  ReadFormatString(format);
  *fVar[0] = var1;
  *fVar[1] = var2;
  *fVar[2] = var3;
  *fVar[3] = var4;
  *fVar[4] = var5;
  *fVar[5] = var6;

  InitializeSubsystems();

  UpdateKinematics();
  
  // Initialize kinematic variables
  fIsVar = new bool[6];
  fLowerLimit = new double[6];
  fUpperLimit = new double[6];
  fStepSize = new double[6];
  fNumberOfSteps = new int[6];
  
  // Assume the independent variables to be constant
  for(int i=0; i<6; ++i) {
    fIsVar[i] = false;
    fLowerLimit[i] = *fVar[i];
    fUpperLimit[i] = *fVar[i];
    fStepSize[i] = 0.0;
    fNumberOfSteps[i] = 1;
  }
}


//_____________________________________________________________________
Kinematics2to3::Kinematics2to3(const Kinematics2to3& toCopy)
  : fKYsystem(0), fKNsystem(0), fYNsystem(0),
    fKvector(0), fYvector(0), fNvector(0), fGvector(0), fDvector(0),
    fIsospin(toCopy.fIsospin), fMd(toCopy.fMd), fMk(toCopy.fMk), 
    fMy(toCopy.fMy), fMn(toCopy.fMn), fIsPhysical(toCopy.fIsPhysical),
    fQsquared(toCopy.fQsquared), fWlab(toCopy.fWlab), fKlab(toCopy.fKlab), 
    fEk(toCopy.fEk), fEy(toCopy.fEy), fEn(toCopy.fEn), fPk(toCopy.fPk), 
    fPy(toCopy.fPy), fPn(toCopy.fPn), fThk(toCopy.fThk), fThy(toCopy.fThy), 
    fThn(toCopy.fThn), fCosthk(toCopy.fCosthk), fCosthy(toCopy.fCosthy), 
    fCosthn(toCopy.fCosthn), fPhik(toCopy.fPhik), fPhiy(toCopy.fPhiy), 
    fPhin(toCopy.fPhin), fWtot(toCopy.fWtot), fWky(toCopy.fWky), fWkn(toCopy.fWkn), 
    fWyn(toCopy.fWyn), fStot(toCopy.fStot), fSky(toCopy.fSky), fSkn(toCopy.fSkn), 
    fSyn(toCopy.fSyn), fTgk(toCopy.fTgk), fTgy(toCopy.fTgy), fTgn(toCopy.fTgn), 
    fKYcosthkcm(toCopy.fKYcosthkcm), fKYcosthycm(toCopy.fKYcosthycm), 
    fKNcosthkcm(toCopy.fKNcosthkcm), fKNcosthncm(toCopy.fKNcosthncm), 
    fYNcosthycm(toCopy.fYNcosthycm), fYNcosthncm(toCopy.fYNcosthncm),
    fVar(0), fE1Var(0), fTh1Var(0), fE2Var(0), fTh2Var(0), fPhi2Var(0),
    fTdk(toCopy.fTdk), fTdy(toCopy.fTdy), fTdn(toCopy.fTdy), 
    fKYphikcm(toCopy.fKYphikcm), fKYphiycm(toCopy.fKYphiycm), 
    fKNphikcm(toCopy.fKNphikcm), fKNphincm(toCopy.fKNphincm), 
    fYNphiycm(toCopy.fYNphiycm), fYNphincm(toCopy.fYNphincm), 
    fKYthkcm(toCopy.fKYthkcm), fKYthycm(toCopy.fKYthycm), 
    fKNthkcm(toCopy.fKNthkcm), fKNthncm(toCopy.fKNthncm), 
    fYNthycm(toCopy.fYNthycm), fYNthncm(toCopy.fYNthncm),
    fEVar(0), fIsVar(0), fLowerLimit(0), fUpperLimit(0),
    fStepSize(0), fNumberOfSteps(0), fStep(toCopy.fStep), fFormat(0),
    fAC(toCopy.fAC), fKinSolution(0), fNrOfPhysicals(toCopy.fNrOfPhysicals)
{
  // Copy constructor

  if(toCopy.fVar) fVar =  new double*[6];
  if(toCopy.fFormat) {
    fFormat = new char[100];
    fKinSolution = new char[4];
    ReadFormatString(toCopy.fFormat);
  }

  fKvector = new nucqcd::LorentzVector(*toCopy.fKvector);
  fYvector = new nucqcd::LorentzVector(*toCopy.fYvector);
  fNvector = new nucqcd::LorentzVector(*toCopy.fNvector);
  fGvector = new nucqcd::LorentzVector(*toCopy.fGvector);
  fDvector = new nucqcd::LorentzVector(*toCopy.fDvector);

  if(toCopy.fKYsystem) fKYsystem = new Kinematics2to2(*toCopy.fKYsystem);
  if(toCopy.fKNsystem) fKNsystem = new Kinematics2to2(*toCopy.fKNsystem);
  if(toCopy.fYNsystem) fYNsystem = new Kinematics2to2(*toCopy.fYNsystem);

  if(toCopy.fIsVar) {  
    fIsVar = new bool[6];
    fLowerLimit = new double[6];
    fUpperLimit = new double[6];
    fStepSize = new double[6];
    fNumberOfSteps = new int[6];
  
    for(int i=0; i<6; ++i) {
      fIsVar[i] = toCopy.fIsVar[i];
      fLowerLimit[i] = toCopy.fLowerLimit[i];
      fUpperLimit[i] = toCopy.fUpperLimit[i];
      fStepSize[i] = toCopy.fStepSize[i];
      fNumberOfSteps[i] = toCopy.fNumberOfSteps[i];
    }
  }
}

//_____________________________________________________________________

//_____________________________________________________________________
Kinematics2to3::~Kinematics2to3()
{
  // Destructor

  delete[] fFormat;
  delete[] fKinSolution;
  delete[] fVar;
  delete[] fIsVar;
  delete[] fLowerLimit;
  delete[] fUpperLimit;
  delete[] fStepSize;
  delete[] fNumberOfSteps;

  delete fKYsystem;
  delete fKNsystem;
  delete fYNsystem;
  delete fKvector;
  delete fYvector;
  delete fNvector;
  delete fGvector;
  delete fDvector;
}

//_____________________________________________________________________
Kinematics2to3& Kinematics2to3::operator=(const Kinematics2to3& toCopy)
{
  // Assignment

  if( this!=&toCopy ) { // avoid self-assignment
    
    ChangeIsospinChannel(toCopy.fIsospin);
    fAC = toCopy.fAC;
    ReadFormatString(toCopy.fFormat);

    fIsPhysical = toCopy.fIsPhysical;
    fQsquared = toCopy.fQsquared;
    fWlab = toCopy.fWlab;
    fKlab = toCopy.fKlab;
    fEk = toCopy.fEk;
    fEy = toCopy.fEy;
    fEn = toCopy.fEn;
    fPk = toCopy.fPk;
    fPy = toCopy.fPy;
    fPn = toCopy.fPn;
    fThk = toCopy.fThk;
    fThy = toCopy.fThy;
    fThn = toCopy.fThn;
    fCosthk = toCopy.fCosthk;
    fCosthy = toCopy.fCosthy;
    fCosthn = toCopy.fCosthn;
    fPhik = toCopy.fPhik;
    fPhiy = toCopy.fPhiy;
    fPhin = toCopy.fPhin;
    fWtot = toCopy.fWtot;
    fWky = toCopy.fWky;
    fWkn = toCopy.fWkn;
    fWyn = toCopy.fWyn;
    fStot = toCopy.fStot;
    fSky = toCopy.fSky;
    fSkn = toCopy.fSkn;
    fSyn = toCopy.fSyn;
    fTgk = toCopy.fTgk;
    fTgy = toCopy.fTgy;
    fTgn = toCopy.fTgn;
    fTdk = toCopy.fTdk;
    fTdy = toCopy.fTdy;
    fTdn = toCopy.fTdn;
    fKYthkcm = toCopy.fKYthkcm;
    fKYthycm = toCopy.fKYthycm;
    fKNthkcm = toCopy.fKNthkcm;
    fKNthncm = toCopy.fKNthncm;
    fYNthycm = toCopy.fYNthycm;
    fYNthncm = toCopy.fYNthncm;
    fKYcosthkcm = toCopy.fKYcosthkcm;
    fKYcosthycm = toCopy.fKYcosthycm;
    fKNcosthkcm = toCopy.fKNcosthkcm;
    fKNcosthncm = toCopy.fKNcosthncm;
    fYNcosthycm = toCopy.fYNcosthycm;
    fYNcosthncm = toCopy.fYNcosthncm;
    fKYphikcm = toCopy.fKYphikcm;
    fKYphiycm = toCopy.fKYphiycm;
    fKNphikcm = toCopy.fKNphikcm;
    fKNphincm = toCopy.fKNphincm;
    fYNphiycm = toCopy.fYNphiycm;
    fYNphincm = toCopy.fYNphincm;

    *fKYsystem = *toCopy.fKYsystem;
    *fKNsystem = *toCopy.fKNsystem;
    *fYNsystem = *toCopy.fYNsystem;

    *fKvector = *toCopy.fKvector;
    *fYvector = *toCopy.fYvector;
    *fNvector = *toCopy.fNvector;
    *fGvector = *toCopy.fGvector;
    *fDvector = *toCopy.fDvector;

    fStep = toCopy.fStep;
    fNrOfPhysicals = toCopy.fNrOfPhysicals;

    if(toCopy.fIsVar) {
      if( !fIsVar) fIsVar = new bool[6];
      if( !fLowerLimit) fLowerLimit = new double[6];
      if( !fUpperLimit) fUpperLimit = new double[6];
      if( !fStepSize) fStepSize = new double[6];
      if( !fNumberOfSteps) fNumberOfSteps = new int[6];
      
      for(int i=0; i<6; ++i) {
	fIsVar[i] = toCopy.fIsVar[i];
	fLowerLimit[i] = toCopy.fLowerLimit[i];
	fUpperLimit[i] = toCopy.fUpperLimit[i];
	fStepSize[i] = toCopy.fStepSize[i];
	fNumberOfSteps[i] = toCopy.fNumberOfSteps[i];
      }
    } else {
      if( fIsVar) { delete[] fIsVar; fIsVar = 0; }
      if( fLowerLimit) { delete[] fLowerLimit; fLowerLimit = 0; }
      if( fUpperLimit) { delete[] fUpperLimit; fUpperLimit = 0; }
      if( fStepSize) { delete[] fStepSize; fStepSize = 0; }
      if( fNumberOfSteps) { delete[] fNumberOfSteps; fNumberOfSteps = 0; }
    }
  }

  return *this;
}

//_____________________________________________________________________
void Kinematics2to3::SetIsospin(int isospin)
{
  // Change isospin channel. Kinematics will immediately be updated.

  ChangeIsospinChannel(isospin);
  UpdateKinematics();
}

//_____________________________________________________________________
void Kinematics2to3::ChangeIsospinChannel(int isospin)
{
  // Set the particle masses according to the isospin channel

  /*   (1)   ===>   g + d   ->  K+ + L0 + n
   *   (2)   ===>   g + d   ->  K+ + S0 + n
   *   (3)   ===>   g + d   ->  K0 + S+ + n
   *   (4)   ===>   g + d   ->  K0 + L0 + p
   *   (5)   ===>   g + d   ->  K0 + S0 + p
   *   (6)   ===>   g + d   ->  K+ + S- + p
   */

  fIsospin = isospin;

  //fMd = M_D;  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  switch( isospin ) 
    {
    case 1:
      fMd = nucqcd::mass_d;
      fMn = nucqcd::mass_n;
      fMk = nucqcd::mass_kp;
      fMy = nucqcd::mass_lambda;
      break;
    case 2:
      fMd = nucqcd::mass_d;
      fMn = nucqcd::mass_n;
      fMk = nucqcd::mass_kp;
      fMy = nucqcd::mass_sigma0;
      break;
    case 3:
      fMd = nucqcd::mass_d;
      fMn = nucqcd::mass_n;
      fMk = nucqcd::mass_k0;
      fMy = nucqcd::mass_sigmap;
      break;
    case 4:
      fMd = nucqcd::mass_d;
      fMn = nucqcd::mass_p;
      fMk = nucqcd::mass_k0;
      fMy = nucqcd::mass_lambda;
      break;
    case 5:
      fMd = nucqcd::mass_d;
      fMn = nucqcd::mass_p;
      fMk = nucqcd::mass_k0;
      fMy = nucqcd::mass_sigma0;
      break;
    case 6:
      fMd = nucqcd::mass_d;
      fMn = nucqcd::mass_p;
      fMk = nucqcd::mass_kp;
      fMy = nucqcd::mass_sigmam;
      break;
    case 7:
      break;
    default:
      std::cerr << "ERROR in Kinematics2to3::ChangeIsospinChannel(int):\n"
		<< "Isospin channel " << isospin << " is not implemented!\n";
      break;
    }  

  // if the subsystems have been created, update their masses
  if(fKYsystem) {   
    fKYsystem->SetMasses(fMd,fMk,fMy);
    fKNsystem->SetMasses(fMd,fMk,fMn);
    fYNsystem->SetMasses(fMd,fMn,fMy);
  }
  
  // Reset the number of physical steps
  fNrOfPhysicals = -1;
}

//_____________________________________________________________________
void Kinematics2to3::SetFormat(const char* format)
{
  /* Begin_Html
     Set the 6 independent variables by providing a ':' separated list of variables.
  
     All variables are in the LAB frame, except when their name explicitly
     suggests differently.
  
     Look <a href="#indepvar">here</a> to find out what your format string should
     look like.
     <p>
     Below is a list of the 'code names' of variables:
     <table>
     <tr><td>Q<sup>2</sup></td><td><tt>qsquared</tt></td></tr>
     <tr><td>photon energy</sub></td><td><tt>wlab</tt></td></tr>
     <tr><td>photon momentum</td><td><tt>klab</tt></td></tr>
     <tr><td>kaon energy</td><td><tt>ek</tt></td></tr>
     <tr><td>hyperon energy</td><td><tt>ey</tt></td></tr>
     <tr><td>nucleon energy</td><td><tt>en</tt></td></tr>
     <tr><td>kaon momentum</td><td><tt>pk</tt></td></tr>
     <tr><td>hyperon momentum</td><td><tt>py</tt></td></tr>
     <tr><td>nucleon momentum</td><td><tt>pn</tt></td></tr>
     <tr><td>kaon theta</td><td><tt>thk</tt></td></tr>
     <tr><td>hyperon theta</td><td><tt>tky</tt></td></tr>
     <tr><td>nucleon theta</td><td><tt>thn</tt></td></tr>
     <tr><td>cos(theta<sub>K</sub>)</td><td><tt>costhk</tt></td></tr>
     <tr><td>cos(theta<sub>Y</sub>)</td><td><tt>costhy</tt></td></tr>
     <tr><td>cos(theta<sub>N</sub>)</td><td><tt>costhn</tt></td></tr>
     <tr><td>kaon phi</td><td><tt>phik</tt></td></tr>
     <tr><td>hyperon phi</td><td><tt>phiy</tt></td></tr>
     <tr><td>nucleon phi</td><td><tt>phin</tt></td></tr>
     <tr><td>total invariant mass</td><td><tt>wtot</tt></td></tr>
     <tr><td>invariant mass KY-system</td><td><tt>wky</tt></td></tr>
     <tr><td>invariant mass KN-system</td><td><tt>wkn</tt></td></tr>
     <tr><td>invariant mass YN-system</td><td><tt>wyn</tt></td></tr>
     <tr><td>total s</td><td><tt>stot</tt></td></tr>
     <tr><td>s of KY-system</td><td><tt>sky</tt></td></tr>
     <tr><td>s of KN-system</td><td><tt>skn</tt></td></tr>
     <tr><td>s of YN-system</td><td><tt>syn</tt></td></tr>
     <tr><td>Mom.transf. photon-kaon</td><td><tt>tgk</tt></td></tr>
     <tr><td>Mom.transf. photon-hyperon</td><td><tt>tgy</tt></td></tr>
     <tr><td>Mom.transf. photon-nucleon</td><td><tt>tgn</tt></td></tr>
     <tr><td>Mom.transf. deut.-kaon</td><td><tt>tdk</tt></td></tr>
     <tr><td>Mom.transf. deut.-hyperon</td><td><tt>tdy</tt></td></tr>
     <tr><td>Mom.transf. deut.-nucleon</td><td><tt>tdn</tt></td></tr>
     <tr><td>kaon theta in KY-system</td><td><tt>kythkcm</tt></td></tr>
     <tr><td>hyperon theta in KY-system</td><td><tt>kythycm</tt></td></tr>
     <tr><td>kaon theta in KN-system</td><td><tt>knthkcm</tt></td></tr>
     <tr><td>nucleon theta in KN-system</td><td><tt>knthncm</tt></td></tr>
     <tr><td>hyperon theta in YN-system</td><td><tt>ynthycm</tt></td></tr>
     <tr><td>nucleon theta in YN-system</td><td><tt>ynthncm</tt></td></tr>
     <tr><td>cos(<tt>kythkcm</tt>)</td><td><tt>kycosthkcm</tt></td></tr>
     <tr><td>cos(<tt>kythycm</tt>)</td><td><tt>kycosthycm</tt></td></tr>
     <tr><td>cos(<tt>knthkcm</tt>)</td><td><tt>kncosthkcm</tt></td></tr>
     <tr><td>cos(<tt>knthncm</tt>)</td><td><tt>kncosthncm</tt></td></tr>
     <tr><td>cos(<tt>ynthycm</tt>)</td><td><tt>yncosthycm</tt></td></tr>
     <tr><td>cos(<tt>ynthncm</tt>)</td><td><tt>yncosthncm</tt></td></tr>
     <tr><td>kaon phi in KY-system</td><td><tt>kyphikcm</tt></td></tr>
     <tr><td>hyperon phi in KY-system</td><td><tt>kyphiycm</tt></td></tr>
     <tr><td>kaon phi in KN-system</td><td><tt>knphikcm</tt></td></tr>
     <tr><td>nucleon phi in KN-system</td><td><tt>knphincm</tt></td></tr>
     <tr><td>hyperon phi in YN-system</td><td><tt>ynphiycm</tt></td></tr>
     <tr><td>nucleon phi in YN-system</td><td><tt>ynphincm</tt></td></tr>
     </table>
     End_Html */ 

  if(!IsFixed()) {
    std::cout << "INFO in Kinematics2to3::SetFormat(const char*): "
		 << "Fixing all variables.\n";
    for(int i=0; i<6; ++i)
      FixVariable(i+1);
  }

  UpdateVariables();

  ReadFormatString(format);
  
  // Finally, set the upper and lower limit
  for(int i=0; i<6; ++i) {
    fLowerLimit[i] = *fVar[i];
    fUpperLimit[i] = *fVar[i];
  }
}

//_____________________________________________________________________
void Kinematics2to3::SetAxisConvention(EAxisConvention ac)
{
  // Change the axis convention. Be sure you know what you're doing before
  // you call this function. This function doesn't check whether what 
  // you are doing makes sense.
  fAC = ac;
}

}  // namespace nucqcd
