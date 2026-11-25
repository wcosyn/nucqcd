/*
 * TKinematics2to3_Variables.cpp
 * TKinematics2to3: setting, stepping and getting the kinematic variables by number or name
 * Part of the implementation of TKinematics2to3 (split from TKinematics2to3.cpp, which holds the class
 * description, constructors and assignment); author of the original: Pieter Vancraeyveld.
 */

#include <Tools/Kinematics/TKinematics2to3.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/TLorentzQuaternion.hpp>
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
void TKinematics2to3::SetVarRange(int varNumber, double lowerLimit,
				  double upperLimit, int numberOfSteps)
{
  // Give independent variable number 'varNumber' a range in the interval [lowerLimit,upperLimit].

  // First we do some checks
  if( varNumber<1 || varNumber>6 )
    std::cerr << "ERROR in TKinematics2to3::SetVarRange(int,double,double,int):\n"
	      << "Index out-of-bounds!\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics2to3::SetVarRange(int,double,double,int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else if( numberOfSteps<2 )
    std::cerr << "ERROR in TKinematics2to3::SetVarRange(int,double,double,int)\n"
	      << "Number of steps need to be >=2.\n";

  else {
    fIsVar[varNumber-1] = true;

    // Set limits
    fLowerLimit[varNumber-1] = lowerLimit;
    fUpperLimit[varNumber-1] = upperLimit;

    // Set number of steps
    fNumberOfSteps[varNumber-1] = numberOfSteps;

    // Calculate stepsize
    fStepSize[varNumber-1] = (upperLimit-lowerLimit) / ((numberOfSteps-1)*1.0);

    // Go to first step
    GoTo(0);

    // Reset the number of physical steps
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics2to3::FixVariable(int varNumber)
{
  // Unset the range of variable number 'varNumber'. Set it to its current value.

  // First we do some checks
  if( varNumber<1 || varNumber>6 )
    std::cerr << "ERROR in TKinematics2to3::FixVariable(int):\n"
	      << "Index out-of-bounds!\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics2to3::FixVariable(int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else if( fIsVar[varNumber-1] ){
    // Determine the current step
    int step1 = GetStep(1);
    int step2 = GetStep(2);
    int step3 = GetStep(3);
    int step4 = GetStep(4);
    int step5 = GetStep(5);
    int step6 = GetStep(6);
    fIsVar[varNumber-1] = false;
    
    // Reset all variables
    fLowerLimit[varNumber-1] = *fVar[varNumber-1];
    fUpperLimit[varNumber-1] = *fVar[varNumber-1];
    fStepSize[varNumber-1] = 0.0;
    fNumberOfSteps[varNumber-1] = 1;

    // Go to first step
    switch(varNumber) {
    case 1: GoTo(0,step2,step3,step4,step5,step6); break;
    case 2: GoTo(step1,0,step3,step4,step5,step6); break;
    case 3: GoTo(step1,step2,0,step4,step5,step6); break;
    case 4: GoTo(step1,step2,step3,0,step5,step6); break;
    case 5: GoTo(step1,step2,step3,step4,0,step6); break;
    case 6: GoTo(step1,step2,step3,step4,step5,0); break;
    }

    // Reset the number of physical steps
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics2to3::FixVariables()
{
  // Calls FixVariable(int) on all variables.
  for(int i=1; i<=6; ++i)
    FixVariable(i);
}

//_____________________________________________________________________
int TKinematics2to3::Next()
{
  // Proceed to the next step in the grid of independent variables.
  // Return value is zero when the end of the grid is reached.

  // Total number of steps
  int totNrOfSteps = 1;
  for(int i=0; i<6; ++i) totNrOfSteps *= fNumberOfSteps[i];

  // Check whether we are at the end of the grid.
  if( fStep == (totNrOfSteps-1) ) return 0;

  GoTo(++fStep);
  return 1;
}

//_____________________________________________________________________
void TKinematics2to3::GoTo(int step)
{
  // Proceed to the 'step'th point in the grid of independent variables.

  // Total number of steps
  int totNrOfSteps = 1;
  for(int i=0; i<6; ++i)
    totNrOfSteps *= fNumberOfSteps[i];

  // First we do some checks
  if( step<0 || step>=totNrOfSteps )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int):\n"
	      << "Invalid step. Should be in range [0,"
	      << totNrOfSteps-1 << "].\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics2to3::GoTo(int)\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else {
    fStep = step;
    
    // Determine step per variable
    int stepVar[6];
    int remainder=step;
    int denum=1;
  
    for(int i=6; i>0; --i)
      {
	denum = 1;
	for(int j=0; j<i-1; ++j)
	  denum *= fNumberOfSteps[j];
	
	stepVar[i-1] = remainder / denum;
	remainder = remainder % denum;
      }
    
    // Set variables
    for(int i=0; i<6; ++i)
      *fVar[i] = fLowerLimit[i] + stepVar[i] * fStepSize[i];

    UpdateKinematics();
  }
}

//_____________________________________________________________________
void TKinematics2to3::GoTo(int stepVar1, int stepVar2, int stepVar3,
			   int stepVar4, int stepVar5, int stepVar6)
{
  // Proceed to point in the grid of independent variables specified by the
  // arguments. 
  // Each stepVari fixes the position in the grid of the 'i'th variable.

  // First we do some checks
  if( stepVar1<0 || stepVar1>=fNumberOfSteps[0] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(1) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[0]-1 << "].\n";

  else if( stepVar2<0 || stepVar2>=fNumberOfSteps[1] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(2) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[1]-1 << "].\n";

  else if( stepVar3<0 || stepVar3>=fNumberOfSteps[2] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(3) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[2]-1 << "].\n";

  else if( stepVar4<0 || stepVar4>=fNumberOfSteps[3] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(4) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[3]-1 << "].\n";

  else if( stepVar5<0 || stepVar5>=fNumberOfSteps[4] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(5) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[4]-1 << "].\n";

  else if( stepVar6<0 || stepVar6>=fNumberOfSteps[5] )
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6)):\n"
	      << "Invalid \"" << GetVarName(6) << "\" step. "
	      << "Should be in range [0," << fNumberOfSteps[5]-1 << "].\n";

  else if( !fVar ) 
    std::cerr << "ERROR in TKinematics2to3::GoTo(int(x6))\n"
	      << "Independent kinematic variables have not been initialized!\n";

  else {
    fStep = stepVar1 
      + stepVar2 * fNumberOfSteps[0]
      + stepVar3 * fNumberOfSteps[0] * fNumberOfSteps[1]
      + stepVar4 * fNumberOfSteps[0] * fNumberOfSteps[1] * fNumberOfSteps[2]
      + stepVar5 * (fNumberOfSteps[0] * fNumberOfSteps[1] 
		    * fNumberOfSteps[2]*fNumberOfSteps[3])
      + stepVar5 * (fNumberOfSteps[0] * fNumberOfSteps[1] * fNumberOfSteps[2] 
		    * fNumberOfSteps[3]*fNumberOfSteps[4]);
    
    // Set variables
    *fVar[0] = fLowerLimit[0] + stepVar1 * fStepSize[0];
    *fVar[1] = fLowerLimit[1] + stepVar2 * fStepSize[1];
    *fVar[2] = fLowerLimit[2] + stepVar3 * fStepSize[2];
    *fVar[3] = fLowerLimit[3] + stepVar4 * fStepSize[3];
    *fVar[4] = fLowerLimit[4] + stepVar5 * fStepSize[4];
    *fVar[5] = fLowerLimit[5] + stepVar6 * fStepSize[5];

    UpdateKinematics();
  }
}

//_____________________________________________________________________
int TKinematics2to3::VariableInfo() const
{
  // Print information on the type of independent variables and their status.

  std::cout << GetName() << ": " << GetTitle()
	       << "\n****************************************************\n";

  for(int var=0; var<6; ++var) {
    std::cout << var+1 << ": " << GetVarName(var+1) << " is ";
    
    if( fIsVar[var] )
      std::cout << "variable.\n"
		   << "   Range:\t\t" << fLowerLimit[var] << " <-> "
		   << fUpperLimit[var]
		   << "\n   Number of steps:\t" << fNumberOfSteps[var]
		   << "\n   Current step:\t" << GetStep(var+1) << "\n";
    
    else
      std::cout << "fixed.\n";
  }
  
  std::cout << "****************************************************\n";

  return GetNumberOfVariables();
}

//_____________________________________________________________________
void TKinematics2to3::SetVar(int varNumber, double value)
{
  if( varNumber<1 || varNumber>6 )
    std::cerr << "ERROR in TKinematics2to3::SetVar(int,double): "
	      << "Index out-of-bounds!\n";

  else if( !fVar )
    std::cerr << "ERROR in TKinematics2to3::SetVar(int,double): "
	      << "Independent kinematic variables have not been initialized!\n";
    
  else if( fIsVar[varNumber-1] )
    std::cerr << "ERROR in TKinematics2to3::SetVar(int,double): "
	      << "This method can not be used because the variable is not fixed.\n";

  else {
    fLowerLimit[varNumber-1] = value;
    fUpperLimit[varNumber-1] = value;
    GoTo(0);
    
    // Reset the number of physical steps
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
void TKinematics2to3::SetVar(double var1, double var2, double var3,
			     double var4, double var5, double var6)
{
  if( !fVar )
    std::cerr << "ERROR in TKinematics2to3::SetVar(double x6): "
	      << "Independent kinematic variables have not been initialized!\n";
    
  else if( !IsFixed() )
    std::cerr << "ERROR in TKinematics2to3::SetVar(double x6): "
	      << "This method can not be used because some variables are not fixed.\n";

  else {
    fLowerLimit[0] = var1;
    fUpperLimit[0] = var1;
    fLowerLimit[1] = var2;
    fUpperLimit[1] = var2;
    fLowerLimit[2] = var3;
    fUpperLimit[2] = var3;
    fLowerLimit[3] = var4;
    fUpperLimit[3] = var4;
    fLowerLimit[4] = var5;
    fUpperLimit[4] = var5;
    fLowerLimit[5] = var6;
    fUpperLimit[5] = var6;
    GoTo(0);

    // Reset the number of physical steps
    fNrOfPhysicals = -1;
  }
}

//_____________________________________________________________________
double TKinematics2to3::GetVar(int varNumber) const
{
  // Returns the value of the independent variable number 'varNumber'
  // as specified in the format string.

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetVar(int):\n"
	      << "Index out-of-bounds!\n";
    return -1e16;
  }

  if( !fVar ) {
    std::cerr << "ERROR in TKinematics2to3::GetVar(int):\n"
	      << "Independent kinematic variables have not been initialized!\n";
    return -1e16;
  }

  return *fVar[varNumber-1];
  
}

//_____________________________________________________________________
double TKinematics2to3::GetVar(const std::string& variable) const
{
  // Returns the value of the variable given by name. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.

  VarGetter ptrToVarGetter; // pointer to a TKinematics2to3 getter
  
  // Look up the getter (table in TKinematics2to3_VarNames.cpp)
  ptrToVarGetter = FindVarGetter(variable, "TKinematics2to3::GetVar(const std::string&)");

  return (this->*ptrToVarGetter)();
}

//_____________________________________________________________________
double* TKinematics2to3::GetVarArray(const std::string& variable) const
{
  // Return an array of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.
  //
  // The returned pointer has a size GetNumberOfSteps() and is owned by the user.

  TKinematics2to3 tk(*this); // temporary copy
  VarGetter ptrToVarGetter; // pointer to a TKinematics2to3 getter
  
  // Look up the getter (table in TKinematics2to3_VarNames.cpp)
  ptrToVarGetter = FindVarGetter(variable, "TKinematics2to3::GetVarArray(std::string)");

  double *array = new double[tk.GetNumberOfSteps()];
  
  // Fill the array
  for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
    tk.GoTo(step);
    array[step] = (tk.*ptrToVarGetter)();
  }    

  return array;
}

//_____________________________________________________________________
double* TKinematics2to3::GetPhysicalVarArray(const std::string& variable) const
{
  // Return an array of the variable of type 'variable'. The variable can effectively
  // by any one for which their exists a getter. It need not be listed in the info
  // for SetFormat(const char*), eg. kyphikcm_deg.
  //
  // The returned pointer has a size GetNumberOfPhysicalSteps() and is owned by the user.

  TKinematics2to3 tk(*this); // temporary copy
  VarGetter ptrToVarGetter; // pointer to a TKinematics2to3 getter

  // Look up the getter (table in TKinematics2to3_VarNames.cpp)
  ptrToVarGetter = FindVarGetter(variable, "TKinematics2to3::GetPhysicalVarArray(std::string)");

  double *array = new double[tk.GetNumberOfPhysicalSteps()];
  
  // Fill the array
  int physicalStep = 0;
  for(int step=0; step<tk.GetNumberOfSteps(); ++step) {
    tk.GoTo(step);
    if(tk.IsPhysical()) array[physicalStep++] = (tk.*ptrToVarGetter)();
  }    

  return array;
}

//_____________________________________________________________________
std::string TKinematics2to3::GetVarName(int varNumber) const
{

  // First we do some checks
  if( varNumber<1 || varNumber>6 ) {
    std::cerr << "ERROR in TKinematics2to3::GetVarName(int):\n"
	      << "Index out-of-bounds!\n";
    return std::string();
  }

  if( !fVar ) {
    std::cerr << "ERROR in TKinematics2to3::GetVarName(int):\n"
	      << "Independent kinematic variables have not been initialized!\n";
    return std::string();
  }

  std::string varName;
  
  // Lookup name
  if( fVar[varNumber-1] == &fQsquared )
    varName = "qsquared";
  else if( fVar[varNumber-1] == &fWlab )
    varName = "wlab";
  else if( fVar[varNumber-1] == &fKlab )
    varName = "klab";
  else if( fVar[varNumber-1] == &fEk )
    varName = "ek";
  else if( fVar[varNumber-1] == &fEy )
    varName = "ey";
  else if( fVar[varNumber-1] == &fEn )
    varName = "en";
  else if( fVar[varNumber-1] == &fPk )
    varName = "pk"; 
  else if( fVar[varNumber-1] == &fPy )
    varName = "py"; 
  else if( fVar[varNumber-1] == &fPn )
    varName = "pn"; 
  else if( fVar[varNumber-1] == &fThk )
    varName = "thk";   
  else if( fVar[varNumber-1] == &fThy )
    varName = "thy";   
  else if( fVar[varNumber-1] == &fThn )
    varName = "thn";   
  else if( fVar[varNumber-1] == &fCosthk )
    varName = "costhk";
  else if( fVar[varNumber-1] == &fCosthy )
    varName = "costhy";
  else if( fVar[varNumber-1] == &fCosthn )
    varName = "costhn";
  else if( fVar[varNumber-1] == &fPhik )
    varName = "phik";
  else if( fVar[varNumber-1] == &fPhiy )
    varName = "phiy";  
  else if( fVar[varNumber-1] == &fPhin )
    varName = "phin";  
  else if( fVar[varNumber-1] == &fWtot )
    varName = "wtot";  
  else if( fVar[varNumber-1] == &fWky )
    varName = "wky";   
  else if( fVar[varNumber-1] == &fWkn )
    varName = "wkn";   
  else if( fVar[varNumber-1] == &fWyn )
    varName = "wyn";   
  else if( fVar[varNumber-1] == &fStot )
    varName = "stot";  
  else if( fVar[varNumber-1] == &fSky )
    varName = "sky";   
  else if( fVar[varNumber-1] == &fSkn )
    varName = "skn";   
  else if( fVar[varNumber-1] == &fSyn )
    varName = "syn";   
  else if( fVar[varNumber-1] == &fTgk )
    varName = "tgk";   
  else if( fVar[varNumber-1] == &fTgy )
    varName = "tgy";   
  else if( fVar[varNumber-1] == &fTgn )
    varName = "tgn";   
  else if( fVar[varNumber-1] == &fTdk )
    varName = "tdk";   
  else if( fVar[varNumber-1] == &fTdy )
    varName = "tdy";   
  else if( fVar[varNumber-1] == &fTdn )
    varName = "tdn";   
  else if( fVar[varNumber-1] == &fKYthkcm )
    varName = "kythkcm";
  else if( fVar[varNumber-1] == &fKYthycm )
    varName = "kythycm";
  else if( fVar[varNumber-1] == &fKYcosthkcm )
    varName = "kycosthkcm";
  else if( fVar[varNumber-1] == &fKYcosthycm )
    varName = "kycosthycm";
  else if( fVar[varNumber-1] == &fKYphikcm )
    varName = "kyphikcm";
  else if( fVar[varNumber-1] == &fKYphiycm )
    varName = "kyphiycm";  
  else if( fVar[varNumber-1] == &fKNthkcm )
    varName = "knthkcm";   
  else if( fVar[varNumber-1] == &fKNthncm )
    varName = "knthncm";   
  else if( fVar[varNumber-1] == &fKNcosthkcm )
    varName = "kncosthkcm";
  else if( fVar[varNumber-1] == &fKNcosthncm )
    varName = "kncosthncm";
  else if( fVar[varNumber-1] == &fKNphikcm )
    varName = "knphikcm";  
  else if( fVar[varNumber-1] == &fKNphincm )
    varName = "knphincm";  
  else if( fVar[varNumber-1] == &fYNthycm )
    varName = "ynthycm";   
  else if( fVar[varNumber-1] == &fYNthncm )
    varName = "ynthncm";   
  else if( fVar[varNumber-1] == &fYNcosthycm )
    varName = "yncosthycm";
  else if( fVar[varNumber-1] == &fYNcosthncm )
    varName = "yncosthncm";
  else if( fVar[varNumber-1] == &fYNphiycm )
    varName = "ynphiycm";  
  else if( fVar[varNumber-1] == &fYNphincm )
    varName = "ynphincm";  
  else {
    std::cerr << "Error in TKinematics2to3::GetVarName(int):\n"
	      << "TKinematics2to3::fVar[" << varNumber-1 << "] "
	      << "points to unknown variable.\n";
  }

  return varName;
}

}  // namespace nucqcd
