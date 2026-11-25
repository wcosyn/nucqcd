/*
 * TKinematics2to3_ReadFormatString.cpp
 * TKinematics2to3::ReadFormatString: parsing of the format string (the choice of independent variables)
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
void TKinematics2to3::ReadFormatString(const char *format)
{
  // Process the format string.

  // Store format in TKinematics2to3::fFormat
  // Read TKinematics2to3::fFormat and let fVar[] point
  // to the 6 independent kinematical points as specified
  // by the user.
  
  // Check that the chosen variables make sense.
  // The user should always provide Q^2 of the photon

  // We identify 3 distinct scenarios:
  // 1) 1+(23)
  // 2) 1+2
  // 3) mandelstam
  // Given the format string this function will determine what
  // scenario is adequate to solve the kinematics. Some specifics
  // will be stored in fKinSolution to assist UpdateKinematics().
  std::strcpy(fKinSolution,"---");

  // Some more info on the different scenarios:
  // 1) The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    One variable should either specify the photon energy
  //    or the invariant mass of the total system or the
  //    (2+3) subsystem (fEVar).
  //    In addition angular info needs to be given about a
  //    second final particle. Both theta (fTh2Var) and phi (fPhi2Var)
  //    should be either in the cm of the (2+3) subsystem or in the
  //    lab frame. This info will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "1YL": scenario 1, 2nd particle = hyperon with angles in LAB
  //    "1K*": scenario 1, 2nd particle = kaon with angles in CM
  //
  // 2) The user provides info about the energy (fE1Var)
  //    and the theta-angle (fTh1Var) of the particle that
  //    fixes the x-axis (fAC).
  //    In addition a second final particle is fully specified,
  //    energy, theta and phi (fE2Var,fTh2Var,fPhi2Var) in the lab frame.
  //    The name of the 2nd particle will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "2N-": scenario 2, 2nd particle = nucleon
  //
  // 3) All variables are mandelstam variables, except for Q^2
  //    obviously. 
  //    S_tot = (P_gamma + P_deut)^2 should always be given.
  //    fAC fixes the 1st particle. This leaves 7 combinations:
  //    (3a) s_12 + s_23 + t_g1 + t_d3
  //    (3b) s_12 + s_23 + t_g1 + t_d2
  //    (3c) s_12 + s_23 + t_d2 + t_d3
  //    (3d) s_12 + t_d2 + t_g1 + t_d3
  //    (3e) s_12 + s_23 + t_g1 + t_g2
  //    (3f) s_12 + s_23 + t_g2 + t_d3
  //    (3g) t_g2 + s_23 + t_g1 + t_d3
  //    The name of the 2nd particle will be stored in fKinSolution.
  //    fKinSolution can look like this:
  //    "3Kb": scenario 3b with 2nd particle = kaon

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!                                                      !!
  // !! ATTENTION                                            !!
  // !! some cases are not valid anymore, until I find a way !!
  // !! to solve them without ambiguities.                   !!
  // !!                                                      !!
  // !! Not valid: Case 1, with angles in the LAB            !!
  // !!            Case 2                                    !!
  // !!                                                      !!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // Set all pointers to independent variables to NULL                 
  fE1Var = fE2Var = fEVar = fTh1Var = fTh2Var = fPhi2Var = 0;

  // Put format string in member and in a temporary string.
  std::strcpy(fFormat,format);
  const std::string formatString = format;

  // Tokenize format string (empty tokens are left out, as with ROOT's Tokenize before 2026-11)
  std::vector<std::string> tokens = nucqcd::split(formatString, ':');

  // Check the number of tokens
  int nrOfTokens = tokens.size();
  if( nrOfTokens != 6 ) {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	 << "Format string should provide 6 independent variables: "
	 << "\"var1:var2:var3:var4:var5:var6\"");
  }

  // Transform tokens to lower case
  for(std::string &tok : tokens)
    tok = nucqcd::to_lower(tok);

  // Check that Q^2 is present
  if( !nucqcd::icontains(formatString,"qsquared") ) {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	 << "Q^2 needs to be one of 6 independent kinematic variables.");
  }

  // Check whether all tokens are mandelstam variables. We do this by
  // counting the s's and t's
  // mandel: the subscenarios a-g still possible ('-' for excluded ones), only in the mandelstam case
  std::optional<std::string> mandel;
  auto count = [&formatString](char c) { return (int)std::count(formatString.begin(), formatString.end(), c); };
  if( count('s')+count('S')-1
      +count('t')+count('T')-2 ==5 
      && nucqcd::icontains(formatString,"stot") ) {

    fKinSolution[0] = '3';

    mandel = std::string("abcdefg"); // list all subscenarios
    
    // find out name of '2nd particle'
    if( fAC == kK ) {
      if( nucqcd::icontains(formatString,"sky") ||
	  nucqcd::icontains(formatString,"syk") ||
	  nucqcd::icontains(formatString,"tgy") ||
	  nucqcd::icontains(formatString,"tyg") )
	fKinSolution[1] = 'Y';
      else
	fKinSolution[1] = 'N';
    } // fAC == kK
    else if( fAC == kY ) {
      if( nucqcd::icontains(formatString,"sky") ||
	  nucqcd::icontains(formatString,"syk") ||
	  nucqcd::icontains(formatString,"tgk") ||
	  nucqcd::icontains(formatString,"tkg") )
	fKinSolution[1] = 'K';
      else
	fKinSolution[1] = 'N';
    } // fAC == kY
    else if( fAC == kN ) {
      if( nucqcd::icontains(formatString,"sny") ||
	  nucqcd::icontains(formatString,"syn") ||
	  nucqcd::icontains(formatString,"tgy") ||
	  nucqcd::icontains(formatString,"tyg") )
	fKinSolution[1] = 'Y';
      else
	fKinSolution[1] = 'K';
    } // fAC == kN
  } // end mandelstam case

  // Determine independent kinematic variables
  for(int i=0; i<nrOfTokens; ++i) 
    {
      const std::string &token = tokens[i];  // (was copied into char[11]: overflow for long tokens)

      //__________________________________
      if( (token=="qsquared") ) {
	fVar[i] = &fQsquared;
      }
      //__________________________________
      else if( (token=="wlab") ) {
	if(fEVar || fKinSolution[0] == '2') {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	fVar[i] = fEVar = &fWlab;
	fKinSolution[0] = '1';
      }
      //__________________________________
      else if( (token=="klab") ) {
	if(fEVar || fKinSolution[0] == '2') {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	fVar[i] = fEVar = &fKlab;
	fKinSolution[0] = '1';
      }
      //__________________________________
      else if( (token=="ek") ) {
	if( fAC == kK ) {
	  if(fE1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE1Var = &fEk;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE2Var = &fEk;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'K';
	}
      }
      //__________________________________
      else if( (token=="ey") ) {
	if( fAC == kY ) {
	  if(fE1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE1Var = &fEy;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE2Var = &fEy;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'Y';
	}	
      }
      //__________________________________
      else if( (token=="en") ) {
	if( fAC == kN ) {
	  if(fE1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE1Var = &fEn;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE2Var = &fEn;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'N';
	}	
      }
      //__________________________________
      else if( (token=="pk") ) {
	if( fAC == kK ) {
	  if(fE1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE1Var = &fPk;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE2Var = &fPk;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'K';
	}
      }
      //__________________________________
      else if( (token=="py") ) {
	if( fAC == kY ) {
	  if(fE1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE1Var = &fPy;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE2Var = &fPy;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'Y';
	}	
      }
      //__________________________________
      else if( (token=="pn") ) {
	if( fAC == kN ) {
	  if(fE1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE1Var = &fPn;
	}
	else {
	  if(fE2Var || fKinSolution[0] == '1' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fE2Var = &fEn;
	  fKinSolution[0] = '2';
	  fKinSolution[1] = 'P';
	}	
      }
      //__________________________________
      else if( (token=="thk") ) {
	if( fAC == kK ) {
	  if(fTh1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh1Var = &fThk;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fThk;
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="thy") ) {
	if( fAC == kY ) {
	  if(fTh1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh1Var = &fThy;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fThy;
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="thn") ) {
	if( fAC == kN ) {
	  if(fTh1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh1Var = &fThn;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fThn;
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="costhk") ) {
	if( fAC == kK ) {
	  if(fTh1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh1Var = &fCosthk;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fCosthk;
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="costhy") ) {
	if( fAC == kY ) {
	  if(fTh1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh1Var = &fCosthy;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fCosthy;
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="costhn") ) {
	if( fAC == kN ) {
	  if(fTh1Var) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh1Var = &fCosthn;
	}
	else {
	  if(fTh2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fCosthn;
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="phik") ) {
	if( fAC == kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane! Because of the axis convention "
	       << "phik is zero by definition.");
	}
	else {
	  if(fPhi2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fPhik;
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="phiy") ) {
	if( fAC == kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane! Because of the axis convention "
	       << "phiy is zero by definition.");
	}
	else {
	  if(fPhi2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fPhiy;
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="phin") ) {
	if( fAC == kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane! Because of the axis convention "
	       << "phin is zero by definition.");
	}
	else {
	  if(fPhi2Var || fKinSolution[2] == '*' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fPhin;
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = 'L';
	}
      }
      //__________________________________
      else if( (token=="wtot") ) {
	if(fEVar || fKinSolution[0] == '2') {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	fVar[i] = fEVar = &fWtot;
	fKinSolution[0] = '1';
      }
      //__________________________________
      else if( (token=="wky") || (token=="wyk") ) {
	if( fAC != kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fEVar || fKinSolution[0] == '2') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fEVar = &fWky;
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( (token=="wkn") || (token=="wnk") ) {
	if( fAC != kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fEVar || fKinSolution[0] == '2') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fEVar = &fWkn;
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( (token=="wyn") || (token=="wny") ) {
	if( fAC != kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fEVar || fKinSolution[0] == '2') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fEVar = &fWyn;
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( (token=="stot") ) {
	if(fEVar || fKinSolution[0] == '2') {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else if(fKinSolution[0] == '3') {
	  fVar[i] = &fStot;
	}
	else {
	  fVar[i] = fEVar = &fStot;
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( (token=="sky") || (token=="syk") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kK && fKinSolution[1]=='Y') ||
	      (fAC==kY && fKinSolution[1]=='K') ) {
	    std::replace(mandel->begin(),mandel->end(),'g','-');
	  }
	  else if( fAC==kN ) {
	    std::replace(mandel->begin(),mandel->end(),'d','-');
	  }
	  else {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = &fSky;
	} // end mandelstam case
	else {
	  if( fAC != kN ) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  else {
	    if(fEVar || fKinSolution[0] == '2') {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fEVar = &fSky;
	    fKinSolution[0] = '1';
	  }
	}
      }
      //__________________________________
      else if( (token=="skn") || (token=="snk") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kK && fKinSolution[1]=='N') ||
	      (fAC==kN && fKinSolution[1]=='K') ) {
	    std::replace(mandel->begin(),mandel->end(),'g','-');
	  }
	  else if( fAC==kY ) {
	    std::replace(mandel->begin(),mandel->end(),'d','-');
	  }
	  else {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = &fSkn;
	} // end mandelstam case
	else {
	  if( fAC != kY ) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  else {
	    if(fEVar || fKinSolution[0] == '2') {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fEVar = &fSkn;
	    fKinSolution[0] = '1';
	  }
	}
      }
      //__________________________________
      else if( (token=="syn") || (token=="sny") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kY && fKinSolution[1]=='N') ||
	      (fAC==kN && fKinSolution[1]=='Y') ) {
	    std::replace(mandel->begin(),mandel->end(),'g','-');
	  }
	  else if( fAC==kK ) {
	    std::replace(mandel->begin(),mandel->end(),'d','-');
	  }
	  else {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = &fSyn;
	} // end mandelstam case
	else {
	  if( fAC != kK ) {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  else {
	    if(fEVar || fKinSolution[0] == '2') {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fEVar = &fSyn;
	    fKinSolution[0] = '1';
	  }
	}
      }
      //__________________________________
      else if( (token=="tgk") || (token=="tkg") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kY || fAC==kN ) && fKinSolution[1]=='K' ) {
	    std::replace(mandel->begin(),mandel->end(),'a','-');
	    std::replace(mandel->begin(),mandel->end(),'b','-');
	    std::replace(mandel->begin(),mandel->end(),'c','-');
	    std::replace(mandel->begin(),mandel->end(),'d','-');
	  }
	  else if( fAC==kK ) {
	    std::replace(mandel->begin(),mandel->end(),'c','-');
	    std::replace(mandel->begin(),mandel->end(),'f','-');
	  }
	  else {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = &fTgk;
	} // end mandelstam case
	else if (fKinSolution[0] == '2') {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if( fAC == kK ) {
	    if(fTh1Var) {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fTh1Var = &fTgk;
	  }
	  else {
	    if(fTh2Var || fKinSolution[2] == 'L' 
	       || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fTh2Var = &fTgk;
	    fKinSolution[1] = 'K';
	    fKinSolution[2] = '*';
	  }
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( (token=="tgy") || (token=="tyg") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kK || fAC==kN ) && fKinSolution[1]=='Y' ) {
	    std::replace(mandel->begin(),mandel->end(),'a','-');
	    std::replace(mandel->begin(),mandel->end(),'b','-');
	    std::replace(mandel->begin(),mandel->end(),'c','-');
	    std::replace(mandel->begin(),mandel->end(),'d','-');
	  }
	  else if( fAC==kY ) {
	    std::replace(mandel->begin(),mandel->end(),'c','-');
	    std::replace(mandel->begin(),mandel->end(),'f','-');
	  }
	  else {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = &fTgy;
	} // end mandelstam case
	else if (fKinSolution[0] == '2') {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if( fAC == kY ) {
	    if(fTh1Var) {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fTh1Var = &fTgy;
	  }
	  else {
	    if(fTh2Var || fKinSolution[2] == 'L' 
	       || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fTh2Var = &fTgy;
	    fKinSolution[1] = 'Y';
	    fKinSolution[2] = '*';
	  }
	  fKinSolution[0] = '1'; 
	}
      }
      //__________________________________
      else if( (token=="tgn") || (token=="tng") ) {
	if(fKinSolution[0] == '3') {
	  if( (fAC==kK || fAC==kY ) && fKinSolution[1]=='N' ) {
	    std::replace(mandel->begin(),mandel->end(),'a','-');
	    std::replace(mandel->begin(),mandel->end(),'b','-');
	    std::replace(mandel->begin(),mandel->end(),'c','-');
	    std::replace(mandel->begin(),mandel->end(),'d','-');
	  }
	  else if( fAC==kN ) {
	    std::replace(mandel->begin(),mandel->end(),'c','-');
	    std::replace(mandel->begin(),mandel->end(),'f','-');
	  }
	  else {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = &fTgn;
	} // end mandelstam case
	else if (fKinSolution[0] == '2') {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if( fAC == kN ) {
	    if(fTh1Var) {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fTh1Var = &fTgn;
	  }
	  else {
	    if(fTh2Var || fKinSolution[2] == 'L' 
	       || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		   << "Format string is not sane!");
	    }
	    fVar[i] = fTh2Var = &fTgn;
	    fKinSolution[1] = 'N';
	    fKinSolution[2] = '*';
	  }
	  fKinSolution[0] = '1';
	}
      }
      //__________________________________
      else if( (token=="tdk") || (token=="tkd") ) {
	if( fKinSolution[0] != '3' || fAC==kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else if( fKinSolution[1] == 'K' ) {
	  std::replace(mandel->begin(),mandel->end(),'a','-');
	  std::replace(mandel->begin(),mandel->end(),'e','-');
	  std::replace(mandel->begin(),mandel->end(),'f','-');
	  std::replace(mandel->begin(),mandel->end(),'g','-');
	}
	else {
	  std::replace(mandel->begin(),mandel->end(),'b','-');
	  std::replace(mandel->begin(),mandel->end(),'e','-');
	}
	fVar[i] = &fTdk;
      }
      //__________________________________
      else if( (token=="tdy") || (token=="tyd") ) {
	if( fKinSolution[0] != '3' || fAC==kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else if( fKinSolution[1] == 'Y' ) {
	  std::replace(mandel->begin(),mandel->end(),'a','-');
	  std::replace(mandel->begin(),mandel->end(),'e','-');
	  std::replace(mandel->begin(),mandel->end(),'f','-');
	  std::replace(mandel->begin(),mandel->end(),'g','-');
	}
	else {
	  std::replace(mandel->begin(),mandel->end(),'b','-');
	  std::replace(mandel->begin(),mandel->end(),'e','-');
	}
	fVar[i] = &fTdy;
      }
      //__________________________________
      else if( (token=="tdn") || (token=="tnd") ) {
	if( fKinSolution[0] != '3' || fAC==kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else if( fKinSolution[1] == 'N' ) {
	  std::replace(mandel->begin(),mandel->end(),'a','-');
	  std::replace(mandel->begin(),mandel->end(),'e','-');
	  std::replace(mandel->begin(),mandel->end(),'f','-');
	  std::replace(mandel->begin(),mandel->end(),'g','-');
	}
	else {
	  std::replace(mandel->begin(),mandel->end(),'b','-');
	  std::replace(mandel->begin(),mandel->end(),'e','-');
	}
	fVar[i] = &fTdn;
      }
      //__________________________________
      else if( (token=="kythkcm") || (token=="ykthkcm") ) {
	if( fAC != kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fKYthkcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="kythycm") || (token=="ykthycm") ) {
	if( fAC != kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fKYthycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="knthkcm") || (token=="nkthkcm") ) {
	if( fAC != kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fKNthkcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="knthncm") || (token=="nkthncm") ) {
	if( fAC != kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fKNthncm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="ynthycm") || (token=="nythycm") ) {
	if( fAC != kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fYNthycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="ynthncm") || (token=="nythncm") ) {
	if( fAC != kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fYNthncm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="kycosthkcm") || (token=="ykcosthkcm") ) {
	if( fAC != kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fKYcosthkcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="kycosthycm") || (token=="ykcosthycm") ) {
	if( fAC != kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fKYcosthycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="kncosthkcm") || (token=="nkcosthkcm") ) {
	if( fAC != kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fKNcosthkcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="kncosthncm") || (token=="nkcosthncm") ) {
	if( fAC != kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fKNcosthncm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="yncosthycm") || (token=="nycosthycm") ) {
	if( fAC != kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fYNcosthycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="yncosthncm") || (token=="nycosthncm") ) {
	if( fAC != kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fTh2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fTh2Var = &fYNcosthncm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="kyphikcm") || (token=="ykphikcm") ) {
	if( fAC != kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fKYphikcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="kyphiycm") || (token=="ykphiycm") ) {
	if( fAC != kN ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fKYphiycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="knphikcm") || (token=="nkphikcm") ) {
	if( fAC != kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fKNphikcm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'K';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="knphincm") || (token=="nkphincm") ) {
	if( fAC != kY ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'K' || fKinSolution[1] == 'Y') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fKNphincm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="ynphiycm") || (token=="nyphiycm") ) {
	if( fAC != kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'N' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fYNphiycm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'Y';
	  fKinSolution[2] = '*';
	}
      }
      //__________________________________
      else if( (token=="ynphincm") || (token=="nyphincm") ) {
	if( fAC != kK ) {
	  throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	       << "Format string is not sane!");
	}
	else {
	  if(fPhi2Var || fKinSolution[0] == 2 || fKinSolution[2] == 'L' 
	     || fKinSolution[1] == 'Y' || fKinSolution[1] == 'K') {
	    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
		 << "Format string is not sane!");
	  }
	  fVar[i] = fPhi2Var = &fYNphincm;
	  fKinSolution[0] = '1';
	  fKinSolution[1] = 'N';
	  fKinSolution[2] = '*';
	}
      }

      //__________________________________
      else {
	fVar[i] = NULL;
	throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in  TKinematics2to3::ReadFormatString(const char*): "
	     << "Kinematic variable \"" << token << "\" unknown!");
      }
    } // end loop over tokens
  
  // Check whether the user provided a meaningful format string
  // in the mandelstam case
  if(mandel) {
    if(std::count(mandel->begin(),mandel->end(),'-') != 6) {
      throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	   << "Format string is not sane!");
    }
    fKinSolution[2] = (*mandel)[mandel->find_first_not_of('-')];  // the one subscenario left
  }

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!                                                      !!
  // !! ATTENTION                                            !!
  // !! some cases are not valid anymore, until I find a way !!
  // !! to solve them without ambiguities.                   !!
  // !!                                                      !!
  // !! Not valid: Case 1, with angles in the LAB            !!
  // !!            Case 2                                    !!
  // !!                                                      !!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if( fKinSolution[0] == '2' ) {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	 << "I cannot solve kinematics with these variables.");
  }
  else if( fKinSolution[0] == '1' && fKinSolution[2] == 'L' )
    cerr << "WARNING in TKinematics2to3::ReadFormatString(const char*): "
	 << "CAVEAT: risk of running into multiple solutions.\n"
	 << "Perhaps TKinematics2to3WithLabAngles is something up your alley.\n";

}

}  // namespace nucqcd
