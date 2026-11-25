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

  // Determine independent kinematic variables.
  // The tokens fall into a few groups that are handled alike; the helpers below hold the checks of a
  // group (the same ones as the former 1000 lines of one branch per token, including their quirks, see
  // the comments there, and the checks of the exceptions happen in the same order). A format string is
  // "not sane" when a variable is given twice, or conflicts with another given one.
  auto insane = [](const char *extra = "") {
    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
                                                       << "Format string is not sane!" << extra);
  };
  // a mandelstam subscenario is excluded by '-'
  auto exclude = [&mandel](const char *scenarios) {
    for(const char *c = scenarios; *c; ++c) std::replace(mandel->begin(), mandel->end(), *c, '-');
  };
  // the variable that fixes the total energy of the system (scenario 1): wlab, klab, wtot, ...
  auto energyVar = [&](int i, double *var) {
    if(fEVar || fKinSolution[0] == '2') insane();
    fVar[i] = fEVar = var;
    fKinSolution[0] = '1';
  };
  // ... which for wky, wkn, wyn is only available with one axis convention
  auto energyVarFor = [&](int i, double *var, EAxisConvention ac) {
    if(fAC != ac) insane();
    energyVar(i, var);
  };
  // energy or momentum of particle `particle` (ac: the axis convention that fixes it): the variable of the
  // particle on the x axis (fE1Var), or of the second particle (fE2Var; `var2` is stored, `solution` is
  // put in fKinSolution[1]; ex1, ex2: the other two particles)
  auto particleEnergy = [&](int i, EAxisConvention ac, double *var1, double *var2, char ex1, char ex2, char solution) {
    if(fAC == ac) {
      if(fE1Var) insane();
      fVar[i] = fE1Var = var1;
    }
    else {
      if(fE2Var || fKinSolution[0] == '1' || fKinSolution[1] == ex1 || fKinSolution[1] == ex2) insane();
      fVar[i] = fE2Var = var2;
      fKinSolution[0] = '2';
      fKinSolution[1] = solution;
    }
  };
  // polar angle (or cosine) in the lab frame of particle `solution` (ac: its axis convention)
  auto labTheta = [&](int i, EAxisConvention ac, double *var, char ex1, char ex2, char solution) {
    if(fAC == ac) {
      if(fTh1Var) insane();
      fVar[i] = fTh1Var = var;
    }
    else {
      if(fTh2Var || fKinSolution[2] == '*' || fKinSolution[1] == ex1 || fKinSolution[1] == ex2) insane();
      fVar[i] = fTh2Var = var;
      fKinSolution[1] = solution;
      fKinSolution[2] = 'L';
    }
  };
  // azimuth in the lab frame of particle `solution`: zero by definition for the axis convention ac
  auto labPhi = [&](int i, EAxisConvention ac, double *var, const char *zeroMessage, char ex1, char ex2, char solution) {
    if(fAC == ac) insane(zeroMessage);
    if(fPhi2Var || fKinSolution[2] == '*' || fKinSolution[1] == ex1 || fKinSolution[1] == ex2) insane();
    fVar[i] = fPhi2Var = var;
    fKinSolution[1] = solution;
    fKinSolution[2] = 'L';
  };
  // angle in the CM of a subsystem (th2 or phi2, `slot`), only with the axis convention ac; scenario 1
  auto cmAngle = [&](int i, EAxisConvention ac, double *&slot, double *var, char ex1, char ex2, char solution) {
    if(fAC != ac) insane();
    // (fKinSolution[0] == 2 compares with the number 2, which never happens: as before)
    if(slot || fKinSolution[0] == 2 || fKinSolution[2] == 'L' || fKinSolution[1] == ex1 || fKinSolution[1] == ex2) insane();
    fVar[i] = slot = var;
    fKinSolution[0] = '1';
    fKinSolution[1] = solution;
    fKinSolution[2] = '*';
  };
  // s_KY, s_KN, s_YN: a mandelstam variable (scenario 3: removes subscenarios, `g` if the axis convention
  // and the second particle are one of the two (ac1, sol1) or (ac2, sol2), else `d` if the axis
  // convention is acRest), otherwise the variable of the total energy of the pair, only with the axis
  // convention acOnly
  auto pairS = [&](int i, double *var, EAxisConvention ac1, char sol1, EAxisConvention ac2, char sol2,
                   EAxisConvention acRest, EAxisConvention acOnly) {
    if(fKinSolution[0] == '3') {
      if( (fAC == ac1 && fKinSolution[1] == sol1) || (fAC == ac2 && fKinSolution[1] == sol2) ) exclude("g");
      else if( fAC == acRest ) exclude("d");
      else insane();
      fVar[i] = var;
    }
    else energyVarFor(i, var, acOnly);
  };
  // t_gK, t_gY, t_gN: a mandelstam variable (scenario 3) or the angle-like variable of particle P
  auto tg = [&](int i, double *var, EAxisConvention ac, char particle, char ex1, char ex2) {
    if(fKinSolution[0] == '3') {
      if( fAC != ac && fKinSolution[1] == particle ) exclude("abcd");
      else if( fAC == ac ) exclude("cf");
      else insane();
      fVar[i] = var;
    }
    else if(fKinSolution[0] == '2') insane();
    else {
      if( fAC == ac ) {
        if(fTh1Var) insane();
        fVar[i] = fTh1Var = var;
      }
      else {
        if(fTh2Var || fKinSolution[2] == 'L' || fKinSolution[1] == ex1 || fKinSolution[1] == ex2) insane();
        fVar[i] = fTh2Var = var;
        fKinSolution[1] = particle;
        fKinSolution[2] = '*';
      }
      fKinSolution[0] = '1';
    }
  };
  // t_dK, t_dY, t_dN: only a mandelstam variable, not for the axis convention ac of its particle
  auto td = [&](int i, double *var, EAxisConvention ac, char particle) {
    if(fKinSolution[0] != '3' || fAC == ac) insane();
    else if(fKinSolution[1] == particle) exclude("aefg");
    else exclude("be");
    fVar[i] = var;
  };

  for(int i=0; i<nrOfTokens; ++i) 
    {
      const std::string &token = tokens[i];  // (was copied into char[11]: overflow for long tokens)
      auto is = [&token](const char *a, const char *b = nullptr) { return token == a || (b && token == b); };

      if( is("qsquared") ) fVar[i] = &fQsquared;
      else if( is("wlab") ) energyVar(i, &fWlab);
      else if( is("klab") ) energyVar(i, &fKlab);
      else if( is("ek") ) particleEnergy(i, kK, &fEk, &fEk, 'N', 'Y', 'K');
      else if( is("ey") ) particleEnergy(i, kY, &fEy, &fEy, 'N', 'K', 'Y');
      else if( is("en") ) particleEnergy(i, kN, &fEn, &fEn, 'Y', 'K', 'N');
      else if( is("pk") ) particleEnergy(i, kK, &fPk, &fPk, 'N', 'Y', 'K');
      else if( is("py") ) particleEnergy(i, kY, &fPy, &fPy, 'N', 'K', 'Y');
      // (pn with another axis convention stores fEn and the name 'P': as before; the others store their own)
      else if( is("pn") ) particleEnergy(i, kN, &fPn, &fEn, 'Y', 'K', 'P');
      else if( is("thk") ) labTheta(i, kK, &fThk, 'N', 'Y', 'K');
      else if( is("thy") ) labTheta(i, kY, &fThy, 'N', 'K', 'Y');
      else if( is("thn") ) labTheta(i, kN, &fThn, 'K', 'Y', 'N');
      else if( is("costhk") ) labTheta(i, kK, &fCosthk, 'N', 'Y', 'K');
      else if( is("costhy") ) labTheta(i, kY, &fCosthy, 'N', 'K', 'Y');
      else if( is("costhn") ) labTheta(i, kN, &fCosthn, 'K', 'Y', 'N');
      else if( is("phik") ) labPhi(i, kK, &fPhik, " Because of the axis convention phik is zero by definition.", 'N', 'Y', 'K');
      else if( is("phiy") ) labPhi(i, kY, &fPhiy, " Because of the axis convention phiy is zero by definition.", 'N', 'K', 'Y');
      else if( is("phin") ) labPhi(i, kN, &fPhin, " Because of the axis convention phin is zero by definition.", 'K', 'Y', 'N');
      else if( is("wtot") ) energyVar(i, &fWtot);
      else if( is("wky", "wyk") ) energyVarFor(i, &fWky, kN);
      else if( is("wkn", "wnk") ) energyVarFor(i, &fWkn, kY);
      else if( is("wyn", "wny") ) energyVarFor(i, &fWyn, kK);
      else if( is("stot") ) {
        if(fEVar || fKinSolution[0] == '2') insane();
        else if(fKinSolution[0] == '3') fVar[i] = &fStot;
        else energyVar(i, &fStot);
      }
      else if( is("sky", "syk") ) pairS(i, &fSky, kK, 'Y', kY, 'K', kN, kN);
      else if( is("skn", "snk") ) pairS(i, &fSkn, kK, 'N', kN, 'K', kY, kY);
      else if( is("syn", "sny") ) pairS(i, &fSyn, kY, 'N', kN, 'Y', kK, kK);
      else if( is("tgk", "tkg") ) tg(i, &fTgk, kK, 'K', 'N', 'Y');
      else if( is("tgy", "tyg") ) tg(i, &fTgy, kY, 'Y', 'N', 'K');
      else if( is("tgn", "tng") ) tg(i, &fTgn, kN, 'N', 'K', 'Y');
      else if( is("tdk", "tkd") ) td(i, &fTdk, kK, 'K');
      else if( is("tdy", "tyd") ) td(i, &fTdy, kY, 'Y');
      else if( is("tdn", "tnd") ) td(i, &fTdn, kN, 'N');
      // CM angles of a subsystem: theta (the second particle), cos theta, phi
      else if( is("kythkcm", "ykthkcm") ) cmAngle(i, kN, fTh2Var, &fKYthkcm, 'N', 'Y', 'K');
      else if( is("kythycm", "ykthycm") ) cmAngle(i, kN, fTh2Var, &fKYthycm, 'N', 'K', 'Y');
      else if( is("knthkcm", "nkthkcm") ) cmAngle(i, kY, fTh2Var, &fKNthkcm, 'N', 'Y', 'K');
      else if( is("knthncm", "nkthncm") ) cmAngle(i, kY, fTh2Var, &fKNthncm, 'K', 'Y', 'N');
      else if( is("ynthycm", "nythycm") ) cmAngle(i, kK, fTh2Var, &fYNthycm, 'N', 'K', 'Y');
      else if( is("ynthncm", "nythncm") ) cmAngle(i, kK, fTh2Var, &fYNthncm, 'Y', 'K', 'N');
      else if( is("kycosthkcm", "ykcosthkcm") ) cmAngle(i, kN, fTh2Var, &fKYcosthkcm, 'N', 'Y', 'K');
      else if( is("kycosthycm", "ykcosthycm") ) cmAngle(i, kN, fTh2Var, &fKYcosthycm, 'N', 'K', 'Y');
      else if( is("kncosthkcm", "nkcosthkcm") ) cmAngle(i, kY, fTh2Var, &fKNcosthkcm, 'N', 'Y', 'K');
      else if( is("kncosthncm", "nkcosthncm") ) cmAngle(i, kY, fTh2Var, &fKNcosthncm, 'K', 'Y', 'N');
      else if( is("yncosthycm", "nycosthycm") ) cmAngle(i, kK, fTh2Var, &fYNcosthycm, 'N', 'K', 'Y');
      else if( is("yncosthncm", "nycosthncm") ) cmAngle(i, kK, fTh2Var, &fYNcosthncm, 'Y', 'K', 'N');
      else if( is("kyphikcm", "ykphikcm") ) cmAngle(i, kN, fPhi2Var, &fKYphikcm, 'N', 'Y', 'K');
      else if( is("kyphiycm", "ykphiycm") ) cmAngle(i, kN, fPhi2Var, &fKYphiycm, 'N', 'K', 'Y');
      else if( is("knphikcm", "nkphikcm") ) cmAngle(i, kY, fPhi2Var, &fKNphikcm, 'N', 'Y', 'K');
      else if( is("knphincm", "nkphincm") ) cmAngle(i, kY, fPhi2Var, &fKNphincm, 'K', 'Y', 'N');
      else if( is("ynphiycm", "nyphiycm") ) cmAngle(i, kK, fPhi2Var, &fYNphiycm, 'N', 'K', 'Y');
      else if( is("ynphincm", "nyphincm") ) cmAngle(i, kK, fPhi2Var, &fYNphincm, 'Y', 'K', 'N');
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
  else if( fKinSolution[0] == '1' && !(fE1Var && fTh1Var && fEVar && fTh2Var && fPhi2Var) ) {
    // Scenario 1 needs all of its variables: the energy of the system (fEVar), the energy and the angle of
    // the particle on the x axis (fE1Var, fTh1Var) and the angles of the second particle (fTh2Var, fPhi2Var).
    // Until 2026-11 a missing one (a variable given twice, or two energies and no angle, ...) was a null
    // pointer that UpdateKinematics dereferenced: a segmentation fault.
    throw std::invalid_argument(nucqcd::ErrorMessage() << "ERROR in TKinematics2to3::ReadFormatString(const char*): "
	 << "Format string is not sane! (Variables are missing: it needs an energy of the system, the energy and the "
	 << "polar angle of the first particle and the polar and azimuthal angle of the second one.)");
  }
  else if( fKinSolution[0] == '1' && fKinSolution[2] == 'L' )
    cerr << "WARNING in TKinematics2to3::ReadFormatString(const char*): "
	 << "CAVEAT: risk of running into multiple solutions.\n"
	 << "Perhaps TKinematics2to3WithLabAngles is something up your alley.\n";

}

}  // namespace nucqcd
