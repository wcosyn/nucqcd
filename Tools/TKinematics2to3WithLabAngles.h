/*! \file TKinematics2to3WithLabAngles.h
 * \brief 2 -> 3 kinematics with lab-frame angles of two final-state particles as input
 * \author Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be), work started on August 20, 2009
 *
 * \addtogroup MePhys
 * @{
 */

#ifndef TKINEMATICS2TO3WITHLABANGLES_H
#define TKINEMATICS2TO3WITHLABANGLES_H

#include "TKinematics2to3.h"
#include <fstream>

/*! \brief TKinematics2to3 with the lab momentum and angles of particle 1 and the lab angles of particle
 * 2 as input; the energy conservation then has 0, 1 or 2 solutions (GetSolution), the object itself
 * holds the first one (e.g. A(e,e'pN) with format "qsquared:wlab:pn:thn:thk:phik")
 */
class TKinematics2to3WithLabAngles : public TKinematics2to3
{
 public:
  /*! constructor, arguments as for TKinematics2to3::TKinematics2to3 */
  TKinematics2to3WithLabAngles(const char *name, const char *title,
			       int isospin, EAxisConvention ac, 
			       const char *format,
			       double var1, double var2, double var3,
			       double var4, double var5, double var6,
			       double Md =0.,double Mn=0., double Mk=0.,double My=0.);
  TKinematics2to3WithLabAngles(const TKinematics2to3WithLabAngles&); ///< copy constructor
  TKinematics2to3WithLabAngles& operator=(const TKinematics2to3WithLabAngles&); ///< assignment
  virtual ~TKinematics2to3WithLabAngles(); ///< destructor

  int GetNumberOfSolutions() const { return fNrOfSolutions; } ///< number of solutions (0, 1 or 2) for the current point
  const TKinematics2to3& GetSolution(int) const; ///< solution 0 or 1
  virtual int     GetNumberOfPhysicalSteps() const; ///< number of physical points in the grid
  virtual double *GetVarArray(const TString& varName) const; ///< values of a variable at all grid points; the caller owns the array
  virtual double *GetPhysicalVarArray(const TString& varName) const; ///< values of a variable at the physical grid points; the caller owns the array
 
 protected:
  virtual void ReadFormatString(const char*); ///< processes the format string
  virtual void UpdateKinematics(); ///< solves the kinematics and builds the solutions
  static void TurnOutputOff(bool=true); ///< silences (or restores) output during the solution

 private:
  /// number of solutions (0, 1 or 2) for the current point
  int fNrOfSolutions;           //! 0, 1 or 2 solutions for current Omega_2
  /// the solutions
  TKinematics2to3 **fSolutions; //! the real solutions
  /// stream to /dev/null
  static const std::ofstream fgDevNull; // stream to /dev/null

  
}; // class TKinematics2to3WithLabAngles

/** @} */
#endif // TKINEMATICS2TO3WITHLABANGLES_H
