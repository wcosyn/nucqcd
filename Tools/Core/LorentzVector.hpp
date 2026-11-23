/*! \file LorentzVector.hpp
 * \brief Four-vectors (nucqcd::LorentzVector), replacing ROOT's TLorentzVector (roadmap step 9,
 * phase 4).
 * \details Same interface (the subset nucqcd uses), conventions and formulas as ROOT 6.40; see
 * Vector3.hpp for how the results are kept bitwise the same as ROOT's. Metric (+,-,-,-): Mag2() =
 * t^2 - |p|^2, and Mag() = -sqrt(-Mag2()) for spacelike vectors. Boost(b) is an active boost with
 * velocity b; BoostVector() = p/E, so v.Boost(-v.BoostVector()) takes v to its rest frame.
 * Differences with ROOT: no TObject base; a bad index throws std::out_of_range.
 * \author Wim Cosyn (ROOT's TLorentzVector by Pasha Murat, Peter Malzacher and others)
 */
#ifndef NUCQCD_LORENTZVECTOR_HPP
#define NUCQCD_LORENTZVECTOR_HPP

#include <Tools/Core/Vector3.hpp>

#include <cmath>
#include <iosfwd>

namespace nucqcd {

class Rotation;
class LorentzRotation;

/*! \brief A four-vector (x, y, z, t), with the interface of ROOT's TLorentzVector that nucqcd uses */
class LorentzVector {
 public:
  LorentzVector() : fP(), fE(0.0) {}                                             ///< null vector
  LorentzVector(double x, double y, double z, double t) : fP(x, y, z), fE(t) {}  ///< vector (x, y, z, t)
  LorentzVector(const Vector3 &p, double e) : fP(p), fE(e) {}                    ///< vector (p, e)
  explicit LorentzVector(const double *a) : fP(a), fE(a[3]) {}                   ///< vector (a[0], a[1], a[2], a[3])

  double operator()(int i) const;  ///< component i (0, 1, 2, 3 = x, y, z, t); throws std::out_of_range otherwise
  double &operator()(int i);       ///< component i (0, 1, 2, 3 = x, y, z, t); throws std::out_of_range otherwise
  double operator[](int i) const { return (*this)(i); }  ///< component i, as operator()
  double &operator[](int i) { return (*this)(i); }       ///< component i, as operator()

  double X() const { return fP.X(); }   ///< x component
  double Y() const { return fP.Y(); }   ///< y component
  double Z() const { return fP.Z(); }   ///< z component
  double T() const { return fE; }       ///< time component
  double Px() const { return X(); }     ///< x component
  double Py() const { return Y(); }     ///< y component
  double Pz() const { return Z(); }     ///< z component
  double E() const { return T(); }      ///< time component
  double Energy() const { return T(); } ///< time component
  double P() const { return fP.Mag(); } ///< length of the three-vector
  Vector3 Vect() const { return fP; }   ///< three-vector part
  void SetX(double a) { fP.SetX(a); }   ///< sets the x component
  void SetY(double a) { fP.SetY(a); }   ///< sets the y component
  void SetZ(double a) { fP.SetZ(a); }   ///< sets the z component
  void SetT(double a) { fE = a; }       ///< sets the time component
  void SetPx(double a) { SetX(a); }     ///< sets the x component
  void SetPy(double a) { SetY(a); }     ///< sets the y component
  void SetPz(double a) { SetZ(a); }     ///< sets the z component
  void SetE(double a) { SetT(a); }      ///< sets the time component
  void SetVect(const Vector3 &p) { fP = p; }  ///< sets the three-vector part
  void SetXYZT(double x, double y, double z, double t) { fP.SetXYZ(x, y, z); SetT(t); }  ///< sets all components
  void SetPxPyPzE(double px, double py, double pz, double e) { SetXYZT(px, py, pz, e); }  ///< sets all components
  /*! sets the three-vector and the energy sqrt(p^2 + m^2) (for m < 0: sqrt(max(p^2 - m^2, 0))) */
  void SetXYZM(double x, double y, double z, double m) {
    if (m >= 0)
      SetXYZT(x, y, z, std::sqrt(x * x + y * y + z * z + m * m));
    else
      {
      double p2m2 = x * x + y * y + z * z - m * m;  // TMath::Max(p2m2, 0.)
      SetXYZT(x, y, z, std::sqrt(p2m2 >= 0. ? p2m2 : 0.));
    }
  }
  void SetVectM(const Vector3 &p, double m) { SetXYZM(p.X(), p.Y(), p.Z(), m); }  ///< three-vector p and mass m
  void SetVectMag(const Vector3 &p, double m) { SetXYZM(p.X(), p.Y(), p.Z(), m); }  ///< three-vector p and mass m
  void GetXYZT(double *a) const { fP.GetXYZ(a); a[3] = fE; }  ///< copies the components into a[0..3]

  double Phi() const { return fP.Phi(); }            ///< azimuth of the three-vector, in [-pi, pi]
  double Theta() const { return fP.Theta(); }        ///< polar angle of the three-vector
  double CosTheta() const { return fP.CosTheta(); }  ///< cos of the polar angle of the three-vector
  double Rho() const { return fP.Mag(); }            ///< length of the three-vector
  void SetTheta(double th) { fP.SetTheta(th); }      ///< sets the polar angle of the three-vector
  void SetPhi(double phi) { fP.SetPhi(phi); }        ///< sets the azimuth of the three-vector
  void SetRho(double rho) { fP.SetMag(rho); }        ///< sets the length of the three-vector
  double Perp2() const { return fP.Perp2(); }        ///< squared transverse component (to z)
  double Perp() const { return fP.Perp(); }          ///< transverse component (to z)
  double Pt() const { return Perp(); }               ///< transverse component (to z)
  void SetPerp(double r) { fP.SetPerp(r); }          ///< sets the transverse component
  double Perp2(const Vector3 &v) const { return fP.Perp2(v); }  ///< squared component perpendicular to v
  double Perp(const Vector3 &v) const { return fP.Perp(v); }    ///< component perpendicular to v
  double Pt(const Vector3 &v) const { return Perp(v); }         ///< component perpendicular to v
  double Angle(const Vector3 &v) const { return fP.Angle(v); }  ///< angle between the three-vector and v

  double Mag2() const { return T() * T() - fP.Mag2(); }  ///< invariant mass squared t^2 - |p|^2
  double Mag() const {  ///< invariant mass (negative for spacelike vectors)
    double mm = Mag2();
    return mm < 0.0 ? -std::sqrt(-mm) : std::sqrt(mm);
  }
  double M2() const { return Mag2(); }  ///< invariant mass squared
  double M() const { return Mag(); }    ///< invariant mass (negative for spacelike vectors)
  double Mt2() const { return E() * E() - Z() * Z(); }  ///< transverse mass squared E^2 - pz^2
  double Mt() const {  ///< transverse mass (negative if Mt2() < 0)
    double mm = Mt2();
    return mm < 0.0 ? -std::sqrt(-mm) : std::sqrt(mm);
  }
  double Beta() const { return fP.Mag() / fE; }  ///< |p|/E
  double Gamma() const {  ///< 1/sqrt(1 - Beta()^2)
    double b = Beta();
    return 1.0 / std::sqrt(1 - b * b);
  }
  double Plus() const { return T() + Z(); }   ///< light-cone component t + z
  double Minus() const { return T() - Z(); }  ///< light-cone component t - z
  double Rapidity() const;                    ///< 0.5 ln((E + pz)/(E - pz))
  double PseudoRapidity() const { return fP.PseudoRapidity(); }  ///< pseudorapidity of the three-vector
  double Eta() const { return PseudoRapidity(); }                 ///< pseudorapidity of the three-vector
  double Dot(const LorentzVector &q) const { return T() * q.T() - Z() * q.Z() - Y() * q.Y() - X() * q.X(); }  ///< Minkowski product
  double operator*(const LorentzVector &q) const { return Dot(q); }  ///< Minkowski product

  Vector3 BoostVector() const { return Vector3(X() / T(), Y() / T(), Z() / T()); }  ///< velocity p/E
  void Boost(double bx, double by, double bz);                       ///< boosts with velocity (bx, by, bz)
  void Boost(const Vector3 &b) { Boost(b.X(), b.Y(), b.Z()); }       ///< boosts with velocity b
  void RotateX(double angle) { fP.RotateX(angle); }                  ///< rotates the three-vector around x
  void RotateY(double angle) { fP.RotateY(angle); }                  ///< rotates the three-vector around y
  void RotateZ(double angle) { fP.RotateZ(angle); }                  ///< rotates the three-vector around z
  void RotateUz(const Vector3 &newUz) { fP.RotateUz(newUz); }        ///< see Vector3::RotateUz
  void Rotate(double a, const Vector3 &axis) { fP.Rotate(a, axis); } ///< rotates the three-vector around axis
  LorentzVector &operator*=(const Rotation &m) { fP *= m; return *this; }     ///< rotates the three-vector
  LorentzVector &Transform(const Rotation &m) { fP.Transform(m); return *this; }  ///< rotates the three-vector
  LorentzVector &operator*=(const LorentzRotation &m);  ///< applies the Lorentz transformation m
  LorentzVector &Transform(const LorentzRotation &m);   ///< applies the Lorentz transformation m

  bool operator==(const LorentzVector &q) const { return Vect() == q.Vect() && T() == q.T(); }  ///< exact equality
  bool operator!=(const LorentzVector &q) const { return Vect() != q.Vect() || T() != q.T(); }  ///< exact inequality
  LorentzVector operator+(const LorentzVector &q) const { return LorentzVector(fP + q.Vect(), fE + q.T()); }  ///< sum
  LorentzVector operator-(const LorentzVector &q) const { return LorentzVector(fP - q.Vect(), fE - q.T()); }  ///< difference
  LorentzVector operator-() const { return LorentzVector(-X(), -Y(), -Z(), -T()); }  ///< opposite vector
  LorentzVector operator*(double a) const { return LorentzVector(a * X(), a * Y(), a * Z(), a * T()); }  ///< scaled vector
  LorentzVector &operator+=(const LorentzVector &q) { fP += q.Vect(); fE += q.T(); return *this; }  ///< adds q
  LorentzVector &operator-=(const LorentzVector &q) { fP -= q.Vect(); fE -= q.T(); return *this; }  ///< subtracts q
  LorentzVector &operator*=(double a) { fP *= a; fE *= a; return *this; }  ///< scales by a

  void Print(std::ostream &out) const;  ///< prints (x,y,z,t) and (P,eta,phi,E)
  void Print() const;                   ///< prints to std::cout

 private:
  [[noreturn]] static void bad_index(int i);  ///< throws std::out_of_range for component index i
  Vector3 fP;  ///< three-vector part
  double fE;   ///< time component
};

/*! scaled four-vector */
inline LorentzVector operator*(double a, const LorentzVector &p) { return LorentzVector(a * p.X(), a * p.Y(), a * p.Z(), a * p.T()); }

inline double LorentzVector::operator()(int i) const {
  switch (i) {
    case 0: return fP.X();
    case 1: return fP.Y();
    case 2: return fP.Z();
    case 3: return fE;
  }
  bad_index(i);
}

inline double &LorentzVector::operator()(int i) {
  switch (i) {
    case 0: return fP(0);
    case 1: return fP(1);
    case 2: return fP(2);
    case 3: return fE;
  }
  bad_index(i);
}

}  // namespace nucqcd

#endif
