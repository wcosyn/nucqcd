/*! \file Vector3.hpp
 * \brief Two- and three-vectors (nucqcd::Vector2, nucqcd::Vector3), replacing ROOT's TVector2 and
 * TVector3 (roadmap step 9, phase 4).
 * \details The interface is the subset of ROOT's that nucqcd uses, with the same names and
 * conventions, so that switching is a rename. The formulas and their order of operations are ROOT's
 * (ROOT 6.40, math/physics), and the functions that ROOT defines in its headers are inline here as
 * well, while the ones ROOT compiles into libPhysics are in Vector3.cpp, compiled without fused
 * multiply-adds (-ffp-contract=off) like ROOT's library. The results are therefore bitwise the same
 * as ROOT's (checked in tests/unit/unit_vectors.cpp). Conventions: Vector3::Phi() is in [-pi, pi]
 * (0 for a vector along z), Vector3::Theta() in [0, pi]; Vector2::Phi() is in [0, 2pi).
 * Differences with ROOT: no TObject base; a bad index throws std::out_of_range instead of printing
 * an error and returning 0.
 * \author Wim Cosyn (ROOT's TVector2/TVector3 by Pasha Murat, Peter Malzacher and others)
 */
#ifndef NUCQCD_VECTOR3_HPP
#define NUCQCD_VECTOR3_HPP

#include <cmath>
#include <iosfwd>

namespace nucqcd {

class Rotation;

/*! \brief A two-vector (x, y), with the interface of ROOT's TVector2 that nucqcd uses */
class Vector2 {
 public:
  Vector2() : fX(0.0), fY(0.0) {}                      ///< null vector
  Vector2(double x, double y) : fX(x), fY(y) {}        ///< vector (x, y)
  explicit Vector2(const double *a) : fX(a[0]), fY(a[1]) {}  ///< vector (a[0], a[1])

  double X() const { return fX; }   ///< x component
  double Y() const { return fY; }   ///< y component
  double Px() const { return fX; }  ///< x component
  double Py() const { return fY; }  ///< y component
  void SetX(double x) { fX = x; }   ///< sets the x component
  void SetY(double y) { fY = y; }   ///< sets the y component
  void Set(double x, double y) { fX = x; fY = y; }      ///< sets both components
  void Set(const Vector2 &v) { fX = v.fX; fY = v.fY; }  ///< copies v

  double Mod2() const { return fX * fX + fY * fY; }  ///< squared length
  double Mod() const;                                ///< length
  double Phi() const;                                ///< azimuth in [0, 2pi) (ROOT: pi + atan2(-y, -x))
  double DeltaPhi(const Vector2 &v) const { return Phi_mpi_pi(Phi() - v.Phi()); }  ///< Phi() - v.Phi() in [-pi, pi)
  void SetMagPhi(double mag, double phi);            ///< sets length |mag| and azimuth phi
  Vector2 Unit() const;                              ///< unit vector (null vector for a null vector)
  Vector2 Ort() const { return Unit(); }             ///< unit vector, as Unit()
  Vector2 Proj(const Vector2 &v) const;              ///< projection onto v
  Vector2 Norm(const Vector2 &v) const;              ///< component perpendicular to v
  Vector2 Rotate(double phi) const;                  ///< this vector rotated by phi

  static double Phi_0_2pi(double x);   ///< x shifted by multiples of 2pi into [0, 2pi) (NaN stays NaN)
  static double Phi_mpi_pi(double x);  ///< x shifted by multiples of 2pi into [-pi, pi) (NaN stays NaN)

  Vector2 &operator+=(const Vector2 &v) { fX += v.fX; fY += v.fY; return *this; }  ///< adds v
  Vector2 &operator-=(const Vector2 &v) { fX -= v.fX; fY -= v.fY; return *this; }  ///< subtracts v
  Vector2 &operator*=(double s) { fX *= s; fY *= s; return *this; }                 ///< scales by s
  Vector2 &operator/=(double s) { fX /= s; fY /= s; return *this; }                 ///< divides by s

 private:
  double fX;  ///< x component
  double fY;  ///< y component
};

inline Vector2 operator+(const Vector2 &a, const Vector2 &b) { return Vector2(a.X() + b.X(), a.Y() + b.Y()); }  ///< sum
inline Vector2 operator-(const Vector2 &a, const Vector2 &b) { return Vector2(a.X() - b.X(), a.Y() - b.Y()); }  ///< difference
inline Vector2 operator*(const Vector2 &v, double s) { return Vector2(v.X() * s, v.Y() * s); }  ///< scaled vector
inline Vector2 operator*(double s, const Vector2 &v) { return Vector2(v.X() * s, v.Y() * s); }  ///< scaled vector
inline Vector2 operator/(const Vector2 &v, double s) { return Vector2(v.X() / s, v.Y() / s); }  ///< v divided by s
inline double operator*(const Vector2 &a, const Vector2 &b) { return a.X() * b.X() + a.Y() * b.Y(); }  ///< scalar product
inline double operator^(const Vector2 &a, const Vector2 &b) { return a.X() * b.Y() - a.Y() * b.X(); }  ///< z component of the cross product

inline Vector2 Vector2::Proj(const Vector2 &v) const { return v * (((*this) * v) / v.Mod2()); }
inline Vector2 Vector2::Norm(const Vector2 &v) const { return *this - Proj(v); }

/*! \brief A three-vector (x, y, z), with the interface of ROOT's TVector3 that nucqcd uses */
class Vector3 {
 public:
  Vector3() : fX(0.0), fY(0.0), fZ(0.0) {}                       ///< null vector
  Vector3(double x, double y, double z) : fX(x), fY(y), fZ(z) {}  ///< vector (x, y, z)
  explicit Vector3(const double *a) : fX(a[0]), fY(a[1]), fZ(a[2]) {}  ///< vector (a[0], a[1], a[2])

  double operator()(int i) const;  ///< component i (0, 1, 2 = x, y, z); throws std::out_of_range otherwise
  double &operator()(int i);       ///< component i (0, 1, 2 = x, y, z); throws std::out_of_range otherwise
  double operator[](int i) const { return operator()(i); }  ///< component i, as operator()
  double &operator[](int i) { return operator()(i); }       ///< component i, as operator()

  double x() const { return fX; }   ///< x component
  double y() const { return fY; }   ///< y component
  double z() const { return fZ; }   ///< z component
  double X() const { return fX; }   ///< x component
  double Y() const { return fY; }   ///< y component
  double Z() const { return fZ; }   ///< z component
  double Px() const { return fX; }  ///< x component
  double Py() const { return fY; }  ///< y component
  double Pz() const { return fZ; }  ///< z component
  void SetX(double x) { fX = x; }   ///< sets the x component
  void SetY(double y) { fY = y; }   ///< sets the y component
  void SetZ(double z) { fZ = z; }   ///< sets the z component
  void SetXYZ(double x, double y, double z) { fX = x; fY = y; fZ = z; }  ///< sets all components
  void GetXYZ(double *a) const { a[0] = fX; a[1] = fY; a[2] = fZ; }      ///< copies the components into a[0..2]

  double Mag2() const { return fX * fX + fY * fY + fZ * fZ; }  ///< squared length
  double Mag() const { return std::sqrt(Mag2()); }           ///< length
  double Perp2() const { return fX * fX + fY * fY; }           ///< squared transverse component (to z)
  double Perp() const;                                        ///< transverse component (to z)
  double Pt() const { return Perp(); }                        ///< transverse component (to z)
  double Perp2(const Vector3 &p) const;                       ///< squared component perpendicular to p
  double Perp(const Vector3 &p) const;                        ///< component perpendicular to p
  double Pt(const Vector3 &p) const { return Perp(p); }       ///< component perpendicular to p
  double Phi() const;       ///< azimuth atan2(y, x) in [-pi, pi], 0 if x = y = 0
  double Theta() const;     ///< polar angle in [0, pi], 0 for the null vector
  double CosTheta() const;  ///< cos of the polar angle, 1 for the null vector
  double Angle(const Vector3 &q) const;  ///< angle with q in [0, pi], 0 if either is null
  double PseudoRapidity() const;         ///< -ln tan(theta/2) (+-1e11 along +-z, 0 for the null vector)
  double Eta() const { return PseudoRapidity(); }  ///< pseudorapidity
  Vector3 Unit() const;        ///< unit vector (the vector itself if null)
  Vector3 Orthogonal() const;  ///< a vector orthogonal to this one
  Vector2 XYvector() const { return Vector2(fX, fY); }  ///< (x, y)

  void SetMag(double mag);      ///< scales to length mag (a null vector stays null)
  void SetPerp(double perp);    ///< scales x and y to transverse component perp (unchanged if it is 0)
  void SetTheta(double theta);  ///< sets the polar angle, keeping length and azimuth
  void SetPhi(double phi);      ///< sets the azimuth, keeping z and the transverse component
  void SetMagThetaPhi(double mag, double theta, double phi);  ///< sets length |mag|, polar angle and azimuth

  double Dot(const Vector3 &p) const { return fX * p.fX + fY * p.fY + fZ * p.fZ; }  ///< scalar product
  Vector3 Cross(const Vector3 &p) const {  ///< cross product this x p
    return Vector3(fY * p.fZ - p.fY * fZ, fZ * p.fX - p.fZ * fX, fX * p.fY - p.fX * fY);
  }

  void RotateX(double angle);                       ///< rotates by angle around the x axis
  void RotateY(double angle);                       ///< rotates by angle around the y axis
  void RotateZ(double angle);                       ///< rotates by angle around the z axis
  void Rotate(double angle, const Vector3 &axis);   ///< rotates by angle around axis
  void RotateUz(const Vector3 &newUz);              ///< transforms from the frame with z axis along the unit vector newUz
  Vector3 &operator*=(const Rotation &m);           ///< applies the rotation m (this = m*this)
  Vector3 &Transform(const Rotation &m);            ///< applies the rotation m (this = m*this)

  bool operator==(const Vector3 &v) const { return v.fX == fX && v.fY == fY && v.fZ == fZ; }  ///< exact equality
  bool operator!=(const Vector3 &v) const { return v.fX != fX || v.fY != fY || v.fZ != fZ; }  ///< exact inequality
  Vector3 &operator+=(const Vector3 &p) { fX += p.fX; fY += p.fY; fZ += p.fZ; return *this; }  ///< adds p
  Vector3 &operator-=(const Vector3 &p) { fX -= p.fX; fY -= p.fY; fZ -= p.fZ; return *this; }  ///< subtracts p
  Vector3 &operator*=(double a) { fX *= a; fY *= a; fZ *= a; return *this; }                   ///< scales by a
  Vector3 operator-() const { return Vector3(-fX, -fY, -fZ); }                                  ///< opposite vector

  void Print(std::ostream &out) const;  ///< prints (x,y,z) and (rho,theta,phi) [degrees]
  void Print() const;                   ///< prints to std::cout

 private:
  [[noreturn]] static void bad_index(int i);  ///< throws std::out_of_range for component index i
  double fX;  ///< x component
  double fY;  ///< y component
  double fZ;  ///< z component
};

Vector3 operator+(const Vector3 &a, const Vector3 &b);  ///< sum
Vector3 operator-(const Vector3 &a, const Vector3 &b);  ///< difference
Vector3 operator*(const Vector3 &p, double a);          ///< scaled vector
Vector3 operator*(double a, const Vector3 &p);          ///< scaled vector
double operator*(const Vector3 &a, const Vector3 &b);   ///< scalar product

inline double Vector3::Perp2(const Vector3 &p) const {
  double tot = p.Mag2();
  double ss = Dot(p);
  double per = Mag2();
  if (tot > 0.0) per -= ss * ss / tot;
  if (per < 0) per = 0;
  return per;
}

inline double Vector3::CosTheta() const {
  double ptot = Mag();
  return ptot == 0.0 ? 1.0 : fZ / ptot;
}

inline Vector3 Vector3::Orthogonal() const {
  double xx = fX < 0.0 ? -fX : fX;
  double yy = fY < 0.0 ? -fY : fY;
  double zz = fZ < 0.0 ? -fZ : fZ;
  if (xx < yy) return xx < zz ? Vector3(0, fZ, -fY) : Vector3(fY, -fX, 0);
  return yy < zz ? Vector3(-fZ, 0, fX) : Vector3(fY, -fX, 0);
}

inline void Vector3::SetMag(double mag) {
  double factor = Mag();
  if (factor != 0) {
    factor = mag / factor;
    SetX(fX * factor);
    SetY(fY * factor);
    SetZ(fZ * factor);
  }
}

inline void Vector3::SetPerp(double perp) {
  double p = Perp();
  if (p != 0.0) {
    fX *= perp / p;
    fY *= perp / p;
  }
}

inline double Vector3::operator()(int i) const {
  switch (i) {
    case 0: return fX;
    case 1: return fY;
    case 2: return fZ;
  }
  bad_index(i);
}

inline double &Vector3::operator()(int i) {
  switch (i) {
    case 0: return fX;
    case 1: return fY;
    case 2: return fZ;
  }
  bad_index(i);
}

}  // namespace nucqcd

#endif
