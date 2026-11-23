// The functions that ROOT 6.40 compiles into libPhysics (TVector2.cxx, TVector3.cxx), with ROOT's
// formulas. Compiled with -ffp-contract=off (Tools/Core/CMakeLists.txt): ROOT's library has no fused
// multiply-adds, so the results are bitwise the same.
#include <Tools/Core/Vector3.hpp>
#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/LorentzAlgebraDetail.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace nucqcd {

// ---- Vector2 ------------------------------------------------------------------------------------

double Vector2::Mod() const { return std::sqrt(fX * fX + fY * fY); }

Vector2 Vector2::Unit() const { return (Mod2()) ? *this / Mod() : Vector2(); }

double Vector2::Phi() const { return detail::pi + detail::root_atan2(-fY, -fX); }

double Vector2::Phi_0_2pi(double x) {
  if (std::isnan(x)) return x;
  while (x >= detail::two_pi) x -= detail::two_pi;
  while (x < 0.) x += detail::two_pi;
  return x;
}

double Vector2::Phi_mpi_pi(double x) {
  if (std::isnan(x)) return x;
  while (x >= detail::pi) x -= detail::two_pi;
  while (x < -detail::pi) x += detail::two_pi;
  return x;
}

Vector2 Vector2::Rotate(double phi) const {
  return Vector2(fX * std::cos(phi) - fY * std::sin(phi), fX * std::sin(phi) + fY * std::cos(phi));
}

void Vector2::SetMagPhi(double mag, double phi) {
  double amag = std::fabs(mag);
  fX = amag * std::cos(phi);
  fY = amag * std::sin(phi);
}

// ---- Vector3 ------------------------------------------------------------------------------------

void Vector3::bad_index(int i) {
  throw std::out_of_range(ErrorMessage() << "nucqcd::Vector3: bad component index " << i);
}

Vector3 &Vector3::operator*=(const Rotation &m) { return *this = m * (*this); }

Vector3 &Vector3::Transform(const Rotation &m) { return *this = m * (*this); }

double Vector3::Angle(const Vector3 &q) const {
  double ptot2 = Mag2() * q.Mag2();
  if (ptot2 <= 0) return 0.0;
  double arg = Dot(q) / std::sqrt(ptot2);
  if (arg > 1.0) arg = 1.0;
  if (arg < -1.0) arg = -1.0;
  return std::acos(arg);
}

double Vector3::Perp() const { return std::sqrt(Perp2()); }

double Vector3::Perp(const Vector3 &p) const { return std::sqrt(Perp2(p)); }

double Vector3::Phi() const { return fX == 0.0 && fY == 0.0 ? 0.0 : detail::root_atan2(fY, fX); }

double Vector3::Theta() const {
  return fX == 0.0 && fY == 0.0 && fZ == 0.0 ? 0.0 : detail::root_atan2(Perp(), fZ);
}

Vector3 Vector3::Unit() const {
  double tot2 = Mag2();
  double tot = (tot2 > 0) ? 1.0 / std::sqrt(tot2) : 1.0;
  return Vector3(fX * tot, fY * tot, fZ * tot);
}

void Vector3::RotateX(double angle) {
  double s = std::sin(angle);
  double c = std::cos(angle);
  double yy = fY;
  fY = c * yy - s * fZ;
  fZ = s * yy + c * fZ;
}

void Vector3::RotateY(double angle) {
  double s = std::sin(angle);
  double c = std::cos(angle);
  double zz = fZ;
  fZ = c * zz - s * fX;
  fX = s * zz + c * fX;
}

void Vector3::RotateZ(double angle) {
  double s = std::sin(angle);
  double c = std::cos(angle);
  double xx = fX;
  fX = c * xx - s * fY;
  fY = s * xx + c * fY;
}

void Vector3::Rotate(double angle, const Vector3 &axis) {
  Rotation trans;
  trans.Rotate(angle, axis);
  operator*=(trans);
}

void Vector3::RotateUz(const Vector3 &newUz) {
  double u1 = newUz.fX;
  double u2 = newUz.fY;
  double u3 = newUz.fZ;
  double up = u1 * u1 + u2 * u2;
  if (up) {
    up = std::sqrt(up);
    double px = fX, py = fY, pz = fZ;
    fX = (u1 * u3 * px - u2 * py + u1 * up * pz) / up;
    fY = (u2 * u3 * px + u1 * py + u2 * up * pz) / up;
    fZ = (u3 * u3 * px - px + u3 * up * pz) / up;
  } else if (u3 < 0.) {
    fX = -fX;
    fZ = -fZ;
  }
}

double Vector3::PseudoRapidity() const {
  double cosTheta = CosTheta();
  if (cosTheta * cosTheta < 1) return -0.5 * std::log((1.0 - cosTheta) / (1.0 + cosTheta));
  if (fZ == 0) return 0;
  return fZ > 0 ? 10e10 : -10e10;
}

void Vector3::SetTheta(double th) {
  double ma = Mag();
  double ph = Phi();
  SetX(ma * std::sin(th) * std::cos(ph));
  SetY(ma * std::sin(th) * std::sin(ph));
  SetZ(ma * std::cos(th));
}

void Vector3::SetPhi(double ph) {
  double xy = Perp();
  SetX(xy * std::cos(ph));
  SetY(xy * std::sin(ph));
}

void Vector3::SetMagThetaPhi(double mag, double theta, double phi) {
  double amag = std::fabs(mag);
  fX = amag * std::sin(theta) * std::cos(phi);
  fY = amag * std::sin(theta) * std::sin(phi);
  fZ = amag * std::cos(theta);
}

void Vector3::Print(std::ostream &out) const {
  out << "(x,y,z)=(" << fX << "," << fY << "," << fZ << ") (rho,theta,phi)=(" << Mag() << ","
      << Theta() * 180. / detail::pi << "," << Phi() * 180. / detail::pi << ")\n";
}

void Vector3::Print() const { Print(std::cout); }

Vector3 operator+(const Vector3 &a, const Vector3 &b) {
  return Vector3(a.X() + b.X(), a.Y() + b.Y(), a.Z() + b.Z());
}

Vector3 operator-(const Vector3 &a, const Vector3 &b) {
  return Vector3(a.X() - b.X(), a.Y() - b.Y(), a.Z() - b.Z());
}

Vector3 operator*(const Vector3 &p, double a) { return Vector3(a * p.X(), a * p.Y(), a * p.Z()); }

Vector3 operator*(double a, const Vector3 &p) { return Vector3(a * p.X(), a * p.Y(), a * p.Z()); }

double operator*(const Vector3 &a, const Vector3 &b) { return a.Dot(b); }

}  // namespace nucqcd
