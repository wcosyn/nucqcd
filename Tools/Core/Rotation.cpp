// The functions that ROOT 6.40 compiles into libPhysics (TRotation.cxx, TLorentzRotation.cxx), with
// ROOT's formulas. Compiled with -ffp-contract=off, see Vector3.cpp.
#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/LorentzAlgebraDetail.hpp>

#include <cmath>
#include <stdexcept>

namespace nucqcd {

// ---- Rotation -----------------------------------------------------------------------------------

double Rotation::operator()(int i, int j) const {
  if (i == 0) {
    if (j == 0) return fxx;
    if (j == 1) return fxy;
    if (j == 2) return fxz;
  } else if (i == 1) {
    if (j == 0) return fyx;
    if (j == 1) return fyy;
    if (j == 2) return fyz;
  } else if (i == 2) {
    if (j == 0) return fzx;
    if (j == 1) return fzy;
    if (j == 2) return fzz;
  }
  throw std::out_of_range(ErrorMessage() << "nucqcd::Rotation: bad indices (" << i << "," << j << ")");
}

Rotation Rotation::operator*(const Rotation &b) const {
  return Rotation(fxx * b.fxx + fxy * b.fyx + fxz * b.fzx, fxx * b.fxy + fxy * b.fyy + fxz * b.fzy,
                  fxx * b.fxz + fxy * b.fyz + fxz * b.fzz, fyx * b.fxx + fyy * b.fyx + fyz * b.fzx,
                  fyx * b.fxy + fyy * b.fyy + fyz * b.fzy, fyx * b.fxz + fyy * b.fyz + fyz * b.fzz,
                  fzx * b.fxx + fzy * b.fyx + fzz * b.fzx, fzx * b.fxy + fzy * b.fyy + fzz * b.fzy,
                  fzx * b.fxz + fzy * b.fyz + fzz * b.fzz);
}

Rotation &Rotation::Rotate(double a, const Vector3 &axis) {
  if (std::fmod(a, 2 * detail::pi) != 0.) {
    double ll = axis.Mag();
    if (ll != 0.0) {
      double sa = std::sin(a), ca = std::cos(a);
      double dx = axis.X() / ll, dy = axis.Y() / ll, dz = axis.Z() / ll;
      Rotation m(ca + (1 - ca) * dx * dx, (1 - ca) * dx * dy - sa * dz, (1 - ca) * dx * dz + sa * dy,
                 (1 - ca) * dy * dx + sa * dz, ca + (1 - ca) * dy * dy, (1 - ca) * dy * dz - sa * dx,
                 (1 - ca) * dz * dx - sa * dy, (1 - ca) * dz * dy + sa * dx, ca + (1 - ca) * dz * dz);
      Transform(m);
    }
  }
  return *this;
}

Rotation &Rotation::RotateX(double a) {
  double c = std::cos(a);
  double s = std::sin(a);
  double x = fyx, y = fyy, z = fyz;
  fyx = c * x - s * fzx;
  fyy = c * y - s * fzy;
  fyz = c * z - s * fzz;
  fzx = s * x + c * fzx;
  fzy = s * y + c * fzy;
  fzz = s * z + c * fzz;
  return *this;
}

Rotation &Rotation::RotateY(double a) {
  double c = std::cos(a);
  double s = std::sin(a);
  double x = fzx, y = fzy, z = fzz;
  fzx = c * x - s * fxx;
  fzy = c * y - s * fxy;
  fzz = c * z - s * fxz;
  fxx = s * x + c * fxx;
  fxy = s * y + c * fxy;
  fxz = s * z + c * fxz;
  return *this;
}

Rotation &Rotation::RotateZ(double a) {
  double c = std::cos(a);
  double s = std::sin(a);
  double x = fxx, y = fxy, z = fxz;
  fxx = c * x - s * fyx;
  fxy = c * y - s * fyy;
  fxz = c * z - s * fyz;
  fyx = s * x + c * fyx;
  fyy = s * y + c * fyy;
  fyz = s * z + c * fyz;
  return *this;
}

Rotation &Rotation::RotateAxes(const Vector3 &newX, const Vector3 &newY, const Vector3 &newZ) {
  double del = 0.001;
  Vector3 w = newX.Cross(newY);
  if (std::fabs(newZ.X() - w.X()) > del || std::fabs(newZ.Y() - w.Y()) > del || std::fabs(newZ.Z() - w.Z()) > del ||
      std::fabs(newX.Mag2() - 1.) > del || std::fabs(newY.Mag2() - 1.) > del || std::fabs(newZ.Mag2() - 1.) > del ||
      std::fabs(newX.Dot(newY)) > del || std::fabs(newY.Dot(newZ)) > del || std::fabs(newZ.Dot(newX)) > del) {
    return *this;
  }
  return Transform(Rotation(newX.X(), newY.X(), newZ.X(), newX.Y(), newY.Y(), newZ.Y(), newX.Z(), newY.Z(), newZ.Z()));
}

void Rotation::AngleAxis(double &angle, Vector3 &axis) const {
  double cosa = 0.5 * (fxx + fyy + fzz - 1);
  double cosa1 = 1 - cosa;
  if (cosa1 <= 0) {
    angle = 0;
    axis = Vector3(0, 0, 1);
  } else {
    double x = 0, y = 0, z = 0;
    if (fxx > cosa) x = std::sqrt((fxx - cosa) / cosa1);
    if (fyy > cosa) y = std::sqrt((fyy - cosa) / cosa1);
    if (fzz > cosa) z = std::sqrt((fzz - cosa) / cosa1);
    if (fzy < fyz) x = -x;
    if (fxz < fzx) y = -y;
    if (fyx < fxy) z = -z;
    angle = std::acos(cosa);
    axis = Vector3(x, y, z);
  }
}

double Rotation::PhiX() const { return (fyx == 0.0 && fxx == 0.0) ? 0.0 : detail::root_atan2(fyx, fxx); }
double Rotation::PhiY() const { return (fyy == 0.0 && fxy == 0.0) ? 0.0 : detail::root_atan2(fyy, fxy); }
double Rotation::PhiZ() const { return (fyz == 0.0 && fxz == 0.0) ? 0.0 : detail::root_atan2(fyz, fxz); }
double Rotation::ThetaX() const { return std::acos(fzx); }
double Rotation::ThetaY() const { return std::acos(fzy); }
double Rotation::ThetaZ() const { return std::acos(fzz); }

// ---- LorentzRotation ----------------------------------------------------------------------------

double LorentzRotation::operator()(int i, int j) const {
  if (i == 0) {
    if (j == 0) return fxx;
    if (j == 1) return fxy;
    if (j == 2) return fxz;
    if (j == 3) return fxt;
  } else if (i == 1) {
    if (j == 0) return fyx;
    if (j == 1) return fyy;
    if (j == 2) return fyz;
    if (j == 3) return fyt;
  } else if (i == 2) {
    if (j == 0) return fzx;
    if (j == 1) return fzy;
    if (j == 2) return fzz;
    if (j == 3) return fzt;
  } else if (i == 3) {
    if (j == 0) return ftx;
    if (j == 1) return fty;
    if (j == 2) return ftz;
    if (j == 3) return ftt;
  }
  throw std::out_of_range(ErrorMessage() << "nucqcd::LorentzRotation: bad indices (" << i << "," << j << ")");
}

void LorentzRotation::SetBoost(double bx, double by, double bz) {
  double bp2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / std::sqrt(1.0 - bp2);
  double bgamma = gamma * gamma / (1.0 + gamma);
  fxx = 1.0 + bgamma * bx * bx;
  fyy = 1.0 + bgamma * by * by;
  fzz = 1.0 + bgamma * bz * bz;
  fxy = fyx = bgamma * bx * by;
  fxz = fzx = bgamma * bx * bz;
  fyz = fzy = bgamma * by * bz;
  fxt = ftx = gamma * bx;
  fyt = fty = gamma * by;
  fzt = ftz = gamma * bz;
  ftt = gamma;
}

LorentzRotation LorentzRotation::MatrixMultiplication(const LorentzRotation &b) const {
  return LorentzRotation(
      fxx * b.fxx + fxy * b.fyx + fxz * b.fzx + fxt * b.ftx, fxx * b.fxy + fxy * b.fyy + fxz * b.fzy + fxt * b.fty,
      fxx * b.fxz + fxy * b.fyz + fxz * b.fzz + fxt * b.ftz, fxx * b.fxt + fxy * b.fyt + fxz * b.fzt + fxt * b.ftt,
      fyx * b.fxx + fyy * b.fyx + fyz * b.fzx + fyt * b.ftx, fyx * b.fxy + fyy * b.fyy + fyz * b.fzy + fyt * b.fty,
      fyx * b.fxz + fyy * b.fyz + fyz * b.fzz + fyt * b.ftz, fyx * b.fxt + fyy * b.fyt + fyz * b.fzt + fyt * b.ftt,
      fzx * b.fxx + fzy * b.fyx + fzz * b.fzx + fzt * b.ftx, fzx * b.fxy + fzy * b.fyy + fzz * b.fzy + fzt * b.fty,
      fzx * b.fxz + fzy * b.fyz + fzz * b.fzz + fzt * b.ftz, fzx * b.fxt + fzy * b.fyt + fzz * b.fzt + fzt * b.ftt,
      ftx * b.fxx + fty * b.fyx + ftz * b.fzx + ftt * b.ftx, ftx * b.fxy + fty * b.fyy + ftz * b.fzy + ftt * b.fty,
      ftx * b.fxz + fty * b.fyz + ftz * b.fzz + ftt * b.ftz, ftx * b.fxt + fty * b.fyt + ftz * b.fzt + ftt * b.ftt);
}

}  // namespace nucqcd
