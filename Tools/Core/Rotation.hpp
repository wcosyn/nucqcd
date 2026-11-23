/*! \file Rotation.hpp
 * \brief Rotations and Lorentz transformations (nucqcd::Rotation, nucqcd::LorentzRotation), replacing
 * ROOT's TRotation and TLorentzRotation (roadmap step 9, phase 4).
 * \details Same interface (the subset nucqcd uses), conventions and formulas as ROOT 6.40; see
 * Vector3.hpp for how the results are kept bitwise the same as ROOT's. The transformations are
 * active (they move vectors, R*v), and R.Transform(M) means R = M*R, i.e. M is applied after R:
 * R.RotateX(a).Boost(b) first rotates, then boosts. R.Inverse() is the transpose for a Rotation and
 * the inverse metric-transpose for a LorentzRotation.
 * Differences with ROOT: no TObject base; bad indices throw std::out_of_range; Rotate with a null
 * axis and RotateAxes with non-orthonormal axes leave the rotation unchanged, as in ROOT, but
 * without ROOT's warning.
 * \author Wim Cosyn (ROOT's TRotation/TLorentzRotation by Peter Malzacher and others)
 */
#ifndef NUCQCD_ROTATION_HPP
#define NUCQCD_ROTATION_HPP

#include <Tools/Core/Vector3.hpp>
#include <Tools/Core/LorentzVector.hpp>

namespace nucqcd {

/*! \brief A 3x3 rotation matrix, with the interface of ROOT's TRotation that nucqcd uses */
class Rotation {
 public:
  Rotation() : fxx(1.0), fxy(0.0), fxz(0.0), fyx(0.0), fyy(1.0), fyz(0.0), fzx(0.0), fzy(0.0), fzz(1.0) {}  ///< identity
  /*! matrix with rows (mxx, mxy, mxz), (myx, myy, myz), (mzx, mzy, mzz); not checked to be a rotation */
  Rotation(double mxx, double mxy, double mxz, double myx, double myy, double myz, double mzx, double mzy, double mzz)
      : fxx(mxx), fxy(mxy), fxz(mxz), fyx(myx), fyy(myy), fyz(myz), fzx(mzx), fzy(mzy), fzz(mzz) {}

  double XX() const { return fxx; }  ///< matrix element xx
  double XY() const { return fxy; }  ///< matrix element xy
  double XZ() const { return fxz; }  ///< matrix element xz
  double YX() const { return fyx; }  ///< matrix element yx
  double YY() const { return fyy; }  ///< matrix element yy
  double YZ() const { return fyz; }  ///< matrix element yz
  double ZX() const { return fzx; }  ///< matrix element zx
  double ZY() const { return fzy; }  ///< matrix element zy
  double ZZ() const { return fzz; }  ///< matrix element zz
  double operator()(int i, int j) const;  ///< matrix element (i, j), 0..2; throws std::out_of_range otherwise

  bool operator==(const Rotation &m) const;  ///< exact equality
  bool operator!=(const Rotation &m) const { return !(*this == m); }  ///< exact inequality
  bool IsIdentity() const;                   ///< exactly the identity?
  Rotation &SetToIdentity();                 ///< sets to the identity

  Vector3 operator*(const Vector3 &p) const {  ///< rotated vector
    return Vector3(fxx * p.X() + fxy * p.Y() + fxz * p.Z(), fyx * p.X() + fyy * p.Y() + fyz * p.Z(),
                   fzx * p.X() + fzy * p.Y() + fzz * p.Z());
  }
  Rotation operator*(const Rotation &b) const;  ///< matrix product this*b
  Rotation &operator*=(const Rotation &m) { return *this = operator*(m); }  ///< this = this*m
  Rotation &Transform(const Rotation &m) { return *this = m.operator*(*this); }  ///< this = m*this
  Rotation Inverse() const { return Rotation(fxx, fyx, fzx, fxy, fyy, fzy, fxz, fyz, fzz); }  ///< inverse (transpose)
  Rotation &Invert() { return *this = Inverse(); }  ///< inverts this rotation

  Rotation &RotateX(double a);                  ///< followed by a rotation by a around the x axis
  Rotation &RotateY(double a);                  ///< followed by a rotation by a around the y axis
  Rotation &RotateZ(double a);                  ///< followed by a rotation by a around the z axis
  Rotation &Rotate(double a, const Vector3 &axis);  ///< followed by a rotation by a around axis
  /*! followed by the rotation that takes the x, y, z axes to the orthonormal right-handed newX,
   * newY, newZ (unchanged if they are not orthonormal within 1e-3) */
  Rotation &RotateAxes(const Vector3 &newX, const Vector3 &newY, const Vector3 &newZ);
  void AngleAxis(double &angle, Vector3 &axis) const;  ///< rotation angle and axis
  double PhiX() const;    ///< azimuth of the rotated x axis
  double PhiY() const;    ///< azimuth of the rotated y axis
  double PhiZ() const;    ///< azimuth of the rotated z axis
  double ThetaX() const;  ///< polar angle of the rotated x axis
  double ThetaY() const;  ///< polar angle of the rotated y axis
  double ThetaZ() const;  ///< polar angle of the rotated z axis

 private:
  double fxx;  ///< matrix element xx
  double fxy;  ///< matrix element xy
  double fxz;  ///< matrix element xz
  double fyx;  ///< matrix element yx
  double fyy;  ///< matrix element yy
  double fyz;  ///< matrix element yz
  double fzx;  ///< matrix element zx
  double fzy;  ///< matrix element zy
  double fzz;  ///< matrix element zz
};

inline bool Rotation::operator==(const Rotation &m) const {
  return fxx == m.fxx && fxy == m.fxy && fxz == m.fxz && fyx == m.fyx && fyy == m.fyy && fyz == m.fyz &&
         fzx == m.fzx && fzy == m.fzy && fzz == m.fzz;
}

inline bool Rotation::IsIdentity() const {
  return fxx == 1.0 && fxy == 0.0 && fxz == 0.0 && fyx == 0.0 && fyy == 1.0 && fyz == 0.0 && fzx == 0.0 &&
         fzy == 0.0 && fzz == 1.0;
}

inline Rotation &Rotation::SetToIdentity() {
  fxx = fyy = fzz = 1.0;
  fxy = fxz = fyx = fyz = fzx = fzy = 0.0;
  return *this;
}

/*! \brief A 4x4 Lorentz transformation (boosts and rotations), with the interface of ROOT's
 * TLorentzRotation that nucqcd uses; components ordered (x, y, z, t) */
class LorentzRotation {
 public:
  /*! identity */
  LorentzRotation()
      : fxx(1.0), fxy(0.0), fxz(0.0), fxt(0.0), fyx(0.0), fyy(1.0), fyz(0.0), fyt(0.0),
        fzx(0.0), fzy(0.0), fzz(1.0), fzt(0.0), ftx(0.0), fty(0.0), ftz(0.0), ftt(1.0) {}
  /*! the rotation r (implicit, as in ROOT) */
  LorentzRotation(const Rotation &r)
      : fxx(r.XX()), fxy(r.XY()), fxz(r.XZ()), fxt(0.0), fyx(r.YX()), fyy(r.YY()), fyz(r.YZ()), fyt(0.0),
        fzx(r.ZX()), fzy(r.ZY()), fzz(r.ZZ()), fzt(0.0), ftx(0.0), fty(0.0), ftz(0.0), ftt(1.0) {}
  LorentzRotation(double bx, double by, double bz) { SetBoost(bx, by, bz); }  ///< boost with velocity (bx, by, bz)
  LorentzRotation(const Vector3 &b) { SetBoost(b.X(), b.Y(), b.Z()); }  ///< boost with velocity b (implicit, as in ROOT)
  /*! matrix with rows (rxx, rxy, rxz, rxt), ..., (rtx, rty, rtz, rtt); not checked */
  LorentzRotation(double rxx, double rxy, double rxz, double rxt, double ryx, double ryy, double ryz, double ryt,
                  double rzx, double rzy, double rzz, double rzt, double rtx, double rty, double rtz, double rtt)
      : fxx(rxx), fxy(rxy), fxz(rxz), fxt(rxt), fyx(ryx), fyy(ryy), fyz(ryz), fyt(ryt),
        fzx(rzx), fzy(rzy), fzz(rzz), fzt(rzt), ftx(rtx), fty(rty), ftz(rtz), ftt(rtt) {}

  double XX() const { return fxx; }  ///< matrix element xx
  double XY() const { return fxy; }  ///< matrix element xy
  double XZ() const { return fxz; }  ///< matrix element xz
  double XT() const { return fxt; }  ///< matrix element xt
  double YX() const { return fyx; }  ///< matrix element yx
  double YY() const { return fyy; }  ///< matrix element yy
  double YZ() const { return fyz; }  ///< matrix element yz
  double YT() const { return fyt; }  ///< matrix element yt
  double ZX() const { return fzx; }  ///< matrix element zx
  double ZY() const { return fzy; }  ///< matrix element zy
  double ZZ() const { return fzz; }  ///< matrix element zz
  double ZT() const { return fzt; }  ///< matrix element zt
  double TX() const { return ftx; }  ///< matrix element tx
  double TY() const { return fty; }  ///< matrix element ty
  double TZ() const { return ftz; }  ///< matrix element tz
  double TT() const { return ftt; }  ///< matrix element tt
  double operator()(int i, int j) const;  ///< matrix element (i, j), 0..3 = x, y, z, t; throws std::out_of_range otherwise

  bool operator==(const LorentzRotation &r) const;  ///< exact equality
  bool operator!=(const LorentzRotation &r) const { return !(*this == r); }  ///< exact inequality
  bool IsIdentity() const;                          ///< exactly the identity?

  LorentzVector VectorMultiplication(const LorentzVector &p) const {  ///< transformed four-vector
    return LorentzVector(fxx * p.X() + fxy * p.Y() + fxz * p.Z() + fxt * p.T(),
                         fyx * p.X() + fyy * p.Y() + fyz * p.Z() + fyt * p.T(),
                         fzx * p.X() + fzy * p.Y() + fzz * p.Z() + fzt * p.T(),
                         ftx * p.X() + fty * p.Y() + ftz * p.Z() + ftt * p.T());
  }
  LorentzRotation MatrixMultiplication(const LorentzRotation &b) const;  ///< matrix product this*b
  LorentzVector operator*(const LorentzVector &p) const { return VectorMultiplication(p); }  ///< transformed four-vector
  LorentzRotation operator*(const LorentzRotation &m) const { return MatrixMultiplication(m); }  ///< matrix product
  LorentzRotation &operator*=(const LorentzRotation &m) { return *this = MatrixMultiplication(m); }  ///< this = this*m
  LorentzRotation &Transform(const LorentzRotation &m) { return *this = m.MatrixMultiplication(*this); }  ///< this = m*this
  LorentzRotation &Transform(const Rotation &m) { return Transform(LorentzRotation(m)); }  ///< this = m*this
  LorentzRotation Inverse() const {  ///< inverse transformation
    return LorentzRotation(fxx, fyx, fzx, -ftx, fxy, fyy, fzy, -fty, fxz, fyz, fzz, -ftz, -fxt, -fyt, -fzt, ftt);
  }
  LorentzRotation &Invert() { return *this = Inverse(); }  ///< inverts this transformation

  LorentzRotation &Boost(double bx, double by, double bz) { return Transform(LorentzRotation(bx, by, bz)); }  ///< followed by a boost
  LorentzRotation &Boost(const Vector3 &b) { return Transform(LorentzRotation(b)); }  ///< followed by a boost with velocity b
  LorentzRotation &RotateX(double angle) { return Transform(Rotation().RotateX(angle)); }  ///< followed by a rotation around x
  LorentzRotation &RotateY(double angle) { return Transform(Rotation().RotateY(angle)); }  ///< followed by a rotation around y
  LorentzRotation &RotateZ(double angle) { return Transform(Rotation().RotateZ(angle)); }  ///< followed by a rotation around z
  LorentzRotation &Rotate(double angle, const Vector3 &axis) { return Transform(Rotation().Rotate(angle, axis)); }  ///< followed by a rotation around axis

 private:
  void SetBoost(double bx, double by, double bz);  ///< sets to the boost with velocity (bx, by, bz) (protected in ROOT)
  double fxx;  ///< matrix element xx
  double fxy;  ///< matrix element xy
  double fxz;  ///< matrix element xz
  double fxt;  ///< matrix element xt
  double fyx;  ///< matrix element yx
  double fyy;  ///< matrix element yy
  double fyz;  ///< matrix element yz
  double fyt;  ///< matrix element yt
  double fzx;  ///< matrix element zx
  double fzy;  ///< matrix element zy
  double fzz;  ///< matrix element zz
  double fzt;  ///< matrix element zt
  double ftx;  ///< matrix element tx
  double fty;  ///< matrix element ty
  double ftz;  ///< matrix element tz
  double ftt;  ///< matrix element tt
};

inline bool LorentzRotation::operator==(const LorentzRotation &r) const {
  return fxx == r.fxx && fxy == r.fxy && fxz == r.fxz && fxt == r.fxt && fyx == r.fyx && fyy == r.fyy &&
         fyz == r.fyz && fyt == r.fyt && fzx == r.fzx && fzy == r.fzy && fzz == r.fzz && fzt == r.fzt &&
         ftx == r.ftx && fty == r.fty && ftz == r.ftz && ftt == r.ftt;
}

inline bool LorentzRotation::IsIdentity() const {
  return fxx == 1.0 && fxy == 0.0 && fxz == 0.0 && fxt == 0.0 && fyx == 0.0 && fyy == 1.0 && fyz == 0.0 &&
         fyt == 0.0 && fzx == 0.0 && fzy == 0.0 && fzz == 1.0 && fzt == 0.0 && ftx == 0.0 && fty == 0.0 &&
         ftz == 0.0 && ftt == 1.0;
}

}  // namespace nucqcd

#endif
