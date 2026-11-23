// Tools/Core: nucqcd::Vector2, Vector3, LorentzVector, Rotation and LorentzRotation, the ROOT-free
// replacements of TVector2, TVector3, TLorentzVector, TRotation and TLorentzRotation (roadmap step 9,
// phase 4). Every method nucqcd uses is compared with ROOT on random and special inputs, and the
// results must be bitwise the same (NaN counts as equal to NaN; +0 and -0 are different). This test
// needs ROOT and is kept as long as ROOT is available.
#include <test_support.hpp>

#include <Tools/Core/LorentzVector.hpp>
#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/Vector3.hpp>

#include <TError.h>
#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TRotation.h>
#include <TVector2.h>
#include <TVector3.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using nucqcd::LorentzRotation;
using nucqcd::LorentzVector;
using nucqcd::Rotation;
using nucqcd::Vector2;
using nucqcd::Vector3;

const int kN = 20000;  // random cases per test

// ROOT warns about zero axes and zero vectors (our classes do the same thing silently)
const bool kQuietRoot = (gErrorIgnoreLevel = kError, true);

bool same(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) return true;
  return std::memcmp(&a, &b, sizeof a) == 0;
}

std::string hex(double a) {
  std::ostringstream s;
  s << std::hexfloat << a;
  return s.str();
}

// Collects mismatches (with the first few described) instead of thousands of CHECKs
struct Comparison {
  std::vector<std::string> failures;
  long checks = 0;
  void operator()(const char *what, double ours, double root) {
    ++checks;
    if (!same(ours, root) && failures.size() < 20)
      failures.push_back(std::string(what) + ": nucqcd " + hex(ours) + " ROOT " + hex(root));
    else if (!same(ours, root))
      failures.push_back("");
  }
  void operator()(const char *what, const Vector3 &o, const TVector3 &r) {
    (*this)(what, o.X(), r.X());
    (*this)(what, o.Y(), r.Y());
    (*this)(what, o.Z(), r.Z());
  }
  void operator()(const char *what, const Vector2 &o, const TVector2 &r) {
    (*this)(what, o.X(), r.X());
    (*this)(what, o.Y(), r.Y());
  }
  void operator()(const char *what, const LorentzVector &o, const TLorentzVector &r) {
    for (int i = 0; i < 4; ++i) (*this)(what, o(i), r(i));
  }
  void operator()(const char *what, const Rotation &o, const TRotation &r) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) (*this)(what, o(i, j), r(i, j));
  }
  void operator()(const char *what, const LorentzRotation &o, const TLorentzRotation &r) {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) (*this)(what, o(i, j), r(i, j));
  }
  void report() {
    INFO(checks << " comparisons, " << failures.size() << " mismatches");
    for (const auto &f : failures)
      if (!f.empty()) UNSCOPED_INFO(f);
    CHECK(failures.empty());
  }
};

// Random numbers: mostly ordinary values, with zeros (both signs), +-1, tiny and huge magnitudes
struct Random {
  std::mt19937_64 gen{20261123};
  double uniform(double a, double b) { return std::uniform_real_distribution<double>(a, b)(gen); }
  int pick(int n) { return std::uniform_int_distribution<int>(0, n - 1)(gen); }
  double number() {
    switch (pick(10)) {
      case 0: return 0.0;
      case 1: return -0.0;
      case 2: return pick(2) ? 1.0 : -1.0;
      case 3: return (pick(2) ? 1 : -1) * std::pow(10., uniform(-12, 12));
      default: return uniform(-10, 10);
    }
  }
  double angle() {
    switch (pick(8)) {
      case 0: return 0.0;
      case 1: return (pick(9) - 4) * M_PI / 2;  // multiples of pi/2, including +-2pi
      case 2: return uniform(-100, 100);
      default: return uniform(-2 * M_PI, 2 * M_PI);
    }
  }
  Vector3 v3() { return Vector3(number(), number(), number()); }
  Vector3 velocity() {  // mostly |b| < 1, sometimes 0, sometimes close to 1, rarely > 1 (NaN)
    Vector3 dir(uniform(-1, 1), uniform(-1, 1), uniform(-1, 1));
    double b;
    switch (pick(10)) {
      case 0: return Vector3();
      case 1: b = 1 - std::pow(10., uniform(-15, -2)); break;
      case 2: b = uniform(1, 1.5); break;
      default: b = uniform(0, 0.99);
    }
    double m = dir.Mag();
    return m > 0 ? Vector3(b * dir.X() / m, b * dir.Y() / m, b * dir.Z() / m) : Vector3(0, 0, b);
  }
};

TVector3 root(const Vector3 &v) { return TVector3(v.X(), v.Y(), v.Z()); }
TVector2 root(const Vector2 &v) { return TVector2(v.X(), v.Y()); }
TLorentzVector root(const LorentzVector &v) { return TLorentzVector(v.X(), v.Y(), v.Z(), v.T()); }

// A random rotation/Lorentz transformation built the same way for both libraries
void random_rotation(Random &rnd, Rotation &o, TRotation &r) {
  for (int k = 0, n = 1 + rnd.pick(4); k < n; ++k) {
    double a = rnd.angle();
    switch (rnd.pick(4)) {
      case 0: o.RotateX(a); r.RotateX(a); break;
      case 1: o.RotateY(a); r.RotateY(a); break;
      case 2: o.RotateZ(a); r.RotateZ(a); break;
      default: {
        Vector3 axis = rnd.v3();
        o.Rotate(a, axis);
        r.Rotate(a, root(axis));
      }
    }
  }
}

void random_lorentz(Random &rnd, LorentzRotation &o, TLorentzRotation &r) {
  for (int k = 0, n = 1 + rnd.pick(5); k < n; ++k) {
    double a = rnd.angle();
    switch (rnd.pick(6)) {
      case 0: o.RotateX(a); r.RotateX(a); break;
      case 1: o.RotateY(a); r.RotateY(a); break;
      case 2: o.RotateZ(a); r.RotateZ(a); break;
      case 3: {
        Vector3 axis = rnd.v3();
        o.Rotate(a, axis);
        r.Rotate(a, root(axis));
        break;
      }
      case 4: {
        Vector3 b = rnd.velocity();
        o.Boost(b);
        r.Boost(root(b));
        break;
      }
      default: {
        Vector3 b = rnd.velocity();
        o.Boost(b.X(), b.Y(), b.Z());
        r.Boost(b.X(), b.Y(), b.Z());
      }
    }
  }
}

}  // namespace

TEST_CASE("Vector3: bitwise the same as TVector3", "[vectors]") {
  Random rnd;
  Comparison cmp;
  for (int n = 0; n < kN; ++n) {
    Vector3 a = rnd.v3(), b = rnd.v3();
    TVector3 ra = root(a), rb = root(b);
    double s = rnd.number(), ang = rnd.angle();
    cmp("Mag2", a.Mag2(), ra.Mag2());
    cmp("Mag", a.Mag(), ra.Mag());
    cmp("Perp2", a.Perp2(), ra.Perp2());
    cmp("Perp", a.Perp(), ra.Perp());
    cmp("Pt", a.Pt(), ra.Pt());
    cmp("Perp2(p)", a.Perp2(b), ra.Perp2(rb));
    cmp("Perp(p)", a.Perp(b), ra.Perp(rb));
    cmp("Phi", a.Phi(), ra.Phi());
    cmp("Theta", a.Theta(), ra.Theta());
    cmp("CosTheta", a.CosTheta(), ra.CosTheta());
    cmp("Angle", a.Angle(b), ra.Angle(rb));
    cmp("PseudoRapidity", a.PseudoRapidity(), ra.PseudoRapidity());
    cmp("Eta", a.Eta(), ra.Eta());
    cmp("Dot", a.Dot(b), ra.Dot(rb));
    cmp("operator*(v,v)", a * b, ra * rb);
    cmp("Unit", a.Unit(), ra.Unit());
    cmp("Orthogonal", a.Orthogonal(), ra.Orthogonal());
    cmp("Cross", a.Cross(b), ra.Cross(rb));
    cmp("XYvector", a.XYvector(), ra.XYvector());
    cmp("XYvector().Phi", a.XYvector().Phi(), ra.XYvector().Phi());
    cmp("a+b", a + b, ra + rb);
    cmp("a-b", a - b, ra - rb);
    cmp("a*s", a * s, ra * s);
    cmp("s*a", s * a, s * ra);
    cmp("-a", -a, -ra);
    for (int i = 0; i < 3; ++i) cmp("operator()", a(i), ra(i));
    for (int i = 0; i < 3; ++i) cmp("operator[]", a[i], ra[i]);
    {
      Vector3 o = a; TVector3 r = ra;
      o += b; r += rb; cmp("+=", o, r);
      o -= a; r -= ra; cmp("-=", o, r);
      o *= s; r *= s; cmp("*=", o, r);
    }
    {
      Vector3 o = a; TVector3 r = ra;
      o.SetMag(s); r.SetMag(s); cmp("SetMag", o, r);
    }
    {
      Vector3 o = a; TVector3 r = ra;
      o.SetPerp(s); r.SetPerp(s); cmp("SetPerp", o, r);
    }
    {
      Vector3 o = a; TVector3 r = ra;
      o.SetTheta(ang); r.SetTheta(ang); cmp("SetTheta", o, r);
    }
    {
      Vector3 o = a; TVector3 r = ra;
      o.SetPhi(ang); r.SetPhi(ang); cmp("SetPhi", o, r);
    }
    {
      Vector3 o; TVector3 r;
      double th = rnd.angle();
      o.SetMagThetaPhi(s, th, ang); r.SetMagThetaPhi(s, th, ang); cmp("SetMagThetaPhi", o, r);
    }
    {
      Vector3 o = a; TVector3 r = ra;
      o.RotateX(ang); r.RotateX(ang); cmp("RotateX", o, r);
      o.RotateY(ang); r.RotateY(ang); cmp("RotateY", o, r);
      o.RotateZ(ang); r.RotateZ(ang); cmp("RotateZ", o, r);
      o.Rotate(ang, b); r.Rotate(ang, rb); cmp("Rotate(angle,axis)", o, r);
    }
    {
      Vector3 o = a; TVector3 r = ra;
      Vector3 u = b.Unit();
      o.RotateUz(u); r.RotateUz(root(u)); cmp("RotateUz", o, r);
    }
    {
      Rotation m; TRotation rm;
      random_rotation(rnd, m, rm);
      Vector3 o = a; TVector3 r = ra;
      o *= m; r *= rm; cmp("*=(Rotation)", o, r);
      o.Transform(m); r.Transform(rm); cmp("Transform(Rotation)", o, r);
    }
  }
  cmp.report();
}

TEST_CASE("Vector2: bitwise the same as TVector2", "[vectors]") {
  Random rnd;
  Comparison cmp;
  for (int n = 0; n < kN; ++n) {
    Vector2 a(rnd.number(), rnd.number()), b(rnd.number(), rnd.number());
    TVector2 ra = root(a), rb = root(b);
    double s = rnd.number(), ang = rnd.angle();
    cmp("Mod2", a.Mod2(), ra.Mod2());
    cmp("Mod", a.Mod(), ra.Mod());
    cmp("Phi", a.Phi(), ra.Phi());
    cmp("DeltaPhi", a.DeltaPhi(b), ra.DeltaPhi(rb));
    cmp("Unit", a.Unit(), ra.Unit());
    cmp("Ort", a.Ort(), ra.Ort());
    cmp("Proj", a.Proj(b), ra.Proj(rb));
    cmp("Norm", a.Norm(b), ra.Norm(rb));
    cmp("Rotate", a.Rotate(ang), ra.Rotate(ang));
    cmp("Phi_0_2pi", Vector2::Phi_0_2pi(ang), TVector2::Phi_0_2pi(ang));
    cmp("Phi_mpi_pi", Vector2::Phi_mpi_pi(ang), TVector2::Phi_mpi_pi(ang));
    cmp("a+b", a + b, ra + rb);
    cmp("a-b", a - b, ra - rb);
    cmp("a*s", a * s, ra * s);
    cmp("s*a", s * a, s * ra);
    cmp("a/s", a / s, ra / s);
    cmp("a*b", a * b, ra * rb);
    cmp("a^b", a ^ b, ra ^ rb);
    {
      Vector2 o; TVector2 r;
      o.SetMagPhi(s, ang); r.SetMagPhi(s, ang); cmp("SetMagPhi", o, r);
    }
  }
  cmp.report();
}

TEST_CASE("LorentzVector: bitwise the same as TLorentzVector", "[vectors]") {
  Random rnd;
  Comparison cmp;
  for (int n = 0; n < kN; ++n) {
    LorentzVector a(rnd.number(), rnd.number(), rnd.number(), rnd.number());
    LorentzVector b(rnd.number(), rnd.number(), rnd.number(), rnd.number());
    TLorentzVector ra = root(a), rb = root(b);
    double s = rnd.number(), ang = rnd.angle(), m = rnd.number();
    Vector3 v = rnd.v3();
    cmp("P", a.P(), ra.P());
    cmp("Vect", a.Vect(), ra.Vect());
    cmp("Phi", a.Phi(), ra.Phi());
    cmp("Theta", a.Theta(), ra.Theta());
    cmp("CosTheta", a.CosTheta(), ra.CosTheta());
    cmp("Rho", a.Rho(), ra.Rho());
    cmp("Perp", a.Perp(), ra.Perp());
    cmp("Perp2", a.Perp2(), ra.Perp2());
    cmp("Perp(v)", a.Perp(v), ra.Perp(root(v)));
    cmp("Angle", a.Angle(v), ra.Angle(root(v)));
    cmp("Mag2", a.Mag2(), ra.Mag2());
    cmp("Mag", a.Mag(), ra.Mag());
    cmp("M2", a.M2(), ra.M2());
    cmp("M", a.M(), ra.M());
    cmp("Mt2", a.Mt2(), ra.Mt2());
    cmp("Mt", a.Mt(), ra.Mt());
    cmp("Beta", a.Beta(), ra.Beta());
    cmp("Gamma", a.Gamma(), ra.Gamma());
    cmp("Plus", a.Plus(), ra.Plus());
    cmp("Minus", a.Minus(), ra.Minus());
    cmp("Rapidity", a.Rapidity(), ra.Rapidity());
    cmp("Eta", a.Eta(), ra.Eta());
    cmp("Dot", a.Dot(b), ra.Dot(rb));
    cmp("operator*(v,v)", a * b, ra * rb);
    cmp("BoostVector", a.BoostVector(), ra.BoostVector());
    cmp("a+b", a + b, ra + rb);
    cmp("a-b", a - b, ra - rb);
    cmp("-a", -a, -ra);
    cmp("a*s", a * s, ra * s);
    cmp("s*a", s * a, s * ra);
    for (int i = 0; i < 4; ++i) cmp("operator()", a(i), ra(i));
    {
      LorentzVector o = a; TLorentzVector r = ra;
      o += b; r += rb; cmp("+=", o, r);
      o -= a; r -= ra; cmp("-=", o, r);
      o *= s; r *= s; cmp("*=", o, r);
    }
    {
      LorentzVector o; TLorentzVector r;
      o.SetXYZM(v.X(), v.Y(), v.Z(), m); r.SetXYZM(v.X(), v.Y(), v.Z(), m); cmp("SetXYZM", o, r);
      o.SetVectM(v, -m); r.SetVectM(root(v), -m); cmp("SetVectM", o, r);
      o.SetPxPyPzE(s, m, ang, v.X()); r.SetPxPyPzE(s, m, ang, v.X()); cmp("SetPxPyPzE", o, r);
      o.SetVect(v); r.SetVect(root(v)); cmp("SetVect", o, r);
      o.SetTheta(ang); r.SetTheta(ang); cmp("SetTheta", o, r);
      o.SetPhi(-ang); r.SetPhi(-ang); cmp("SetPhi", o, r);
      o.SetRho(s); r.SetRho(s); cmp("SetRho", o, r);
      o.SetE(m); r.SetE(m); cmp("SetE", o, r);
    }
    {
      Vector3 bv = rnd.velocity();
      LorentzVector o = a; TLorentzVector r = ra;
      o.Boost(bv); r.Boost(root(bv)); cmp("Boost(b)", o, r);
      o.Boost(-bv.X(), bv.Y(), -bv.Z()); r.Boost(-bv.X(), bv.Y(), -bv.Z()); cmp("Boost(bx,by,bz)", o, r);
      o.Boost(-o.BoostVector()); r.Boost(-r.BoostVector()); cmp("Boost(-BoostVector)", o, r);
    }
    {
      LorentzVector o = a; TLorentzVector r = ra;
      o.RotateX(ang); r.RotateX(ang); cmp("RotateX", o, r);
      o.RotateY(ang); r.RotateY(ang); cmp("RotateY", o, r);
      o.RotateZ(ang); r.RotateZ(ang); cmp("RotateZ", o, r);
      o.Rotate(ang, v); r.Rotate(ang, root(v)); cmp("Rotate", o, r);
    }
    {
      Rotation rot; TRotation rrot;
      random_rotation(rnd, rot, rrot);
      LorentzVector o = a; TLorentzVector r = ra;
      o *= rot; r *= rrot; cmp("*=(Rotation)", o, r);
      o.Transform(rot); r.Transform(rrot); cmp("Transform(Rotation)", o, r);
      LorentzRotation lr; TLorentzRotation rlr;
      random_lorentz(rnd, lr, rlr);
      o *= lr; r *= rlr; cmp("*=(LorentzRotation)", o, r);
      o.Transform(lr); r.Transform(rlr); cmp("Transform(LorentzRotation)", o, r);
      cmp("LorentzRotation*v", lr * a, rlr * ra);
    }
  }
  cmp.report();
}

TEST_CASE("Rotation: bitwise the same as TRotation", "[vectors]") {
  Random rnd;
  Comparison cmp;
  for (int n = 0; n < kN; ++n) {
    Rotation a, b;
    TRotation ra, rb;
    random_rotation(rnd, a, ra);
    random_rotation(rnd, b, rb);
    Vector3 v = rnd.v3();
    cmp("random rotation", a, ra);
    cmp("a*b", a * b, ra * rb);
    cmp("a*v", a * v, ra * root(v));
    cmp("Inverse", a.Inverse(), ra.Inverse());
    cmp("PhiX", a.PhiX(), ra.PhiX());
    cmp("PhiY", a.PhiY(), ra.PhiY());
    cmp("PhiZ", a.PhiZ(), ra.PhiZ());
    cmp("ThetaX", a.ThetaX(), ra.ThetaX());
    cmp("ThetaY", a.ThetaY(), ra.ThetaY());
    cmp("ThetaZ", a.ThetaZ(), ra.ThetaZ());
    {
      Rotation o = a; TRotation r = ra;
      o *= b; r *= rb; cmp("*=", o, r);
      o.Transform(b); r.Transform(rb); cmp("Transform", o, r);
      o.Invert(); r.Invert(); cmp("Invert", o, r);
    }
    {
      double angle, rangle;
      Vector3 axis;
      TVector3 raxis;
      a.AngleAxis(angle, axis);
      ra.AngleAxis(rangle, raxis);
      cmp("AngleAxis angle", angle, rangle);
      cmp("AngleAxis axis", axis, raxis);
    }
    {
      // RotateAxes with the (orthonormal) columns of b, and with bad axes (unchanged)
      Vector3 nx(b.XX(), b.YX(), b.ZX()), ny(b.XY(), b.YY(), b.ZY()), nz(b.XZ(), b.YZ(), b.ZZ());
      Rotation o = a; TRotation r = ra;
      o.RotateAxes(nx, ny, nz); r.RotateAxes(root(nx), root(ny), root(nz)); cmp("RotateAxes", o, r);
      o.RotateAxes(nx, nx, nz); r.RotateAxes(root(nx), root(nx), root(nz)); cmp("RotateAxes (bad)", o, r);
    }
    CHECK(a.IsIdentity() == ra.IsIdentity());
    CHECK((a == b) == (ra == rb));
  }
  cmp.report();
  CHECK(Rotation().IsIdentity());
  CHECK_THROWS_AS(Rotation()(3, 0), std::out_of_range);
}

TEST_CASE("LorentzRotation: bitwise the same as TLorentzRotation", "[vectors]") {
  Random rnd;
  Comparison cmp;
  for (int n = 0; n < kN; ++n) {
    LorentzRotation a, b;
    TLorentzRotation ra, rb;
    random_lorentz(rnd, a, ra);
    random_lorentz(rnd, b, rb);
    Vector3 bv = rnd.velocity();
    LorentzVector v(rnd.number(), rnd.number(), rnd.number(), rnd.number());
    cmp("random transformation", a, ra);
    cmp("boost (bx,by,bz)", LorentzRotation(bv.X(), bv.Y(), bv.Z()), TLorentzRotation(bv.X(), bv.Y(), bv.Z()));
    cmp("boost (Vector3)", LorentzRotation(bv), TLorentzRotation(root(bv)));
    cmp("a*b", a * b, ra * rb);
    cmp("MatrixMultiplication", a.MatrixMultiplication(b), ra.MatrixMultiplication(rb));
    cmp("a*v", a * v, ra * root(v));
    cmp("VectorMultiplication", a.VectorMultiplication(v), ra.VectorMultiplication(root(v)));
    cmp("Inverse", a.Inverse(), ra.Inverse());
    {
      Rotation rot; TRotation rrot;
      random_rotation(rnd, rot, rrot);
      cmp("from Rotation", LorentzRotation(rot), TLorentzRotation(rrot));
      LorentzRotation o = a; TLorentzRotation r = ra;
      o.Transform(rot); r.Transform(rrot); cmp("Transform(Rotation)", o, r);
      o.Transform(b); r.Transform(rb); cmp("Transform(LorentzRotation)", o, r);
      o *= b; r *= rb; cmp("*=", o, r);
      o.Invert(); r.Invert(); cmp("Invert", o, r);
    }
    CHECK(a.IsIdentity() == ra.IsIdentity());
    CHECK((a == b) == (ra == rb));
  }
  cmp.report();
  CHECK(LorentzRotation().IsIdentity());
  CHECK_THROWS_AS(LorentzRotation()(0, 4), std::out_of_range);
}

TEST_CASE("Vectors: conventions", "[vectors]") {
  // Phi in [-pi, pi] for Vector3, [0, 2pi) for Vector2; Theta of the null vector is 0
  CHECK(Vector3(-1, 0, 0).Phi() == M_PI);
  CHECK(Vector3(-1, -0.0, 0).Phi() == -M_PI);
  CHECK(Vector3(0, -1, 0).XYvector().Phi() == 1.5 * M_PI);
  CHECK(Vector3().Theta() == 0.0);
  CHECK(Vector3().CosTheta() == 1.0);
  // active transformations; Transform(m) applies m after the existing transformation
  LorentzRotation r;
  r.RotateZ(M_PI / 2).Boost(0, 0, 0.6);
  LorentzVector p = r * LorentzVector(1, 0, 0, 2);
  CHECK_THAT(p.X(), Catch::Matchers::WithinAbs(0, 1e-15));
  CHECK_THAT(p.Y(), Catch::Matchers::WithinRel(1., 1e-15));
  CHECK_THAT(p.Z(), Catch::Matchers::WithinRel(0.6 * 1.25 * 2, 1e-15));
  CHECK_THAT(p.T(), Catch::Matchers::WithinRel(1.25 * 2, 1e-15));
  // the boost to the rest frame
  LorentzVector q(0.3, -0.4, 1.2, 2.5);
  q.Boost(-q.BoostVector());
  CHECK_THAT(q.Vect().Mag(), Catch::Matchers::WithinAbs(0, 1e-15));
  CHECK_THROWS_AS(Vector3()(3), std::out_of_range);
  CHECK_THROWS_AS(LorentzVector()(-1), std::out_of_range);
}
