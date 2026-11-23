// The functions that ROOT 6.40 compiles into libPhysics (TLorentzVector.cxx), with ROOT's formulas.
// Compiled with -ffp-contract=off, see Vector3.cpp.
#include <Tools/Core/LorentzVector.hpp>
#include <Tools/Core/Rotation.hpp>
#include <Tools/Core/ErrorMessage.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace nucqcd {

void LorentzVector::bad_index(int i) {
  throw std::out_of_range(ErrorMessage() << "nucqcd::LorentzVector: bad component index " << i);
}

void LorentzVector::Boost(double bx, double by, double bz) {
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / std::sqrt(1.0 - b2);
  double bp = bx * X() + by * Y() + bz * Z();
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;
  SetX(X() + gamma2 * bp * bx + gamma * bx * T());
  SetY(Y() + gamma2 * bp * by + gamma * by * T());
  SetZ(Z() + gamma2 * bp * bz + gamma * bz * T());
  SetT(gamma * (T() + bp));
}

double LorentzVector::Rapidity() const { return 0.5 * std::log((E() + Pz()) / (E() - Pz())); }

LorentzVector &LorentzVector::operator*=(const LorentzRotation &m) { return *this = m.VectorMultiplication(*this); }

LorentzVector &LorentzVector::Transform(const LorentzRotation &m) { return *this = m.VectorMultiplication(*this); }

void LorentzVector::Print(std::ostream &out) const {
  out << "(x,y,z,t)=(" << X() << "," << Y() << "," << Z() << "," << T() << ") (P,eta,phi,E)=(" << P() << ","
      << Eta() << "," << Phi() << "," << T() << ")\n";
}

void LorentzVector::Print() const { Print(std::cout); }

}  // namespace nucqcd
