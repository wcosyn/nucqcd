#include <Tools/StructureFunctions/NuclStructure.hpp>
#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/constants.hpp>

#include <stdexcept>

namespace nucqcd {


namespace {
double nucleon_mass(bool proton) { return proton ? nucqcd::mass_p : nucqcd::mass_n; }

/* the kinematics (Q^2, x, W^2) from two of them, with the on-shell relation W^2 = Q^2(1/x-1) + M^2 */
struct OnShellPoint { double Q2, x, Wsq; };
OnShellPoint on_shell(bool proton, double var1, double var2, int switchvar) {
  const double mass = nucleon_mass(proton);
  switch (switchvar) {
    case 0: return {var1, var2, var1 * (1. / var2 - 1.) + mass * mass};
    case 1: return {var1, var1 / (var2 - mass * mass + var1), var2};
    case 2: return {(var2 - mass * mass) / (1. / var1 - 1.), var1, var2};
    default: throw std::invalid_argument(nucqcd::ErrorMessage() << "non supported format in constructor");
  }
}

/* "CB" here stands for Christy-Bosted up to W^2 = 9 GeV^2 and SLAC above */
std::string parametrization(const std::string &name) {
  return !name.compare("CB") ? std::string("CB+SLAC@9") : name;
}
}  // namespace

NuclStructure::NuclStructure(bool pr, double var1, double var2, int switchvar, std::string nm)
  : NuclStructure(pr, on_shell(pr, var1, var2, switchvar).Q2, on_shell(pr, var1, var2, switchvar).x,
                  on_shell(pr, var1, var2, switchvar).Wsq, nm) {}

NuclStructure::NuclStructure(bool pr, double Q2in, double xin, double Wsqin, std::string nm)
  : name(nm), proton(pr), mass(nucleon_mass(pr)), x(xin), Q2(Q2in), Wsq(Wsqin),
    structure(parametrization(nm)) {}

void NuclStructure::getF(double &F1, double &F2){
  structure.getF_xQW(F1,F2,proton,x,Q2,Wsq);
}

double NuclStructure::getF1(){
  double F1,F2;
  getF(F1,F2);
  return F1;
}

double NuclStructure::getF2(){
  double F1,F2;
  getF(F1,F2);
  return F2;
}

}  // namespace nucqcd

