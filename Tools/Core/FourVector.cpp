#include <Tools/Core/FourVector.hpp>
#include <Tools/Core/LorentzVector.hpp>


//_____________________________________________________________________
FourVector<double> ToFourVector(const nucqcd::LorentzVector& v)
{
  return FourVector<double>(v.E(),v.X(),v.Y(),v.Z());
}

//_____________________________________________________________________
nucqcd::LorentzVector ToLorentzVector(const FourVector<double>& v)
{
  return nucqcd::LorentzVector(v[1],v[2],v[3],v[0]);
}
