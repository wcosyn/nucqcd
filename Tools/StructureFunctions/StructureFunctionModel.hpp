/*! \file StructureFunctionModel.hpp
 * \brief The named nucleon structure function parametrizations of NucleonStructure, as an enumeration
 *
 * Besides these names NucleonStructure takes a grammar ("LHAPDF:\<set\>", "\<low\>+\<high\>@\<W2\>"), which
 * cannot be an enumeration; the plain names are parsed here, in one place, when a NucleonStructure is
 * constructed, so a misspelled name fails there (std::invalid_argument, listing the valid names) and not
 * when a structure function is first asked for.
 */
#ifndef STRUCTUREFUNCTIONMODEL_HPP
#define STRUCTUREFUNCTIONMODEL_HPP

#include <Tools/Core/ErrorMessage.hpp>

#include <stdexcept>
#include <string>

namespace nucqcd {

/*! \brief Plain parametrizations of NucleonStructure (see its constructor for the references) */
enum class StructureFunctionModel {
  SLAC,     ///< "SLAC": Bodek's SLAC parametrization
  Alekhin,  ///< "Alekhin": leading twist parametrization by Alekhin
  CB,       ///< "CB": Christy & Bosted
  CTEQ,     ///< "CTEQ": F2 from the CTEQ pdfs
  HMRS,     ///< "HMRS": F2 from the HMRS pdfs
  MSTW      ///< "MSTW": the MSTW 2008 LO pdfs through LHAPDF (alias of "LHAPDF:MSTW2008lo68cl")
};

/*! \return the name of a model as NucleonStructure takes it (case as in StructureFunctionModel) */
inline std::string structure_function_model_name(StructureFunctionModel model) {
  switch (model) {
    case StructureFunctionModel::SLAC: return "SLAC";
    case StructureFunctionModel::Alekhin: return "Alekhin";
    case StructureFunctionModel::CB: return "CB";
    case StructureFunctionModel::CTEQ: return "CTEQ";
    case StructureFunctionModel::HMRS: return "HMRS";
    case StructureFunctionModel::MSTW: return "MSTW";
  }
  throw std::invalid_argument("structure_function_model_name: invalid StructureFunctionModel");
}

/*! \param name one of the names of StructureFunctionModel, case sensitive (as NucleonStructure always was)
 * \return the model
 * \throw std::invalid_argument for any other name; the message lists the valid ones */
inline StructureFunctionModel parse_structure_function_model(const std::string &name) {
  for (StructureFunctionModel m : {StructureFunctionModel::SLAC, StructureFunctionModel::Alekhin, StructureFunctionModel::CB,
                                   StructureFunctionModel::CTEQ, StructureFunctionModel::HMRS, StructureFunctionModel::MSTW})
    if (name == structure_function_model_name(m)) return m;
  throw std::invalid_argument(nucqcd::ErrorMessage() << "NucleonStructure: invalid name for structure function \"" << name
                                                     << "\" (known: SLAC, Alekhin, CB, CTEQ, HMRS, MSTW, "
                                                        "LHAPDF:<set>[:<polarized set>], <low>+<high>[@<W^2 in GeV^2>])");
}

}  // namespace nucqcd

#endif
