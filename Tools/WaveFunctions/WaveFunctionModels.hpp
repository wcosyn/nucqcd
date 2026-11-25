/*! \file WaveFunctionModels.hpp
 * \brief The wave function parametrizations that can be chosen by name, as enumerations
 *
 * Deuteron (Deuteron, OldDeuteron) and He3 (He3wf) wave functions are chosen by a name, which programs
 * read from their command line. The names are parsed in one place, here, into the enumerations; the
 * classes then switch on the enumeration (the compiler warns about a missing case). An unknown name is an
 * std::invalid_argument that lists the valid names, thrown where the name is parsed. The classes keep
 * their constructors and factories that take a name, and have overloads that take the enumeration, so
 * library code can name a parametrization without a string that could be misspelled.
 */
#ifndef WAVEFUNCTIONMODELS_HPP
#define WAVEFUNCTIONMODELS_HPP

#include <string>
#include <vector>

namespace nucqcd {

/*! \brief Deuteron wave function parametrizations (Deuteron::Wavefunction::CreateWavefunction) */
enum class DeuteronModel {
  CDBonn,     ///< "cdbonn": CD-Bonn
  CDBonn87,   ///< "cdbonn87"
  Paris,      ///< "paris": Paris
  AV18,       ///< "av18": Argonne v18
  AV18b,      ///< "av18b"
  Nijmegen3,  ///< "nijm3": Nijmegen III
  GrossIIb,   ///< "grossiib": Gross IIb
  GrossL4,    ///< "grossl4": Gross lambda = 0.4
  GrossWJC1,  ///< "grosswjc1": Gross WJC-1
  GrossWJC2,  ///< "grosswjc2": Gross WJC-2
  Test        ///< "test": test wave function
};

/*! \return the name of a deuteron model as used on command lines and in messages (lower case, see DeuteronModel) */
std::string deuteron_model_name(DeuteronModel model);
/*! \return all deuteron models, in the order of the enumeration */
const std::vector<DeuteronModel> &deuteron_models();
/*! \param name name of a deuteron wave function, case insensitive
 * \return the model
 * \throw std::invalid_argument for an unknown name; the message lists the known ones */
DeuteronModel parse_deuteron_model(const std::string &name);

/*! \brief He3 wave function parametrizations (He3wf) */
enum class He3Model {
  AV18,      ///< "AV18"
  AV18_UIX,  ///< "AV18-UIX": AV18 with the Urbana IX three-nucleon force
  AV18_TML,  ///< "AV18-TML": AV18 with the Tucson-Melbourne force
  CDBonn,    ///< "CDBonn"
  CDBonn_TML ///< "CDBonn-TML": CD-Bonn with the Tucson-Melbourne force
};

/*! \return the name of a He3 model (case as in He3Model) */
std::string he3_model_name(He3Model model);
/*! \return all He3 models, in the order of the enumeration */
const std::vector<He3Model> &he3_models();
/*! \param name name of a He3 wave function, case sensitive as before
 * \return the model
 * \throw std::invalid_argument for an unknown name; the message lists the known ones */
He3Model parse_he3_model(const std::string &name);

}  // namespace nucqcd

#endif
