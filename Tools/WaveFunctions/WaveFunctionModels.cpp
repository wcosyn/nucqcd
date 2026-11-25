#include <Tools/WaveFunctions/WaveFunctionModels.hpp>

#include <Tools/Core/ErrorMessage.hpp>
#include <Tools/Core/StringUtil.hpp>

#include <stdexcept>

namespace nucqcd {

namespace {
struct DeuteronEntry {
  DeuteronModel model;
  const char *name;
};
// the order is that of the list in the error message ("known: ...")
const DeuteronEntry deuteron_table[] = {
    {DeuteronModel::CDBonn, "cdbonn"},       {DeuteronModel::CDBonn87, "cdbonn87"},   {DeuteronModel::Paris, "paris"},
    {DeuteronModel::AV18, "av18"},           {DeuteronModel::AV18b, "av18b"},         {DeuteronModel::Nijmegen3, "nijm3"},
    {DeuteronModel::GrossIIb, "grossiib"},   {DeuteronModel::GrossL4, "grossl4"},     {DeuteronModel::GrossWJC1, "grosswjc1"},
    {DeuteronModel::GrossWJC2, "grosswjc2"}, {DeuteronModel::Test, "test"}};

struct He3Entry {
  He3Model model;
  const char *name;
};
const He3Entry he3_table[] = {{He3Model::AV18, "AV18"},
                              {He3Model::AV18_UIX, "AV18-UIX"},
                              {He3Model::AV18_TML, "AV18-TML"},
                              {He3Model::CDBonn, "CDBonn"},
                              {He3Model::CDBonn_TML, "CDBonn-TML"}};
}  // namespace

std::string deuteron_model_name(DeuteronModel model) {
  for (const auto &e : deuteron_table)
    if (e.model == model) return e.name;
  throw std::invalid_argument("deuteron_model_name: invalid DeuteronModel");
}

const std::vector<DeuteronModel> &deuteron_models() {
  static const std::vector<DeuteronModel> all = [] {
    std::vector<DeuteronModel> v;
    for (const auto &e : deuteron_table) v.push_back(e.model);
    return v;
  }();
  return all;
}

DeuteronModel parse_deuteron_model(const std::string &name) {
  for (const auto &e : deuteron_table)
    if (nucqcd::iequals(name, e.name)) return e.model;
  std::string known;
  for (const auto &e : deuteron_table) known += (known.empty() ? "" : ", ") + std::string(e.name);
  throw std::invalid_argument(nucqcd::ErrorMessage() << "TDeuteron::Wavefunction::CreateWavefunction: unknown wave function \"" << name
                                                     << "\" (known: " << known << ")");
}

std::string he3_model_name(He3Model model) {
  for (const auto &e : he3_table)
    if (e.model == model) return e.name;
  throw std::invalid_argument("he3_model_name: invalid He3Model");
}

const std::vector<He3Model> &he3_models() {
  static const std::vector<He3Model> all = [] {
    std::vector<He3Model> v;
    for (const auto &e : he3_table) v.push_back(e.model);
    return v;
  }();
  return all;
}

He3Model parse_he3_model(const std::string &name) {
  for (const auto &e : he3_table)
    if (name == e.name) return e.model;
  std::string known;
  for (const auto &e : he3_table) known += (known.empty() ? "" : ", ") + std::string(e.name);
  throw std::invalid_argument(nucqcd::ErrorMessage() << "He3wf: invalid name for the he3 wave function: " << name
                                                     << " (known: " << known << ")");
}

}  // namespace nucqcd
