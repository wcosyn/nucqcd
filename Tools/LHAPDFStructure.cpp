#include "LHAPDFStructure.hpp"
#include "ErrorMessage.hpp"
#include "NucleonStructure.hpp"
#include "constants.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>

#ifdef NUCQCD_HAVE_LHAPDF
#include <LHAPDF/LHAPDF.h>
#else
namespace LHAPDF {
class PDF {};  // never created, but std::unique_ptr<PDF> needs a complete type
}
#endif

namespace nucqcd {

std::string LHAPDFOptions::key() const {
  std::ostringstream s;
  s.precision(17);
  s << unpolarized << '|' << polarized << '|' << nf_unpol << '|' << nf_pol << '|' << R << '|' << target_mass << '|'
    << n_ww_x << '|' << n_ww_Q2;
  return s.str();
}

std::shared_ptr<const LHAPDFStructure> LHAPDFStructure::get(const LHAPDFOptions &options) {
  static std::mutex mutex;
  static std::map<std::string, std::shared_ptr<const LHAPDFStructure>> cache;
  std::lock_guard<std::mutex> lock(mutex);
  auto &entry = cache[options.key()];
  if (!entry) entry = std::make_shared<const LHAPDFStructure>(options);
  return entry;
}

namespace {
// e_q^2 for d, u, s, c, b, t (PDG ids 1..6)
constexpr double charge2[7] = {0., 1. / 9., 4. / 9., 1. / 9., 4. / 9., 1. / 9., 4. / 9.};

// bilinear interpolation in table[i*nq+j] on the nodes lx, lq, clamped to the edges
double bilinear(const std::vector<double> &lx, const std::vector<double> &lq, const std::vector<double> &table,
                double x, double q) {
  x = std::clamp(x, lx.front(), lx.back());
  q = std::clamp(q, lq.front(), lq.back());
  const std::size_t nq = lq.size();
  const std::size_t i = std::min<std::size_t>(std::upper_bound(lx.begin(), lx.end(), x) - lx.begin(), lx.size() - 1) - 1;
  const std::size_t j = std::min<std::size_t>(std::upper_bound(lq.begin(), lq.end(), q) - lq.begin(), nq - 1) - 1;
  const double tx = (x - lx[i]) / (lx[i + 1] - lx[i]);
  const double tq = (q - lq[j]) / (lq[j + 1] - lq[j]);
  return (1 - tx) * (1 - tq) * table[i * nq + j] + tx * (1 - tq) * table[(i + 1) * nq + j] +
         (1 - tx) * tq * table[i * nq + j + 1] + tx * tq * table[(i + 1) * nq + j + 1];
}

std::vector<double> log_nodes(double lo, double hi, int n) {
  std::vector<double> v(n);
  for (int i = 0; i < n; i++) v[i] = std::log(lo) + (std::log(hi) - std::log(lo)) * i / (n - 1);
  return v;
}
}  // namespace

#ifdef NUCQCD_HAVE_LHAPDF

bool LHAPDFStructure::available() { return true; }

bool LHAPDFStructure::installedSet(const std::string &set) {
  const std::string name = set.substr(0, set.find('/'));
  const auto &sets = LHAPDF::availablePDFSets();
  return std::find(sets.begin(), sets.end(), name) != sets.end();
}

LHAPDFStructure::LHAPDFStructure(const LHAPDFOptions &options) : opts(options) {
  if (opts.unpolarized.empty() && opts.polarized.empty())
    throw std::invalid_argument(ErrorMessage() << "LHAPDFStructure: no LHAPDF set given");
  if (opts.n_ww_x < 3 || opts.n_ww_Q2 < 2)
    throw std::invalid_argument(ErrorMessage() << "LHAPDFStructure: too few nodes for the WW table");
  LHAPDF::setVerbosity(0);
  for (const std::string *set : {&opts.unpolarized, &opts.polarized})
    if (!set->empty() && !installedSet(*set))
      throw std::invalid_argument(ErrorMessage() << "LHAPDFStructure: LHAPDF set " << *set
                                                 << " is not installed (lhapdf install <set>)");
  if (!opts.unpolarized.empty()) unpol.reset(LHAPDF::mkPDF(opts.unpolarized));
  if (!opts.polarized.empty()) {
    pol.reset(LHAPDF::mkPDF(opts.polarized));
    pol_xmin = pol->xMin();
    pol_q2min = pol->q2Min();
    pol_q2max = pol->q2Max();
  }
}

double LHAPDFStructure::partonSum(bool polarized, bool proton, double x, double Q2, int nf) const {
  const LHAPDF::PDF *pdf = polarized ? pol.get() : unpol.get();
  if (!pdf)
    throw std::logic_error(ErrorMessage() << "LHAPDFStructure: no " << (polarized ? "polarized" : "unpolarized")
                                          << " LHAPDF set given");
  if (!(x > 0. && x < 1.)) return 0.;
  std::vector<double> xf;  // flavours -6..6 at index id+6
  pdf->xfxQ2(x, Q2, xf);
  double sum = 0.;
  for (int f = 1; f <= nf; f++) {
    const int fl = (!proton && f <= 2) ? 3 - f : f;  // isospin: u <-> d for the neutron
    sum += charge2[f] * (xf[fl + 6] + xf[-fl + 6]);
  }
  return sum;
}

double LHAPDFStructure::xfxQ(bool polarized, int id, double x, double Q) const {
  const LHAPDF::PDF *pdf = polarized ? pol.get() : unpol.get();
  if (!pdf)
    throw std::logic_error(ErrorMessage() << "LHAPDFStructure: no " << (polarized ? "polarized" : "unpolarized")
                                          << " LHAPDF set given");
  return pdf->xfxQ(id, x, Q);
}

#else  // built without LHAPDF

bool LHAPDFStructure::available() { return false; }
bool LHAPDFStructure::installedSet(const std::string &) { return false; }

LHAPDFStructure::LHAPDFStructure(const LHAPDFOptions &options) : opts(options) {
  throw std::runtime_error(ErrorMessage() << "LHAPDFStructure: nucqcd was built without LHAPDF");
}

double LHAPDFStructure::partonSum(bool, bool, double, double, int) const {
  throw std::runtime_error(ErrorMessage() << "LHAPDFStructure: nucqcd was built without LHAPDF");
}

double LHAPDFStructure::xfxQ(bool, int, double, double) const {
  throw std::runtime_error(ErrorMessage() << "LHAPDFStructure: nucqcd was built without LHAPDF");
}

#endif

LHAPDFStructure::~LHAPDFStructure() = default;

void LHAPDFStructure::getF(double &F1, double &F2, bool proton, double x, double Q2) const {
  F2 = partonSum(false, proton, x, Q2 * 1.E-06, opts.nf_unpol);
  if (F2 == 0.) {
    F1 = 0.;
    return;
  }
  const double R = opts.R < 0. ? NucleonStructure::getr1998(x, Q2) : opts.R;
  const double mass = proton ? mass_p : mass_n;
  F1 = F2 / (2. * x * (1. + R)) * (opts.target_mass ? 1. + 4. * mass * mass * x * x / Q2 : 1.);
}

double LHAPDFStructure::getG1(bool proton, double x, double Q2) const {
  const double xg1 = partonSum(true, proton, x, Q2 * 1.E-06, opts.nf_pol);
  return xg1 == 0. ? 0. : 0.5 * xg1 / x;
}

void LHAPDFStructure::buildWWTable() const {
  // x nodes log-spaced up to 0.95, then x = 1 (where the integral vanishes)
  const int nx = opts.n_ww_x, nq = opts.n_ww_Q2;
  ww_lx = log_nodes(std::max(pol_xmin, 1e-5), 0.95, nx - 1);
  ww_lx.push_back(0.);
  ww_lq = log_nodes(std::max(pol_q2min, 1.), std::min(pol_q2max, 1e5), nq);
  ww_p.assign(nx * nq, 0.);
  ww_n.assign(nx * nq, 0.);
  // Gauss-Legendre in t = ln y on each interval between neighbouring x nodes, accumulated from x = 1 down
  constexpr int n_sub = 40;
  gsl_integration_glfixed_table *gl = gsl_integration_glfixed_table_alloc(n_sub);
  for (int j = 0; j < nq; j++) {
    const double Q2 = std::exp(ww_lq[j]);  // GeV^2
    double cum_p = 0., cum_n = 0.;
    for (int i = nx - 2; i >= 0; i--) {
      double seg_p = 0., seg_n = 0.;
      for (int k = 0; k < n_sub; k++) {
        double t, w;
        gsl_integration_glfixed_point(ww_lx[i], ww_lx[i + 1], k, &t, &w, gl);
        const double y = std::exp(t);
        seg_p += w * 0.5 * partonSum(true, true, y, Q2, opts.nf_pol) / y;
        seg_n += w * 0.5 * partonSum(true, false, y, Q2, opts.nf_pol) / y;
      }
      cum_p += seg_p;
      cum_n += seg_n;
      ww_p[i * nq + j] = cum_p;
      ww_n[i * nq + j] = cum_n;
    }
  }
  gsl_integration_glfixed_table_free(gl);
}

double LHAPDFStructure::getG1plusg2(bool proton, double x, double Q2) const {
  if (!pol) throw std::logic_error(ErrorMessage() << "LHAPDFStructure: no polarized LHAPDF set given");
  if (!(x > 0. && x < 1.)) return 0.;
  std::call_once(ww_once, [this] { buildWWTable(); });
  return bilinear(ww_lx, ww_lq, proton ? ww_p : ww_n, std::log(x), std::log(Q2 * 1.E-06));
}

double LHAPDFStructure::getG1plusg2_direct(bool proton, double x, double Q2) const {
  if (!pol) throw std::logic_error(ErrorMessage() << "LHAPDFStructure: no polarized LHAPDF set given");
  if (!(x > 0. && x < 1.)) return 0.;
  struct Params { const LHAPDFStructure *self; bool proton; double Q2; };
  Params params{this, proton, Q2};
  gsl_function F;
  F.function = [](double t, void *p) {  // integrand in t = ln y: g1(e^t)
    const auto *par = static_cast<const Params *>(p);
    return par->self->getG1(par->proton, std::exp(t), par->Q2);
  };
  F.params = &params;
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc(1000);
  double result, error;
  gsl_error_handler_t *old = gsl_set_error_handler_off();
  const int status = gsl_integration_qag(&F, std::log(x), 0., 1e-9, 1e-8, 1000, GSL_INTEG_GAUSS61, ws, &result, &error);
  gsl_set_error_handler(old);
  gsl_integration_workspace_free(ws);
  if (status && status != GSL_EROUND)
    throw std::runtime_error(ErrorMessage() << "LHAPDFStructure: WW integration failed: " << gsl_strerror(status));
  return result;
}

}  // namespace nucqcd
