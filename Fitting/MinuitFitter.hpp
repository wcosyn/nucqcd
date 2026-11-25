/*! \file MinuitFitter.hpp
 * \brief nucqcd::MinuitFitter: the part of ROOT's TFitter/TMinuit interface that the fit programs use,
 * on top of the standalone Minuit2 library
 *
 * The fit programs (`progs/F2_largex_fit`, `progs/deepsdis`) used ROOT's `TVirtualFitter` (TMinuit) and
 * its commands ("MIGRAD", "SIMPLEX", "MINOS", ...), and peeked into TMinuit's status strings. This class
 * keeps their method names and conventions (external parameter numbering, `limits 0,0 = none`, the legacy
 * FCN signature `fcn(npar, gin, f, par, iflag)`), so the programs changed little, but needs no ROOT: it uses
 * the standalone Minuit2 (https://github.com/GooFit/Minuit2, built from the ROOT sources).
 *
 * Differences from TMinuit that matter (Minuit2 is a rewrite, not the same algorithm; results agree
 * within the fit errors, not digit by digit):
 * - "IMPROVE" does nothing: Minuit2 has no equivalent (TMinuit's IMPROVE looks for a better minimum by
 *   random sampling around the current one).
 * - "SET RAN" and "SET LIM" are accepted and ignored (nothing to seed; limits are never removed).
 * - Converged() is Minuit2's FunctionMinimum::IsValid (edm below the goal, call limit not reached),
 *   the counterpart of TMinuit's status "CONVERGED ".
 * - The parameter errors are Minuit2's (Hesse-based, strategy dependent), as TMinuit's were MIGRAD's.
 */
#ifndef NUCQCD_MINUITFITTER_HPP
#define NUCQCD_MINUITFITTER_HPP

#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinosError.h>
#include <Minuit2/MnGlobalCorrelationCoeff.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameterState.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace nucqcd {

/*! \brief Minuit2 behind the TFitter-like interface described in the file comment */
class MinuitFitter {
public:
  /*! the legacy TMinuit FCN: npar (number of parameters), gin (gradient, unused), f (the function value to
   * set), par (all parameters, fixed ones included), iflag (always 4: a normal function evaluation) */
  typedef void (*Fcn)(int &npar, double *gin, double &f, double *par, int iflag);

  /*! \param ndim number of parameters (as TVirtualFitter::Fitter(0, ndim); the parameters are added with SetParameter) */
  explicit MinuitFitter(int ndim) : ndim_(ndim), parameters_(ndim), set_(ndim, false) {}

  /*! \param fcn the function to minimize (a chi-square: the error definition is 1) */
  void SetFCN(Fcn fcn) { fcn_.reset(new Function(fcn)); }

  /*! Defines parameter \p i (0 <= i < ndim), as TFitter::SetParameter
   * \param name name \param value start value \param step initial step size
   * \param lo lower limit \param hi upper limit; lo == hi == 0 means no limits */
  void SetParameter(int i, const char *name, double value, double step, double lo, double hi) {
    check(i);
    parameters_[i] = Def{name, value, step, lo, hi, false};
    set_[i] = true;
    built_ = false;
  }

  /*! Fixes parameter \p i at its current value (also at its current value in the running fit) */
  void FixParameter(int i) {
    check(i);
    parameters_[i].fixed = true;
    if (built_) state_.Fix(i);
  }

  int GetNumberTotalParameters() const { return ndim_; }
  int GetNumberFreeParameters() const {
    int n = 0;
    for (const auto &p : parameters_) n += p.fixed ? 0 : 1;
    return n;
  }
  /*! TFitter returned TMinuit's call limit; Minuit2 takes the limit with each command, this returns the
   * one of the last MIGRAD/SIMPLEX (0 before the first: Minuit2's own default) */
  int GetMaxIterations() const { return static_cast<int>(max_calls_); }

  /*! Executes a TMinuit-style command, as TFitter::ExecuteCommand: "SET PRINT" (level), "SET STR"
   * (strategy 0, 1 or 2), "MIGRAD" (max calls, tolerance), "SIMPLEX" (max calls, tolerance), "MINOS" (max
   * calls), "IMPROVE" (does nothing), "SHO COV" (prints the covariance matrix), "SET RAN" and "SET LIM"
   * (ignored). \param command the command (case insensitive, one or more spaces are one space)
   * \param args its arguments \param nargs the number of arguments in \p args
   * \return 0 on success \throw std::invalid_argument for an unknown command */
  int ExecuteCommand(const char *command, double *args, int nargs) {
    const std::string cmd = normalize(command);
    auto arg = [&](int k, double fallback) { return (args && k < nargs) ? args[k] : fallback; };
    if (cmd == "SET PRINT") {
      // TMinuit's -1..3 (0: minimal) to Minuit2's 0..3 (0: errors only)
      ROOT::Minuit2::MnPrint::SetGlobalLevel(std::min(3, std::max(0, static_cast<int>(arg(0, 0.)))));
    } else if (cmd == "SET STR") {
      strategy_ = static_cast<unsigned>(arg(0, 1.));
    } else if (cmd == "MIGRAD" || cmd == "SIMPLEX") {
      build();
      max_calls_ = static_cast<unsigned>(arg(0, 0.));
      const double tolerance = arg(1, 0.1);
      ROOT::Minuit2::MnStrategy strategy(strategy_);
      if (cmd == "MIGRAD") {
        ROOT::Minuit2::MnMigrad migrad(*fcn_, state_, strategy);
        minimum_ = std::make_unique<ROOT::Minuit2::FunctionMinimum>(migrad(max_calls_, tolerance));
      } else {
        ROOT::Minuit2::MnSimplex simplex(*fcn_, state_, strategy);
        minimum_ = std::make_unique<ROOT::Minuit2::FunctionMinimum>(simplex(max_calls_, tolerance));
      }
      state_ = minimum_->UserState();
      minos_.clear();
    } else if (cmd == "MINOS") {
      if (!minimum_) throw std::logic_error("MinuitFitter: MINOS before MIGRAD");
      ROOT::Minuit2::MnMinos minos(*fcn_, *minimum_, ROOT::Minuit2::MnStrategy(strategy_));
      minos_.assign(ndim_, std::pair<double, double>(0., 0.));
      for (int i = 0; i < ndim_; ++i) {
        if (parameters_[i].fixed) continue;
        const ROOT::Minuit2::MinosError e = minos.Minos(i, static_cast<unsigned>(arg(0, 0.)));
        minos_[i] = std::pair<double, double>(e.Upper(), e.Lower());  // TMinuit: eplus > 0, eminus < 0
      }
    } else if (cmd == "IMPROVE" || cmd == "SET RAN" || cmd == "SET LIM") {
      // no Minuit2 equivalent, see the file comment
    } else if (cmd == "SHO COV" || cmd == "SHOW COV") {
      std::cout << "Covariance matrix of the free parameters:\n";
      const int n = GetNumberFreeParameters();
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) std::cout << GetCovarianceMatrixElement(i, j) << " ";
        std::cout << "\n";
      }
    } else {
      throw std::invalid_argument("MinuitFitter::ExecuteCommand: unknown command \"" + std::string(command) + "\"");
    }
    return 0;
  }

  /*! value of parameter \p i (external number, as TFitter) */
  double GetParameter(int i) const {
    check(i);
    return built_ ? state_.Value(static_cast<unsigned>(i)) : parameters_[i].value;
  }
  /*! parabolic error of parameter \p i (the step size before a minimization); 0 for a fixed parameter, as
   * in TMinuit (Minuit2 reports the step size, which made the programs' "fix parameters at their limits"
   * loop endless: a fixed parameter at its limit was fixed again and again) */
  double GetParError(int i) const {
    check(i);
    if (parameters_[i].fixed) return 0.;
    return built_ ? state_.Error(static_cast<unsigned>(i)) : parameters_[i].step;
  }
  /*! covariance matrix element of two free parameters, numbered as the parameters (0 for fixed ones or
   * without a minimization) */
  double GetCovarianceMatrixElement(int i, int j) const {
    check(i);
    check(j);
    if (!built_ || !state_.HasCovariance() || parameters_[i].fixed || parameters_[j].fixed) return 0.;
    return state_.Covariance()(state_.IntOfExt(i), state_.IntOfExt(j));
  }
  /*! as TFitter::GetErrors: MINOS errors of parameter \p i (eplus > 0, eminus < 0; 0 if MINOS was not run or the
   * parameter is fixed), the parabolic error and the global correlation coefficient (0 without covariance) */
  int GetErrors(int i, double &eplus, double &eminus, double &eparab, double &globcc) const {
    check(i);
    eplus = eminus = globcc = 0.;
    eparab = GetParError(i);
    if (!minos_.empty()) {
      eplus = minos_[i].first;
      eminus = minos_[i].second;
    }
    if (built_ && state_.HasCovariance() && !parameters_[i].fixed) {
      const ROOT::Minuit2::MnGlobalCorrelationCoeff cc = state_.GlobalCC();
      if (cc.IsValid()) globcc = cc.GlobalCC()[state_.IntOfExt(i)];
    }
    return 0;
  }

  /*! the last MIGRAD or SIMPLEX converged (TMinuit: status "CONVERGED ") */
  bool Converged() const { return minimum_ && minimum_->IsValid(); }
  /*! the error matrix of the last MIGRAD or SIMPLEX is accurate (TMinuit: "ERROR MATRIX ACCURATE ") */
  bool CovarianceAccurate() const { return minimum_ && minimum_->HasAccurateCovar(); }
  /*! the error matrix had to be forced positive definite (TMinuit: "ERR MATRIX NOT POS-DEF") */
  bool CovarianceForcedPosDef() const { return minimum_ && minimum_->HasMadePosDefCovar(); }

private:
  struct Def {
    std::string name;
    double value = 0., step = 0., lo = 0., hi = 0.;
    bool fixed = false;
  };

  class Function : public ROOT::Minuit2::FCNBase {
  public:
    explicit Function(Fcn fcn) : fcn_(fcn) {}
    double operator()(const std::vector<double> &par) const override {
      std::vector<double> p(par);
      int npar = static_cast<int>(p.size());
      double gin[1] = {0.}, f = 0.;
      fcn_(npar, gin, f, p.data(), 4);
      return f;
    }
    double Up() const override { return 1.; }  // chi-square

  private:
    Fcn fcn_;
  };

  void check(int i) const {
    if (i < 0 || i >= ndim_) throw std::out_of_range("MinuitFitter: parameter number out of range");
  }

  static std::string normalize(const char *command) {
    std::string s;
    for (const char *c = command; *c; ++c) {
      const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(*c)));
      if (u == ' ' && (s.empty() || s.back() == ' ')) continue;
      s += u;
    }
    while (!s.empty() && s.back() == ' ') s.pop_back();
    return s;
  }

  // builds the Minuit2 parameter state from the definitions at the first MIGRAD/SIMPLEX
  void build() {
    if (!fcn_) throw std::logic_error("MinuitFitter: no FCN set");
    if (built_) return;
    for (int i = 0; i < ndim_; ++i)
      if (!set_[i]) throw std::logic_error("MinuitFitter: parameter " + std::to_string(i) + " was not defined");
    ROOT::Minuit2::MnUserParameters pars;
    for (const Def &d : parameters_) {
      if (d.lo == 0. && d.hi == 0.) pars.Add(d.name, d.value, d.step);
      else pars.Add(d.name, d.value, d.step, d.lo, d.hi);
    }
    state_ = ROOT::Minuit2::MnUserParameterState(pars);
    for (int i = 0; i < ndim_; ++i)
      if (parameters_[i].fixed) state_.Fix(i);
    built_ = true;
  }

  int ndim_;
  std::vector<Def> parameters_;
  std::vector<bool> set_;
  std::unique_ptr<Function> fcn_;
  ROOT::Minuit2::MnUserParameterState state_;
  std::unique_ptr<ROOT::Minuit2::FunctionMinimum> minimum_;
  std::vector<std::pair<double, double>> minos_;  // (eplus, eminus) per parameter after MINOS
  bool built_ = false;
  unsigned strategy_ = 1;
  unsigned max_calls_ = 0;
};

}  // namespace nucqcd

#endif
