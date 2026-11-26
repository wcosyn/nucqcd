// Benchmark of the integrators used in nucqcd (roadmap step 6, "Remove the integrators rombergerN and
// cube_romb?"): rombergerN (nested 1D Romberg, now progs/attic/rombergerN.hpp), numint::cube_romb,
// numint::cube_adaptive (hcubature of cubature, see Integrator/CMakeLists.txt) and Cuba's Cuhre (numint::cuhre), at
// requested relative precisions 1e-3, 1e-5, 1e-7.
//  Part A: test integrands with known values in 3D: smooth, a step (sphere boundary inside the box), a
//          kink, a narrow peak, a non-rectangular domain mapped onto the box, and a Glauber-like
//          density integral over an impact parameter, azimuth and z with a cut. Reported: relative error
//          with respect to the exact (or GSL high-precision) value, number of evaluations, time.
//  Part B: the real integrands, through the integrator flag of the classes (1 = cube_romb,
//          2 = cube_adaptive): Glauber grid points (GlauberGridThick) for a proton along z and off axis,
//          and the plane-wave A(e,e'p) cross section (EMQECross). Reference: cube_adaptive at 1e-6 (1e-8
//          took hours for one grid). Every run has a time limit (60 s in part A, 600 s in part B).
// Usage: integrators [parts] with parts any of A, B, C (default AB; C: the evaluation cap of the A(e,e'p)
// cross section, see part_c). Cuba runs with CUBACORES=0, i.e. without worker processes
// (set CUBACORES to compare timings with more cores). Grids are computed in empty temporary directories (NUCQCD_GRID_DIR), never read
// from a cache.
#include <Integrator/numint/numint.hpp>
#include <Tools/Core/Utilfunctions.hpp>
#include <progs/attic/rombergerN.hpp>  // the archived routine, method 0 of the benchmark
#include <Tools/Core/constants.hpp>
#include <Tools/Core/RunLog.hpp>
#include <Tools/Kinematics/ElectronKinematics.hpp>
#include <Tools/Kinematics/Kinematics2to2.hpp>
#include <Glauber/FastParticle.hpp>
#include <Glauber/GlauberGridThick.hpp>
#include <Glauber/MeanFieldNucleusThick.hpp>
#include <Knockout/EMQECross.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_erf.h>

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>
using namespace nucqcd;

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

// ---- Part A: test integrands on the unit cube [0,1]^3 (some rescaled inside) -----------------------

struct Integrand {
  std::string name;
  std::function<double(double, double, double)> f;
  double exact;
};

long g_count = 0;  // number of evaluations
const Integrand *g_current = nullptr;

double eval(double x, double y, double z) {
  ++g_count;
  return g_current->f(x, y, z);
}

// Woods-Saxon-like density, for the Glauber-like case
double ws(double r) { return r < 8. ? 1. / (1. + std::exp((r - 4.) / 0.5)) : 0.; }
// Glauber-like: b in [0,8], phi in [0,2pi], z in [1,8] mapped from the unit cube; distance to the point
// (b0 = 1.5, phi0 = 0) in the transverse plane; cut at r = 8
double glauber_like(double u, double v, double w) {
  double b = 8. * u, phi = 2 * nucqcd::pi * v, z = 1. + 7. * w;
  double r = std::sqrt(b * b + 1.5 * 1.5 - 2 * b * 1.5 * std::cos(phi) + z * z);
  return 8. * 2 * nucqcd::pi * 7. * b * ws(r);
}

// GSL reference for the Glauber-like integral: nested adaptive Gauss-Kronrod (qag, 21 points, relative
// 1e-11) over [0,1]^3; the error handler is off (roundoff warnings at the cut are expected)
struct GslCtx {
  std::function<double(double, double, double)> *f;
  gsl_integration_workspace *w[3];
  double x, y;
};
double gsl_z(double z, void *p) { auto *c = static_cast<GslCtx *>(p); return (*c->f)(c->x, c->y, z); }
double gsl_y(double y, void *p) {
  auto *c = static_cast<GslCtx *>(p);
  c->y = y;
  gsl_function F{gsl_z, c};
  double r, e;
  gsl_integration_qag(&F, 0, 1, 0, 1e-11, 2000, GSL_INTEG_GAUSS21, c->w[2], &r, &e);
  return r;
}
double gsl_x(double x, void *p) {
  auto *c = static_cast<GslCtx *>(p);
  c->x = x;
  gsl_function F{gsl_y, c};
  double r, e;
  gsl_integration_qag(&F, 0, 1, 0, 1e-11, 2000, GSL_INTEG_GAUSS21, c->w[1], &r, &e);
  return r;
}
double gsl_nested(std::function<double(double, double, double)> f) {
  GslCtx c{&f, {gsl_integration_workspace_alloc(2000), gsl_integration_workspace_alloc(2000), gsl_integration_workspace_alloc(2000)}, 0, 0};
  gsl_function F{gsl_x, &c};
  double r, e;
  gsl_integration_qag(&F, 0, 1, 0, 1e-11, 2000, GSL_INTEG_GAUSS21, c.w[0], &r, &e);
  for (auto *w : c.w) gsl_integration_workspace_free(w);
  return r;
}

std::vector<Integrand> integrands() {
  std::vector<Integrand> v;
  const double g1 = std::sqrt(nucqcd::pi) / 2 * std::erf(1.);
  v.push_back({"smooth exp(-r^2)", [](double x, double y, double z) { return std::exp(-(x * x + y * y + z * z)); }, g1 * g1 * g1});
  v.push_back({"step: inside sphere", [](double x, double y, double z) { return x * x + y * y + z * z < 1. ? 1. : 0.; }, nucqcd::pi / 6});
  v.push_back({"kink: (1-r^2) inside sphere", [](double x, double y, double z) { double r2 = x * x + y * y + z * z; return r2 < 1. ? 1. - r2 : 0.; }, nucqcd::pi / 15});
  const double s = 0.02, c[3] = {0.3, 0.6, 0.45};
  double peak = 1;
  for (double ci : c) peak *= s * std::sqrt(2 * nucqcd::pi) * 0.5 * (std::erf((1 - ci) / (s * std::sqrt(2.))) + std::erf(ci / (s * std::sqrt(2.))));
  v.push_back({"narrow peak (sigma 0.02)", [=](double x, double y, double z) {
                 return std::exp(-((x - c[0]) * (x - c[0]) + (y - c[1]) * (y - c[1]) + (z - c[2]) * (z - c[2])) / (2 * s * s)); }, peak});
  v.push_back({"simplex mapped onto cube", [](double x, double u, double w) { return std::exp(x + x * u + x * u * w) * x * x * u; },
               std::pow(std::exp(1.) - 1., 3) / 6});
  v.push_back({"Glauber-like WS with cut", glauber_like, gsl_nested(glauber_like)});
  return v;
}

// rombergerN nested three times
struct RombState { double acc; int maxsteps; double est_y, est_z; double x, y; };
void romb_z(double z, double *res, va_list ap) { RombState *s = va_arg(ap, RombState *); res[0] = eval(s->x, s->y, z); }
void romb_y(double y, double *res, va_list ap) {
  RombState *s = va_arg(ap, RombState *); s->y = y;
  rombergerN(romb_z, 0., 1., 1, res, s->acc, 3, s->maxsteps, &s->est_z, s);
}
void romb_x(double x, double *res, va_list ap) {
  RombState *s = va_arg(ap, RombState *); s->x = x;
  rombergerN(romb_y, 0., 1., 1, res, s->acc, 3, s->maxsteps, &s->est_y, s);
}

void mdf_exec(const numint::array<double, 3> &x, void *, double &ret) { ret = eval(x[0], x[1], x[2]); }

struct Result { double value; long nevals; double time; int flag; };  // flag: integrator's return/fail code (0 = ok)

Result run(int method, double prec) {
  g_count = 0;
  auto t0 = Clock::now();
  double value = 0;
  int flag = 0;
  numint::mdfunction<double, 3> mdf;
  mdf.func = mdf_exec;
  mdf.param = nullptr;
  numint::array<double, 3> lo = {{0., 0., 0.}}, hi = {{1., 1., 1.}};
  unsigned neval = 0;
  switch (method) {
    case 0: {  // rombergerN, as in the code: min 3, max 8 Romberg steps per dimension
      RombState s{prec, 8, 0., 0., 0., 0.};
      double est_x = 0.;
      rombergerN(romb_x, 0., 1., 1, &value, prec, 3, s.maxsteps, &est_x, &s);
      break;
    }
    case 1: flag = numint::cube_romb(mdf, lo, hi, 1e-12, prec, value, neval); break;
    case 2: flag = numint::cube_adaptive(mdf, lo, hi, 1e-12, prec, 100, 20000000, value, neval, 0); break;
    case 3: {
      int nregions, nev, fail;
      double err, prob;
      numint::cuhre(mdf, lo, hi, 1, prec, 1e-12, 0, 0, 20000000, 0, nullptr, nregions, nev, fail, value, err, prob);
      flag = fail;
      g_count = nev;  // as counted by Cuba (its own count, in case it evaluates in worker processes)
      break;
    }
  }
  return {value, g_count, seconds_since(t0), flag};
}

const char *method_name[] = {"rombergerN", "cube_romb", "cube_adaptive", "Cuhre"};

// cube_romb has no limit on the number of evaluations (up to 12 halvings per dimension), so every run
// is done in a child process that is killed after a time limit
const int kTimeLimit = 60;  // seconds
bool run_limited(int method, double prec, Result &r) {
  int fd[2];
  if (pipe(fd) != 0) return false;
  pid_t pid = fork();
  if (pid == 0) {
    close(fd[0]);
    alarm(kTimeLimit);
    Result res = run(method, prec);
    if (write(fd[1], &res, sizeof res) != sizeof res) _exit(1);
    _exit(0);
  }
  close(fd[1]);
  ssize_t n = read(fd[0], &r, sizeof r);
  close(fd[0]);
  int status = 0;
  waitpid(pid, &status, 0);
  return n == sizeof r && WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

void part_a() {
  std::printf("# Part A: test integrands on [0,1]^3 (rel.err = |value/exact - 1|)\n");
  std::printf("# flag: return code of cube_romb/cube_adaptive or Cuba's fail (0 = converged; rombergerN reports nothing)\n");
  std::printf("# %-28s %-14s %8s %12s %12s %10s %5s\n", "integrand", "method", "req.prec", "rel.err", "evaluations", "time [s]", "flag");
  for (const Integrand &in : integrands()) {
    g_current = &in;
    for (double prec : {1e-3, 1e-5, 1e-7})
      for (int m = 0; m < 4; ++m) {
        Result r;
        if (run_limited(m, prec, r))
          std::printf("  %-28s %-14s %8.0e %12.3e %12ld %10.3f %5d\n", in.name.c_str(), method_name[m], prec,
                      std::abs(r.value / in.exact - 1), r.nevals, r.time, r.flag);
        else
          std::printf("  %-28s %-14s %8.0e %12s %12s %10s\n", in.name.c_str(), method_name[m], prec, "timeout", "-",
                      (">" + std::to_string(kTimeLimit)).c_str());
        std::fflush(stdout);
      }
  }
}

// ---- Part B: real integrands through the classes' integrator flag ---------------------------------

std::string fresh_grid_dir(const std::string &tag) {
  auto dir = std::filesystem::temp_directory_path() / ("nucqcd_integrator_benchmark_" + tag);
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  setenv("NUCQCD_GRID_DIR", dir.c_str(), 1);
  return dir.string();
}

// Runs f() in a child process, killed after `limit` seconds; f returns the values to compare. On success
// out holds the values and the last element is the time [s].
const int kTimeLimitB = 600;  // seconds
bool limited(const std::function<std::vector<double>()> &f, std::vector<double> &out) {
  int fd[2];
  if (pipe(fd) != 0) return false;
  pid_t pid = fork();
  if (pid == 0) {
    close(fd[0]);
    alarm(kTimeLimitB);
    auto t0 = Clock::now();
    std::vector<double> v = f();
    v.push_back(seconds_since(t0));
    size_t n = v.size();
    if (write(fd[1], &n, sizeof n) != sizeof n) _exit(1);
    if (write(fd[1], v.data(), n * sizeof(double)) != ssize_t(n * sizeof(double))) _exit(1);
    _exit(0);
  }
  close(fd[1]);
  size_t n = 0;
  bool ok = read(fd[0], &n, sizeof n) == sizeof n;
  if (ok) {
    out.resize(n);
    size_t got = 0;
    while (got < n * sizeof(double)) {
      ssize_t r = read(fd[0], reinterpret_cast<char *>(out.data()) + got, n * sizeof(double) - got);
      if (r <= 0) break;
      got += r;
    }
    ok = got == n * sizeof(double);
  }
  close(fd[0]);
  int status = 0;
  waitpid(pid, &status, 0);
  return ok && WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

const double kRefPrec = 1e-6;  // reference: cube_adaptive at this precision

void part_b() {
  const std::string share = nucqcd::share_dir;
  const double pts[][3] = {{0.5, -0.8, 0.4}, {2., 0., 3.}, {4., 0.9, 1.}, {1., 0.5, 5.5}};
  std::printf("\n# Part B1: Glauber grid (carbon, 12x6x3 grid, 1500 MeV proton): max rel. deviation of 4 grid points\n");
  std::printf("# (FSI factor with and without CT) from cube_adaptive at %.0e; runs over %d s are stopped\n", kRefPrec, kTimeLimitB);
  std::printf("# %-14s %-14s %8s %12s %10s\n", "proton", "method", "req.prec", "max.rel.dev", "time [s]");
  for (double theta : {0.0, 0.5}) {
    auto grid_values = [&](int integrator, double prec) {
      return [=, &share, &pts]() {
        fresh_grid_dir("B1_" + std::to_string(integrator) + "_" + std::to_string(prec) + "_" + std::to_string(theta));
        MeanFieldNucleusThick carbon(MeanFieldNucleus::C, share);
        FastParticle proton(0, 0, 1500., theta, 0., 8., 0., 1., 1., share);
        GlauberGridThick grid(12, 6, 3, &carbon, prec, integration_method(integrator), share);
        grid.addParticle(proton);
        grid.fillGrids();
        std::vector<double> v;
        for (auto &p : pts)
          for (auto c : {grid.getFsiGridFull_interp3(p[0], p[1], p[2]), grid.getFsiCtGridFull_interp3(p[0], p[1], p[2])}) {
            v.push_back(c.real());
            v.push_back(c.imag());
          }
        return v;
      };
    };
    const char *label = theta == 0 ? "along z" : "theta=0.5 rad";
    std::vector<double> ref, vals;
    if (!limited(grid_values(2, kRefPrec), ref)) {
      std::printf("  %-14s reference did not finish within %d s\n", label, kTimeLimitB);
      continue;
    }
    std::printf("  %-14s %-14s %8.0e %12s %10.2f   (reference)\n", label, method_name[2], kRefPrec, "-", ref.back());
    for (double prec : {1e-3, 1e-5})
      for (int integrator : {1, 2}) {
        if (!limited(grid_values(integrator, prec), vals)) {
          std::printf("  %-14s %-14s %8.0e %12s %10s\n", label, method_name[integrator], prec, "timeout", (">" + std::to_string(kTimeLimitB)).c_str());
        } else {
          double dev = 0;
          for (size_t i = 0; i + 1 < ref.size(); i += 2)
            dev = std::max(dev, std::hypot(vals[i] - ref[i], vals[i + 1] - ref[i + 1]) / std::hypot(ref[i], ref[i + 1]));
          std::printf("  %-14s %-14s %8.0e %12.3e %10.2f\n", label, method_name[integrator], prec, dev, vals.back());
        }
        std::fflush(stdout);
      }
  }

  std::printf("\n# Part B2: plane-wave A(e,e'p) cross section (EMQECross, carbon, Q2=1.8 GeV2, pm ~ 150 MeV/c): rel.\n");
  std::printf("# deviation from cube_adaptive at %.0e, per shell; runs over %d s are stopped\n", kRefPrec, kTimeLimitB);
  std::printf("# %-6s %-14s %8s %12s %10s\n", "shell", "method", "req.prec", "rel.dev", "time [s]");
  int nshells = MeanFieldNucleusThick(MeanFieldNucleus::C, share).getPLevels();
  for (int shell = 0; shell < nshells; ++shell) {
    auto cross = [&](int integrator, double prec) {
      return [=, &share]() {
        fresh_grid_dir("B2");
        MeanFieldNucleusThick carbon(MeanFieldNucleus::C, share);
        std::unique_ptr<ElectronKinematics> elec(ElectronKinematics::CreateWithBeamEnergy(3245.));
        Kinematics2to2 kin(carbon.getMassA(), carbon.getMassA_min_proton(), nucqcd::mass_p, "qsquared:wlab:pklab", 1.8E06, 1000., 150.);
        EMQECross obs(*elec, &carbon, prec, integration_method(integrator), share, false, 1.);
        return std::vector<double>{obs.getDiffCross(kin, 2, 0, 0, 0, 1, shell, 0., 2000000, true, false, 1., 1.)};
      };
    };
    std::vector<double> ref, v;
    if (!limited(cross(2, kRefPrec), ref)) {
      std::printf("  %-6d reference did not finish within %d s\n", shell, kTimeLimitB);
      continue;
    }
    std::printf("  %-6d %-14s %8.0e %12s %10.3f   (reference)\n", shell, method_name[2], kRefPrec, "-", ref.back());
    for (double prec : {1e-3, 1e-5})
      for (int integrator : {1, 2}) {
        if (!limited(cross(integrator, prec), v))
          std::printf("  %-6d %-14s %8.0e %12s %10s\n", shell, method_name[integrator], prec, "timeout", (">" + std::to_string(kTimeLimitB)).c_str());
        else
          std::printf("  %-6d %-14s %8.0e %12.3e %10.3f\n", shell, method_name[integrator], prec, std::abs(v[0] / ref[0] - 1), v.back());
        std::fflush(stdout);
      }
  }
}

// Part C: the accuracy of the A(e,e'p) cross sections is set by the evaluation cap maxEval that the caller
// passes (cube_adaptive stops at maxEval for the 24-component amplitude integral, whatever the requested
// precision): scan maxEval at requested precision 1e-5, plane wave and RMSGA (one Glauber grid, computed once).
void part_c() {
  const std::string share = nucqcd::share_dir;
  std::printf("\n# Part C: A(e,e'p) cross section (EMQECross, carbon, Q2=1.8 GeV2, pm ~ 150 MeV/c) with cube_adaptive at\n");
  std::printf("# requested precision 1e-5 versus the evaluation cap maxEval: rel. deviation from maxEval = 2e7\n");
  std::printf("# %-6s %-6s %10s %12s %10s\n", "shell", "FSI", "maxEval", "rel.dev", "time [s]");
  fresh_grid_dir("C");  // one grid directory: the RMSGA grid is computed once and then read from the cache
  MeanFieldNucleusThick carbon(MeanFieldNucleus::C, share);
  std::unique_ptr<ElectronKinematics> elec(ElectronKinematics::CreateWithBeamEnergy(3245.));
  Kinematics2to2 kin(carbon.getMassA(), carbon.getMassA_min_proton(), nucqcd::mass_p, "qsquared:wlab:pklab", 1.8E06, 1000., 150.);
  EMQECross obs(*elec, &carbon, 1e-5, IntegrationMethod::Adaptive, share, false, 1.);
  for (int shell = 0; shell < carbon.getPLevels(); ++shell)
    for (int pw : {1, 0}) {
      obs.getDiffCross(kin, 2, pw ? 0 : 1, 0, 0, pw, shell, 0., 20000, true, false, 1., 1.);  // fills the grid once
      double ref = obs.getDiffCross(kin, 2, pw ? 0 : 1, 0, 0, pw, shell, 0., 20000000, true, false, 1., 1.);
      for (int maxeval : {2000, 20000, 200000, 2000000}) {
        auto t0 = Clock::now();
        double v = obs.getDiffCross(kin, 2, pw ? 0 : 1, 0, 0, pw, shell, 0., maxeval, true, false, 1., 1.);
        std::printf("  %-6d %-6s %10d %12.3e %10.3f\n", shell, pw ? "no" : "RMSGA", maxeval, std::abs(v / ref - 1), seconds_since(t0));
        std::fflush(stdout);
      }
    }
}

}  // namespace

int main(int argc, char *argv[]) {
  nucqcd::RunLog runlog(argc, argv, {"parts (any of A, B, C)"});
  setenv("CUBACORES", "0", 0);  // Cuba evaluates in this process (no workers), for comparable timings
  gsl_set_error_handler_off();
  const std::string parts = argc > 1 ? argv[1] : "AB";  // C on request
  if (parts.find('A') != std::string::npos) part_a();
  if (parts.find('B') != std::string::npos) part_b();
  if (parts.find('C') != std::string::npos) part_c();
  return 0;
}
