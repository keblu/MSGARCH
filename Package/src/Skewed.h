#ifndef SKEWED_H  // include guard
#define SKEWED_H

#include <RcppArmadillo.h>
#include "Utils.h"
using namespace Rcpp;

template <typename underlying>
class Skewed {
  underlying f1;  // symmetric distribution
  double xi;      // degree of assymetry
  double xi_lb;   // lower bound for "xi"
  double xi2;     // := xi^2
  double num;     // := 1 / (xi + 1/xi)
  double mu_xi, sig_xi, cutoff;
  double pcut;
  double lncst;  // constant that accounts for asymmetry for kernel calculation
  double intgrl_1, intgrl_2;
  int Nsi;  // number of subintervals for composite Simpson rule

  // Composite Simpson rule
  double compositeSimpsons(const double& a, const double& b, const double& c,
                           const int& p) {
    double x = a, dx = (b - a) / (2 * Nsi), dx2 = 2 * dx;
    double I1, I2, I3 = pow(c - x, p) * f1.pdf(x);
    double out = 0;
    for (int i = 0; i < Nsi; i++) {
      I1 = I3;
      I2 = pow(c - x - dx, p) * f1.pdf(x + dx);
      I3 = pow(c - x - dx2, p) * f1.pdf(x + dx2);
      out += dx / 3 * (I1 + 4 * I2 + I3);
      x += dx2;
    }
    return out;
  }

 public:
  double Eabsz;    // := E[|z|]
  double EzIpos;   // := E[z * I(z>=0)]  = Prob(z >= 0) * E[z | z >= 0]
  double EzIneg;   // := E[z * I(z<0)]   = Prob(z < 0)  * E[z | z < 0]
  double Ez2Ineg;  // := E[z^2 * I(z<0)] = Prob(z < 0)  * E[z^2 | z < 0]

  // constructor
  Skewed() {
    xi_lb = 0.01;
    Nsi = 5;
  }

  // constructor function called by higher-level classes (e.g. Garch).
  // arguments are passed by reference and modified according to the
  // distribution
  void constructor(std::string& name, int& nb_coeffs,
                   NumericVector& coeffs_mean, NumericVector& coeffs_sd,
                   NumericVector& Sigma0, CharacterVector& label,
                   NumericVector& lower, NumericVector& upper) {
    f1.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label,
                   lower, upper);
    name.append("skew");
    nb_coeffs++;
    label.push_back("xi");
    coeffs_mean.push_back(1);
    coeffs_sd.push_back(10);
    Sigma0.push_back(1);
    lower.push_back(xi_lb);
    upper.push_back(100);
  }

  // Load parameters and setup (this function should always be called first)
  void loadparam(const NumericVector& theta, int& Ind) {
    f1.loadparam(theta, Ind);
    f1.set_M1();
    xi = theta[Ind];
    xi2 = pow(xi, 2);
    num = 1 / (xi + 1 / xi);
    mu_xi = f1.M1 * (xi - 1 / xi);
    sig_xi =
        sqrt((1 - pow(f1.M1, 2)) * (xi2 + 1 / xi2) + 2 * pow(f1.M1, 2) - 1);
    pcut = num / xi;
    cutoff = -mu_xi / sig_xi;
    prep_moments1();
    prep_moments2();
  }

  // check prior
  bool calc_r1() { return f1.calc_r1() && xi > xi_lb; }

  // setup for "calc_kernel"
  void prep_kernel() {
    f1.prep_kernel();
    lncst = log(2 * sig_xi * num);
  }

  // returns CDF evaluated at "x"
  double calc_cdf(const double& x) {
    double tmp = sig_xi * x + mu_xi;
    return ((x < cutoff) ? 2 / xi * num * f1.cdf(tmp * xi)
                         : 2 * num * (xi * f1.cdf(tmp / xi) + 1 / xi) - 1);
  }

  // returns kernel of a single observation (must call "prep_kernel" first)
  double calc_kernel(const volatility& vol, const double& yi) {
    double sig = sqrt(vol.h);
    double yi_xi =
        ((yi >= sig * cutoff) ? 1 / xi : xi) * (sig_xi * yi + sig * mu_xi);
    return lncst + f1.kernel(vol, yi_xi);
  }

  double calc_pdf(const double& x) {
    prep_kernel();
    volatility unit;
    unit.h = 1;
    unit.lnh = 0;
    double LLd = calc_kernel(unit, x);
    LLd = ((LLd < LND_MIN) ? LND_MIN : LLd);
    return (exp(LLd));
  }

  // setup for moments of order 1
  void prep_moments1() {
    intgrl_1 = ((xi >= 1) ? compositeSimpsons(0, mu_xi / xi, mu_xi / xi, 1)
                          : compositeSimpsons(xi * mu_xi, 0, xi * mu_xi, 1));
  }

  // setup for moments of order 2
  void prep_moments2() {
    intgrl_2 = ((xi >= 1) ? compositeSimpsons(0, mu_xi / xi, mu_xi / xi, 2)
                          : compositeSimpsons(xi * mu_xi, 0, xi * mu_xi, 2));
  }

  // set Eabsz := E[|z|]
  void set_Eabsz() {
    Eabsz = 2 / sig_xi * num *
            (f1.M1 + 2 * ((xi >= 1) ? xi2 : -1 / xi2) * intgrl_1);
  }

  // set EzIpos := E[z * I(z>=0)] = Prob(z >= 0) * E[z | z >= 0]
  void set_EzIpos() {
    EzIpos = 2 / sig_xi * num *
             (0.5 * f1.M1 + ((xi >= 1) ? xi2 : -1 / xi2) * intgrl_1);
  }

  // set EzIneg := E[z * I(z<0)] = Prob(z < 0) * E[z | z < 0]
  void set_EzIneg() {
    EzIneg = -2 / sig_xi * num *
             (0.5 * f1.M1 + ((xi >= 1) ? xi2 : -1 / xi2) * intgrl_1);
  }

  // set Ez2Ineg := E[z^2 * I(z<0)] = Prob(z < 0) * E[z^2 | z < 0]
  void set_Ez2Ineg() {
    double xi3 = xi2 * xi, xi4 = xi3 * xi, sig2_xi = pow(sig_xi, 2),
           M1_2 = pow(f1.M1, 2);
    Ez2Ineg = ((xi >= 1)
                   ? 2 / sig2_xi * num *
                         (0.5 / xi3 * (1 + M1_2 * (xi4 - 1)) + xi3 * intgrl_2)
                   : 2 / (xi3 * sig2_xi) * num *
                         (0.5 - 0.5 * M1_2 * (1 - xi4) - intgrl_2));
  }

  // returns a random vector of length "n"
  NumericVector rndgen(const int& n) {
    NumericVector out(n);
    NumericVector u = runif(n, 0, 1);
    for (int i = 0; i < n; i++)
      out[i] =
          ((u[i] < pcut)
               ? (f1.invsample(0.5 * u[i] * (xi2 + 1)) / xi - mu_xi) / sig_xi
               : (f1.invsample(0.5 * u[i] * (1 + 1 / xi2) - 0.5 / xi2 + 0.5) *
                      xi -
                  mu_xi) /
                     sig_xi);
    return out;
  }
  double calc_invsample(const double& x) {
    double out =
        ((x < pcut)
             ? (f1.invsample(0.5 * x * (xi2 + 1)) / xi - mu_xi) / sig_xi
             : (f1.invsample(0.5 * x * (1 + 1 / xi2) - 0.5 / xi2 + 0.5) * xi -
                mu_xi) /
                   sig_xi);
    return (out);
  }
};

#endif  // Skewed.h