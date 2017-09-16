#ifndef SYMMETRIC_H  // include guard
#define SYMMETRIC_H

#include <RcppArmadillo.h>
#include "Utils.h"
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename underlying>
class Symmetric {
  underlying f1;  // symmetric distribution
 public:
  double Eabsz;    // := E[|z|]
  double EzIpos;   // := E[z * I(z>=0)]  = Prob(z >= 0) * E[z | z >= 0]
  double EzIneg;   // := E[z * I(z<0)]   = Prob(z < 0)  * E[z | z < 0]
  double Ez2Ineg;  // := E[z^2 * I(z<0)] = Prob(z < 0)  * E[z^2 | z < 0]

  // constructor function called by higher-level classes (e.g. GARCH).
  // arguments are passed by reference and modified according to the
  // distribution
  void constructor(std::string& name, int& nb_coeffs,
                   NumericVector& coeffs_mean, NumericVector& coeffs_sd,
                   NumericVector& Sigma0, CharacterVector& label,
                   NumericVector& lower, NumericVector& upper) {
    f1.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label,
                   lower, upper);
    name.append("sym");
  }

  // nested call (this function should always be called first)
  void loadparam(const NumericVector& theta, int& Ind) {
    f1.loadparam(theta, Ind);
    prep_moments1();
  }

  // check prior (nested call)
  bool calc_r1() { return f1.calc_r1(); }

  // setup for "calc_kernel" (nested call)
  void prep_kernel() { f1.prep_kernel(); }

  double calc_cdf(const double& x) { return f1.cdf(x); }

  // returns kernel of a single observation (must call "prep_kernel" first)
  double calc_kernel(const volatility& vol, const double& yi) {
    return f1.kernel(vol, yi);  // if in A (log density);  // if not
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

  void prep_moments1() { f1.set_M1(); }  // prep-function for moments of order 1
  void set_Eabsz() { Eabsz = f1.M1; }    // = E[|z|]
  void set_EzIpos() {
    EzIpos = 0.5 * f1.M1;
  }  // = E[z * I(z>=0)] = Prob(z >= 0) * E[z | z >= 0]
  void set_EzIneg() {
    EzIneg = -0.5 * f1.M1;
  }  // = E[z * I(z<0)] = Prob(z < 0) * E[z | z < 0]
  void set_Ez2Ineg() {
    Ez2Ineg = 0.5;
  }  // = E[z^2 * I(z<0)] = Prob(z < 0) * E[z^2 | z < 0]

  // returns a random vector of length "n"
  NumericVector rndgen(const int& n) {
    NumericVector out(n);
    NumericVector u = runif(n, 0, 1);
    for (int i = 0; i < n; i++) out[i] = f1.invsample(u[i]);
    return out;
  }

  double calc_invsample(const double& x) { return (f1.invsample(x)); }
};

#endif  // Symmetric.h