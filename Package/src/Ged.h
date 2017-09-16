#ifndef GED_H  // include guard
#define GED_H

#include <RcppArmadillo.h>
#include "Utils.h"
using namespace Rcpp;

//---------------------- Generalized error distribution ----------------------//
class Ged {
  double nu;      // shape parameter
  double nu_lb;   // lower bound for "nu"
  double lncst;   // constant term in "kernel"
  double cst;     // constant term in "PDF"
  double lambda;  // lambda
 public:
  double M1;  // E[|z|]

  // constructor
  Ged() { nu_lb = 0.7; }

  // constructor function called by higher-level classes (e.g. Garch).
  // arguments are passed by reference and are modified to include "nu"
  void constructor(std::string& name, int& nb_coeffs,
                   NumericVector& coeffs_mean, NumericVector& coeffs_sd,
                   NumericVector& Sigma0, CharacterVector& label,
                   NumericVector& lower, NumericVector& upper) {
    name.append("ged_");
    nb_coeffs++;
    label.push_back("nu");  // nu
    coeffs_mean.push_back(2),
    coeffs_sd.push_back(1e4);            // mean and standard deviation of prior distribution
    Sigma0.push_back(10);  // Sigma0
    lower.push_back(nu_lb), upper.push_back(20);  // lower and upper bounds
  }

  // set "nu" (this function should always be called first)
  void loadparam(const NumericVector& theta, int& Ind) {
    nu = theta[Ind];
    lambda = sqrt(pow(2, -2 / nu) * exp(lgammal(1 / nu) - lgammal(3 / nu)));
    cst = nu / (lambda * pow(2, 1 + 1 / nu) * exp(lgammal(1 / nu)));
    Ind++;
  }

  // check prior
  bool calc_r1() { return nu > nu_lb; }

  // set M1 := E[|z|]
  void set_M1() {
    M1 =
        0.5 * lambda * pow(8, 1 / nu) * exp(lgammal(1 / nu + 0.5)) / sqrt(M_PI);
  }

  // set constant term of "kernel"
  void prep_kernel() { lncst = log(cst); }

  // returns loglikelihood of a single observation (must call "prep_kernel"
  // first)
  double kernel(const volatility& vol, const double& yi) {
    return lncst - 0.5 * vol.lnh -
           0.5 * pow(fabs(yi / (sqrt(vol.h) * lambda)), nu);
  }

  // returns PDF evaluated at "x"
  double pdf(const double& x) {
    prep_kernel();
    volatility unit;
    unit.h = 1;
    unit.lnh = 0;
    double LLd = kernel(unit, x);
    LLd = ((LLd < LND_MIN) ? LND_MIN : LLd);
    return (exp(LLd));
  }

  // returns CDF evaluated at "x"
  double cdf(const double& x) {
    return (
        (x < 0)
            ? 0.5 * (1 - R::pgamma(0.5 * pow(-x / lambda, nu), 1 / nu, 1, 1, 0))
            : 0.5 *
                  (1 + R::pgamma(0.5 * pow(x / lambda, nu), 1 / nu, 1, 1, 0)));
  }

  // applies inverse transform sampling on a Uniform(0,1) draw
  double invsample(const double& u) {
    return (
        (u < 0.5)
            ? -lambda * pow(2 * R::qgamma(1 - 2 * u, 1 / nu, 1, 1, 0), 1 / nu)
            : lambda * pow(2 * R::qgamma(2 * u - 1, 1 / nu, 1, 1, 0), 1 / nu));
  }
};

#endif  // Ged.h