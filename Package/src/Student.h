#ifndef STUDENT_H  // include guard
#define STUDENT_H

#include <RcppArmadillo.h>
#include "Utils.h"
using namespace Rcpp;

//---------------------- Student-t distribution ----------------------//
class Student {
  double nu;     // degrees of freedom
  double nu_lb;  // lower bound for "nu"
  double lncst;  // constant term in "kernel"
  double cst;    // constant term in "PDF"
  double P;  // factor to standardize R's definition of Student-t distribution
 public:
  double M1;  // E[|z|]

  // constructor
  Student() { nu_lb = 2.1; }

  // constructor function called by higher-level classes (e.g. Garch).
  // arguments are passed by reference and are modified to include "nu"
  void constructor(std::string& name, int& nb_coeffs,
                   NumericVector& coeffs_mean, NumericVector& coeffs_sd,
                   NumericVector& Sigma0, CharacterVector& label,
                   NumericVector& lower, NumericVector& upper) {
    name.append("std_");
    nb_coeffs++;
    label.push_back("nu");  // nu
    coeffs_mean.push_back(10),
    coeffs_sd.push_back(1e4);           // mean and standard deviation of prior distribution
    Sigma0.push_back(10);  // Sigma0
    lower.push_back(nu_lb), upper.push_back(1000);  // lower and upper bounds
  }

  // set "nu" (this function should always be called first)
  void loadparam(const NumericVector& theta, int& Ind) {
    nu = theta[Ind];
    P = sqrt(nu / (nu - 2));
    cst =
        P * exp(lgammal(0.5 * (nu + 1)) - lgammal(0.5 * nu)) / sqrt(nu * M_PI);
    Ind++;
  }

  // check prior
  bool calc_r1() { return nu > nu_lb; }

  // set M1 := E[|z|]
  void set_M1() {
    M1 = sqrt((nu - 2) / M_PI) *
         exp(lgammal(0.5 * (nu - 1)) - lgammal(0.5 * nu));
  }

  // set constant term of "kernel"
  void prep_kernel() {
    lncst = lgammal(0.5 * (nu + 1)) - lgammal(0.5 * nu) - 0.5 * log(M_PI) +
            0.5 * nu * log(nu - 2);
  }

  // returns loglikelihood of a single observation (must call "prep_kernel"
  // first)
  double kernel(const volatility& vol, const double& yi) {
    return lncst + 0.5 * nu * vol.lnh -
           0.5 * (nu + 1) * log(vol.h * (nu - 2) + pow(yi, 2));
  }

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
  double cdf(const double& x) { return R::pt(x * P, nu, 1, 0); }

  // applies inverse transform sampling on a Uniform(0,1) draw
  double invsample(const double& u) { return R::qt(u, nu, 1, 0) / P; }
};

#endif  // Student.h