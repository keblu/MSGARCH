#ifndef tGARCH_H  // include guard
#define tGARCH_H

#include <RcppArmadillo.h>
using namespace Rcpp;

template <typename distribution>
class tGARCH {
  distribution fz;  // distribution of innovations (e.g. Symmetric<Normal>)
  double alpha0, alpha1, alpha2, beta;  // coefficients
 public:
  std::string name;  // name of the model
  int nb_coeffs;     // total number of coefficients (including those of "fz")
  int nb_coeffs_model;
  CharacterVector label;  // labels of coefficients
  NumericVector coeffs_mean,
      coeffs_sd;  // means and standard deviations for prior distributions
  NumericVector Sigma0;        // diagonal of matrix Sigma0
  NumericVector lower, upper;  // lower and upper bounds for coefficients
  double ineq_lb, ineq_ub;  // lower and upper bounds for inequality constraint

  // constructor
  tGARCH() {
    ineq_lb         = 1e-6;
    ineq_ub         = 0.99999999;
    label           = CharacterVector::create("alpha0", "alpha1", "alpha2", "beta");
    coeffs_mean     = NumericVector::create(0.125, 0.05, 0.01, 0.8);
    coeffs_sd       = NumericVector::create(1e4, 1e4, 1e4, 1e4);
    Sigma0          = NumericVector::create(1, 1, 1, 1);
    lower           = NumericVector::create(1e-7, 1e-6, 1e-4, 0.0);
    upper           = NumericVector::create(100, 10, 10, 10);
    nb_coeffs       = label.size();
    nb_coeffs_model = 4;
    name            = "tGARCH_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label,
                   lower, upper);
  }

  // set all parameters (including those of the distribution).
  // this function should always be called first
  void loadparam(const NumericVector& theta) {
    alpha0 = theta[0], alpha1 = theta[1], alpha2 = theta[2], beta = theta[3];
    int Ind = 4;
    fz.loadparam(theta, Ind);
  }

  void set_sd(const NumericVector& new_sd) { coeffs_sd = new_sd; }

  void set_mean(const NumericVector& new_mean) { coeffs_mean = new_mean; }

  void prep_ineq_vol() {
    fz.set_EzIneg();
    fz.set_Ez2Ineg();
  };

  // inequality constraint function
  double ineq_func() {
    return pow(alpha1, 2) + pow(beta, 2) -
           2 * (alpha1 + alpha2) * beta * fz.EzIneg -
           (pow(alpha1, 2) - pow(alpha2, 2)) * fz.Ez2Ineg;
  }

  // computes r1
  bool calc_r1() {
    return fz.calc_r1() && alpha0 >= lower[0] && alpha1 >= lower[1] &&
           alpha2 >= lower[2] && beta >= lower[3] && (ineq_func() < ineq_ub);
  }

  // initialize volatility to its undonditional expected value
  volatility set_vol() {
    volatility out;
    out.fh = alpha0 / (1 + (alpha1 + alpha2) * fz.EzIneg - beta);
    out.h = pow(out.fh, 2);
    out.lnh = log(out.h);
    return out;
  }

  // increment volatility
  void increment_vol(volatility& vol, const double& yim1) {
    vol.fh = alpha0 + beta * vol.fh + ((yim1 >= 0) ? alpha1 : -alpha2) * yim1;
    vol.h = pow(vol.fh, 2);
    vol.lnh = log(vol.h);
  }

  // some nested-call functions
  void prep_kernel() { fz.prep_kernel(); }
  NumericVector rndgen(const int& n) { return fz.rndgen(n); }
  double calc_pdf(const double& x) { return fz.calc_pdf(x); }
  double calc_cdf(const double& x) { return fz.calc_cdf(x); }
  double calc_kernel(const volatility& vol, const double& yi) {
    return fz.calc_kernel(vol, yi);
  }
};

#endif  // tGARCH.h
