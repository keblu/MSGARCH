#ifndef sGARCH_H  // include guard
#define sGARCH_H

#include <RcppArmadillo.h>
using namespace Rcpp;

template <typename distribution>
class sGARCH {
  distribution fz;  // distribution of innovations (e.g. Symmetric<Normal>)
  double alpha0, alpha1, beta;  // coefficients
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
  sGARCH() {
    ineq_lb         = 1e-6;
    ineq_ub         = 0.99999999;
    label           = CharacterVector::create("alpha0", "alpha1", "beta");
    coeffs_mean     = NumericVector::create(0.1, 0.1, 0.8);
    coeffs_sd       = NumericVector::create(1e4, 1e4, 1e4);
    Sigma0          = NumericVector::create(1, 1, 1);
    lower           = NumericVector::create(1e-7, 1e-6, 0.0);
    upper           = NumericVector::create(100, 0.9999, 0.9999);
    nb_coeffs       = label.size();
    nb_coeffs_model = 3;
    name        = "sGARCH_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label,
                   lower, upper);
  }

  // set all parameters (including those of the distribution).
  // this function should always be called first
  void loadparam(const NumericVector& theta) {
    alpha0 = theta[0], alpha1 = theta[1], beta = theta[2];
    int Ind = 3;
    fz.loadparam(theta, Ind);
  }

  NumericVector get_sd() { return (coeffs_sd); }

  void set_sd(const NumericVector& new_sd) { coeffs_sd = new_sd; }

  void set_mean(const NumericVector& new_mean) { coeffs_mean = new_mean; }

  // empty
  void prep_ineq_vol(){};

  // inequality constraint function
  double ineq_func() { return alpha1 + beta; }

  // check prior
  bool calc_r1() {
    return fz.calc_r1() && alpha0 >= lower[0] && alpha1 >= lower[1] &&
           beta >= lower[2] && (ineq_func() < ineq_ub);
  }

  // initialize volatility
  volatility set_vol(const double& y0) {
    volatility out;
    out.h = alpha0 / (1 - alpha1 - beta);
    out.lnh = log(out.h);
    return out;
  }

  // increment volatility
  void increment_vol(volatility& vol, const double& yim1) {
    vol.h = alpha0 + alpha1 * pow(yim1, 2) + beta * vol.h;
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

#endif  // sGARCH.h
