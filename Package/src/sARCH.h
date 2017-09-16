#ifndef sARCH_H  // include guard
#define sARCH_H

#include <RcppArmadillo.h>
using namespace Rcpp;

template <typename distribution>
class sARCH {
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
  sARCH() {
    ineq_lb         = 1e-6;
    ineq_ub         = 0.99999999;
    label           = CharacterVector::create("alpha0", "alpha1");
    coeffs_mean     = NumericVector::create(0.1, 0.1);
    coeffs_sd       = NumericVector::create(1e4, 1e4);
    Sigma0          = NumericVector::create(1, 1);
    lower           = NumericVector::create(1e-6, 1e-6);
    upper           = NumericVector::create(100, 0.9999);
    nb_coeffs       = label.size();
    nb_coeffs_model = 2;
    name            = "sARCH_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label,
                   lower, upper);
  }

  // set all parameters (including those of the distribution).
  // this function should always be called first
  void loadparam(const NumericVector& theta) {
    alpha0 = theta[0], alpha1 = theta[1];
    int Ind = 2;
    fz.loadparam(theta, Ind);
  }

  void set_sd(const NumericVector& new_sd) { coeffs_sd = new_sd; }

  void set_mean(const NumericVector& new_mean) { coeffs_mean = new_mean; }

  // empty
  void prep_ineq_vol(){};

  // inequality constraint function
  double ineq_func() { return alpha1; }

  // check prior
  bool calc_r1() {
    return fz.calc_r1() && alpha0 >= lower[0] && alpha1 >= lower[1] &&
           (ineq_func() < ineq_ub);
  }

  // initialize volatility
  volatility set_vol(const double& y0) {
    volatility out;
    out.h = alpha0 / (1 - alpha1);
    out.lnh = log(out.h);
    return out;
  }

  // increment volatility
  void increment_vol(volatility& vol, const double& yim1) {
    vol.h = alpha0 + alpha1 * pow(yim1, 2);
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

#endif  // sARCH.h
