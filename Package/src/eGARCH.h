#ifndef eGARCH_H  // include guard
#define eGARCH_H

#include <RcppArmadillo.h>
using namespace Rcpp;

template <typename distribution>
class eGARCH {
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
  eGARCH() {
    ineq_lb         = -0.99999999;
    ineq_ub         = 0.99999999;
    label           = CharacterVector::create("alpha0", "alpha1", "alpha2", "beta");
    coeffs_mean     = NumericVector::create(0.00, 0.2, -0.1, 0.8);
    coeffs_sd       = NumericVector::create(1e4, 1e4, 1e4, 1e4);
    Sigma0          = NumericVector::create(1, 1, 1, 1);
    lower           = NumericVector::create(-50, -5, -5, -0.9999);
    upper           = NumericVector::create(50, 5, 5, 0.9999);
    nb_coeffs       = label.size();
    nb_coeffs_model = 4;
    name            = "eGARCH_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label,
                   lower, upper);
  }

  // set all parameters (including those of the distribution).
  // this function should always be called first
  void loadparam(const NumericVector& theta) {
    alpha0 = theta[0], alpha1 = theta[1], alpha2 = theta[2], beta = theta[3];
    int Ind = 4;
    fz.loadparam(theta, Ind);
    fz.prep_moments1();
    fz.set_Eabsz();
  }

  void set_sd(const NumericVector& new_sd) { coeffs_sd = new_sd; }

  void set_mean(const NumericVector& new_mean) { coeffs_mean = new_mean; }

  // empty
  void prep_ineq_vol(){};

  // inequality constraint function
  double ineq_func() { return beta; }

  // computes r1
  bool calc_r1() {
    double ineq_theta = ineq_func();
    return fz.calc_r1() && (ineq_theta > ineq_lb && ineq_theta < ineq_ub);
  }

  // initialize volatility to its undonditional expected value
  volatility set_vol() {
    volatility out;
    out.lnh = alpha0 / (1 - beta);
    out.h = exp(out.lnh);
    return out;
  }

  // increment volatility
  void increment_vol(volatility& vol, const double& yim1) {
    double z = yim1 / sqrt(vol.h);
    vol.lnh =
        alpha0 + alpha1 * (fabs(z) - fz.Eabsz) + alpha2 * z + beta * vol.lnh;
    vol.h = exp(vol.lnh);
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

#endif  // eGARCH.h
