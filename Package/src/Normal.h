#ifndef NORMAL_H  // include guard
#define NORMAL_H

#include <RcppArmadillo.h>
#include "Utils.h"
using namespace Rcpp;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Normal {
  double lncst;  // constant term for "kernel"
 public:
  double M1;
  // constructor
  Normal() { lncst = -0.5 * log(2 * M_PI); }

  // constructor function called by higher-level classes (e.g. Garch).
  // here, the argument "name" is passed by reference and modified
  void constructor(std::string& name, int&, NumericVector&, NumericVector&,
                   NumericVector&, CharacterVector&, NumericVector&,
                   NumericVector&) {
    name.append("norm_");
  }

  // empty
  void loadparam(const NumericVector&, int&){};

  // check prior (always TRUE)
  bool calc_r1() { return 1; }

  void set_M1() { M1 = sqrt(2 / M_PI); }

  // empty
  void prep_kernel(){};

  // returns loglikelihood of a single observation
  double kernel(const volatility& vol, const double& yi) {
    return lncst - 0.5 * pow(yi, 2) / vol.h - 0.5 * vol.lnh;
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
  double cdf(const double& x) { return R::pnorm(x, 0, 1, 1, 0); }

  // applies inverse transform sampling on a Uniform(0,1) draw
  double invsample(const double& u) { return R::qnorm(u, 0, 1, 1, 0); }
};

#endif  // Normal.h