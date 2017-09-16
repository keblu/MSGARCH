#ifndef UTILS_H  // include guard
#define UTILS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

typedef std::pair<double, double> pair;
const double LND_MIN = log(DBL_MIN) + 1;
struct volatility {
  double h;    // variance
  double lnh;  // log(h)
  double fh;   // some arbitraty function of "h" (e.g. sqrt(h) as in the Tgarch
               // model)
};

struct prior {
  bool r1;    // TRUE if the coefficients respect the prior, FALSE if not
  double r2;  // loglikelihood of the coefficients
  double r3;
};

// signum function
inline double signum(const double& x) { return (0 < x) - (x < 0); }

// check if x is infinite OR a NaN
inline bool IsInfNan(const double& x) {
  return traits::is_infinite<REALSXP>(x) || ISNAN(x);
}

// Auxiliary function for adaptive Simpson's Rule:
// https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
template <typename C>
double adaptiveSimpsonsAux(double (C::*f)(const double&), C* pf, double a,
                           double b, double epsilon, double S, double fa,
                           double fb, double fc, int bottom) {
  double c = (a + b) / 2, h = b - a;
  double d = (a + c) / 2, e = (c + b) / 2;
  double fd = (pf->*f)(d), fe = (pf->*f)(e);
  double Sleft = (h / 12) * (fa + 4 * fd + fc);
  double Sright = (h / 12) * (fc + 4 * fe + fb);
  double S2 = Sleft + Sright;
  if (bottom <= 100 || fabs(S2 - S) / (1e-10 + fabs(S2)) <= epsilon ||
      IsInfNan(S2))
    return S2 + (S2 - S) / 8;
  return adaptiveSimpsonsAux(f, pf, a, c, epsilon, Sleft, fa, fc, fd,
                             bottom - 1) +
         adaptiveSimpsonsAux(f, pf, c, b, epsilon, Sright, fc, fb, fe,
                             bottom - 1);
}

// Adaptive Simpson's Rule in integrate in interval [a, b]
// https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
template <typename C>
double adaptiveSimpsons(
    double (C::*f)(const double&), C* pf, double a, double b, double epsilon,
    int maxRecursion) {  // error tolerance and recursion cap
  double c = (a + b) / 2, h = b - a;
  double fa = (pf->*f)(a), fb = (pf->*f)(b), fc = (pf->*f)(c);
  double S = (h / 6) * (fa + 4 * fc + fb);
  return adaptiveSimpsonsAux(f, pf, a, b, epsilon, S, fa, fb, fc, maxRecursion);
}

template <typename T>
void MyConcatenate(T& x, T y) {
  int n = y.size();
  for (int i = 0; i < n; i++) x.push_back(y[i]);
}

// computes cumulative sum of vector of integer
inline int MyCumsum(const IntegerVector& x, const int& n) {
  int out = 0;
  for (int i = 0; i < n; i++) out += x[i];
  return out;
}

// samples the state given a probability vector. the output is in [0,
// P.size()-1]
inline int sampleState(const NumericVector& P) {
  double u = runif(1, 0, 1)[0], cumP = P[0];
  int ct = 1, ct_max = P.size() - 1;
  while (u > cumP && ct <= ct_max) cumP += P[ct], ct++;
  return ct - 1;
}

// Rcpp implementation of v %*% M
inline NumericVector matrixProd(const NumericVector& v,
                                const NumericMatrix& M) {
  int n = v.size();
  NumericVector out(n);
  for (int i = 0; i < n; i++) out[i] = sum(v * M(_, i));
  return out;
}

#endif  // Utils.h
