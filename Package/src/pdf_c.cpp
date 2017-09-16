/*################################################################################
## This code is part of the rugarch package of Ghalanos (2016)
## see https://cran.r-project.org/web/packages/rugarch/index.html
#################################################################################*/
// #include <R.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

/*
* -----------------------------------------
* Key Functions
* -----------------------------------------
*/

double dnormstd(const double& x) {
  double pdf;
  pdf = exp(-0.5 * x * x) / sqrt(2.0 * PI);
  if (pdf == 0.0) pdf = 0.0 + 2.22507e-24;
  return pdf;
}
double signum(const double& x) {
  double res = -(x < 0) + (x > 0);
  return res;
}
double dhuge(void) { return HUGE_VAL; }
double Heaviside(const double& x, const double& a) {
  return ((signum(x - a) + 1.0) / 2.0);
}
double depsilon(void) {
  double r;
  r = 1.0;
  while (1.0 < (double)(1.0 + r)) {
    r = r / 2.0;
  }
  return (2.0 * r);
}

/*
* -----------------------------------------
* Student Distribution
* -----------------------------------------
*/

double xdt(const double& x, const double& nu) {
  double a, b, pdf;
  a = Rf_gammafn((nu + 1.0) / 2.0) / sqrt(PI * nu);
  b = Rf_gammafn(nu / 2.0) * pow((1.0 + (x * x) / nu), ((nu + 1.0) / 2.0));
  pdf = a / b;
  return pdf;
}

double dstdstd(const double& x, const double& nu) {
  double pdf, s;
  if (nu <= 2) {
    pdf = 999;
  } else {
    s = sqrt(nu / (nu - 2.0));
    pdf = s * xdt(x * s, nu);
  }
  return pdf;
}

/*
* -----------------------------------------
* Skew Student Distribution (Fernandez & Steel)
* -----------------------------------------
*/
double dsstdstd(const double& x, const double& xi, const double& nu) {
  double mu, m1, beta, sigma, z, g, pdf, a, b, xxi;
  xxi = xi;
  a = 1.0 / 2.0;
  b = nu / 2.0;
  beta = (Rf_gammafn(a) / Rf_gammafn(a + b)) * Rf_gammafn(b);
  m1 = 2.0 * sqrt(nu - 2.0) / (nu - 1.0) / beta;
  mu = m1 * (xi - 1.0 / xi);
  sigma = sqrt((1.0 - pow(m1, 2)) * (pow(xi, 2) + 1.0 / (pow(xi, 2))) +
               2.0 * pow(m1, 2) - 1.0);
  z = x * sigma + mu;
  if (z == 0) {
    xxi = 1;
  }
  if (z < 0) {
    xxi = 1 / xi;
  }
  g = 2.0 / (xi + 1.0 / xi);
  pdf = g * dstdstd(z / xxi, nu) * sigma;
  return pdf;
}

/*
* -----------------------------------------
* Skew Normal Distribution
* -----------------------------------------
*/

double dsnormstd(const double& x, const double& xi) {
  double pdf;
  double mu, sigma, z, xxi, g;
  double m1 = 2.0 / sqrt(2.0 * PI);
  double m12 = m1 * m1;
  double xi2 = xi * xi;
  mu = m1 * (xi - 1.0 / xi);
  sigma = sqrt((1 - m12) * (xi2 + 1.0 / xi2) + 2 * m12 - 1);
  z = x * sigma + mu;
  xxi = (z < 0) ? 1.0 / xi : xi;
  g = 2.0 / (xi + 1.0 / xi);
  pdf = g * dnormstd(z / xxi) * sigma;
  return pdf;
}

/*
* -----------------------------------------
* Generalized Error Distribution
* -----------------------------------------
*/

double dgedstd(const double& x, const double& nu) {
  double lambda, g, pdf;
  lambda = sqrt(pow(1.0 / 2.0, 2.0 / nu) * Rf_gammafn(1.0 / nu) /
                Rf_gammafn(3.0 / nu));
  g = nu / (lambda * (pow(2.0, 1.0 + (1.0 / nu))) * Rf_gammafn(1.0 / nu));
  pdf = g * exp(-0.5 * pow(fabs(x / lambda), nu));
  return pdf;
}

/*
* -----------------------------------------
* Skew Generalized Error Distribution (Fernandez & Steel)
* -----------------------------------------
*/

double dsgedstd(const double& x, const double& xi, const double& nu) {
  double lambda, m1, mu, sigma, z, g, pdf, xxi;
  xxi = xi;
  lambda = sqrt(pow(1.0 / 2.0, (2 / nu)) * Rf_gammafn(1.0 / nu) /
                Rf_gammafn(3.0 / nu));
  g = nu / (lambda * (pow(2.0, 1.0 + (1.0 / nu))) * Rf_gammafn(1.0 / nu));
  m1 = pow(2.0, (1.0 / nu)) * lambda * Rf_gammafn(2.0 / nu) /
       Rf_gammafn(1.0 / nu);
  mu = m1 * (xi - 1.0 / xi);
  sigma = (1 - pow(m1, 2.0)) * (pow(xi, 2.0) + 1.0 / (pow(xi, 2.0))) +
          2.0 * (pow(m1, 2)) - 1.0;
  sigma = sqrt(sigma);
  z = x * sigma + mu;
  if (z == 0) {
    xxi = 1;
  }
  if (z < 0) {
    xxi = 1 / xi;
  }
  g = 2.0 / (xi + 1.0 / xi);
  pdf = g * dgedstd(z / xxi, nu) * sigma;
  return pdf;
}

/*
* dZ is a standardized variable
* xi, nu
* skewness, shape
*/

double ddist_theta_param(const double& dZ, const std::string& sDist,
                         const bool& bSkew, const bool& bLog,
                         const double& dXi = 1.0, const double& dNu = 7.0) {
  double dPdf = 0.0;

  if (bSkew) {
    if (sDist == "norm") {
      dPdf = dsnormstd(dZ, dXi);
    }
    if (sDist == "std") {
      dPdf = dsstdstd(dZ, dXi, dNu);
    }
    if (sDist == "ged") {
      dPdf = dsgedstd(dZ, dXi, dNu);
    }
  } else {
    if (sDist == "norm") {
      dPdf = dnormstd(dZ);
    }
    if (sDist == "std") {
      dPdf = dstdstd(dZ, dNu);
    }
    if (sDist == "ged") {
      dPdf = dgedstd(dZ, dNu);
    }
  }

  if (dPdf < 1e-50) {
    dPdf = 1e-50;
  }

  if (bLog) {
    dPdf = log(dPdf);
  }

  return dPdf;
}

//[[Rcpp::export]]
double dUnivLike(const arma::vec& vZ, const std::string& sDist,
                 const bool& bSkew, const double& dXi = 1.0,const double& dNu = 7.0) {
  double dLLK = 0.0;

  int iT = vZ.size();
  int i;

  for (i = 0; i < iT; i++) {
    dLLK += ddist_theta_param(vZ(i), sDist, bSkew, true, dXi, dNu);
  }

  return dLLK;
}
