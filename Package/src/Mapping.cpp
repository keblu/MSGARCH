#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

const double dLowerShape = 4.0;
const double dUpperShape = 50.0;

const double dLowerSkewFS = 0.50;
const double dUpperSkewFS = 1.50;

const double dNumericalUpperLimit = 1e50;
const double dNumericalLowerLimit = -1e50;
const double dNumericalUnderflow = 1e-50;

double Map(const double& dX, const double& dL, const double& dU) {
  double dMap = dL + (dU - dL) / (1.0 + exp(-dX));
  return dMap;
}

double Unmap(const double& dG, const double& dL, const double& dU) {
  double dUnmap = log((dG - dL) / (dU - dG));
  return dUnmap;
}

double CheckScale(double dScale) {
  if (dScale > dNumericalUpperLimit) {
    dScale = dNumericalLowerLimit;
  }
  if (dScale < dNumericalUnderflow) {
    dScale = dNumericalUnderflow;
  }
  return dScale;
}

double CheckLocation(double dLocation) {
  if (dLocation > dNumericalUpperLimit) {
    dLocation = dNumericalLowerLimit;
  }
  if (dLocation < dNumericalLowerLimit) {
    dLocation = dNumericalUnderflow;
  }
  return dLocation;
}

//[[Rcpp::export]]
arma::vec MapParameters_univ(const arma::vec& vTheta_tilde,
                             const std::string& Dist, const bool& bSkew) {
  int iL = vTheta_tilde.size();
  arma::vec vTheta(iL);

  if (!bSkew) {
    if ((Dist == "std") || (Dist == "ged")) {
      double dNu_tilde = vTheta_tilde(0);

      double dNu = Map(dNu_tilde, dLowerShape, dUpperShape);

      vTheta(0) = dNu;
    }
  } else {
    if ((Dist == "std") || (Dist == "ged")) {
      double dNu_tilde = vTheta_tilde(0);
      double dXi_tilde = vTheta_tilde(1);

      double dNu = Map(dNu_tilde, dLowerShape, dUpperShape);
      double dXi = Map(dXi_tilde, dLowerSkewFS, dUpperSkewFS);

      vTheta(0) = dNu;
      vTheta(1) = dXi;
    }

    if (Dist == "norm") {
      double dXi_tilde = vTheta_tilde(0);

      double dXi = Map(dXi_tilde, dLowerSkewFS, dUpperSkewFS);

      vTheta(0) = dXi;
    }
  }

  return vTheta;
}

//[[Rcpp::export]]
arma::vec UnmapParameters_univ(const arma::vec& vTheta, const std::string& Dist,
                               const bool& bSkew) {
  int iL = vTheta.size();
  arma::vec vTheta_tilde(iL);

  if (!bSkew) {
    if ((Dist == "std") || (Dist == "ged")) {
      double dNu = vTheta(0);

      double dNu_tilde = Unmap(dNu, dLowerShape, dUpperShape);

      vTheta_tilde(0) = dNu_tilde;
    }
  } else {
    if ((Dist == "std") || (Dist == "ged")) {
      double dNu = vTheta(0);
      double dXi = vTheta(1);

      double dNu_tilde = Unmap(dNu, dLowerShape, dUpperShape);
      double dXi_tilde = Unmap(dXi, dLowerSkewFS, dUpperSkewFS);

      vTheta_tilde(0) = dNu_tilde;
      vTheta_tilde(1) = dXi_tilde;
    }

    if (Dist == "norm") {
      double dXi = vTheta(0);

      double dXi_tilde = Unmap(dXi, dLowerSkewFS, dUpperSkewFS);

      vTheta_tilde(0) = dXi_tilde;
    }
  }

  return vTheta_tilde;
}

//////////////////////////// SIMPLEX MAPPING ////////////////////

double Logit(double dP) {
  if (dP < 1e-10) {
    dP = 1e-10;
  }
  if (dP > 1.0 - 1e-10) {
    dP = 1.0 - 1e-10;
  }
  double dLogit = log(dP) - log(1.0 - dP);
  return dLogit;
}
double LogitInv(const double& dLogit) {
  double logx = 0.0;
  double logy = dLogit;
  double dFoo = 0.0;
  if (logx > logy) {
    dFoo = logx + log(1.0 + exp(logy - logx));
  } else {
    dFoo = logy + log(1.0 + exp(logx - logy));
  }
  double dP = exp(dLogit - dFoo);
  return dP;
}
//[[Rcpp::export]]
arma::vec SimplexUnmapping(const arma::vec& vOmega, const int& iK) {
  int k;
  double dFoo = 1.0;
  arma::vec vPhi(iK - 1);
  for (k = 0; k < iK - 1; k++) {
    if (k == 0) {
      vPhi(k) = Logit(vOmega(k));
    } else {
      vPhi(k) = Logit(vOmega(k) / dFoo);
    }
    dFoo = dFoo * (1.0 - LogitInv(vPhi(k)));
  }
  return vPhi;
}
//[[Rcpp::export]]
arma::vec SimplexMapping(const arma::vec& vPhi, const int& iK) {
  int k;
  // arma::vec vOmega(iK);
  arma::vec vOmega(iK - 1);
  double dLogit_foo = LogitInv(vPhi(0));
  vOmega(0) = dLogit_foo;
  double dFoo = log(1.0 - vOmega(0));
  if (iK > 2) {
    for (k = 1; k < iK - 1; k++) {
      dLogit_foo = LogitInv(vPhi(k));
      vOmega(k) = exp(vPhi(k) - log(1.0 + exp(vPhi(k))) + dFoo);
      dFoo += log(1.0 - dLogit_foo);
    }
  }
  // vOmega(iK - 1) = exp(dFoo);
  return vOmega;
}
