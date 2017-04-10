#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

const double dLowerShape = 4.0;
const double dUpperShape = 50.0;

const double dLowerSkewFS = 0.50;
const double dUpperSkewFS = 1.50;

const double dNumericalUpperLimit = 1e50;
const double dNumericalLowerLimit = -1e50;
const double dNumericalUnderflow  = 1e-50;

double Map(double dX, double dL,double dU) {
  double dMap =  dL + ( dU - dL ) / (1.0 + exp( - dX ));
  return dMap;
}

double Unmap(double dG, double dL,double dU) {
  double dUnmap = log((dG-dL)/(dU-dG));
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
arma::vec MapParameters_univ(arma::vec vTheta_tilde, std::string Dist, bool bSkew){
  
  int iL = vTheta_tilde.size();
  arma::vec vTheta(iL);
  
  if (!bSkew) {
    
    if((Dist == "std") | (Dist == "ged")){
      
      double dNu_tilde  = vTheta_tilde(0);
      
      double dNu  = Map(dNu_tilde,dLowerShape,dUpperShape);
      
      vTheta(0) = dNu;
    }
  } else {
    
    if((Dist == "std") | (Dist == "ged")){
      
      double dNu_tilde    = vTheta_tilde(0);
      double dXi_tilde    = vTheta_tilde(1);
      
      double dNu  = Map(dNu_tilde,dLowerShape,dUpperShape);
      double dXi    = Map(dXi_tilde, dLowerSkewFS, dUpperSkewFS);
      
      vTheta(0) = dNu;
      vTheta(1) = dXi;
      
    }
    
    if(Dist == "norm"){
      
      double dXi_tilde  = vTheta_tilde(0);
      
      double dXi  = Map(dXi_tilde, dLowerSkewFS, dUpperSkewFS);
      
      vTheta(0) = dXi;
      
    }
  }
  
  return vTheta;
  
}

//[[Rcpp::export]]
arma::vec UnmapParameters_univ(arma::vec vTheta, std::string Dist, bool bSkew){
  
  int iL = vTheta.size();
  arma::vec vTheta_tilde(iL);
  
  if (!bSkew) {
    
    if((Dist == "std") | (Dist == "ged")){
      
      double dNu  = vTheta(0);
      
      double dNu_tilde  = Unmap(dNu,dLowerShape,dUpperShape);
      
      vTheta_tilde(0) = dNu_tilde;
    }
  } else {
    
    if((Dist == "std") | (Dist == "ged")){
      
      double dNu    = vTheta(0);
      double dXi    = vTheta(1);
      
      double dNu_tilde  = Unmap(dNu,dLowerShape,dUpperShape);
      double dXi_tilde    = Unmap(dXi, dLowerSkewFS, dUpperSkewFS);
      
      vTheta_tilde(0) = dNu_tilde;
      vTheta_tilde(1) = dXi_tilde;
      
    }
    
    if(Dist == "norm"){
      
      double dXi  = vTheta(0);
      
      double dXi_tilde  = Unmap(dXi, dLowerSkewFS, dUpperSkewFS);
      
      vTheta_tilde(0) = dXi_tilde;
      
    }
  }
  
  return vTheta_tilde;
  
}