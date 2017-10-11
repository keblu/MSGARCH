#include <string.h>
#include <RcppArmadillo.h>
#include <math.h> 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

arma::mat makePositiveDefinite(arma::mat& semiDefMat) {
  try{
    // make positive definite, by first taking an eigenvalue decomposition
    arma::vec eigVal;
    arma::mat eigVec;
    arma::eig_sym(eigVal,eigVec,semiDefMat);
    // what should the minimal eigenvalue be? Note that the last element of eigVal is the largest. If all eigenvalues are negative, then return the zero matrix.
    if (eigVal.max() < 0) {
      semiDefMat *= 0;
    } else {
      double minEig = 2.23e-16;
      for (int uu=0; uu < eigVal.n_elem - 1; uu++) {
        if (eigVal(uu) < minEig) {
          eigVal(uu) = minEig;
        }
      }
      semiDefMat = eigVec * diagmat(eigVal) * trans(eigVec);
    }
    return semiDefMat;
  } catch( std::exception& __ex__) {
    forward_exception_to_r(__ex__);
  } catch(...) {
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return 0;
}

// [[Rcpp::export]]
arma::mat f_RCPP_adaptMCMC(const arma::rowvec theta0,
                      Rcpp::Function func,
					            double acc_rate,
                      arma::mat sigma,					  
                      int n_mcmc){
  
  //intialization
  int l_param = theta0.size();
  arma::mat chain(n_mcmc, l_param);
  arma::mat matrixRandomParameters(n_mcmc,l_param); // random variable scale with xtra$sigma
  for (int i = 0; i < l_param; i++) {
    matrixRandomParameters.col(i) = as<arma::vec>(rnorm(n_mcmc));
  }
  arma::vec eigVal;
  arma::mat eigVec;
  arma::vec UnifRandomVariable = as<arma::vec>(runif(n_mcmc - 1));
  arma::rowvec vectorProposalParameter(l_param);
  arma::rowvec currentParams = theta0;
  double currentLikelihood = Rcpp::as<double>(func(theta0));
  double proposalLikelihood = 0;
  double probab = 0;
  double score = 0;
  double gamma = 0.5;
  double B;
  double alpha;
  double adaptrate;
  chain.row(0) = currentParams;
  arma::rowvec draw(l_param);
  for (int position = 1; position < n_mcmc; position++) {
    currentParams = chain.row(position-1);
    draw = matrixRandomParameters.row(position-1);
  	vectorProposalParameter = currentParams + draw * arma::chol(sigma);
    proposalLikelihood = Rcpp::as<double>(func(vectorProposalParameter)); 
    probab = proposalLikelihood - currentLikelihood;
    score = std::min(1.0,std::exp(probab));
    if(UnifRandomVariable(position-1) < score) { 
      chain.row(position) = vectorProposalParameter; 
      currentLikelihood = proposalLikelihood;
    } else {
      chain.row(position) = chain.row(position-1);
    }
    arma::mat S = arma::chol(sigma).t();
    adaptrate = std::min(5.0, l_param * pow(position + 1,-gamma));
    alpha = std::min(1.0, score);
    B = arma::sum(draw % draw);
    sigma = S * (sigma.eye(l_param,l_param) + adaptrate * (alpha - acc_rate) * draw.t() * draw / B) * S.t();
    arma::eig_sym(eigVal,eigVec,sigma);
    sigma = makePositiveDefinite(sigma);
  }
  return(chain);
}


