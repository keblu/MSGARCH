#include <RcppArmadillo.h>
#include "EM.h"

using namespace arma;
using namespace Rcpp;

List FFBS(const arma::mat& allprobs, const arma::vec& delta,
          const arma::mat& mGamma, const int& K, const int& T) {
  arma::mat lalpha = zeros(K, T);
  arma::mat lbeta = zeros(K, T);

  arma::vec foo(K);
  double sumfoo, lscale;
  int i;

  foo = delta % allprobs.row(0).t();
  sumfoo = sum(foo);
  lscale = log(sumfoo);
  foo = foo / sumfoo;

  lalpha.col(0) = log(foo) + lscale;
  for (i = 1; i < T; i++) {
    foo = (foo.t() * mGamma).t() % allprobs.row(i).t();
    sumfoo = sum(foo);
    lscale = lscale + log(sumfoo);
    foo = foo / sumfoo;
    lalpha.col(i) = log(foo) + lscale;
  }
  for (i = 0; i < K; i++) {
    foo(i) = 1.0 / K;
  }
  lscale = log(K);
  for (i = T - 2; i >= 0; i--) {
    foo = mGamma * (allprobs.row(i + 1).t() % foo);
    lbeta.col(i) = log(foo) + lscale;
    sumfoo = sum(foo);
    foo = foo / sumfoo;
    lscale = lscale + log(sumfoo);
  }

  List FS;
  FS["lalpha"] = lalpha;
  FS["lbeta"] = lbeta;

  return FS;
}

//[[Rcpp::export]]
List Decoding_HMM(const arma::mat& allprobs, const arma::mat& mGamma,
                  const int& T, const int& K) {
  arma::vec vDelta = getDelta(mGamma, K);

  List FilterSmoother = FFBS(allprobs, vDelta, mGamma, K, T);

  arma::mat lalpha = FilterSmoother["lalpha"];
  arma::mat lbeta = FilterSmoother["lbeta"];

  arma::mat SmoothProb(T, K);
  arma::mat PredictedProb(T + 1, K);
  arma::mat FilteredProb(T, K);

  double c = max(lalpha.col(T - 1));
  double llk = c + log(sum(exp(lalpha.col(T - 1) - c)));
  double llk_foo = 0.0;

  int j;
  int t;

  for (j = 0; j < K; j++) {
    SmoothProb.col(j) = exp(lalpha.row(j) + lbeta.row(j) - llk).t();
  }

  for (t = 0; t < T; t++) {
    c = max(lalpha.col(t));
    llk_foo = c + log(sum(exp(lalpha.col(t) - c)));

    FilteredProb.row(t) = exp(lalpha.col(t).t() - llk_foo);

    // if(t < T-1){
    PredictedProb.row(t + 1) = FilteredProb.row(t) * mGamma;
    // }
  }

  PredictedProb.row(0) = vDelta.t();

  List lOut;

  lOut["PredictedProb"] = PredictedProb;
  lOut["FilteredProb"] = FilteredProb;
  lOut["SmoothProb"] = SmoothProb;

  return lOut;
}
