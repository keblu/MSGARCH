#include <RcppArmadillo.h>
#include "Decoding.h"

using namespace Rcpp;
using namespace arma;

double MixtDensityScale(const arma::vec& vOmega, const arma::vec& vD_log,
                        const int& M) {
  arma::vec wp_log = log(vOmega) + vD_log;

  double dK = max(wp_log);

  arma::vec wp_log_scaled = wp_log - dK;

  double dLK = 0;
  for (int i = 0; i < M; i++) {
    dLK += exp(wp_log_scaled(i));
  }

  double dLLK = dK + log(dLK);
  // double dLK = as_scalar(vOmega.t() * exp(vD_log ));

  if (dLLK < -1e150) {
    dLLK = -1e50;
  }

  return dLLK;
}

double abs3(const double& x) {
  double abs_x = x;
  if (abs_x < 0) abs_x = -abs_x;
  return abs_x;
}

arma::vec AccessListVectors_vec(const List& list,
                                const std::string& element_name) {
  SEXP foo = wrap(as<NumericVector>(list[element_name]));
  arma::vec vec_out = Rcpp::as<arma::vec>(foo);
  return vec_out;
}

arma::mat AccessListVectors_mat(const List& list,
                                const std::string& element_name) {
  SEXP foo = wrap(as<NumericMatrix>(list[element_name]));
  arma::mat mat_out = Rcpp::as<arma::mat>(foo);
  return mat_out;
}

arma::cube array2cube_2(const SEXP& myArray) {
  Rcpp::NumericVector vecArray(myArray);
  Rcpp::IntegerVector arrayDims = vecArray.attr("dim");

  arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1],
                       arrayDims[2], false);

  return (cubeArray);
}

//[[Rcpp::export]]
arma::vec getDelta(const arma::mat& gamma, const int& m) {
  arma::mat I = eye(m, m);
  arma::mat Umat = ones(m, m);

  arma::vec Uvec(m);
  Uvec.fill(1);

  arma::mat foo = (I - gamma + Umat).t();

  arma::mat delta = (foo).i() * Uvec;

  return delta;
}

List StartingValueEM_HMM(const arma::vec& vY, const int& K) {
  double dMu = mean(vY);
  double dSigma2 = var(vY);

  double start = 0.8;
  double end = 1.2;
  double by = (end - start) / (K * 1.0);

  double foo = start - by;

  arma::vec vMu(K);
  arma::vec vSigma2(K);

  arma::mat mGamma(K, K);

  mGamma.fill(0.1 / (K - 1.0));

  for (int j = 0; j < K; j++) {
    vMu(j) = dMu * foo;
    vSigma2(j) = dSigma2 * foo;
    mGamma(j, j) = 0.9;
    foo += by;
  }

  List out;
  out["vMu"] = vMu;
  out["vSigma2"] = vSigma2;
  out["mGamma"] = mGamma;

  return out;
}

List StartingValueEM_MM(const arma::vec& vY, const int& K) {
  int iT = vY.size();
  int j;
  int t;

  double dMu = mean(vY);
  double dSigma2 = var(vY);

  double start = 0.8;
  double end = 1.2;
  double by = (end - start) / (K * 1.0);

  double foo = start - by;

  arma::vec vMu(K);
  arma::vec vSigma2(K);

  arma::vec vP(K);
  vP.fill(1.0 / (K * 1.0));

  for (j = 0; j < K; j++) {
    vMu(j) = dMu * foo;
    vSigma2(j) = dSigma2 * foo;
    foo += by;
  }

  // initialize weights
  arma::mat mLLK(K, iT);
  arma::mat mW(K, iT);
  arma::vec vLLK(iT);
  for (t = 0; t < iT; t++) {
    for (j = 0; j < K; j++) {
      mLLK(j, t) = Rf_dnorm4(vY(t), vMu(j), pow(vSigma2(j), 0.5), 1);
    }
    vLLK(t) = MixtDensityScale(vP, mLLK.col(t), K);
    for (j = 0; j < K; j++) {
      mW(j, t) = exp(log(vP(j)) + mLLK(j, t) - vLLK(t));
    }
  }

  for (j = 0; j < K; j++) {
    vP(j) = accu(mW.row(j));
  }

  vP = vP / (iT * 1.0);

  List out;
  out["vMu"] = vMu;
  out["vSigma2"] = vSigma2;
  out["vP"] = vP;

  return out;
}

arma::mat GaussianLk(const arma::vec& vY, const arma::vec& vMu,
                     const arma::vec& vSigma2, const int& K, const int& T,
                     const int& lg) {
  arma::mat lk(T, K);

  int i, j;

  for (i = 0; i < T; i++) {
    for (j = 0; j < K; j++) {
      lk(i, j) = R::dnorm4(vY(i), vMu(j), sqrt(vSigma2(j)), lg);
      if (lk(i, j) < 1e-250 && !lg) lk(i, j) = 1e-250;
    }
  }

  return lk;
}

List HMMlalphabeta(const arma::vec vY, const arma::mat mGamma,
                   const arma::vec vMu, const arma::vec vSigma2, const int T,
                   const int K) {
  arma::vec vDelta = getDelta(mGamma, K);

  arma::mat allprobs = GaussianLk(vY, vMu, vSigma2, K, T, 0);

  List FB = FFBS(allprobs, vDelta, mGamma, K, T);

  FB["allprobs"] = allprobs;

  return FB;
}

int WhichMax(arma::vec vX) {
  int iK = vX.size();
  int k;

  int iMax = 0;
  double dMax = vX(0);

  for (k = 1; k < iK; k++) {
    if (vX(k) > dMax) {
      iMax = k;
      dMax = vX(k);
    }
  }

  return iMax;
}

// mLLK is K x T

//[[Rcpp::export]]
arma::vec Viterbi(const arma::mat& mLLK, const arma::mat& mGamma,
                  const int& iK) {
  int iT = mLLK.n_cols;
  int t;
  int k;

  arma::vec vStates(iT);
  arma::vec vDelta = getDelta(mGamma, iK);
  arma::mat mXi(iK, iT);

  arma::mat mLK(iK, iT);
  for (t = 0; t < iT; t++) {
    for (k = 0; k < iK; k++) {
      mLK(k, t) = exp(mLLK(k, t));
    }
  }

  arma::vec vFoo = vDelta % mLK.col(0);
  arma::vec vMax(iK);
  arma::vec vDecoded(iT);

  mXi.col(0) = vFoo / accu(vFoo);

  for (t = 1; t < iT; t++) {
    for (k = 0; k < iK; k++) {
      vMax(k) = max(mXi.col(t - 1) % mGamma.col(k));
    }
    vFoo = vMax % mLK.col(t);
    mXi.col(t) = vFoo / accu(vFoo);
  }

  vDecoded(iT - 1) = WhichMax(mXi.col(iT - 1));

  for (t = iT - 2; t >= 0; t--) {
    vFoo = mGamma.col(vDecoded(t + 1)) % mXi.col(t);
    vDecoded(t) = WhichMax(vFoo);
  }

  return vDecoded;
}

//[[Rcpp::export]]
List EM_HMM(const arma::vec& vY, const int& K, const int& maxIter = 1e3,
            const double& tol = 1e-8, const bool& constraintZero = true) {
  List lStarting = StartingValueEM_HMM(vY, K);
  arma::vec vMu = AccessListVectors_vec(lStarting, "vMu");
  arma::vec vSigma2 = AccessListVectors_vec(lStarting, "vSigma2");
  arma::mat mGamma = AccessListVectors_mat(lStarting, "mGamma");

  if (constraintZero) {
    vMu.zeros();
  }
  arma::vec vMu_Next = vMu;
  arma::vec vSigma2_Next = vSigma2;
  arma::mat mGamma_Next = mGamma;

  int T = vY.size();

  List fb;
  arma::mat lalpha(K, T);
  arma::mat lbeta(K, T);
  arma::mat allprobs(T, K);
  arma::mat SmoothProb(T, K);
  arma::mat PredictedProb(T, K);
  arma::mat FilteredProb(T, K);

  int iter = 0;
  int i, j, b;

  arma::vec LLKSeries(maxIter + 1);

  double eps = 1.0;

  fb = HMMlalphabeta(vY, mGamma, vMu, vSigma2, T, K);
  lalpha = AccessListVectors_mat(fb, "lalpha");

  double c = max(lalpha.col(T - 1));
  double llk = c + log(sum(exp(lalpha.col(T - 1) - c)));

  LLKSeries(0) = llk;

  while (eps > tol && iter < maxIter) {
    fb = HMMlalphabeta(vY, mGamma, vMu, vSigma2, T, K);
    lalpha = AccessListVectors_mat(fb, "lalpha");
    lbeta = AccessListVectors_mat(fb, "lbeta");
    allprobs = AccessListVectors_mat(fb, "allprobs");

    c = max(lalpha.col(T - 1));
    llk = c + log(sum(exp(lalpha.col(T - 1) - c)));

    for (j = 0; j < K; j++) {
      for (b = 0; b < K; b++) {
        mGamma_Next(j, b) =
            mGamma(j, b) * sum(exp(lalpha.row(j).subvec(0, T - 2).t() +
                                   log(allprobs.col(b).subvec(1, T - 1)) +
                                   lbeta.row(b).subvec(1, T - 1).t() - llk));
      }
    }

    for (j = 0; j < K; j++) {
      mGamma_Next.row(j) = mGamma_Next.row(j) / sum(mGamma_Next.row(j));
      // Update Mu and Sigma
      SmoothProb.col(j) = exp(lalpha.row(j) + lbeta.row(j) - llk).t();

      if (!constraintZero) {
        vMu_Next(j) = sum(SmoothProb.col(j) % vY) / sum(SmoothProb.col(j));
      }

      vSigma2_Next(j) = sum(SmoothProb.col(j) % pow(vY - vMu_Next(j), 2.0)) /
                        sum(SmoothProb.col(j));
    }
    // Store the llk
    LLKSeries(iter) = llk;
    iter += 1;

    if (iter > 10)
      eps = abs3((llk - LLKSeries(iter - 2)) / (LLKSeries(iter - 2) + 1.0));

    // Update Parameters

    vMu = vMu_Next;
    vSigma2 = vSigma2_Next;
    mGamma = mGamma_Next;
  }
  double llk_foo = 0;

  for (i = 0; i < T; i++) {
    c = max(lalpha.col(i));
    llk_foo = c + log(sum(exp(lalpha.col(i) - c)));
    FilteredProb.row(i) = exp(lalpha.col(i).t() - llk_foo);
    if (i < T - 1) {
      PredictedProb.row(i + 1) = FilteredProb.row(i) * mGamma;
    }
  }
  PredictedProb.row(0) = SmoothProb.row(0);

  arma::vec vDelta = getDelta(mGamma_Next, K);

  // Decoding
  arma::mat mLLK = log(allprobs);
  arma::vec vDecoding = Viterbi(mLLK.t(), mGamma, K);

  LLKSeries = LLKSeries.subvec(0, iter - 2);

  List EMOut;

  EMOut["SmoothProb"] = SmoothProb;
  EMOut["PredictedProb"] = PredictedProb;
  EMOut["FilteredProb"] = FilteredProb;
  EMOut["vDecoding"] = vDecoding;

  EMOut["lalpha"] = lalpha;
  EMOut["lbeta"] = lbeta;
  EMOut["LLKSeries"] = LLKSeries;
  EMOut["mLLK"] = mLLK;
  EMOut["vMu"] = vMu_Next;
  EMOut["vSigma2"] = vSigma2_Next;
  EMOut["vY"] = vY;
  EMOut["mGamma"] = mGamma_Next;
  EMOut["vDelta"] = vDelta;
  EMOut["eps"] = eps;
  EMOut["iter"] = iter;

  return EMOut;
}

//[[Rcpp::export]]
List EM_MM(const arma::vec& vY, const int& K, const int& maxIter = 1e3,
           const double& tol = 1e-8, const bool& constraintZero = true) {
  List lStarting = StartingValueEM_MM(vY, K);
  arma::vec vMu = AccessListVectors_vec(lStarting, "vMu");
  arma::vec vSigma2 = AccessListVectors_vec(lStarting, "vSigma2");
  arma::vec vP = AccessListVectors_vec(lStarting, "vP");

  if (constraintZero) {
    vMu.zeros();
  }
  arma::vec vMu_Next = vMu;
  arma::vec vSigma2_Next = vSigma2;
  arma::vec vP_Next = vP;

  int T = vY.size();

  int iter = 0;
  int j, t;

  arma::vec LLKSeries(maxIter + 1);

  double eps = 1.0;

  arma::mat mLLK(K, T);
  arma::mat mW(K, T);
  arma::vec vLLK(T);
  for (t = 0; t < T; t++) {
    for (j = 0; j < K; j++) {
      mLLK(j, t) = Rf_dnorm4(vY(t), vMu(j), pow(vSigma2(j), 0.5), 1);
    }
    vLLK(t) = MixtDensityScale(vP, mLLK.col(t), K);
    for (j = 0; j < K; j++) {
      mW(j, t) = exp(log(vP(j)) + mLLK(j, t) - vLLK(t));
    }
  }

  LLKSeries(0) = accu(vLLK);

  while (eps > tol && iter < maxIter) {
    mLLK.zeros();
    mW.zeros();
    vLLK.zeros();
    vMu_Next.zeros();
    vSigma2_Next.zeros();

    for (t = 0; t < T; t++) {
      for (j = 0; j < K; j++) {
        mLLK(j, t) = Rf_dnorm4(vY(t), vMu(j), pow(vSigma2(j), 0.5), 1);
      }
      vLLK(t) = MixtDensityScale(vP, mLLK.col(t), K);
      for (j = 0; j < K; j++) {
        mW(j, t) = exp(log(vP(j)) + mLLK(j, t) - vLLK(t));
        if (!constraintZero) {
          vMu_Next(j) += mW(j, t) * vY(t);
        }
      }
    }

    for (j = 0; j < K; j++) {
      vP_Next(j) = accu(mW.row(j));
      vMu_Next(j) = vMu_Next(j) / vP_Next(j);
      for (t = 0; t < T; t++) {
        vSigma2_Next(j) += mW(j, t) * pow(vY(t) - vMu_Next(j), 2.0);
      }
      vSigma2_Next(j) = vSigma2_Next(j) / vP_Next(j);
    }

    vP_Next = vP_Next / (T * 1.0);

    // Store the llk
    LLKSeries(iter) = accu(vLLK);
    iter += 1;

    if (iter > 10)
      eps = abs3((LLKSeries(iter - 1) - LLKSeries(iter - 2)) /
                 (LLKSeries(iter - 2) + 1.0));

    // Update Parameters

    vMu = vMu_Next;
    vSigma2 = vSigma2_Next;
    vP = vP_Next;
  }

  // Decoding
  arma::vec vDecoding(T);
  for (t = 0; t < T; t++) {
    vDecoding(t) = WhichMax(log(vP) + mLLK.col(t) - vLLK(t));
  }

  List EMOut;

  EMOut["mW"] = mW;

  LLKSeries = LLKSeries.subvec(0, iter - 2);

  EMOut["LLKSeries"] = LLKSeries;
  EMOut["mLLK"] = mLLK;
  EMOut["vDecoding"] = vDecoding;
  EMOut["vMu"] = vMu;
  EMOut["vSigma2"] = vSigma2;
  EMOut["vY"] = vY;
  EMOut["vP"] = vP;
  EMOut["eps"] = eps;
  EMOut["iter"] = iter;

  return EMOut;
}
