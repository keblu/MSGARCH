#ifndef DECODING_H
#define DECODING_H

Rcpp::List FFBS(const arma::mat& allprobs,const arma::vec& delta,const arma::mat& mGamma,const int& K,
                const int& T);
Rcpp::List Decoding_HMM(const arma::mat& allprobs,const arma::mat& mGamma,const int& T,const int& K);

#endif
