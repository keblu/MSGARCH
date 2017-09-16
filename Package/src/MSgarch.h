#ifndef MSGARCH_H  // include guard
#define MSGARCH_H

#include "SingleRegime.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//======================================== AUXILIARY FUNCTIONS
//============================================//

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//========================================== MS-GARCH class
//===============================================//

// type definition (vector of pointers to Base)
typedef std::vector<Base*> many;
typedef std::vector<volatility> volatilityVector;

// MS-GARCH class
class MSgarch {
  many specs;           // vector of pointers to Base objects
  int K;                // number of models
  NumericMatrix P;      // transition-probability matrix
  NumericVector PLast;  // all transition-probability matrix at each step
  NumericVector P0;     // initial distribution of states
  double P_mean;        // mean for the prior on transition-probabilities
  double P_sd;          // sd for the prior on transition-probabilities
  double LND_MIN;       // minimum loglikelihood allowed
public:
  std::vector<std::string> name;
  NumericVector theta0;
  NumericVector Sigma0;
  CharacterVector label;
  NumericVector lower;
  NumericVector upper;
  NumericVector ineq_lb;
  NumericVector ineq_ub;
  IntegerVector NbParams;  // number of parameters for each model (excluding the
  // transition probabilities)
  IntegerVector NbParamsModel;
  // constructor
  MSgarch(List L) {
    // extract pointers to models
    K = L.length();  // number of models
    Environment env;
    for (List::iterator it = L.begin(); it != L.end();
    ++it) {  // loop over models
      env = *it;
      specs.push_back(
        static_cast<Base*>(R_ExternalPtrAddr(env.get(".pointer"))));
    }
    
    // loop over models
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      name.push_back((*it)->spec_name());
      MyConcatenate(theta0, (*it)->spec_theta0());
      MyConcatenate(Sigma0, (*it)->spec_Sigma0());
      MyConcatenate(label, (*it)->spec_label());
      MyConcatenate(lower, (*it)->spec_lower());
      MyConcatenate(upper, (*it)->spec_upper());
      MyConcatenate(ineq_lb, NumericVector::create((*it)->spec_ineq_lb()));
      MyConcatenate(ineq_ub, NumericVector::create((*it)->spec_ineq_ub()));
      NbParams.push_back((*it)->spec_nb_coeffs());
      NbParamsModel.push_back((*it)->spec_nb_coeffs_model());
    }
    P0 = rep(1.0 / K, K);
    PLast = rep(1.0 / K, K);
    P_mean = 1 / K;
    P_sd = 100;
    LND_MIN = log(DBL_MIN) + 1;
    
    // add transition-probabilities
    if (K > 1) {
      int NbP = K * (K - 1);
      NumericVector P_theta0  = rep(1.0 / K, NbP);  // theta0
      NumericVector P_Sigma0  = rep(1.0, NbP);      // Sigma0
      NumericVector P_lower   = rep(0.0, NbP);      // lower
      NumericVector P_upper   = rep(1.0, NbP);      // upper
      NumericVector P_ineq_lb = rep(0.0, K);       // ineq_lb
      NumericVector P_ineq_ub = rep(1.0, K);       // ineq_ub
      CharacterVector P_label(NbP, "P");           // label
      
      // concatenates all stuff above
      MyConcatenate(theta0, P_theta0);
      MyConcatenate(Sigma0, P_Sigma0);
      MyConcatenate(label, P_label);
      MyConcatenate(lower, P_lower);
      MyConcatenate(upper, P_upper);
      MyConcatenate(ineq_lb, P_ineq_lb);
      MyConcatenate(ineq_ub, P_ineq_ub);
    }
  }
  
  // set the parameters (including those of the distribution) of all models
  // the last elements of theta should be those of the transition-probability
  // matrix
  // this function should always be called first
  void loadparam(const NumericVector&);
  
  // to be called before 'calc_prior', 'ineq_func' or 'set_vol'
  void prep_ineq_vol() {
    for (many::iterator it = specs.begin(); it != specs.end(); ++it)
      (*it)->spec_prep_ineq_vol();
  }
  
  // to be called before 'calc_kernel'
  void prep_kernel() {
    for (many::iterator it = specs.begin(); it != specs.end(); ++it)
      (*it)->spec_prep_kernel();
  }
  
  // loglikelihood of a single observation for all models
  NumericVector calc_kernel(const volatilityVector& vol, const double& yi) {
    NumericVector lnd(K);
    int k = 0;
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      lnd[k] = (*it)->spec_calc_kernel(vol[k], yi);
      k++;
    }
    return lnd;
  }
  
  // initialize all volatilities to their undonditional expected value
  volatilityVector set_vol(const double& y0) {
    volatilityVector vol(K);
    int k = 0;
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      vol[k] = (*it)->spec_set_vol(y0);
      k++;
    }
    return vol;
  }
  
  // increment all volatilities
  void increment_vol(volatilityVector& vol, const double& yim1) {
    int k = 0;
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      (*it)->spec_increment_vol(vol[k], yim1);
      k++;
    }
  }
  
  NumericVector get_mean() {
    NumericVector out;
    for (many::iterator it = specs.begin(); it != specs.end();
    ++it) {  // loop over models
      MyConcatenate(out, (*it)->get_mean());
    }
    return (out);
  }
  
  void set_mean(NumericVector& new_mean) {
    int k = 0;
    NumericVector out(sum(NbParams));
    for (many::iterator it = specs.begin(); it != specs.end();
    ++it) {  // loop over models
      (*it)->set_mean(extract_theta_it(new_mean, k));
      k++;
    }
  }
  
  NumericVector get_sd() {
    NumericVector out;
    for (many::iterator it = specs.begin(); it != specs.end();
    ++it) {  // loop over models
      MyConcatenate(out, (*it)->get_sd());
    }
    return (out);
  }
  
  void set_sd(NumericVector& new_sd) {
    int k = 0;
    NumericVector out(sum(NbParams));
    for (many::iterator it = specs.begin(); it != specs.end();
    ++it) {  // loop over models
      (*it)->set_sd(extract_theta_it(new_sd, k));
      k++;
    }
  }
  
  // extract parameter vector of model 'k', where k is in [0, K-1]
  NumericVector extract_theta_it(const NumericVector& theta, const int& k) {
    int start = MyCumsum(NbParams, k);
    NumericVector theta_it(theta.begin() + start,
                           theta.begin() + start + NbParams[k]);
    return theta_it;
  }
  
  // extract transition-probability from state 'k', where k is in [0, K-1]
  NumericVector extract_P_it(const NumericVector& theta, const int& k) {
    int Tot_NbParams = sum(NbParams);
    NumericVector P_it(theta.begin() + Tot_NbParams + k * (K - 1),
                       theta.begin() + Tot_NbParams + (k + 1) * (K - 1));
    P_it.push_back(1 - sum(P_it));
    return P_it;
  }
  
  // simulate a random inovation from model 'k', where k is in [0, K-1]
  double rndgen(const int& k) {
    many::iterator it = specs.begin() + k;
    return (*it)->spec_rndgen(1)[0];
  }
  
  // inequality function
  NumericVector ineq_func(const NumericVector& theta) {
    NumericVector out;
    loadparam(theta);
    prep_ineq_vol();
    for (many::iterator it = specs.begin(); it != specs.end(); ++it)
      out.push_back((*it)->spec_ineq_func());
    if (K > 1) {
      NumericMatrix Psub = P(Range(0, K - 1), Range(0, K - 2));
      for (int i = 0; i < K; i++) out.push_back(sum(Psub(i, _)));
    }
    return out;
  }
  
  NumericVector get_p_last() { return PLast; }
  
  int get_K() { return K; }
  
  List f_get_Pstate(const NumericVector&, const NumericVector&);
  
  // check prior
  prior calc_prior(const NumericVector&);
  
  arma::cube calc_ht(NumericMatrix&, const NumericVector&);
  
  NumericVector f_pdf(const NumericVector&, const NumericVector&,
                      const NumericVector&, const bool&);
  
  arma::cube f_pdf_its(const NumericVector&, const NumericVector&,
                       const NumericMatrix&, const bool&);
  
  NumericVector f_cdf(const NumericVector&, const NumericVector&,
                      const NumericVector&, const bool&);
  
  arma::cube f_cdf_its(const NumericVector&, const NumericVector&,
                       const NumericMatrix&, const bool&);
  
  // model simulation
  Rcpp::List f_sim(const int&, const int&, const NumericVector&);
  
  Rcpp::List f_simAhead(const NumericVector&, const int&,  const int&, const NumericVector&,
                        const NumericVector&);
  
  Rcpp::List f_rnd(const int&, const NumericVector&, const NumericVector&);
  // compute loglikelihood matrix
  NumericMatrix calc_lndMat(const NumericVector&);
  
  NumericMatrix f_unc_vol(NumericMatrix&, const NumericVector&);
  
  // apply Hamilton filter
  double HamiltonFilter(const NumericMatrix&);
  
  // Model evaluation
  NumericVector eval_model(NumericMatrix&, const NumericVector&, const bool&);
};

//---------------------- load parameters of all models  ----------------------//
inline void MSgarch::loadparam(const NumericVector& theta) {
  // load the parameters of each model and the transition-probability matrix
  NumericMatrix P_mat(K, K);
  int k = 0;
  for (many::iterator it = specs.begin(); it != specs.end();
  ++it) {  // loop over models
    NumericVector theta_it =
    extract_theta_it(theta, k);  // parameters of model 'it'
    NumericVector P_it =
      extract_P_it(theta, k);  // transition probabilities from model 'it'
    (*it)->spec_loadparam(theta_it);
    P_mat(k, _) = P_it;
    k++;
  }
  P = P_mat;
  arma::mat I = arma::eye(K,K);
  arma::mat Umat = arma::ones(K,K);
    arma::vec Uvec(K);
  Uvec.fill(1);
  arma::mat foo = (I - as<arma::mat>(P_mat) + Umat).t();
  
  arma::vec delta = (foo).i() * Uvec;
  for(int i = 0; i < K; i++){
    P0(i) = delta(i);
  }
}

//------------------------------ Prior calculation
//------------------------------//
inline prior MSgarch::calc_prior(const NumericVector& theta) {
  // compute prior of individual models
  bool r1_joint = 1;    // joint r1 of the models
  double r2_joint = 0;  // joint r2 of the models
  double r3_joint = 0;
  int k = 0;
  prior pr;
  for (many::iterator it = specs.begin(); it != specs.end();
  ++it) {  // loop over models
    NumericVector theta_it =
    extract_theta_it(theta, k);  // parameters of model 'it'
    NumericVector P_it =
      extract_P_it(theta, k);  // transition probabilities from model 'it'
    pr = (*it)->spec_calc_prior(theta_it);  // prior of model 'it'
    r1_joint = r1_joint && pr.r1 && is_true(all((0 < P_it) & (P_it < 1)));
    r2_joint += pr.r2;
    r3_joint += sum(dunif(P_it, 0.0, 1.0, 1)) + pr.r3;
    k++;
  }
  // return result
  prior out;
  out.r1 = r1_joint;
  out.r2 = ((out.r1) ? r2_joint : -1e10);
  out.r3 = r3_joint;
  return out;
}

inline NumericMatrix MSgarch::f_unc_vol(NumericMatrix& all_thetas,
                                        const NumericVector& y) {
  int K = get_K();
  int nb_thetas = all_thetas.nrow();
  volatilityVector vol;
  NumericVector theta_j;
  NumericMatrix ht(nb_thetas, K);
  
  for (int j = 0; j < nb_thetas; j++) {  // loop over vectors of parameters
    theta_j = all_thetas(j, _);
    loadparam(theta_j);
    prep_ineq_vol();
    vol = set_vol(y[0]);
    // initialize volatility
    for (int s = 0; s < K; s++) {
      ht(j, s) = vol[s].h;
    }
  }
  return ht;
}

inline arma::cube MSgarch::calc_ht(NumericMatrix& all_thetas,
                                   const NumericVector& y) {
  int K = get_K();
  int nb_obs = y.size();
  int nb_thetas = all_thetas.nrow();
  volatilityVector vol;
  NumericVector theta_j;
  arma::cube ht(nb_obs + 1, nb_thetas, K);
  
  for (int j = 0; j < nb_thetas; j++) {  // loop over vectors of parameters
    theta_j = all_thetas(j, _);
    loadparam(theta_j);
    prep_ineq_vol();
    vol = set_vol(y[0]);
    // initialize volatility
    for (int s = 0; s < K; s++) {
      ht(0, j, s) = vol[s].h;
    }
    
    for (int i = 1; i <= nb_obs; i++) {  // loop over observations
      increment_vol(vol, y[i - 1]);      // increment all volatilities
      
      for (int s = 0; s < K; s++) {
        ht(i, j, s) = vol[s].h;
      }
    }
  }
  return ht;
}

inline NumericVector MSgarch::f_pdf(const NumericVector& x,
                                    const NumericVector& theta,
                                    const NumericVector& y,
                                    const bool& is_log) {
  // computes volatility
  int s = 0;
  int nx = x.size();
  int ny = y.size();
  double sig;
  NumericVector tmp(nx);
  NumericVector out(nx);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol(y[0]);  // initialize volatility
  
  for (int t = 0; t < ny; t++) increment_vol(vol, y[t]);
  
  HamiltonFilter(calc_lndMat(y));
  
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
    sig = sqrt(vol[s].h);
    // computes PDF
    for (int i = 0; i < nx; i++) {
      tmp[i] = (*it)->spec_calc_pdf(x[i] / sig) / sig;  //
      out[i] = out[i] + tmp[i] * PLast[s];
    }
    s++;
  }
  
  if (is_log) {
    for (int i = 0; i < nx; i++) {
      out[i] = log(tmp[i]);
    }
  }
  
  return out;
}

inline arma::cube MSgarch::f_pdf_its(const NumericVector& theta,
                                     const NumericVector& y,
                                     const NumericMatrix& x,
                                     const bool& is_log) {
  // computes volatility
  int s = 0;
  int ny = y.size();
  int nx = x.nrow();
  double sig;
  arma::cube tmp(nx, ny, K);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol(y[0]);  // initialize volatility
  
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
    sig = sqrt(vol[s].h);
    for (int ix = 0; ix < nx; ix++) {
      tmp(ix, 0, s) = (*it)->spec_calc_pdf(x(ix, 0) / sig) / sig;  //
    }
    s++;
  }
  
  for (int i = 1; i < ny; i++) {
    s = 0;
    increment_vol(vol, y[i - 1]);
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      sig = sqrt(vol[s].h);
      for (int ix = 0; ix < nx; ix++) {
        tmp(ix, i, s) = (*it)->spec_calc_pdf(x(ix, i) / sig) / sig;  //
      }
      s++;
    }
  }
  
  return tmp;
}

inline NumericVector MSgarch::f_cdf(const NumericVector& x,
                                    const NumericVector& theta,
                                    const NumericVector& y,
                                    const bool& is_log) {
  // computes volatility
  int s = 0;
  int nx = x.size();
  int ny = y.size();
  double sig;
  NumericVector tmp(nx);
  NumericVector out(nx);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol(y[0]);  // initialize volatility
  
  for (int t = 0; t < ny; t++) increment_vol(vol, y[t]);
  
  HamiltonFilter(calc_lndMat(y));
  
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
    sig = sqrt(vol[s].h);
    // computes CDF
    for (int i = 0; i < nx; i++) {
      tmp[i] = (*it)->spec_calc_cdf(x[i] / sig);
      out[i] = out[i] + tmp[i] * PLast[s];
    }
    s++;
  }
  
  if (is_log) {
    for (int i = 0; i < nx; i++) {
      out[i] = log(tmp[i]);
    }
  }
  
  return out;
}

inline arma::cube MSgarch::f_cdf_its(const NumericVector& theta,
                                     const NumericVector& y,
                                     const NumericMatrix& x,
                                     const bool& is_log) {
  // computes volatility
  int s = 0;
  int ny = y.size();
  double sig;
  int nx = x.nrow();
  arma::cube tmp(nx, ny, K);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol(y[0]);  // initialize volatility
  for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
    sig = sqrt(vol[s].h);
    for (int ix = 0; ix < nx; ix++) {
      tmp(ix, 0, s) = (*it)->spec_calc_cdf(x(ix, 0) / sig);  //
    }
    s++;
  }
  for (int i = 1; i < ny; i++) {
    s = 0;
    increment_vol(vol, y[i - 1]);
    for (many::iterator it = specs.begin(); it != specs.end(); ++it) {
      sig = sqrt(vol[s].h);
      for (int ix = 0; ix < nx; ix++) {
        tmp(ix, i, s) = (*it)->spec_calc_cdf(x(ix, i) / sig);  //
      }
      s++;
    }
  }
  
  return tmp;
}

//------------------------------ Model simulation
//------------------------------//
inline List MSgarch::f_sim(const int& n, const int& m, const NumericVector& theta) {
  // setup
  K = get_K();
  NumericMatrix y(m,n);  // observations
  NumericMatrix S(m,n);  // states of the Markov chain
  arma::cube CondVol(m,n,K);
  loadparam(theta);       // load parameters
  
  // generate first draw
  double z;            // random innovation from initial state
  prep_ineq_vol();                    // prep for 'set_vol'
  volatilityVector vol;  // initialize all volatilities
  for (int i = 0; i < m; i++) {
    S(i,0) = sampleState(P0);             // sample initial state
    z = rndgen(S(i,0));
    vol = set_vol(z);
    for (int s = 0; s < K; s++) {
      CondVol(i, 0, s) = sqrt(vol[s].h);
    }
    y(i,0) = z * sqrt(vol[S(i,0)].h);
    for (int t = 1; t < n; t++) {
      S(i,t) = sampleState(P(S(i,t - 1), _));  // sample new state
      z = rndgen(S(i,t));                    // sample new innovation
      increment_vol(vol, y(i,t - 1));        // increment all volatilities
      y(i,t) = z * sqrt(vol[S(i,t)].h);        // new draw
      for (int s = 0; s < K; s++) {
        CondVol(i, t, s) = sqrt(vol[s].h);
      }
    }
  }
  return (List::create(Rcpp::Named("draws") = y, Rcpp::Named("state") = S,  Rcpp::Named("CondVol") = CondVol));
}

inline List MSgarch::f_simAhead(const NumericVector& y, const int& n, const int& m,
                                const NumericVector& theta,
                                const NumericVector& P0_) {
  // setup
  K = get_K();
  int nb_obs = y.size();  // total number of observations to simulate
  NumericMatrix y_sim(m, n);
  NumericMatrix S(m, n);
  arma::cube CondVol(m,n,K);
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prep for 'set_vol'
  volatilityVector vol0 = set_vol(y[0]);
  double z;
  for (int t = 1; t <= nb_obs; t++) {
    increment_vol(vol0, y[t - 1]);  // increment all volatilities
  }
  for (int i = 0; i < m; i++) {
    S(i,0) = sampleState(P0_);           // sample initial state
    z = rndgen(S(i,0));
    y_sim(i,0) = z * sqrt(vol0[S(i,0)].h);  // first draw
  }
  volatilityVector vol = vol0;
  for (int i = 0; i < m; i++) {
    for (int s = 0; s < K; s++) {
      CondVol(i, 0, s) = sqrt(vol[s].h);
    }
    for (int t = 1; t < n; t++) {
      S(i,t) = sampleState(P(S(i,t - 1), _));  // sample new state
      z = rndgen(S(i,t));                    // sample new innovation
      increment_vol(vol, y_sim(i,t - 1));    // increment all volatilities
      y_sim(i,t) = z * sqrt(vol[S(i,t)].h);
      for (int s = 0; s < K; s++) {
        CondVol(i, t, s) = sqrt(vol[s].h);
      }
    }  // new draw
    vol = vol0;
  }
  return (List::create(Rcpp::Named("draws") = y_sim, Rcpp::Named("state") = S,  Rcpp::Named("CondVol") = CondVol));
}

inline List MSgarch::f_rnd(const int& n, const NumericVector& theta,
                           const NumericVector& y) {
  // setup
  int nb_obs = y.size();
  NumericVector draw(n);  // draw
  IntegerVector S(n);     // states of the Markov chain
  loadparam(theta);       // load parameters
  double z;
  prep_ineq_vol();  // prep for 'set_vol'
  volatilityVector vol = set_vol(y[0]);
  
  for (int i = 1; i <= nb_obs; i++) {  // loop over observations
    increment_vol(vol, y[i - 1]);      // increment all volatilities
  }
  HamiltonFilter(calc_lndMat(y));
  
  // increment over time
  for (int i = 0; i < n; i++) {
    S[i] = sampleState(PLast);        // sample new state
    z = rndgen(S[i]);                 // sample new innovation
    draw[i] = z * sqrt(vol[S[i]].h);  // new draw
  }
  NumericVector yy(draw.begin(), draw.end());
  NumericVector SS(S.begin(), S.end());
  return (
      Rcpp::List::create(Rcpp::Named("draws") = yy, Rcpp::Named("state") = SS));
}

//------------------------------ Compute loglikelihood matrix
//------------------------------//
inline NumericMatrix MSgarch::calc_lndMat(const NumericVector& y) {
  // set up
  int nb_obs = y.size();
  NumericMatrix lndMat(K, nb_obs - 1);
  
  // initialize
  volatilityVector vol = set_vol(y[0]);
  prep_kernel();
  
  // loop over observations
  for (int t = 1; t < nb_obs; t++) {
    increment_vol(vol, y[t - 1]);               // increment all volatilities
    lndMat(_, t - 1) = calc_kernel(vol, y[t]);  // calc all kernels
  }
  return lndMat;
}

//-------------------------------------  Hamilton filter
//-------------------------------------//
inline double MSgarch::HamiltonFilter(const NumericMatrix& lndMat) {
  int n_step = lndMat.ncol();
  double lnd = 0, min_lnd, delta, sum_tmp;
  NumericVector Pspot, Ppred, lndCol, tmp;
  
  // first step
  Pspot = clone(P0);             // Prob(St | I(t)
  Ppred = matrixProd(Pspot, P);  // one-step-ahead Prob(St | I(t-1))
  lndCol = lndMat(_, 0);
  min_lnd = min(lndCol),
    delta =
      ((min_lnd < LND_MIN) ? LND_MIN - min_lnd : 0);  // handle over/under-flows
  tmp = Ppred *
    exp(lndCol + delta);  // unormalized one-step-ahead Prob(St | I(t))
  
  // remaining steps
  for (int t = 1; t < n_step; t++) {
    sum_tmp = sum(tmp);
    lnd += -delta + log(sum_tmp);  // increment loglikelihood
    Pspot = tmp / sum_tmp;
    Ppred = matrixProd(Pspot, P);
    lndCol = lndMat(_, t);
    min_lnd = min(lndCol),
      delta = ((min_lnd < LND_MIN) ? LND_MIN - min_lnd
                 : 0);  // handle over/under-flows
    tmp = Ppred * exp(lndCol + delta);
  }
  sum_tmp = sum(tmp);
  lnd += -delta + log(sum_tmp);  // increment loglikelihood
  Pspot = tmp / sum_tmp;
  PLast = matrixProd(Pspot, P);
  return lnd;
}

inline List MSgarch::f_get_Pstate(const NumericVector& theta,
                                  const NumericVector& y) {
  // init
  loadparam(theta);  // load parameters
  prep_ineq_vol();   // prepare functions related to volatility
  volatilityVector vol = set_vol(y[0]);   // initialize volatility
  NumericMatrix lndMat = calc_lndMat(y);  // likelihood in each state
  
  int n_step = lndMat.ncol();
  double lnd = 0, min_lnd, delta, sum_tmp;
  NumericVector Pspot, Ppred, lndCol, tmp;
  arma::mat PtmpSpot(n_step + 1, K);
  arma::mat PtmpPred(n_step + 2, K);
  arma::mat PtmpSmooth(n_step + 2, K);
  
  // first step
  Pspot = clone(P0);             // Prob(St | I(t))
  Ppred = matrixProd(Pspot, P);  // one-step-ahead Prob(St | I(t-1))
  for (int i = 0; i < K; i++) {
    PtmpSpot(0, i) = Pspot(i);
  }
  
  for (int i = 0; i < K; i++) {
    PtmpPred(0, i) = Pspot(i);
  }
  for (int i = 0; i < K; i++) {
    PtmpPred(1, i) = Ppred(i);
  }
  lndCol = lndMat(_, 0);
  min_lnd = min(lndCol),
    delta =
      ((min_lnd < LND_MIN) ? LND_MIN - min_lnd : 0);  // handle over/under-flows
  tmp = Ppred *
    exp(lndCol + delta);  // unormalized one-step-ahead Prob(St | I(t))
  
  // remaining steps
  for (int t = 1; t < n_step; t++) {
    sum_tmp = sum(tmp);
    lnd += -delta + log(sum_tmp);  // increment loglikelihood
    Pspot = tmp / sum_tmp;
    Ppred = matrixProd(Pspot, P);
    for (int i = 0; i < K; i++) {
      PtmpSpot(t, i) = Pspot(i);
    }
    for (int i = 0; i < K; i++) {
      PtmpPred(t + 1, i) = Ppred(i);
    }
    lndCol = lndMat(_, t);
    min_lnd = min(lndCol),
      delta = ((min_lnd < LND_MIN) ? LND_MIN - min_lnd
                 : 0);  // handle over/under-flows
    tmp = Ppred * exp(lndCol + delta);
  }
  sum_tmp = sum(tmp);
  lnd += -delta + log(sum_tmp);  // increment loglikelihood
  Pspot = tmp / sum_tmp;
  PLast = matrixProd(Pspot, P);
  
  for (int i = 0; i < K; i++) {
    PtmpSpot(n_step, i) = Pspot(i);
  }
  for (int i = 0; i < K; i++) {
    PtmpPred(n_step + 1, i) = PLast(i);
  }
  for (int i = 0; i < K; i++) {
    PtmpSmooth(n_step + 1, i) = PLast(i);
  }
  arma::mat tmpMat(1, 2);
  arma::mat tmpMat2(2, 1);
  for (int t = n_step; t >= 0; t--) {
    tmpMat = (PtmpSmooth.row(t + 1) / PtmpPred.row(t + 1));
    tmpMat2 = (as<arma::mat>(P) * tmpMat.t());
    PtmpSmooth.row(t) = PtmpSpot.row(t) % tmpMat2.t();
  }
  
  return List::create(
    Rcpp::Named("FiltProb") = PtmpSpot, Rcpp::Named("PredProb") = PtmpPred,
    Rcpp::Named("SmoothProb") = PtmpSmooth, Rcpp::Named("LL") = lndMat);
}

//------------------------------------- Model evaluation
//-------------------------------------//
inline NumericVector MSgarch::eval_model(NumericMatrix& all_thetas,
                                         const NumericVector& y,
                                         const bool& do_prior) {
  // set up
  int nb_thetas = all_thetas.nrow();
  NumericVector lnd(nb_thetas), theta_j(all_thetas.ncol());
  prior pr;
  double tmp;
  
  // loop over each vector of parameters
  for (int j = 0; j < nb_thetas; j++) {
    theta_j = all_thetas(j, _);  // extract parameters
    loadparam(theta_j);          // load parameters
    prep_ineq_vol();
    pr = calc_prior(theta_j);
    if (do_prior == true) {
      lnd[j] = pr.r2 + pr.r3;
    } else {
      lnd[j] = pr.r2;
    }
    pr = calc_prior(theta_j);
    tmp = 0;
    if (pr.r1) tmp += HamiltonFilter(calc_lndMat(y));
    lnd[j] += tmp;
  }
  return lnd;
}

#endif  // MSgarch.h
