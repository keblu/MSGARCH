#ifndef SingleRegime_H // include guard
#define SingleRegime_H

#include "spec.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

// The class "SingleRegime" below is templated in terms of the model to use (e.g. Garch<Symmetric<Normal> >)

/// Base virtual class with polymorphic abilities
class Base {
public:

    // returns data members
    virtual std::string spec_name() = 0;
    virtual NumericVector spec_theta0() = 0;
    virtual NumericVector spec_Sigma0() = 0;
    virtual CharacterVector spec_label() = 0;
    virtual NumericVector spec_lower() = 0;
    virtual NumericVector spec_upper() = 0;
    virtual double spec_ineq_lb() = 0;
    virtual double spec_ineq_ub() = 0;
    virtual int spec_nb_coeffs() = 0;
    virtual int spec_nb_coeffs_model() = 0;
    // emulate function members
    virtual void spec_loadparam(const NumericVector&) = 0;
    virtual prior spec_calc_prior(const NumericVector&) = 0;
    virtual void spec_prep_ineq_vol() = 0;
    virtual double spec_ineq_func() = 0;
    virtual bool spec_calc_r1() = 0;
    virtual volatility spec_set_vol(const double&) = 0;
    virtual void spec_increment_vol(volatility&, const double&) = 0;
    virtual void spec_prep_kernel() = 0 ;
    virtual NumericVector spec_rndgen(const int&) = 0;
    virtual double spec_calc_pdf(const double&) = 0 ;
    virtual double spec_calc_cdf(const double&) = 0;
    virtual double spec_calc_kernel(const volatility&, const double&)  = 0;

    virtual ~Base() = 0;
};
inline Base::~Base() {}

//------------------------ Master class ------------------------//
template <typename Model>
class SingleRegime : public Base
{
    Model spec;
public:
    std::string name;
    NumericVector theta0;
    NumericVector Sigma0;
    CharacterVector label;
    NumericVector lower;
    NumericVector upper;
    double ineq_lb;
    double ineq_ub;
    IntegerVector NbParams;
    IntegerVector NbParamsModel;
    // constructor
    SingleRegime() {
        name    = spec.name;
        theta0  = spec.coeffs_mean;
        Sigma0  = spec.Sigma0;
        label   = spec.label;
        lower   = spec.lower;
        upper   = spec.upper;
        ineq_lb = spec.ineq_lb;
        ineq_ub = spec.ineq_ub;
        NbParams.push_back(spec.nb_coeffs);
        NbParamsModel.push_back(spec.nb_coeffs_model);
    }

    // functions accessible from R
    double ineq_func(const NumericVector& theta) {
        spec.loadparam(theta);
        spec.prep_ineq_vol();
        return spec.ineq_func();
    }
    prior calc_prior(const NumericVector&);
    NumericVector f_sim(const int&, const NumericVector&, const int&);
    NumericVector f_pdf(const NumericVector&, const NumericVector&, const NumericVector&, const bool&);
    NumericVector f_pdf_its( const NumericVector&, const NumericVector&, const bool&);
    NumericVector f_cdf(const NumericVector&, const NumericVector&, const NumericVector&, const bool&);
    NumericVector f_cdf_its( const NumericVector&, const NumericVector&, const bool&);
    NumericVector f_rnd(const int&, const NumericVector&, const NumericVector&);
    NumericVector f_unc_vol(NumericMatrix&, const NumericVector&);
    NumericMatrix calc_ht(NumericMatrix&, const NumericVector&);
    NumericVector eval_model(NumericMatrix&, const NumericVector&);

    // Handles to 'spec' data members
    std::string spec_name() {
        return spec.name;
    }
    NumericVector spec_theta0() {
        return spec.coeffs_mean;
    }
    NumericVector spec_Sigma0() {
        return spec.Sigma0;
    }
    CharacterVector spec_label() {
        return spec.label;
    }
    NumericVector spec_lower() {
        return spec.lower;
    }
    NumericVector spec_upper() {
        return spec.upper;
    }
    double spec_ineq_lb() {
        return spec.ineq_lb;
    }
    double spec_ineq_ub() {
        return spec.ineq_ub;
    }
    int spec_nb_coeffs() {
        return spec.nb_coeffs;
    }
    int spec_nb_coeffs_model() {
      return spec.nb_coeffs_model;
    }

    // Handles to 'spec' funtion members
    void spec_loadparam(const NumericVector& theta) {
        spec.loadparam(theta);
    }
    prior spec_calc_prior(const NumericVector& theta) {
        return calc_prior(theta);
    }
    void spec_prep_ineq_vol() {
        spec.prep_ineq_vol();
    }
    double spec_ineq_func() {
        return spec.ineq_func();
    }
    bool spec_calc_r1() {
        return spec.calc_r1();
    }
    volatility spec_set_vol(const double& y0) {
        return spec.set_vol(y0);
    }
    
    void spec_increment_vol(volatility& vol, const double& yim1) {
        spec.increment_vol(vol, yim1);
    }
    void spec_prep_kernel() {
        spec.prep_kernel();
    }
    NumericVector spec_rndgen(const int& n) {
        return spec.rndgen(n);
    }
    double spec_calc_pdf(const double& x) {
        return spec.calc_pdf(x);
    }
    double spec_calc_cdf(const double& x) {
        return spec.calc_cdf(x);
    }
    double spec_calc_kernel(const volatility& vol, const double& yi) {
        return spec.calc_kernel(vol, yi);
    }
};

//---------------------- Prior calculation ----------------------//
template <typename Model>
prior SingleRegime<Model>::calc_prior(const NumericVector& theta) {
    bool   r1 = spec.calc_r1();
    double r2 = -1e10;
    if (r1) {
        r2 = 0;
        for (int i = 0; i < spec.nb_coeffs; i++) {
            r2 += R::dnorm(theta[i], spec.coeffs_mean[i], spec.coeffs_sd[i], 1);
        }
    }
    prior out;
    out.r1 = r1;
    out.r2 = r2;
    return out;
}

//---------------------- Model simulation ----------------------//
template <typename Model>
NumericVector SingleRegime<Model>::f_sim(const int& n, const NumericVector& theta, const int& burnin) {
    int ntot = n + burnin;                 // total number of observations to simulate
    spec.loadparam(theta);                 // load parameters
    NumericVector z = spec.rndgen(ntot);   // generate vector of random innovations
    spec.prep_ineq_vol();                  // prepare functions related to volatility
    volatility vol  = spec.set_vol(z[0]);  // initialize volatility
    NumericVector y(ntot);
    y[0] = z[0] * sqrt(vol.h);
    for (int t = 1; t < ntot; t++) {
        spec.increment_vol(vol, y[t-1]);
        y[t] = z[t] * sqrt(vol.h);
    }
    NumericVector yy(y.begin() + burnin, y.end());
    return yy;
}

//---------------------- Calculates PDF ----------------------//
template <typename Model>
NumericVector SingleRegime<Model>::f_pdf(const NumericVector& x, const NumericVector& theta,
                                     const NumericVector& y, const bool& is_log) {
    // computes volatility
    spec.loadparam(theta);                  // load parameters
    spec.prep_ineq_vol();                   // prepare functions related to volatility
    volatility vol  = spec.set_vol(y[0]);   // initialize volatility
    int ny = y.size();
    for (int t = 0; t < ny; t++)
        spec.increment_vol(vol, y[t]);
    double sig = sqrt(vol.h);

    // computes PDF
    int nx = x.size();
    double tmp;
    NumericVector out(nx);
    for (int i = 0; i < nx; i++) {
        tmp    = spec.calc_pdf(x[i]/sig) / sig; //must divide by sig because of variable transformation
        out[i] = ((is_log)?  log(tmp) : tmp);
    }
    return out;
}

template <typename Model>
NumericVector SingleRegime<Model>::f_pdf_its(const NumericVector& theta,
                                         const NumericVector& y, const bool& is_log) {
  // computes volatility
  double sig;
  spec.loadparam(theta);                  // load parameters
  spec.prep_ineq_vol();                   // prepare functions related to volatility
  int ny = y.size();
  
  volatility vol  = spec.set_vol(y[0]);   // initialize volatility
  NumericVector out(ny - 1);
  for (int i = 1; i < ny; i++) { 
    spec.increment_vol(vol, y[i-1]);
      sig = sqrt(vol.h);
      out[i-1] = spec.calc_pdf(y[i]/sig) / sig; //
      out[i-1] = ((is_log)?  log(out[i-1]) : out[i-1]);
    }
  
  return out;
}

//---------------------- Calculates CDF ----------------------//
template <typename Model>
NumericVector SingleRegime<Model>::f_cdf(const NumericVector& x, const NumericVector& theta,
                                     const NumericVector& y, const bool& is_log) {
    // computes volatility
    spec.loadparam(theta);                  // load parameters
    spec.prep_ineq_vol();                   // prepare functions related to volatility
    volatility vol  = spec.set_vol(y[0]);   // initialize volatility
    int ny = y.size();
    for (int t = 0; t < ny; t++)
        spec.increment_vol(vol, y[t]);
    double sig = sqrt(vol.h);

    // computes CDF
    int nx = x.size();
    double tmp;
    NumericVector out(nx);
    for (int i = 0; i < nx; i++) {
        tmp    = spec.calc_cdf(x[i]/sig);
        out[i] = ((is_log)?  log(tmp) : tmp);
    }
    return out;
}

template <typename Model>
NumericVector SingleRegime<Model>::f_cdf_its(const NumericVector& theta,
                                             const NumericVector& y, const bool& is_log) {
  // computes volatility
  double sig;
  spec.loadparam(theta);                  // load parameters
  spec.prep_ineq_vol();                   // prepare functions related to volatility
  int ny = y.size();
  
  volatility vol  = spec.set_vol(y[0]);   // initialize volatility
  NumericVector out(ny-1);
  
  for (int i = 1; i < ny; i++) { 
    spec.increment_vol(vol, y[i-1]);
    sig = sqrt(vol.h);
    out[i-1] = spec.calc_cdf(y[i]/sig); //
    out[i-1] = ((is_log)?  log(out[i-1]) : out[i-1]);
  }
  
  return out;
}

//---------------------- Generates 1-day-ahead simulations ----------------------//
template <typename Model>
NumericVector SingleRegime<Model>::f_rnd(const int& n, const NumericVector& theta, const NumericVector& y) {

    spec.loadparam(theta);                  // load parameters
    spec.prep_ineq_vol();                   // prepare functions related to volatility
    volatility vol  = spec.set_vol(y[0]);   // initialize volatility
    int ny = y.size();
    for (int t = 1; t <= ny; t++)
        spec.increment_vol(vol, y[t-1]);
    return sqrt(vol.h) * spec.rndgen(n);
}

//---------------------- Conditional variance calculation ----------------------//
template <typename Model>
NumericMatrix SingleRegime<Model>::calc_ht(NumericMatrix& all_thetas, const NumericVector& y) {
    int nb_obs    = y.size();
    int nb_thetas = all_thetas.nrow();
    volatility vol;
    NumericVector theta_j;
    NumericMatrix ht(nb_obs+1, nb_thetas);
    for (int j = 0; j < nb_thetas; j++) {    // loop over vectors of parameters
        theta_j = all_thetas(j,_);
        spec.loadparam(theta_j);
        spec.prep_ineq_vol();
        vol     = spec.set_vol(y[0]);          // initialize volatility
        ht(0,j) = vol.h;
        for (int i = 1; i <= nb_obs; i++) {    // loop over observations
            spec.increment_vol(vol, y[i-1]);        // increment volatility
            ht(i,j) = vol.h;
        }
    }
    return ht;
}

//---------------------- Conditional variance calculation ----------------------//
template <typename Model>
NumericVector SingleRegime<Model>::f_unc_vol(NumericMatrix& all_thetas, const NumericVector& y) {
  
  int nb_thetas = all_thetas.nrow();
  volatility vol;
  NumericVector theta_j;
  NumericVector ht(nb_thetas);
  for (int j = 0; j < nb_thetas; j++) {    // loop over vectors of parameters
    theta_j = all_thetas(j,_);
    spec.loadparam(theta_j);
    spec.prep_ineq_vol();
    vol     = spec.set_vol(y[0]);          // initialize volatility
    ht(j) = vol.h;
    }
  return ht;
}

//---------------------- Model evaluation ----------------------//
template <typename Model>
NumericVector SingleRegime<Model>::eval_model(NumericMatrix& all_thetas, const NumericVector& y) {
    double tmp;
    int nb_obs    = y.size();
    int nb_thetas = all_thetas.nrow();
    prior pr;
    volatility vol;
    NumericVector lnd(nb_thetas);
    NumericVector theta_j;
    for (int j = 0; j < nb_thetas; j++) {           // loop over vectors of parameters
        theta_j = all_thetas(j,_);
        spec.loadparam(theta_j);
        spec.prep_ineq_vol();
        pr      = calc_prior(theta_j);
        lnd[j]  = pr.r2;
        if (pr.r1) {                                         // if prior satisfied
            tmp = 0;
            vol = spec.set_vol(y[0]);                          // initialize volatility
            spec.prep_kernel();
            for (int i = 1; i < nb_obs; i++) {          // loop over observations
                spec.increment_vol(vol, y[i-1]);                 // increment volatility
                tmp += spec.calc_kernel(vol, y[i]);      // increment kernel
            }
            lnd[j]  += tmp;
        }
    }
    return lnd;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//==================================== CLASS DEFINITIONS ========================================//

typedef SingleRegime <sGARCH <Symmetric <Normal>  > >   sGARCH_norm_sym;
typedef SingleRegime <gjrGARCH   <Symmetric <Normal>  > >   gjrGARCH_norm_sym;
typedef SingleRegime <eGARCH <Symmetric <Normal>  > >   eGARCH_norm_sym;
typedef SingleRegime <tGARCH <Symmetric <Normal>  > >   tGARCH_norm_sym;
typedef SingleRegime <GAS    <Symmetric <Normal>  > >   GAS_norm_sym;

typedef SingleRegime <sGARCH <Skewed <Normal>  > >   sGARCH_norm_skew;
typedef SingleRegime <gjrGARCH   <Skewed <Normal>  > >   gjrGARCH_norm_skew;
typedef SingleRegime <eGARCH <Skewed <Normal>  > >   eGARCH_norm_skew;
typedef SingleRegime <tGARCH <Skewed <Normal>  > >   tGARCH_norm_skew;
typedef SingleRegime <GAS    <Skewed <Normal>  > >   GAS_norm_skew;

typedef SingleRegime <sGARCH <Symmetric <Student>  > >   sGARCH_std_sym;
typedef SingleRegime <gjrGARCH   <Symmetric <Student>  > >   gjrGARCH_std_sym;
typedef SingleRegime <eGARCH <Symmetric <Student>  > >   eGARCH_std_sym;
typedef SingleRegime <tGARCH <Symmetric <Student>  > >   tGARCH_std_sym;
typedef SingleRegime <GAS    <Symmetric <Student>  > >   GAS_std_sym;

typedef SingleRegime <sGARCH <Skewed <Student>  > >   sGARCH_std_skew;
typedef SingleRegime <gjrGARCH   <Skewed <Student>  > >   gjrGARCH_std_skew;
typedef SingleRegime <eGARCH <Skewed <Student>  > >   eGARCH_std_skew;
typedef SingleRegime <tGARCH <Skewed <Student>  > >   tGARCH_std_skew;
typedef SingleRegime <GAS    <Skewed <Student>  > >   GAS_std_skew;

typedef SingleRegime <sGARCH <Symmetric <Ged>  > >   sGARCH_ged_sym;
typedef SingleRegime <gjrGARCH   <Symmetric <Ged>  > >   gjrGARCH_ged_sym;
typedef SingleRegime <eGARCH <Symmetric <Ged>  > >   eGARCH_ged_sym;
typedef SingleRegime <tGARCH <Symmetric <Ged>  > >   tGARCH_ged_sym;
typedef SingleRegime <GAS    <Symmetric <Ged>  > >   GAS_ged_sym;

typedef SingleRegime <sGARCH <Skewed <Ged>  > >   sGARCH_ged_skew;
typedef SingleRegime <gjrGARCH   <Skewed <Ged>  > >   gjrGARCH_ged_skew;
typedef SingleRegime <eGARCH <Skewed <Ged>  > >   eGARCH_ged_skew;
typedef SingleRegime <tGARCH <Skewed <Ged>  > >   tGARCH_ged_skew;
typedef SingleRegime <GAS    <Skewed <Ged>  > >   GAS_ged_skew;

#endif // SingleRegime.h
