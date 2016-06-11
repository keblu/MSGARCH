#ifndef SPEC_H // include guard
#define SPEC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

struct volatility 
{
  double h;        // variance
  double lnh;      // log(h)
  double fh;       // some arbitraty function of "h" (e.g. sqrt(h) as in the Tgarch model)
};

struct prior
{
  bool r1;       // TRUE if the coefficients respect the prior, FALSE if not
  double r2;     // loglikelihood of the coefficients
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Normal 
{
  double lncst;    // constant term for "kernel" 
public:    
<<<<<<< HEAD
  double M1;    
=======
  
>>>>>>> origin/master
  // constructor
  Normal() {
    lncst = - 0.5*log(2*M_PI);
  }
  
  // constructor function called by higher-level classes (e.g. Garch). 
  // here, the argument "name" is passed by reference and modified
  void constructor(std::string& name, int&, NumericVector&, NumericVector&, NumericVector&,
                   CharacterVector&, NumericVector&, NumericVector&) {
    name.append("normal_");
  }
  
  // empty
  void loadparam(const NumericVector&, int&) { }; 
  
  // check prior (always TRUE)
  bool calc_r1() {return 1;}                      
  
<<<<<<< HEAD
  void set_M1() {M1 = sqrt(2 / M_PI);}
  
=======
>>>>>>> origin/master
  // empty
  void prep_kernel() { };                         
  
  // returns loglikelihood of a single observation
  double kernel(const volatility& vol, const double& yi) {      
    return lncst - 0.5*pow(yi, 2)/vol.h - 0.5*vol.lnh; 
  }
  
  // returns PDF evaluated at "x"
  double pdf(const double& x) {return R::dnorm(x, 0, 1, 0);}    
  
  // returns CDF evaluated at "x"
  double cdf(const double& x) {return R::pnorm(x, 0, 1, 1, 0);}     
  
  // applies inverse transform sampling on a Uniform(0,1) draw
  double invsample(const double& u) {return R::qnorm(u, 0, 1, 1, 0);} 
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename underlying>
class Symmetric
{  
  underlying f1;    // symmetric distribution
public:
<<<<<<< HEAD
  double Eabsz;     // := E[|z|] 
  double EzIpos;    // := E[z * I(z>=0)]  = Prob(z >= 0) * E[z | z >= 0]
  double EzIneg;    // := E[z * I(z<0)]   = Prob(z < 0)  * E[z | z < 0]
  double Ez2Ineg;   // := E[z^2 * I(z<0)] = Prob(z < 0)  * E[z^2 | z < 0]
=======
>>>>>>> origin/master
  
  // constructor function called by higher-level classes (e.g. GARCH). 
  // arguments are passed by reference and modified according to the distribution
  void constructor(std::string& name, int& nb_coeffs, NumericVector& coeffs_mean, NumericVector& coeffs_sd,
                   NumericVector& Sigma0, CharacterVector& label, NumericVector& lower, NumericVector& upper) {
    f1.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label, lower, upper);
    name.append("sym");
  }
  
  // nested call (this function should always be called first)  
  void loadparam(const NumericVector& theta, int& Ind) {f1.loadparam(theta, Ind);}
  
  // check prior (nested call)
  bool calc_r1() {return f1.calc_r1();}
  
  // setup for "calc_kernel" (nested call)
  void prep_kernel() {f1.prep_kernel();}
  
  // returns PDF evaluated at "x" (nested call)
  double calc_pdf(const double& x) {return f1.pdf(x);}
  
  // returns CDF evaluated at "x" (nested call)
  double calc_cdf(const double& x) {return f1.cdf(x);}
  
  // returns kernel of a single observation (must call "prep_kernel" first)
  double calc_kernel(const volatility& vol, const double& yi) {
    return f1.kernel(vol, yi);                     // if in A (log density);  // if not  
  }
  
<<<<<<< HEAD
  void prep_moments1() {f1.set_M1();}         // prep-function for moments of order 1
  void set_Eabsz() {Eabsz = f1.M1;}           // = E[|z|] 
  void set_EzIpos() {EzIpos =  0.5 * f1.M1;}  // = E[z * I(z>=0)] = Prob(z >= 0) * E[z | z >= 0]
  void set_EzIneg() {EzIneg = - 0.5 * f1.M1;} // = E[z * I(z<0)] = Prob(z < 0) * E[z | z < 0]
  void set_Ez2Ineg() {Ez2Ineg = 0.5;}         // = E[z^2 * I(z<0)] = Prob(z < 0) * E[z^2 | z < 0]
  
=======
>>>>>>> origin/master
  
  // returns a random vector of length "n"
  NumericVector rndgen(const int& n) {
    NumericVector out(n);
    NumericVector u = runif(n, 0, 1);
    for (int i = 0; i < n; i++) 
      out[i] = f1.invsample(u[i]);
    return out;
  }
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename distribution>
class Garch 
{
  distribution fz;                // distribution of innovations (e.g. Symmetric<Normal>)
  double alpha0, alpha1, beta;    // coefficients
public:
  std::string name;                      // name of the model
  int nb_coeffs;                         // total number of coefficients (including those of "fz")
  int nb_coeffs_model;  
  CharacterVector label;                 // labels of coefficients
  NumericVector coeffs_mean, coeffs_sd;  // means and standard deviations for prior distributions 
  NumericVector Sigma0;                  // diagonal of matrix Sigma0 
  NumericVector lower, upper;            // lower and upper bounds for coefficients
  
  double ineq_lb, ineq_ub;               // lower and upper bounds for inequality constraint
  
  // constructor
  Garch() {
    ineq_lb     = 1e-4;
    ineq_ub     = 0.9999;
    label       = CharacterVector::create("alpha0", "alpha1", "beta" );
    coeffs_mean = NumericVector::create(   0.1,      0.1,      0.8   );
    coeffs_sd   = NumericVector::create(   2,        2,        2     );
    Sigma0      = NumericVector::create(   1,        1,        1     );
    lower       = NumericVector::create(   1e-4,     1e-4,     1e-4  );
    upper       = NumericVector::create(   100,      0.9999,   0.9999);
    nb_coeffs   = label.size();
    nb_coeffs_model = 3;
    name        = "Garch_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label, lower, upper);  
  }
  
  // set all parameters (including those of the distribution).
  // this function should always be called first 
  void loadparam(const NumericVector& theta) {
    alpha0  = theta[0], alpha1 = theta[1], beta = theta[2];
    int Ind = 3;
    fz.loadparam(theta, Ind);
  }
  
  // empty
  void prep_ineq_vol() { };
  
  // inequality constraint function
  double ineq_func() {return alpha1 + beta;}
  
  // check prior 
  bool calc_r1() {
    return fz.calc_r1() 
    && alpha0 > lower[0] && alpha1 > lower[1] && beta > lower[2]
    && (ineq_func() < ineq_ub);
  }
  
  // initialize volatility 
  volatility set_vol(const double& y0) {
    volatility out;
    out.h   = alpha0 / (1 - alpha1 - beta); 
    out.lnh = log(out.h);
    return out;
  }
  
  // increment volatility 
  void increment_vol(volatility& vol, const double& yim1) { 
    vol.h   = alpha0 + alpha1 * pow(yim1, 2) + beta * vol.h; 
    vol.lnh = log(vol.h);
  }
  
  // some nested-call functions
  void prep_kernel() {fz.prep_kernel();}
  NumericVector rndgen(const int& n) {return fz.rndgen(n);}
  double calc_pdf(const double& x) {return fz.calc_pdf(x);}
  double calc_cdf(const double& x) {return fz.calc_cdf(x);}
  double calc_kernel(const volatility& vol, const double& yi) {
    return fz.calc_kernel(vol, yi);
  }
};

<<<<<<< HEAD
template <typename distribution>
class Gjr 
{
  distribution fz;                      // distribution of innovations (e.g. Symmetric<Normal>)
  double alpha0, alpha1, alpha2, beta;  // coefficients
public:
  std::string name;                      // name of the model
  int nb_coeffs;                         // total number of coefficients (including those of "fz")
  int nb_coeffs_model;  
  CharacterVector label;                 // labels of coefficients
  NumericVector coeffs_mean, coeffs_sd;  // means and standard deviations for prior distributions 
  NumericVector Sigma0;                  // diagonal of matrix Sigma0 
  NumericVector lower, upper;            // lower and upper bounds for coefficients
  double ineq_lb, ineq_ub;               // lower and upper bounds for inequality constraint
  
  // constructor
  Gjr() {
    ineq_lb     = 1e-4;
    ineq_ub     = 0.9999;
    label       = CharacterVector::create("alpha0", "alpha1", "alpha2", "beta" );
    coeffs_mean = NumericVector::create(   0.1,      0.05,     0.1,      0.8   );
    coeffs_sd   = NumericVector::create(   2,        2,        2,        2     );
    Sigma0      = NumericVector::create(   1,        1,        1,        1     );
    lower       = NumericVector::create(   1e-4,     1e-6,     1e-4,     1e-4  );
    upper       = NumericVector::create(   100,      0.9999,   10,      0.9999);
    nb_coeffs   = label.size();
    nb_coeffs_model = 4;
    name        = "Gjr_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label, lower, upper); 
  }
  
  // set all parameters (including those of the distribution).
  // this function should always be called first 
  void loadparam(const NumericVector& theta) {
    alpha0  = theta[0], alpha1 = theta[1], alpha2 = theta[2], beta = theta[3];
    int Ind = 4;
    fz.loadparam(theta, Ind);
  }
  
  void prep_ineq_vol() {
    fz.set_Ez2Ineg(); 
  }
  
  // inequality constraint function   
  double ineq_func() {return alpha1 + fz.Ez2Ineg * alpha2 + beta;}
  
  // computes r1 
  bool calc_r1() {
    return    fz.calc_r1()  
    && alpha0 > lower[0] && alpha1 > lower[1] && alpha2 > lower[2] && beta > lower[3] 
    && (ineq_func() < ineq_ub);
  }
  
  // initialize volatility to its undonditional expected value  
  volatility set_vol(const double& y0) {
    volatility out;
    out.h   = alpha0 / (1 - alpha1 - fz.Ez2Ineg * alpha2 - beta); 
    out.lnh = log(out.h);
    return out;
  }
  
  // increment volatility
  void increment_vol(volatility& vol, const double& yim1) { 
    vol.h   = alpha0 + alpha1 * pow(yim1, 2) + beta * vol.h + ((yim1 < 0)? alpha2*pow(yim1, 2) : 0);  
    vol.lnh = log(vol.h);
  }
  
  // some nested-call functions
  void prep_kernel() {fz.prep_kernel();}
  NumericVector rndgen(const int& n) {return fz.rndgen(n);}
  double calc_pdf(const double& x) {return fz.calc_pdf(x);}
  double calc_cdf(const double& x) {return fz.calc_cdf(x);}
  double calc_kernel(const volatility& vol, const double& yi) {
    return fz.calc_kernel(vol, yi);
  }
};

template <typename distribution>
class Egarch 
{
  distribution fz;                      // distribution of innovations (e.g. Symmetric<Normal>)
  double alpha0, alpha1, alpha2, beta;  // coefficients
public:
  std::string name;                      // name of the model
  int nb_coeffs;                         // total number of coefficients (including those of "fz")
  int nb_coeffs_model;  
  CharacterVector label;                 // labels of coefficients
  NumericVector coeffs_mean, coeffs_sd;  // means and standard deviations for prior distributions 
  NumericVector Sigma0;                  // diagonal of matrix Sigma0 
  NumericVector lower, upper;            // lower and upper bounds for coefficients
  double ineq_lb, ineq_ub;               // lower and upper bounds for inequality constraint
  
  // constructor
  Egarch() {
    ineq_lb     = -0.9999; 
    ineq_ub     = 0.9999; 
    label       = CharacterVector::create("alpha0", "alpha1", "alpha2", "beta" );
    coeffs_mean = NumericVector::create(   0.01,     0.2,     -0.1,      0.8   ); 
    coeffs_sd   = NumericVector::create(   2,        2,        2,        2     );
    Sigma0      = NumericVector::create(   1,        1,        1,        1     );
    lower       = NumericVector::create(  -50,      -5,       -5,       -0.9999);
    upper       = NumericVector::create(   50,       5,        5,        0.9999);  
    nb_coeffs   = label.size();
    nb_coeffs_model = 4;
    name        = "Egarch_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label, lower, upper); 
  }
  
  // set all parameters (including those of the distribution).
  // this function should always be called first 
  void loadparam(const NumericVector& theta) {
    alpha0  = theta[0], alpha1 = theta[1], alpha2 = theta[2], beta = theta[3];
    int Ind = 4;
    fz.loadparam(theta, Ind);
    fz.prep_moments1();
    fz.set_Eabsz();
  }
  
  // empty
  void prep_ineq_vol() { };
  
  // inequality constraint function   
  double ineq_func() {return beta;}
  
  // computes r1 
  bool calc_r1() {
    double ineq_theta = ineq_func();
    return fz.calc_r1() && (ineq_theta > ineq_lb && ineq_theta < ineq_ub);
  }
  
  // initialize volatility to its undonditional expected value  
  volatility set_vol(const double& y0) {
    volatility out;
    out.lnh = alpha0 / (1 - beta);
    out.h   = exp(out.lnh);
    return out;
  }
  
  // increment volatility
  void increment_vol(volatility& vol, const double& yim1) { 
    double z = yim1 / sqrt(vol.h);
    vol.lnh = alpha0 + alpha1 * (fabs(z) - fz.Eabsz) + alpha2 * z + beta * vol.lnh ;
    vol.h   = exp(vol.lnh);
  }
  
  // some nested-call functions
  void prep_kernel() {fz.prep_kernel();}
  NumericVector rndgen(const int& n) {return fz.rndgen(n);}
  double calc_pdf(const double& x) {return fz.calc_pdf(x);}
  double calc_cdf(const double& x) {return fz.calc_cdf(x);}
  double calc_kernel(const volatility& vol, const double& yi) {
    return fz.calc_kernel(vol, yi);
  }
};

template <typename distribution>
class Tgarch 
{
  distribution fz;                      // distribution of innovations (e.g. Symmetric<Normal>)
  double alpha0, alpha1, alpha2, beta;  // coefficients
public:
  std::string name;                      // name of the model
  int nb_coeffs;                         // total number of coefficients (including those of "fz")
  int nb_coeffs_model; 
  CharacterVector label;                 // labels of coefficients
  NumericVector coeffs_mean, coeffs_sd;  // means and standard deviations for prior distributions 
  NumericVector Sigma0;                  // diagonal of matrix Sigma0 
  NumericVector lower, upper;            // lower and upper bounds for coefficients
  double ineq_lb, ineq_ub;               // lower and upper bounds for inequality constraint
  
  // constructor
  Tgarch() {
    ineq_lb     = 1e-4;
    ineq_ub     = 0.9999;
    label       = CharacterVector::create("alpha0", "alpha1", "alpha2", "beta" );
    coeffs_mean = NumericVector::create(   0.125,     0.05,     0.1,     0.8   );
    coeffs_sd   = NumericVector::create(   2,        2,        2,        2     );
    Sigma0      = NumericVector::create(   1,        1,        1,        1     );
    lower       = NumericVector::create(   1e-4,     1e-6,     1e-4,     1e-4  );
    upper       = NumericVector::create(   100,      10,      10,        10    );
    nb_coeffs   = label.size();
    nb_coeffs_model = 4;
    name        = "Tgarch_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label, lower, upper);
  }
  
  // set all parameters (including those of the distribution).
  // this function should always be called first 
  void loadparam(const NumericVector& theta) {
    alpha0  = theta[0], alpha1 = theta[1], alpha2 = theta[2], beta = theta[3];
    int Ind = 4;
    fz.loadparam(theta, Ind);
  }
  
  void prep_ineq_vol() {
    fz.prep_moments1();
    fz.set_EzIneg();
    fz.set_Ez2Ineg();
  };
  
  // inequality constraint function   
  double ineq_func() {
    return pow(alpha1, 2) + pow(beta, 2) - 2 * (alpha1 + alpha2) * beta * fz.EzIneg 
    - (pow(alpha1,2) - pow(alpha2,2)) * fz.Ez2Ineg;
  }
  
  // computes r1 
  bool calc_r1() {
    return fz.calc_r1()  
    && alpha0 > lower[0] && alpha1 > lower[1] && alpha2 > lower[2] && beta > lower[3]    
    && (ineq_func() < ineq_ub);
  }
  
  // initialize volatility to its undonditional expected value 
  volatility set_vol(const double& y0) {
    volatility out;
    out.fh  = alpha0 / (1 + (alpha1 + alpha2) * fz.EzIneg - beta);  
    out.h   = pow(out.fh, 2);
    out.lnh = log(out.h);
    return out;
  }
  
  // increment volatility
  void increment_vol(volatility& vol, const double& yim1) { 
    vol.fh  = alpha0 + beta * vol.fh + ((yim1 >= 0)? alpha1 : -alpha2) * yim1;
    vol.h   = pow(vol.fh, 2);
    vol.lnh = log(vol.h);
  }
  
  // some nested-call functions
  void prep_kernel() {fz.prep_kernel();}
  NumericVector rndgen(const int& n) {return fz.rndgen(n);}
  double calc_pdf(const double& x) {return fz.calc_pdf(x);}
  double calc_cdf(const double& x) {return fz.calc_cdf(x);}
  double calc_kernel(const volatility& vol, const double& yi) {
    return fz.calc_kernel(vol, yi);
  }
};


=======
>>>>>>> origin/master
#endif // spec.h