#ifndef SPEC_H // include guard
#define SPEC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

typedef std::pair<double, double> pair;

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

// signum function
inline double signum(const double& x) {return (0 < x) - (x < 0);}

// check if x is infinite OR a NaN
inline bool IsInfNan(const double& x) {return traits::is_infinite<REALSXP>(x) || ISNAN(x);}

// Auxiliary function for adaptive Simpson's Rule: https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
template <typename C>  
double adaptiveSimpsonsAux(double (C::*f)(const double&), C* pf, double a, double b, 
                           double epsilon, double S, double fa, double fb, double fc, int bottom) {                 
  double c  = (a + b)/2, h = b - a;                                                                  
  double d  = (a + c)/2, e = (c + b)/2;                                                              
  double fd = (pf->*f)(d), fe = (pf->*f)(e);                                                                      
  double Sleft  = (h/12) * (fa + 4*fd + fc);                                                           
  double Sright = (h/12) * (fc + 4*fe + fb);                                                          
  double S2 = Sleft + Sright; 
  if (bottom <= 100 || fabs(S2 - S)/(1e-10 + fabs(S2)) <= epsilon || IsInfNan(S2))                                                     
    return S2 + (S2 - S)/8;                                                                        
  return   adaptiveSimpsonsAux(f, pf, a, c, epsilon, Sleft,  fa, fc, fd, bottom-1)                     
    + adaptiveSimpsonsAux(f, pf, c, b, epsilon, Sright, fc, fb, fe, bottom-1);                     
}

// Adaptive Simpson's Rule in integrate in interval [a, b] https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
template <typename C>  
double adaptiveSimpsons(double (C::*f)(const double&), C* pf, double a, double b, 
                        double epsilon, int maxRecursion) {   // error tolerance and recursion cap
  double c  = (a + b)/2, h = b - a;                                                                  
  double fa = (pf->*f)(a), fb = (pf->*f)(b), fc = (pf->*f)(c);                                                           
  double S  = (h/6) * (fa + 4*fc + fb);                                                                
  return adaptiveSimpsonsAux(f, pf, a, b, epsilon, S, fa, fb, fc, maxRecursion);                   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Normal 
{
  double lncst;    // constant term for "kernel" 
public:    
  double M1;    
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
  
  void set_M1() {M1 = sqrt(2 / M_PI);}
  
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
  
  // returns the ratio involved in the GAS(1,1) model
  double Sratio(const double& x) {return -x;}     
  
};

//---------------------- Student-t distribution ----------------------//
class Student
{
  double nu;          // degrees of freedom
  double nu_lb;       // lower bound for "nu"
  double lncst;       // constant term in "kernel"
  double cst;         // constant term in "PDF"
  double P;           // factor to standardize R's definition of Student-t distribution
public:
  double M1;                   // E[|z|] 
  
  // constructor
  Student() {
    nu_lb = 2.1;
  }
  
  // constructor function called by higher-level classes (e.g. Garch). 
  // arguments are passed by reference and are modified to include "nu" 
  void constructor(std::string& name, int& nb_coeffs, NumericVector& coeffs_mean, NumericVector& coeffs_sd,
                   NumericVector& Sigma0, CharacterVector& label, NumericVector& lower, NumericVector& upper) {
    name.append("student_");
    nb_coeffs++;
    label.push_back("nu");                              // nu   
    coeffs_mean.push_back(10), coeffs_sd.push_back(10); // mean and standard deviation of prior distribution 
    Sigma0.push_back(10);                               // Sigma0 
    lower.push_back(nu_lb), upper.push_back(100);       // lower and upper bounds
  }
  
  // set "nu" (this function should always be called first)  
  void loadparam(const NumericVector& theta, int& Ind) { 
    nu  = theta[Ind];
    P   = sqrt(nu / (nu-2));
    cst = P * exp(lgamma(0.5*(nu+1)) - lgamma(0.5*nu)) / sqrt(nu * M_PI);
    Ind++;
  }
  
  // check prior
  bool calc_r1() {return nu > nu_lb;}   
  
  // set M1 := E[|z|]
  void set_M1() {                     
    M1 = sqrt((nu - 2) / M_PI) * exp(lgamma(0.5 * (nu - 1)) - lgamma(0.5 * nu));
  } 
  
  // set constant term of "kernel" 
  void prep_kernel() {
    lncst = lgamma(0.5 * (nu + 1)) - lgamma(0.5 * nu) - 0.5 * log(M_PI) + 0.5 * nu * log(nu - 2);
  }
  
  // returns loglikelihood of a single observation (must call "prep_kernel" first)
  double kernel(const volatility& vol, const double& yi) {
    return lncst + 0.5*nu*vol.lnh - 0.5*(nu+1)*log(vol.h*(nu-2) + pow(yi, 2));
  }
  
  // returns PDF evaluated at "x"
  double pdf(const double& x) {return cst * pow(1 + pow(P*x, 2)/nu, -0.5*(nu+1));} 
  
  // returns CDF evaluated at "x"
  double cdf(const double& x) {return R::pt(x * P, nu, 1, 0);} 
  
  
  // applies inverse transform sampling on a Uniform(0,1) draw
  double invsample(const double& u) {return R::qt(u, nu, 1, 0) / P;} 
  
  // returns the ratio involved in the GAS(1,1) model
  double Sratio(const double& x) {return - (nu + 1)*x / (nu - 2 + pow(x, 2));}
  
};

//---------------------- Generalized error distribution ----------------------//
class Ged
{
  double nu;        // shape parameter
  double nu_lb;     // lower bound for "nu"
  double lncst;     // constant term in "kernel"
  double cst;       // constant term in "PDF"
  double lambda;    // lambda
public:
  double M1;                   // E[|z|] 
  
  // constructor
  Ged() {
    nu_lb = 0.7;
  }
  
  // constructor function called by higher-level classes (e.g. Garch). 
  // arguments are passed by reference and are modified to include "nu" 
  void constructor(std::string& name, int& nb_coeffs, NumericVector& coeffs_mean, NumericVector& coeffs_sd,
                   NumericVector& Sigma0, CharacterVector& label, NumericVector& lower, NumericVector& upper) {
    name.append("ged_");
    nb_coeffs++;
    label.push_back("nu");                                // nu
    coeffs_mean.push_back(2), coeffs_sd.push_back(1);     // mean and standard deviation of prior distribution 
    Sigma0.push_back(10);                                 // Sigma0
    lower.push_back(nu_lb), upper.push_back(20);          // lower and upper bounds
  }
  
  // set "nu" (this function should always be called first)     
  void loadparam(const NumericVector& theta, int& Ind) {
    nu     = theta[Ind];
    lambda = sqrt(pow(2, -2/nu) * exp(lgamma(1/nu)-lgamma(3/nu)));
    cst    = nu / (lambda * pow(2, 1 + 1/nu) *exp(lgamma(1/nu)));
    Ind++;
  }
  
  // check prior 
  bool calc_r1() {return nu > nu_lb;}
  
  // set M1 := E[|z|]
  void set_M1() {M1 = 0.5 * lambda * pow(8, 1/nu) * exp(lgamma(1/nu + 0.5)) / sqrt(M_PI);} 
  
  // set constant term of "kernel" 
  void prep_kernel() {lncst = log(cst);}
  
  // returns loglikelihood of a single observation (must call "prep_kernel" first)
  double kernel(const volatility& vol, const double& yi) {
    return lncst - 0.5 * vol.lnh - 0.5 * pow(fabs(yi / (sqrt(vol.h) * lambda)), nu);
  }
  
  // returns PDF evaluated at "x"
  double pdf(const double& x) {return cst * exp(-0.5 * pow(fabs(x / lambda), nu));} 
  
  // returns CDF evaluated at "x"
  double cdf(const double& x) {
    return ((x < 0)?  0.5 * (1 - R::pgamma(0.5 * pow(-x / lambda, nu), 1/nu, 1, 1, 0))
              :  0.5 * (1 + R::pgamma(0.5 * pow( x / lambda, nu), 1/nu, 1, 1, 0)));
  } 
  
  
  // applies inverse transform sampling on a Uniform(0,1) draw
  double invsample(const double& u) {
    return ((u < 0.5)?  -lambda * pow(2 * R::qgamma(1-2*u, 1/nu, 1, 1, 0), 1/nu)
              :   lambda * pow(2 * R::qgamma(2*u-1, 1/nu, 1, 1, 0), 1/nu) );
  } 
  
  double Sratio(const double& x) {
    return ((x == 0)?  0 
              :  - 0.5 * nu * signum(x) * exp((nu - 1) * log(fabs(x)) - nu * log(lambda)));
  }
  
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename underlying>
class Symmetric
{  
  underlying f1;    // symmetric distribution
public:
  double Eabsz;     // := E[|z|] 
  double EzIpos;    // := E[z * I(z>=0)]  = Prob(z >= 0) * E[z | z >= 0]
  double EzIneg;    // := E[z * I(z<0)]   = Prob(z < 0)  * E[z | z < 0]
  double Ez2Ineg;   // := E[z^2 * I(z<0)] = Prob(z < 0)  * E[z^2 | z < 0]
  
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
  
  void prep_moments1() {f1.set_M1();}         // prep-function for moments of order 1
  void set_Eabsz() {Eabsz = f1.M1;}           // = E[|z|] 
  void set_EzIpos() {EzIpos =  0.5 * f1.M1;}  // = E[z * I(z>=0)] = Prob(z >= 0) * E[z | z >= 0]
  void set_EzIneg() {EzIneg = - 0.5 * f1.M1;} // = E[z * I(z<0)] = Prob(z < 0) * E[z | z < 0]
  void set_Ez2Ineg() {Ez2Ineg = 0.5;}         // = E[z^2 * I(z<0)] = Prob(z < 0) * E[z^2 | z < 0]
  
  
  // returns a random vector of length "n"
  NumericVector rndgen(const int& n) {
    NumericVector out(n);
    NumericVector u = runif(n, 0, 1);
    for (int i = 0; i < n; i++) 
      out[i] = f1.invsample(u[i]);
    return out;
  }
  
  // returns the ratio involved in the GAS(1,1) model
  double calc_S(const double& x) {return 1 + x * f1.Sratio(x);}
  
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename underlying>
class Skewed
{  
  underlying f1;               // symmetric distribution
  double xi;                   // degree of assymetry
  double xi_lb;                // lower bound for "xi"
  double xi2;                  // := xi^2
  double num;                  // := 1 / (xi + 1/xi)
  double mu_xi, sig_xi, cutoff; 
  double pcut;
  double lncst;                // constant that accounts for asymmetry for kernel calculation
  double intgrl_1, intgrl_2;
  int Nsi;                     // number of subintervals for composite Simpson rule
  
  // Composite Simpson rule 
  double compositeSimpsons(const double& a, const double& b, const double& c, const int& p) {
    double x = a,  dx = (b-a)/(2*Nsi),  dx2 = 2*dx;
    double I1,  I2,  I3 = pow(c-x, p)*f1.pdf(x);
    double out = 0;
    for(int i = 0; i < Nsi; i++) {
      I1 = I3;
      I2 = pow(c-x-dx, p) * f1.pdf(x + dx);
      I3 = pow(c-x-dx2,p) * f1.pdf(x + dx2);
      out += dx/3 * (I1 + 4*I2 + I3);
      x += dx2;
    }
    return out;
  }
  
public:
  double Eabsz;                // := E[|z|] 
  double EzIpos;               // := E[z * I(z>=0)]  = Prob(z >= 0) * E[z | z >= 0]
  double EzIneg;               // := E[z * I(z<0)]   = Prob(z < 0)  * E[z | z < 0]
  double Ez2Ineg;              // := E[z^2 * I(z<0)] = Prob(z < 0)  * E[z^2 | z < 0]
  
  // constructor
  Skewed() {
    xi_lb = 0.01;
    Nsi = 5;
  }
  
  // constructor function called by higher-level classes (e.g. Garch). 
  // arguments are passed by reference and modified according to the distribution
  void constructor(std::string& name, int& nb_coeffs, NumericVector& coeffs_mean, NumericVector& coeffs_sd,
                   NumericVector& Sigma0, CharacterVector& label, NumericVector& lower, NumericVector& upper) {
    f1.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label, lower, upper);
    name.append("skew");
    nb_coeffs++;
    label.push_back("xi"); 
    coeffs_mean.push_back(1);     
    coeffs_sd.push_back(1);      
    Sigma0.push_back(1);   
    lower.push_back(xi_lb);
    upper.push_back(100);
  }
  
  // Load parameters and setup (this function should always be called first)  
  void loadparam(const NumericVector& theta, int& Ind) {
    f1.loadparam(theta, Ind);
    f1.set_M1();
    xi     = theta[Ind];
    xi2    = pow(xi, 2);
    num    = 1 / (xi + 1/xi);
    mu_xi  = f1.M1 * (xi - 1/xi);
    sig_xi = sqrt( (1 - pow(f1.M1, 2))*(xi2 + 1/xi2) + 2*pow(f1.M1, 2)  - 1);
    pcut   = num / xi;
    cutoff = - mu_xi / sig_xi;
  }
  
  // check prior 
  bool calc_r1() {return f1.calc_r1() && xi > xi_lb;}
  
  // setup for "calc_kernel" 
  void prep_kernel() {
    f1.prep_kernel();
    lncst = log(2 * sig_xi * num);
  }
  
  // returns PDF evaluated at "x"
  double calc_pdf(const double& x) {
    return 2*sig_xi*num * f1.pdf((sig_xi * x + mu_xi) * ((x < cutoff)? xi : 1/xi));
  }
  
  // returns CDF evaluated at "x"
  double calc_cdf(const double& x) {
    double tmp = sig_xi * x + mu_xi;
    return ((x < cutoff)?  2/xi * num * f1.cdf(tmp * xi)
              :  2 * num * (xi * f1.cdf(tmp / xi) + 1/xi) - 1);
  }
  
  // returns kernel of a single observation (must call "prep_kernel" first)
  double calc_kernel(const volatility& vol, const double& yi) {
    double sig = sqrt(vol.h);        
    double yi_xi = ((yi >= sig*cutoff) ? 1/xi : xi) * (sig_xi * yi + sig * mu_xi);
    return lncst + f1.kernel(vol, yi_xi);           
  }
  
  // returns the S variable involved in the GAS(1,1) model 
  double calc_S(const double& x) {
    double tmp = sig_xi * x + mu_xi;
    return 1 + x * ((x < cutoff)?  sig_xi * xi * f1.Sratio(tmp * xi)
                      :  sig_xi / xi * f1.Sratio(tmp / xi));
  }
  
  // maximum value of S
  double find_Smax() {
    double lb, ub, out;
    if (cutoff <= 0) 
      lb = cutoff, ub = 0;
    else
      lb = 0, ub = cutoff;
    if (lb != ub) {
      pair tmp0 = findMax(&Skewed::calc_S, this, lb, ub, 10, 1e-4, 10);
      out       = tmp0.first;
      if (out > 2) {
        double tmp1 = calc_S(tmp0.second - 1e-3),  tmp2 = calc_S(tmp0.second + 1e-3);
        out         = ((tmp1 < tmp2)? tmp2 : tmp1);
      }
    } 
    else
      out = 1;
    return (out);
  }
  
  // setup for moments of order 1
  void prep_moments1() {
    intgrl_1 = ((xi >= 1)? compositeSimpsons(0, mu_xi/xi, mu_xi/xi, 1) : compositeSimpsons(xi*mu_xi, 0, xi*mu_xi, 1));
  }
  
  // setup for moments of order 2
  void prep_moments2() {
    intgrl_2 = ((xi >= 1)? compositeSimpsons(0, mu_xi/xi, mu_xi/xi, 2) : compositeSimpsons(xi*mu_xi, 0, xi*mu_xi, 2));
  }
  
  // set Eabsz := E[|z|] 
  void set_Eabsz() {
    Eabsz = 2/sig_xi * num * (f1.M1 + 2 * ((xi >= 1)? xi2 : -1/xi2) * intgrl_1);
  } 
  
  // set EzIpos := E[z * I(z>=0)] = Prob(z >= 0) * E[z | z >= 0]
  void set_EzIpos() {
    EzIpos = 2/sig_xi * num * (0.5*f1.M1 + ((xi >= 1)? xi2 : -1/xi2) * intgrl_1);   
  }
  
  // set EzIneg := E[z * I(z<0)] = Prob(z < 0) * E[z | z < 0]
  void set_EzIneg() {
    EzIneg = - 2/sig_xi * num * (0.5*f1.M1 + ((xi >= 1)? xi2 : -1/xi2) * intgrl_1);
  }
  
  // set Ez2Ineg := E[z^2 * I(z<0)] = Prob(z < 0) * E[z^2 | z < 0]
  void set_Ez2Ineg() {
    double xi3 = xi2*xi,  xi4 = xi3*xi,  sig2_xi = pow(sig_xi, 2),  M1_2 = pow(f1.M1,  2);
    Ez2Ineg = ((xi >= 1)?  2/sig2_xi * num * ( 0.5/xi3 * (1 + M1_2*(xi4-1)) + xi3 * intgrl_2 )
                 :  2/(xi3 * sig2_xi) * num * (0.5 - 0.5*M1_2*(1-xi4) - intgrl_2 ));
  } 
  
  // returns a random vector of length "n"   
  NumericVector rndgen(const int& n) {
    NumericVector out(n);
    NumericVector u = runif(n, 0, 1);
    for (int i = 0; i < n; i++) 
      out[i] = ((u[i] < pcut)? (f1.invsample(0.5 * u[i] * (xi2 + 1)) / xi - mu_xi) / sig_xi
                  : (f1.invsample(0.5 * u[i] * (1 + 1/xi2) - 0.5/xi2 + 0.5) * xi - mu_xi  ) / sig_xi);
    return out;
  }
};


template <typename distribution>
class sGARCH 
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
  sGARCH() {
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
    name        = "sGARCH_";
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

template <typename distribution>
class gjrGARCH
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
  gjrGARCH() {
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
    name        = "gjrGARCH_";
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
class eGARCH
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
  eGARCH() {
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
    name        = "eGARCH_";
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
class tGARCH 
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
  tGARCH() {
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
    name        = "tGARCH_";
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

//---------------------- Gas ----------------------// 
template <typename distribution>
class GAS 
{
  distribution fz;                // distribution of innovations (e.g. Symmetric<Normal>)
  double alpha0, alpha1, beta;    // coefficients
  double invES2;                  // = 1 / E[S^2]
  bool fz_r1;                     // indicates if distribution "fz" has sensible parameters
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
  GAS() {
    ineq_lb     = 1e-4; 
    ineq_ub     = 0.9999; 
    label       = CharacterVector::create("alpha0", "alpha1", "beta");
    coeffs_mean = NumericVector::create(   0.1,      0.1,      0.9  );
    coeffs_sd   = NumericVector::create(   2,        2,        2    );
    Sigma0      = NumericVector::create(   1,        1,        1    );
    lower       = NumericVector::create(   1e-4,     1e-4,     1e-4 );
    upper       = NumericVector::create(   100,      10,     0.9999);
    nb_coeffs   = label.size();
    nb_coeffs_model = 3;
    name        = "GAS_";
    fz.constructor(name, nb_coeffs, coeffs_mean, coeffs_sd, Sigma0, label, lower, upper);
  }
  
  // returns the S variable as defined in the "GAS.pdf" document (nested call)
  double calc_S(const double& z) {return fz.calc_S(z);}
  
  // integrated function for computing invES2
  double integrand_ES2(const double& z) {
    double out = pow(calc_S(z), 2) * fz.calc_pdf(z);
    // returns "0" in case of numerical indetermination
    return (IsInfNan(out)? 0 : out); 
  }
  
  // computes ES2 if the distribution has sensible parameters
  void set_invES2() { 
    fz_r1  = fz.calc_r1();
    invES2 = (fz_r1? 1 / adaptiveSimpsons(&GAS::integrand_ES2, this, -1000, 1000, 1e-5, 160)  :  NAN);
  }
  
  // set all parameters (including those of the distribution).
  // this function should always be called first 
  void loadparam(const NumericVector& theta) {
    alpha0 = theta[0],  alpha1 = theta[1],  beta = theta[2];
    int Ind = 3;
    fz.loadparam(theta, Ind);
    set_invES2();           
  }
  
  // empty
  void prep_ineq_vol() { };
  
  // TODO  
  double ineq_func() {
    return beta - 2 * alpha1 * invES2;
  } 
  
  // computes r1 
  bool calc_r1() {
    return fz_r1
    && alpha0 > lower[0] && alpha1 > lower[1] && beta > lower[2] && beta < upper[2] 
    && (ineq_func() > ineq_lb);
  }
  
  // initialize volatility to its undonditional expected value  
  volatility set_vol(const double& y0) {
    volatility out;
    out.h   = alpha0 / (1 - beta); 
    out.lnh = log(out.h);
    return out;
  }
  
  // increment volatility
  void increment_vol(volatility& vol, const double& yim1) { 
    double S = calc_S(yim1 / sqrt(vol.h));
    vol.h   = alpha0 + vol.h * (beta - 2 * alpha1 * invES2 * S);  
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


#endif // spec.h
