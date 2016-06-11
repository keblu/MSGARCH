#include "MSgarch.h"

//============================================================================//
//============================================================================//
//================================== Egarch ==================================//
//============================================================================//
//============================================================================//
RCPP_MODULE(Egarch){  
  
  // Egarch-normal-symmetric
  class_<Egarch_normal_sym>("Egarch_normal_sym")  
  .constructor()
  .field( "name",        &Egarch_normal_sym ::name )
  .field( "theta0",      &Egarch_normal_sym ::theta0 )
  .field( "Sigma0",      &Egarch_normal_sym ::Sigma0 )
  .field( "label",       &Egarch_normal_sym ::label )
  .field( "lower",       &Egarch_normal_sym ::lower )
  .field( "upper",       &Egarch_normal_sym ::upper )
  .field( "ineq_lb",     &Egarch_normal_sym ::ineq_lb )
  .field( "ineq_ub",     &Egarch_normal_sym ::ineq_ub )
  .field( "NbParams",    &Egarch_normal_sym ::NbParams )
  .field( "NbParamsModel",&Egarch_normal_sym::NbParamsModel)
  .method( "f_sim",      &Egarch_normal_sym ::f_sim )
  .method( "f_pdf",      &Egarch_normal_sym ::f_pdf )
  .method( "f_cdf",      &Egarch_normal_sym ::f_cdf )
  .method( "f_rnd",      &Egarch_normal_sym ::f_rnd )
  .method( "calc_ht",    &Egarch_normal_sym ::calc_ht )
  .method( "eval_model", &Egarch_normal_sym ::eval_model )
  .method( "ineq_func",  &Egarch_normal_sym ::ineq_func )
  .method( "f_unc_vol",  &Egarch_normal_sym ::f_unc_vol) 
  ;
}
