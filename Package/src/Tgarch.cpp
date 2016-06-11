#include "MSgarch.h"

//============================================================================//
//============================================================================//
//================================== Tgarch ==================================//
//============================================================================//
//============================================================================//
RCPP_MODULE(Tgarch){  
  
  // Tgarch-normal-symmetric
  class_<Tgarch_normal_sym>("Tgarch_normal_sym")  
  .constructor()
  .field( "name",        &Tgarch_normal_sym ::name )
  .field( "theta0",      &Tgarch_normal_sym ::theta0 )
  .field( "Sigma0",      &Tgarch_normal_sym ::Sigma0 )
  .field( "label",       &Tgarch_normal_sym ::label )
  .field( "lower",       &Tgarch_normal_sym ::lower )
  .field( "upper",       &Tgarch_normal_sym ::upper )
  .field( "ineq_lb",     &Tgarch_normal_sym ::ineq_lb )
  .field( "ineq_ub",     &Tgarch_normal_sym ::ineq_ub )
  .field( "NbParams",    &Tgarch_normal_sym ::NbParams )
  .field( "NbParamsModel",&Tgarch_normal_sym ::NbParamsModel)
  .method( "f_sim",      &Tgarch_normal_sym ::f_sim )
  .method( "f_pdf",      &Tgarch_normal_sym ::f_pdf )
  .method( "f_cdf",      &Tgarch_normal_sym ::f_cdf )
  .method( "f_rnd",      &Tgarch_normal_sym ::f_rnd )
  .method( "calc_ht",    &Tgarch_normal_sym ::calc_ht )
  .method( "eval_model", &Tgarch_normal_sym ::eval_model )
  .method( "ineq_func",  &Tgarch_normal_sym ::ineq_func )
  .method( "f_unc_vol",  &Tgarch_normal_sym ::f_unc_vol)
  ;
}
