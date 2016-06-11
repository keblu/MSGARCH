#include "MSgarch.h"

//===========================================================================//
//===========================================================================//
//================================== Garch ==================================//
//===========================================================================//
//===========================================================================//
RCPP_MODULE(Garch){   
  // Garch-normal-symmetric
  class_<Garch_normal_sym>("Garch_normal_sym")  
  .constructor()
  .field( "name",        &Garch_normal_sym ::name )
  .field( "theta0",      &Garch_normal_sym ::theta0 )
  .field( "Sigma0",      &Garch_normal_sym ::Sigma0 )
  .field( "label",       &Garch_normal_sym ::label )
  .field( "lower",       &Garch_normal_sym ::lower )
  .field( "upper",       &Garch_normal_sym ::upper )
  .field( "ineq_lb",     &Garch_normal_sym ::ineq_lb )
  .field( "ineq_ub",     &Garch_normal_sym ::ineq_ub )
  .field( "NbParams",    &Garch_normal_sym ::NbParams )
  .field( "NbParamsModel",&Garch_normal_sym::NbParamsModel)
  .method( "f_sim",      &Garch_normal_sym ::f_sim )
  .method( "f_pdf",      &Garch_normal_sym ::f_pdf )
  .method( "f_cdf",      &Garch_normal_sym ::f_cdf )
  .method( "f_rnd",      &Garch_normal_sym ::f_rnd )
  .method( "calc_ht",    &Garch_normal_sym ::calc_ht )
  .method( "eval_model", &Garch_normal_sym ::eval_model )
  .method( "ineq_func",  &Garch_normal_sym ::ineq_func )
  .method( "f_unc_vol",  &Garch_normal_sym ::f_unc_vol)
  ;
}