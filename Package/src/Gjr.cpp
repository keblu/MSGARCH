#include "MSgarch.h"

//=========================================================================//
//=========================================================================//
//================================== Gjr ==================================//
//=========================================================================//
//=========================================================================//
RCPP_MODULE(Gjr){   
  
  // Gjr-normal-symmetric
  class_<Gjr_normal_sym>("Gjr_normal_sym")  
  .constructor()
  .field( "name",        &Gjr_normal_sym ::name )
  .field( "theta0",      &Gjr_normal_sym ::theta0 )
  .field( "Sigma0",      &Gjr_normal_sym ::Sigma0 )
  .field( "label",       &Gjr_normal_sym ::label )
  .field( "lower",       &Gjr_normal_sym ::lower )
  .field( "upper",       &Gjr_normal_sym ::upper )
  .field( "ineq_lb",     &Gjr_normal_sym ::ineq_lb )
  .field( "ineq_ub",     &Gjr_normal_sym ::ineq_ub )
  .field( "NbParams",    &Gjr_normal_sym ::NbParams )
  .field( "NbParamsModel",&Gjr_normal_sym ::NbParamsModel)
  .method( "f_sim",      &Gjr_normal_sym ::f_sim )
  .method( "f_pdf",      &Gjr_normal_sym ::f_pdf )
  .method( "f_cdf",      &Gjr_normal_sym ::f_cdf )
  .method( "f_rnd",      &Gjr_normal_sym ::f_rnd )
  .method( "calc_ht",    &Gjr_normal_sym ::calc_ht )
  .method( "eval_model", &Gjr_normal_sym ::eval_model )
  .method( "ineq_func",  &Gjr_normal_sym ::ineq_func )
  .method( "f_unc_vol",  &Gjr_normal_sym ::f_unc_vol)
  ;
}
