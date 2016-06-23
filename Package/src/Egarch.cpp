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
  // Egarch-student-symmetric
  class_<Egarch_student_sym>("Egarch_student_sym")  
    .constructor()
    .field( "name",        &Egarch_student_sym ::name )
    .field( "theta0",      &Egarch_student_sym ::theta0 )
    .field( "Sigma0",      &Egarch_student_sym ::Sigma0 )
    .field( "label",       &Egarch_student_sym ::label )
    .field( "lower",       &Egarch_student_sym ::lower )
    .field( "upper",       &Egarch_student_sym ::upper )
    .field( "ineq_lb",     &Egarch_student_sym ::ineq_lb )
    .field( "ineq_ub",     &Egarch_student_sym ::ineq_ub )
    .field( "NbParams",    &Egarch_student_sym ::NbParams )
	.field( "NbParamsModel",&Egarch_student_sym::NbParamsModel)
    .method( "f_sim",      &Egarch_student_sym ::f_sim )
    .method( "f_pdf",      &Egarch_student_sym ::f_pdf )
    .method( "f_cdf",      &Egarch_student_sym ::f_cdf )
    .method( "f_rnd",      &Egarch_student_sym ::f_rnd )
    .method( "calc_ht",    &Egarch_student_sym ::calc_ht )
    .method( "eval_model", &Egarch_student_sym ::eval_model )
    .method( "ineq_func",  &Egarch_student_sym ::ineq_func )
    .method( "f_unc_vol",  &Egarch_student_sym ::f_unc_vol)
  ;
  // Egarch-ged-symmetric
  class_<Egarch_ged_sym>("Egarch_ged_sym")  
    .constructor()
    .field( "name",        &Egarch_ged_sym ::name )
    .field( "theta0",      &Egarch_ged_sym ::theta0 )
    .field( "Sigma0",      &Egarch_ged_sym ::Sigma0 )
    .field( "label",       &Egarch_ged_sym ::label )
    .field( "lower",       &Egarch_ged_sym ::lower )
    .field( "upper",       &Egarch_ged_sym ::upper )
    .field( "ineq_lb",     &Egarch_ged_sym ::ineq_lb )
    .field( "ineq_ub",     &Egarch_ged_sym ::ineq_ub )
    .field( "NbParams",    &Egarch_ged_sym ::NbParams )
	.field( "NbParamsModel",&Egarch_ged_sym::NbParamsModel)
    .method( "f_sim",      &Egarch_ged_sym ::f_sim )
    .method( "f_pdf",      &Egarch_ged_sym ::f_pdf )
    .method( "f_cdf",      &Egarch_ged_sym ::f_cdf )
    .method( "f_rnd",      &Egarch_ged_sym ::f_rnd )
    .method( "calc_ht",    &Egarch_ged_sym ::calc_ht )
    .method( "eval_model", &Egarch_ged_sym ::eval_model )
    .method( "ineq_func",  &Egarch_ged_sym ::ineq_func )
    .method( "f_unc_vol",  &Egarch_ged_sym ::f_unc_vol)
  ;
}
