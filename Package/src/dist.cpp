#include "MSgarch.h"


RCPP_MODULE(dist){  
  
  class_<norm_sym>("norm_sym")  
	.constructor()
	.method( "f_pdf",      &norm_sym ::calc_pdf )
	.method( "f_cdf",      &norm_sym ::calc_cdf )
	.method( "f_invsample",&norm_sym ::calc_invsample )
	.method( "f_rnd",      &norm_sym ::rndgen )
	.method ("load_theta", &norm_sym ::loadparam)
  ;
    class_<norm_skew>("norm_skew")  
	.constructor()
	.method( "f_pdf",      &norm_skew ::calc_pdf )
	.method( "f_cdf",      &norm_skew ::calc_cdf )
	.method( "f_invsample",&norm_skew ::calc_invsample )
	.method( "f_rnd",      &norm_skew ::rndgen )
	.method ("load_theta", &norm_skew ::loadparam)
  ;
   class_<std_sym>("std_sym")  
	.constructor()
	.method( "f_pdf",      &std_sym ::calc_pdf )
	.method( "f_cdf",      &std_sym ::calc_cdf )
	.method( "f_invsample",&std_sym ::calc_invsample )
	.method( "f_rnd",      &std_sym ::rndgen )
	.method ("load_theta", &std_sym ::loadparam)
  ;
    class_<std_skew>("std_skew")  
	.constructor()
	.method( "f_pdf",      &std_skew ::calc_pdf )
	.method( "f_cdf",      &std_skew ::calc_cdf )
	.method( "f_invsample", &std_skew ::calc_invsample )
	.method( "f_rnd",      &std_skew ::rndgen )
	.method ("load_theta", &std_skew ::loadparam)
  ;
    class_<ged_sym>("ged_sym")  
	.constructor()
	.method( "f_pdf",      &ged_sym ::calc_pdf )
	.method( "f_cdf",      &ged_sym ::calc_cdf )
	.method( "f_invsample",&ged_sym ::calc_invsample )
	.method( "f_rnd",      &ged_sym ::rndgen )
	.method ("load_theta", &ged_sym ::loadparam)
  ;
    class_<ged_skew>("ged_skew")  
	.constructor()
	.method( "f_pdf",      &ged_skew ::calc_pdf )
	.method( "f_cdf",      &ged_skew ::calc_cdf )
	.method( "f_invsample",&ged_skew ::calc_invsample )
	.method( "f_rnd",      &ged_skew ::rndgen )
	.method ("load_theta", &ged_skew ::loadparam)
  ;
}
