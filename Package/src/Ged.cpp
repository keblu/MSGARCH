#include "Symmetric.h"
#include "Skewed.h"
#include "Ged.h"

typedef Symmetric<Ged> ged_sym;
typedef Skewed<Ged> ged_skew;

RCPP_MODULE(Ged) {
  class_<ged_sym>("ged_sym")
      .constructor()
      .method("f_pdf", &ged_sym::calc_pdf)
      .method("f_cdf", &ged_sym::calc_cdf)
      .method("f_invsample", &ged_sym::calc_invsample)
      .method("f_rnd", &ged_sym::rndgen)
      .method("set_Eabsz", &ged_sym::set_Eabsz)
      .field("Eabsz", &ged_sym::Eabsz)
      .method("set_EzIpos", &ged_sym::set_EzIpos)
      .field("EzIpos", &ged_sym::EzIpos)
      .method("set_EzIneg", &ged_sym::set_EzIneg)
      .field("EzIneg", &ged_sym::EzIneg)
      .method("set_Ez2Ineg", &ged_sym::set_Ez2Ineg)
      .field("Ez2Ineg", &ged_sym::Ez2Ineg)
      .method("load_theta", &ged_sym::loadparam);
  class_<ged_skew>("ged_skew")
      .constructor()
      .method("f_pdf", &ged_skew::calc_pdf)
      .method("f_cdf", &ged_skew::calc_cdf)
      .method("f_invsample", &ged_skew::calc_invsample)
      .method("f_rnd", &ged_skew::rndgen)
      .method("set_Eabsz", &ged_skew::set_Eabsz)
      .field("Eabsz", &ged_skew::Eabsz)
      .method("set_EzIpos", &ged_skew::set_EzIpos)
      .field("EzIpos", &ged_skew::EzIpos)
      .method("set_EzIneg", &ged_skew::set_EzIneg)
      .field("EzIneg", &ged_skew::EzIneg)
      .method("set_Ez2Ineg", &ged_skew::set_Ez2Ineg)
      .field("Ez2Ineg", &ged_skew::Ez2Ineg)
      .method("load_theta", &ged_skew::loadparam);
}