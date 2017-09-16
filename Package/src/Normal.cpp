#include "Symmetric.h"
#include "Skewed.h"
#include "Normal.h"

typedef Symmetric<Normal> norm_sym;
typedef Skewed<Normal> norm_skew;

RCPP_MODULE(Normal) {
  class_<norm_sym>("norm_sym")
      .constructor()
      .method("f_pdf", &norm_sym::calc_pdf)
      .method("f_cdf", &norm_sym::calc_cdf)
      .method("f_invsample", &norm_sym::calc_invsample)
      .method("f_rnd", &norm_sym::rndgen)
      .method("set_Eabsz", &norm_sym::set_Eabsz)
      .field("Eabsz", &norm_sym::Eabsz)
      .method("set_EzIpos", &norm_sym::set_EzIpos)
      .field("EzIpos", &norm_sym::EzIpos)
      .method("set_EzIneg", &norm_sym::set_EzIneg)
      .field("EzIneg", &norm_sym::EzIneg)
      .method("set_Ez2Ineg", &norm_sym::set_Ez2Ineg)
      .field("Ez2Ineg", &norm_sym::Ez2Ineg)
      .method("load_theta", &norm_sym::loadparam);
  class_<norm_skew>("norm_skew")
      .constructor()
      .method("f_pdf", &norm_skew::calc_pdf)
      .method("f_cdf", &norm_skew::calc_cdf)
      .method("f_invsample", &norm_skew::calc_invsample)
      .method("f_rnd", &norm_skew::rndgen)
      .method("set_Eabsz", &norm_skew::set_Eabsz)
      .field("Eabsz", &norm_skew::Eabsz)
      .method("set_EzIpos", &norm_skew::set_EzIpos)
      .field("EzIpos", &norm_skew::EzIpos)
      .method("set_EzIneg", &norm_skew::set_EzIneg)
      .field("EzIneg", &norm_skew::EzIneg)
      .method("set_Ez2Ineg", &norm_skew::set_Ez2Ineg)
      .field("Ez2Ineg", &norm_skew::Ez2Ineg)
      .method("load_theta", &norm_skew::loadparam);
}