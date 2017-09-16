#include "Symmetric.h"
#include "Skewed.h"
#include "Student.h"

typedef Symmetric<Student> std_sym;
typedef Skewed<Student> std_skew;

RCPP_MODULE(Student) {
  class_<std_sym>("std_sym")
      .constructor()
      .method("f_pdf", &std_sym::calc_pdf)
      .method("f_cdf", &std_sym::calc_cdf)
      .method("f_invsample", &std_sym::calc_invsample)
      .method("f_rnd", &std_sym::rndgen)
      .method("set_Eabsz", &std_sym::set_Eabsz)
      .field("Eabsz", &std_sym::Eabsz)
      .method("set_EzIpos", &std_sym::set_EzIpos)
      .field("EzIpos", &std_sym::EzIpos)
      .method("set_EzIneg", &std_sym::set_EzIneg)
      .field("EzIneg", &std_sym::EzIneg)
      .method("set_Ez2Ineg", &std_sym::set_Ez2Ineg)
      .field("Ez2Ineg", &std_sym::Ez2Ineg)
      .method("load_theta", &std_sym::loadparam);
  class_<std_skew>("std_skew")
      .constructor()
      .method("f_pdf", &std_skew::calc_pdf)
      .method("f_cdf", &std_skew::calc_cdf)
      .method("f_invsample", &std_skew::calc_invsample)
      .method("f_rnd", &std_skew::rndgen)
      .method("set_Eabsz", &std_skew::set_Eabsz)
      .field("Eabsz", &std_skew::Eabsz)
      .method("set_EzIpos", &std_skew::set_EzIpos)
      .field("EzIpos", &std_skew::EzIpos)
      .method("set_EzIneg", &std_skew::set_EzIneg)
      .field("EzIneg", &std_skew::EzIneg)
      .method("set_Ez2Ineg", &std_skew::set_Ez2Ineg)
      .field("Ez2Ineg", &std_skew::Ez2Ineg)
      .method("load_theta", &std_skew::loadparam);
}