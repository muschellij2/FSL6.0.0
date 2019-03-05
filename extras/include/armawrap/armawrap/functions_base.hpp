#ifndef __FUNCTIONS_BASE_HPP__
#define __FUNCTIONS_BASE_HPP__

#include <math.h>

namespace armawrap {
  /**
   * Standalone versions of AWBase::methods, which accept armawrap:: types.
   */
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Minimum(             wT &a) { return Minimum(             a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Maximum(             wT &a) { return Maximum(             a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type MinimumAbsoluteValue(wT &a) { return MinimumAbsoluteValue(a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type MaximumAbsoluteValue(wT &a) { return MaximumAbsoluteValue(a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type SumSquare(           wT &a) { return SumSquare(           a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type SumAbsoluteValue(    wT &a) { return SumAbsoluteValue(    a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Sum(                 wT &a) { return Sum(                 a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Norm1(               wT &a) { return Norm1(               a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type NormInfinity(        wT &a) { return NormInfinity(        a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type NormFrobenius(       wT &a) { return NormFrobenius(       a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Trace(               wT &a) { return Trace(               a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Determinant(         wT &a) { return Determinant(         a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value, AWLogAndSign<typename wT::elem_type> >::type LogDeterminant(      wT &a) { return LogDeterminant(      a.get_at_ref()); }
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type IsZero(              wT &a) { return IsZero(              a.get_at_ref()); }

  /*
   * Those same functions which accept arma:: types so, e.g. 
   * expressions and subviews may be passed to them.
   */

  template<typename elem_type, typename T>
  inline elem_type 
  Minimum(arma::Base<elem_type, T> const &a) { 
    return a.get_ref().min(); 
  }

  template<typename elem_type, typename T>
  inline elem_type 
  Maximum(arma::Base<elem_type, T> const &a) { 
    return a.get_ref().max(); 
  }

  template<typename elem_type, typename T>
  inline elem_type 
  MinimumAbsoluteValue(arma::Base<elem_type, T> const &a) { 
    return arma::abs(a.get_ref()).eval().min(); 
  }

  template<typename elem_type, typename T>
  inline elem_type 
  MaximumAbsoluteValue(arma::Base<elem_type, T> const &a) { 
    return arma::abs(a.get_ref()).eval().max(); 
  }

  template<typename elem_type, typename T>
  inline elem_type 
  SumSquare(arma::Base<elem_type, T> const &a) { 
    return arma::accu(arma::square(a.get_ref())); 
  }

  template<typename elem_type, typename T>
  inline elem_type 
  SumAbsoluteValue(arma::Base<elem_type, T> const &a) { 
    return arma::accu(arma::abs(a.get_ref())); 
  }

  template<typename elem_type, typename T>
  inline elem_type 
  Sum(arma::Base<elem_type, T> const &a) { 
    return arma::accu(a.get_ref()); 
  }

  template<typename elem_type, typename T>
  inline elem_type 
  Trace(arma::Base<elem_type, T> const &a) {
    return arma::trace(a.get_ref());
  }

  template<typename elem_type, typename T>
  inline elem_type 
  Norm1(arma::Base<elem_type, T> const &a) {
    return arma::as_scalar(
      arma::max(
        arma::vectorise(
          arma::sum(
            arma::abs(a.get_ref()), 0))));
  }

  template<typename elem_type, typename T>
  inline elem_type 
  NormInfinity(arma::Base<elem_type, T> const &a) {
    return arma::as_scalar(
      arma::max(
        arma::vectorise(
          arma::sum(
            arma::abs(a.get_ref()), 1))));
  }

  template<typename elem_type, typename T>
  inline elem_type 
  NormFrobenius(arma::Base<elem_type, T> const &a) {
    return sqrt(SumSquare(a));
  }

  template<typename elem_type, typename T>
  inline elem_type 
  Determinant(arma::Base<elem_type, T> const &a) {
    return arma::det(a.get_ref());
  }

  template<typename elem_type, typename T>
  inline AWLogAndSign<elem_type> 
  LogDeterminant(arma::Base<elem_type, T> const &a) {
    double val;
    double sign;
    arma::log_det(val, sign, a.get_ref());
    return AWLogAndSign<elem_type>(val, sign);
  }

  template<typename elem_type, typename T>
  inline bool 
  IsZero(arma::Base<elem_type, T> const &a) {
    return arma::all(arma::vectorise(a.get_ref()) == 0);
  }
}

#endif /* __FUNCTIONS_BASE_HPP__ */
