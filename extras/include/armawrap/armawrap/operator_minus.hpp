/**
 * Elementwise minus operator.
 */
#ifndef __OPERATOR_MINUS_HPP__
#define __OPERATOR_MINUS_HPP__


namespace armawrap {

  // Object - Object (same type)
  template<typename T1, typename T2>
  inline
  typename enable_if<
    either_are_armawrap_type<T1, T2>::value &&
  arma::is_same_type<
    typename T1::elem_type,
    typename T2::elem_type>::value,
    AWEGlue<
      typename armawrap_type_map<T1>::type,
      typename armawrap_type_map<T2>::type,
      arma::eglue_minus> >::type
    operator-(const T1 &X, const T2 &Y) {
    return AWEGlue<typename armawrap_type_map<T1>::type,
                   typename armawrap_type_map<T2>::type, 
                   arma::eglue_minus>(X, Y);
  }

  // Object - Scalar
  template<typename T>
  inline
  typename enable_if<
    is_armawrap_type<T>::value,
    AWEOp<
      typename armawrap_type_map<T>::type, 
      arma::eop_scalar_minus_post> >::type
  operator-(const T &X, const typename T::elem_type k) {
    return AWEOp<typename armawrap_type_map<T>::type, 
                 arma::eop_scalar_minus_post>(X, k);
  }

  // Scalar - Object
  template<typename T>
  inline
  typename enable_if<
    is_armawrap_type<T>::value,
    AWEOp<
      typename armawrap_type_map<T>::type, 
      arma::eop_scalar_minus_pre> >::type
  operator-(const typename T::elem_type k, const T &X) {
    return AWEOp<typename armawrap_type_map<T>::type, 
                 arma::eop_scalar_minus_pre>(X, k);
  }

  // -Object (unary negation)
  template<typename T>
  inline 
  typename enable_if<
    is_armawrap_type<T>::value,
    AWEOp<
      typename armawrap_type_map<T>::type, 
      arma::eop_neg> >::type
  operator-(const T &X) {
    return AWEOp<typename armawrap_type_map<T>::type, 
                 arma::eop_neg>(X);
  }
}

#endif
