/**
 * Schur product (element-wise multiplication).
 */
#ifndef __FUNCTION_SCHUR_HPP__
#define __FUNCTION_SCHUR_HPP__

namespace armawrap {

  // SP (arma type, arma type)
  template<typename T1, typename T2>
  inline
  typename enable_if<
    both_are_arma_type<T1, T2>::value &&
  arma::is_same_type<
    typename T1::elem_type, 
    typename T2::elem_type>::value,
    AWEGlue<T1, T2, arma::eglue_schur> >::type
    SP(const T1 &X, const T2 &Y) {
    return AWEGlue<T1, T2, arma::eglue_schur>(X, Y);
  }

  // SP(arma/armawrap type, arma/armawrap type) (same data type)
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
      arma::eglue_schur> >::type
    SP(const T1 &X, const T2 &Y) {
    return AWEGlue<typename armawrap_type_map<T1>::type, 
                   typename armawrap_type_map<T2>::type,
                   arma::eglue_schur>(X, Y);
  }
}

#endif
