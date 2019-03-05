/**
 * Horizontal concatenation operator.
 */
#ifndef __OPERATOR_PIPE_HPP__
#define __OPERATOR_PIPE_HPP__

namespace armawrap {
  template<typename T1, typename T2>
  inline
  typename enable_if<
    either_are_armawrap_type<T1, T2>::value &&
  arma::is_same_type<
    typename T1::elem_type,
    typename T2::elem_type>::value,
    AWGlue<
      typename armawrap_type_map<T1>::type,
      typename armawrap_type_map<T2>::type,
      arma::glue_join> >::type
    operator|(const T1 &X, const T2 &Y) {

    // The third parameter (1) indiciates 
    // horizontal concatenation
    return AWGlue<typename armawrap_type_map<T1>::type,
                  typename armawrap_type_map<T2>::type,
                  arma::glue_join>(X, Y, 1);
  }
}

#endif
