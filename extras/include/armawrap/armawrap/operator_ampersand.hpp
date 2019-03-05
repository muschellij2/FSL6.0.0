/**
 * Vertical concatenation operator.
 */
#ifndef __OPERATOR_AMPERSAND_HPP__
#define __OPERATOR_AMPERSAND_HPP__

namespace armawrap {

  // Object & Object (same type)
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
    operator&(const T1 &X, const T2 &Y) {

    // The third parameter is 0 for horizontal (column-wise) 
    // concatenation, and 1 for horizontal concatenation
    return AWGlue<typename armawrap_type_map<T1>::type,
                  typename armawrap_type_map<T2>::type,
                  arma::glue_join>(X, Y, 0);
  }
}

#endif
