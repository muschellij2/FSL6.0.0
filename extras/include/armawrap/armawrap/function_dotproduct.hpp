#ifndef __FUNCTION_DOTPRODUCT_HPP__
#define __FUNCTION_DOTPRODUCT_HPP__

namespace armawrap {

  template<typename T1, typename T2>
  inline
  typename enable_if<
    either_are_armawrap_type<T1, T2>::value &&
  arma::is_same_type<
    typename T1::elem_type, 
    typename T2::elem_type>::value,
    typename T1::elem_type>::type
    DotProduct(const T1 &X, const T2 &Y) {
    return arma::dot<typename armawrap_type_map<T1>::type, 
                     typename armawrap_type_map<T2>::type>(X, Y);
  }
}

#endif /* __FUNCTION_DOTPRODUCT_HPP__ */
