/**
 * Object level equality/inequality operators. Newmat does not support 
 * element-wise equality/inequality testing, nor does it support 
 * scalar <-> object testing.
 */
#ifndef __OPERATOR_EQUALITY_HPP__
#define __OPERATOR_EQUALITY_HPP__

namespace armawrap {

  // Object == Object
  template<typename T1, typename T2>
  inline 
  typename enable_if<
    either_are_armawrap_type<T1, T2>::value,
    bool>::type
  operator==(const T1 &X, const T2 &Y) {

    return arma::all(arma::vectorise(
      arma::operator==<typename armawrap_type_map<T1>::type,
                       typename armawrap_type_map<T2>::type>(X, Y)));
  }

  // Object != Object
  template<typename T1, typename T2>
  inline 
  typename enable_if<
    either_are_armawrap_type<T1, T2>::value,
    bool>::type
  operator!=(const T1 &X, const T2&Y) {

    return !arma::all(arma::vectorise(
      arma::operator==<typename armawrap_type_map<T1>::type,
                       typename armawrap_type_map<T2>::type>(X, Y)));
  }
}

#endif
