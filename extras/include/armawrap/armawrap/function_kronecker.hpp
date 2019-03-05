/**
 * Kronecker product.
 */
#ifndef __FUNCTION_KRONECKER_HPP__
#define __FUNCTION_KRONECKER_HPP__


namespace armawrap {
  template<typename T1, typename T2>
  inline
  const arma::Glue<
    typename armawrap_type_map<T1>::type, 
    typename armawrap_type_map<T2>::type,
    arma::glue_kron>
  KP(const T1 &X, const T2 &Y) {

    /*
     * I don't know why, but calling arma::kron will not compile for
     * heterogeneous types (e.g. Mat, Col). The compiler seems to
     * get confused with the types.
     */
    return arma::Glue<typename armawrap_type_map<T1>::type, 
                      typename armawrap_type_map<T2>::type, 
                      arma::glue_kron>(X.get_ref(), Y.get_ref());
  }
}

#endif
