#ifndef __FUNCTION_CHOLESKY_HPP__
#define __FUNCTION_CHOLESKY_HPP__

/*
 * Cholesky decomposition.
 */

namespace armawrap {
  template<typename T>
  inline
  arma::Op<typename armawrap_type_map<T>::type, arma::op_chol>
  Cholesky(const T &X) {
    return arma::chol(X);
  }
}

#endif /* __FUNCTION_CHOLESKY_HPP__ */
