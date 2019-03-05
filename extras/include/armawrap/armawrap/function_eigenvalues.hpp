#ifndef __FUNCTION_EIGENVALUES_HPP__
#define __FUNCTION_EIGENVALUES_HPP__

namespace armawrap {
  template<typename aT, typename dT>
  inline void EigenValues(const aT &A, dT &D) {

    arma::Col<typename aT::elem_type> vals;
    arma::eig_sym(vals, A);
    D = vals;
  }

  template<typename aT, typename dT, typename vT>
  inline void EigenValues(const aT &A, dT &D, vT &V) {
    arma::Col<typename aT::elem_type> vals;
    arma::eig_sym(vals, V, A);
    D = vals;
  }


  template<typename xT, typename dT, typename aT, typename vT>
  inline void Jacobi(const xT &X, dT &D, aT &A, vT &V, bool eivec) {
    EigenValues(A, D, V);
  }

  template<typename xT, typename dT>
  inline void Jacobi(const xT &X, dT &D) {
    EigenValues(X, D);
  }
  
  template<typename xT, typename dT, typename vT>
  inline void Jacobi(const xT &X, dT &D, vT &V) {
    EigenValues(X, D, V);
  }
}

#endif /* __FUNCTION_EIGENVALUES_HPP__ */
