#ifndef __FUNCTION_SVD_HPP__
#define __FUNCTION_SVD_HPP__

/*
 * Singular value decomposition.
 */

namespace armawrap {
  
  template<typename T1, typename T2, typename T3, typename T4>
  inline  void SVD(const T1  &A, 
                   T2        &D,
                   T3        &U, 
                   T4        &V,
                   bool withU=true,
                   bool withV=true) {

    arma::Col<typename T1::elem_type> s;
    arma::Mat<typename T1::elem_type> Ut;
    arma::Mat<typename T1::elem_type> At(A);

    arma::svd_econ(Ut, s, V, At);

    if (!withV) V = 0;
    if ( withU) U = Ut.submat(0, 0, At.n_rows-1, At.n_cols-1);

    D = s;
  }

  template<typename T1, typename T2>
  inline void SVD(const T1 &A, T2 &D) {

    AWMatrix<typename T1::elem_type> U;
    AWMatrix<typename T1::elem_type> V;

    SVD(A, D, U, V, false, false);
  }

  template<typename T1, typename T2, typename T3>
  inline void SVD(const T1 &A, 
                  T2       &D,
                  T3       &U,
                  bool withU=true) {

    AWMatrix<typename T1::elem_type> V(D.Nrows(), D.Ncols());
    SVD(A, D, U, V, withU, false);
  }
}

#endif /* __FUNCTION_SVD_HPP__ */
