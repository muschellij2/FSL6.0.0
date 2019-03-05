#ifndef __FUNCTION_SORT_HPP__
#define __FUNCTION_SORT_HPP__

namespace armawrap {
  template<typename elem_type>
  inline void _sort(AWMatrix<elem_type> &X, const char *dir) {

    arma::Col<elem_type> Xc = arma::vectorise<typename armawrap_type_map<AWMatrix<elem_type> >::type>(X);
    arma::Mat<elem_type> Xs = arma::sort(Xc, dir);

    Xs.reshape(X.n_cols, X.n_rows);
    X = Xs.t();
  }

  template<typename elem_type>
  inline void _sort(AWDiagonalMatrix<elem_type> &X, const char *dir) {
    X = arma::sort(X.diag(0), dir).eval();
  }

  template<typename elem_type>
  inline void _sort(AWRowVector<elem_type> &X, const char *dir) {
    X = arma::sort<arma::Row<elem_type> >(X, dir);
  }
  

  template<typename elem_type>
  inline void _sort(AWColVector<elem_type> &X, const char *dir) {
    X = arma::sort<arma::Col<elem_type> >(X, dir);
  }

  template<typename T> inline void SortAscending( T &X) { return _sort(X, "ascend");  }
  template<typename T> inline void SortDescending(T &X) { return _sort(X, "descend"); }
}

#endif /* __FUNCTION_SORT_HPP__ */
