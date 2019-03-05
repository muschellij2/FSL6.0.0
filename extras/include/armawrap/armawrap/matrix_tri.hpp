#ifndef __MATRIX_TRI_HPP__
#define __MATRIX_TRI_HPP__

#include <cmath>

/*
 * Upper/Lower triangular and symmetric matrix types.
 */

namespace armawrap {

  template<typename elem_type>
  class AWUpperTriangularMatrix : public arma::Mat<elem_type>, 
                                  public armawrap::AWBase<elem_type, 
                                                          AWUpperTriangularMatrix<elem_type>,
                                                          arma::Mat<elem_type> >{ 
  public:

    using armawrap::AWBase<elem_type, 
                           AWUpperTriangularMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator<<;

    using armawrap::AWBase<elem_type, 
                           AWUpperTriangularMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();

    using armawrap::AWBase<elem_type, 
                           AWUpperTriangularMatrix<elem_type>, 
                           arma::Mat<elem_type> >::t;

    using armawrap::AWBase<elem_type, 
                           AWUpperTriangularMatrix<elem_type>, 
                           arma::Mat<elem_type> >::ReSize; 

    /*
     * AWUpperTriangularMatrix i;
     */
    inline AWUpperTriangularMatrix() : arma::Mat<elem_type>() {}

    inline virtual ~AWUpperTriangularMatrix() { }

    /*
     * AWUpperTriangularMatrix i(10);
     */
    inline AWUpperTriangularMatrix(const int nelems) : arma::Mat<elem_type>(nelems, nelems) {
      arma::Mat<elem_type>::zeros();
    }

    /*
     * AWUpperTriangularMatrix i(some expression);
     */
    template<typename T>
    inline AWUpperTriangularMatrix(const arma::Base<elem_type, T> &X) : arma::Mat<elem_type>() {
      operator=(X);
    }

    /*
     * AWUpperTriangularMatrix i;
     * i = (some expression);
     */  
    template<typename T> 
    inline const AWUpperTriangularMatrix & 
    operator=(const arma::Base<elem_type, T> &X) {
      arma::Mat<elem_type>::operator=(arma::trimatu(X.get_ref()));
      return *this;
    }

    /*
     * AWUpperTriangularMatrix i(100);
     * i = 99; // set all elements to 99.
     */
    inline const AWUpperTriangularMatrix &
    operator=(const elem_type k) {
      arma::Mat<elem_type>::zeros();
      operator=(arma::trimatu(arma::Mat<elem_type>(this->Nrows(), this->Nrows()).fill(k)));
      return *this;
    }

    /* 
     * Newmat only stores a half matrix for triangular
     * matrices, but armadillo stores a full matrix.
     */
    inline virtual int Storage() const { 
      return arma::Mat<elem_type>::n_rows * (arma::Mat<elem_type>::n_rows + 1) / 2.0;
    }

    inline void ReSize(int nelem) {
      ReSize(nelem, nelem);
    }
  };


  template<typename elem_type>
  class AWLowerTriangularMatrix : public arma::Mat<elem_type>, 
                                  public armawrap::AWBase<elem_type, 
                                                          AWLowerTriangularMatrix<elem_type>,
                                                          arma::Mat<elem_type> >{ 
  public:

    using armawrap::AWBase<elem_type, 
                           AWLowerTriangularMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator<<;

    using armawrap::AWBase<elem_type, 
                           AWLowerTriangularMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();

    using armawrap::AWBase<elem_type, 
                           AWLowerTriangularMatrix<elem_type>, 
                           arma::Mat<elem_type> >::t;

    using armawrap::AWBase<elem_type, 
                           AWLowerTriangularMatrix<elem_type>, 
                           arma::Mat<elem_type> >::ReSize; 

    /*
     * AWLowerTriangularMatrix i;
     */
    inline AWLowerTriangularMatrix() : arma::Mat<elem_type>() {}

    inline virtual ~AWLowerTriangularMatrix() { }

    /*
     * AWLowerTriangularMatrix i(10);
     */
    inline AWLowerTriangularMatrix(const int nelems) : arma::Mat<elem_type>(nelems, nelems) { }

    /*
     * AWLowerTriangularMatrix i(some expression);
     */
    template<typename T>
    inline AWLowerTriangularMatrix(const arma::Base<elem_type, T> &X) : arma::Mat<elem_type>() {
      operator=(X);
    }

    /*
     * AWLowerTriangularMatrix i;
     * i = (some expression);
     */  
    template<typename T> 
    inline const AWLowerTriangularMatrix & 
    operator=(const arma::Base<elem_type, T> &X) {
      arma::Mat<elem_type>::operator=(arma::trimatl(X.get_ref()));
      return *this;
    }

    /*
     * AWLowerTriangularMatrix i(100);
     * i = 99; // set all elements to 99.
     */
    inline const AWLowerTriangularMatrix &
    operator=(const elem_type k) {
      arma::Mat<elem_type>::zeros();
      operator=(arma::trimatl(arma::Mat<elem_type>(this->Nrows(), this->Nrows()).fill(k)));
      return *this;
    }

    inline virtual int Storage() const { 
      return arma::Mat<elem_type>::n_rows * (arma::Mat<elem_type>::n_rows + 1) / 2.0;
    }
 
    inline void ReSize(int nelem) {
      ReSize(nelem, nelem);
    }
  };


  template<typename elem_type>
  class AWSymmetricMatrix : public arma::Mat<elem_type>, 
                            public armawrap::AWBase<elem_type, 
                                                    AWSymmetricMatrix<elem_type>,
                                                    arma::Mat<elem_type> >{ 
  public:

    using armawrap::AWBase<elem_type, 
                           AWSymmetricMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator<<;

    using armawrap::AWBase<elem_type, 
                           AWSymmetricMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();

    using armawrap::AWBase<elem_type, 
                           AWSymmetricMatrix<elem_type>, 
                           arma::Mat<elem_type> >::t;

    using armawrap::AWBase<elem_type,
                           AWSymmetricMatrix<elem_type>,
                           arma::Mat<elem_type> >::ReSize; 
 
    /*
     * AWSymmetricMatrix i;
     */
    inline AWSymmetricMatrix() : arma::Mat<elem_type>() {}

    inline virtual ~AWSymmetricMatrix() { }

    /*
     * AWSymmetricMatrix i(10);
     */
    inline AWSymmetricMatrix(const int nelems) : arma::Mat<elem_type>(nelems, nelems) { }

    /* AWSymmetricMatrix a(some expression); */
    template<typename T>
    inline AWSymmetricMatrix(const arma::Base<elem_type, T> &X) : arma::Mat<elem_type>() {
      operator=(X);
    }

    /*
     * AWSymmetricMatrix i;
     * i = (some expression);
     */  
    template<typename T> 
    inline const AWSymmetricMatrix & 
    operator=(const arma::Base<elem_type, T> &X) {
      arma::Mat<elem_type>::operator=(arma::symmatl(X.get_ref()));
      return *this;
    }

    /*
     * AWSymmetricMatrix i(100);
     * i = 99; // set all elements to 99.
     */
    inline const AWSymmetricMatrix &
    operator=(const elem_type k) {
      arma::Mat<elem_type>::fill(k);
      return *this;
    }

    inline virtual int Storage() const { 
      return arma::Mat<elem_type>::n_rows * (arma::Mat<elem_type>::n_rows + 1) / 2.0;
    }

    inline void ReSize(int nelem) {
      ReSize(nelem, nelem);
    }
  };
}

#endif /* __MATRIX_TRI_HPP__ */
