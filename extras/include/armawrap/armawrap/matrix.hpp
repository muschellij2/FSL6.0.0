#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

namespace armawrap {
  template<typename elem_type>
  class AWMatrix : public arma::Mat<elem_type>, 
                   public armawrap::AWBase<elem_type, 
                                           AWMatrix<elem_type>,
                                           arma::Mat<elem_type> > {

  public:
  
    /*
     * Prevent multiple inheritance ambiguity:
     * use armawrap::AWBase::methods, and 
     * not arma::Mat::methods.
     */
    using armawrap::AWBase<elem_type, 
                           AWMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator<<;

    using armawrap::AWBase<elem_type, 
                           AWMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();

    using armawrap::AWBase<elem_type, 
                           AWMatrix<elem_type>, 
                           arma::Mat<elem_type> >::t;

    /* 
     * Destructor - marked virtual to suppress compiler warnings.
     */
    inline virtual ~AWMatrix() { }

    /*
     * Constructors
     */

    /* AWMatrix a; */
    inline AWMatrix() : arma::Mat<elem_type>() {}

    /* AWMatrix a(10, 20); */
    inline AWMatrix(const int nrows, const int ncols) : arma::Mat<elem_type>(nrows, ncols) {}

    /* AWMatrix a(some expression); */
    template<typename T>
    inline AWMatrix(const arma::Base<elem_type, T>& X) : arma::Mat<elem_type>() {
      operator=(X);
    }

    /*
     * AWMatrix a;
     * a = (some expression); 
     */
    template<typename T> 
    inline const AWMatrix & 
    operator=(const arma::Base<elem_type, T> &X) {
      arma::Mat<elem_type>::operator=(X.get_ref());
      return *this;
    }

    /*
     * AWMatrix a;
     * a = (some scalar);
     * In armadillo, this resizes the matrix.
     * In NewMat, this sets all values in the matrix to the scalar.
     */
    inline const AWMatrix &
    operator=(const elem_type k) {
      arma::Mat<elem_type>::fill(k);
      return *this;
    }
  };
}

#endif /* __MATRIX_HPP__ */
