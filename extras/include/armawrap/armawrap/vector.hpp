#ifndef __VECTOR_HPP__
#define __VECTOR_HPP__


namespace armawrap {

  /*
   * Row and column vectors. The AWColVector and AWRowVector classes 
   * inherit from AWMatrix, rather than from arma::Row/Col. This is 
   * because, in Newmat, the Row and ColumnVectors inherit from Matrix. 
   * So a Row/ColumnVector must be able to be accessed as if it is a 
   * Matrix.
   */
  
  template<typename elem_type>
  class AWColVector : public AWMatrix<elem_type> {

  public:

    // Make sure that the AWBase::ReSize methods are
    // also available (otherwise they will be hidden
    // by the ReSize declaration at the bottom)
    using armawrap::AWBase<elem_type,
                           AWMatrix<elem_type>,
                           arma::Mat<elem_type> >::ReSize; 
 

    // AWColVector a;
    inline AWColVector() : AWMatrix<elem_type>() {}

    // AWColVector a(100);
    // This is marked explicit, to prevent implicit conversions from int to vector
    // If 0 is passed in for nrows, we force the underlying arma::mat to size (0x0),
    // instead of being size (0x1)
    explicit inline AWColVector(const int nrows) :
      AWMatrix<elem_type>(nrows, nrows > 0 ? 1 : 0) {}

    // AWColVector a(some expression);
    template<typename T>
    inline AWColVector(const arma::Base<elem_type, T> &X) : AWMatrix<elem_type>(X) { }

    // AWColVector a;
    // a = (some expression); 
    template<typename T> 
    inline const AWColVector & 
    operator=(const arma::Base<elem_type, T> &X) {
      AWMatrix<elem_type>::operator=(X);
      return *this;
    }

    // AWColVector a;
    // a = (some scalar);
    inline const AWColVector &
    operator=(const elem_type k) {
      AWMatrix<elem_type>::operator=(k);
      return *this;
    }

    inline void ReSize(int nelem) {
      arma::Mat<elem_type>::set_size(nelem, 1);
    } 
  };

  template<typename elem_type>
  class AWRowVector : public AWMatrix<elem_type> {

  public:

    // Make sure that the AWBase::ReSize methods are
    // also available (otherwise they will be hidden
    // by the ReSize declaration at the bottom)
    using armawrap::AWBase<elem_type,
                           AWMatrix<elem_type>,
                           arma::Mat<elem_type> >::ReSize; 

    // AWRowVector a;
    inline AWRowVector() : AWMatrix<elem_type>() {}

    // AWRowVector a(100);
    // This is marked explicit, to prevent implicit conversions from int to vector
    // If 0 is passed in for ncols, we force the underlying arma::mat to size (0x0),
    // instead of being size (1x0) 
    explicit inline AWRowVector(const int ncols) :
      AWMatrix<elem_type>(ncols > 0 ? 1 : 0, ncols) {}

    // AWRowVector a(some expression);
    template<typename T>
    inline AWRowVector(const arma::Base<elem_type, T> &X) : AWMatrix<elem_type>(X) { }

    // AWRowVector a;
    // a = (some expression); 
    template<typename T> 
    inline const AWRowVector & 
    operator=(const arma::Base<elem_type, T> &X) {
      AWMatrix<elem_type>::operator=(X);
      return *this;
    }

    // AWRowVector a;
    // a = (some scalar);
    inline const AWRowVector &
    operator=(const elem_type k) {
      AWMatrix<elem_type>::operator=(k);
      return *this;
    }

    inline void ReSize(int nelem) {
      arma::Mat<elem_type>::set_size(1, nelem);
    }
  };  
}

#endif /* __VECTOR_HPP__*/
