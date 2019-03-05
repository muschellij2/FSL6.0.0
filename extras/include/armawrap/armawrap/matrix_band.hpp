#ifndef __MATRIX_BAND_HPP__
#define __MATRIX_BAND_HPP__

/*
 * BandMatrix, UpperBandMatrix, LowerBandMatrix, and SymmetricBandMatrix
 * classes.
 */

/*
 * The AWMatrixBandWidth class encapsulates the (lower, upper) band widths
 * of a band matrix. Copied more or less verbatim from the newmat source.
 */
namespace armawrap {

  class AWMatrixBandWidth {
  public:
    int lower;
    int upper;

    inline AWMatrixBandWidth(const int l, const int u) : lower(l), upper(u) { }
    inline AWMatrixBandWidth(const int i)              : lower(i), upper(i) { }

    inline virtual ~AWMatrixBandWidth() { }

    inline AWMatrixBandWidth operator+(const AWMatrixBandWidth &bw) const {
      int l = bw.lower; 
      int u = bw.upper;
      l = (lower < 0 || l < 0) ? -1 : (lower > l) ? lower : l;
      u = (upper < 0 || u < 0) ? -1 : (upper > u) ? upper : u;
      return AWMatrixBandWidth(l, u);
    }

    inline AWMatrixBandWidth operator*(const AWMatrixBandWidth &bw) const {
      int l = bw.lower; 
      int u = bw.upper;
      l = (lower < 0 || l < 0) ? -1 : lower + l;
      u = (upper < 0 || u < 0) ? -1 : upper + u;
      return AWMatrixBandWidth(l, u);
    }

    inline AWMatrixBandWidth minimum(const AWMatrixBandWidth &bw) const {
      int l = bw.lower; 
      int u = bw.upper;
      if ((lower >= 0) && ( (l < 0) || (l > lower) )) l = lower;
      if ((upper >= 0) && ( (u < 0) || (u > upper) )) u = upper;
      return AWMatrixBandWidth(l, u);
    }

    inline AWMatrixBandWidth t() const { 
      return AWMatrixBandWidth(upper,lower); 
    }

    inline bool operator==(const AWMatrixBandWidth &bw) const { 
      return (lower == bw.lower) && (upper == bw.upper); 
    }

    inline bool operator!=(const AWMatrixBandWidth &bw) const { 
      return !operator==(bw); 
    }

    inline int Upper() const { return upper; }
    inline int Lower() const { return lower; }
  };


  template<typename elem_type> 
  class AWBandMatrix : public arma::Mat<elem_type>, 
                       public armawrap::AWBase<elem_type, 
                                               AWBandMatrix<elem_type>,
                                               arma::Mat<elem_type> >{ 

  protected:
    int lower;
    int upper;

  public:

    using armawrap::AWBase<elem_type, 
                           AWBandMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator<<;

    using armawrap::AWBase<elem_type, 
                           AWBandMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();

    using armawrap::AWBase<elem_type, 
                           AWBandMatrix<elem_type>, 
                           arma::Mat<elem_type> >::t;

    using armawrap::AWBase<elem_type, 
                           AWBandMatrix<elem_type>, 
                           arma::Mat<elem_type> >::ReSize; 

    /* AWBandMatrix i; */
    inline AWBandMatrix() : arma::Mat<elem_type>() {
      arma::Mat<elem_type>::zeros();
    }

    inline virtual ~AWBandMatrix() { }

    /*
     * AWBandMatrix i(100, 10, 20);
     */
    inline AWBandMatrix(const int nelems, 
                        const int lower, 
                        const int upper) : 
      arma::Mat<elem_type>(nelems, nelems),
      lower(lower), upper(upper) {
      arma::Mat<elem_type>::zeros();
    }

    /*
     * AWBandMatrix i(some expression);
     */
    template<typename T>
    inline AWBandMatrix(const arma::Base<elem_type, T> &X) : arma::Mat<elem_type>() {
      operator=(X);
    }

    /*
     * AWBandMatrix i(100, 10, 20);
     * i = (some expression);
     */
    template<typename T> 
    inline const AWBandMatrix & 
    operator=(const arma::Base<elem_type, T> &X) {

      arma::Mat<elem_type> exp = X.get_ref();

      // start at 0 so the main diagonal is copied
      for (int i = 0; i <= lower; i++) { arma::Mat<elem_type>::diag(-i) = exp.diag(-i); };
      for (int i = 1; i <= upper; i++) { arma::Mat<elem_type>::diag( i) = exp.diag( i); };

      return *this;
    }

    /*
     * AWBandMatrix i(100, 10, 20);
     * i = 99; // set all non-zero elements to 99.
     */
    inline const AWBandMatrix &
    operator=(const elem_type k) {

      for (int i = 0; i <= lower; i++) { arma::Mat<elem_type>::diag(-i).fill(k); };
      for (int i = 1; i <= upper; i++) { arma::Mat<elem_type>::diag( i).fill(k); }; 

      return *this;
    }

    inline AWMatrixBandWidth BandWidth() const {
      return AWMatrixBandWidth(lower, upper);
    }

    /* 
     * Newmat only stores the values on the bands for band 
     * matrices, but armadillo stores a full matrix.
     */
    inline virtual int Storage() const { 

      // However, newmat seems to over-allocate 
      // the storage required for a band matrix. 
      int nrows = arma::Mat<elem_type>::n_rows; 
      return nrows * (1 + lower + upper);
    }

    inline void ReSize(const int nelem) {
      ReSize(nelem, nelem);
    }

    inline void ReSize(const int nelem, const int lb, const int ub) {
      ReSize(nelem);
      lower = lb;
      upper = ub;
    }
  };

  template<typename elem_type>
  class AWUpperBandMatrix : public AWBandMatrix<elem_type> {
  public:

    /* AWUpperBandMatrix i; */
    inline AWUpperBandMatrix() : AWBandMatrix<elem_type>() { }

    /*
     * AWUpperBandMatrix i(100, 10);
     */
    inline AWUpperBandMatrix(const int nelems, const int upper) : 
      AWBandMatrix<elem_type>(nelems, 0, upper) { }
  };

  template<typename elem_type>
  class AWLowerBandMatrix : public AWBandMatrix<elem_type> {
  public:

    /* AWLowerBandMatrix i; */
    inline AWLowerBandMatrix() : AWBandMatrix<elem_type>() { }

    /*
     * AWLowerBandMatrix i(100, 10);
     */
    inline AWLowerBandMatrix(const int nelems, const int lower) : 
      AWBandMatrix<elem_type>(nelems, lower, 0) { } 
  };

  template<typename elem_type>
  class AWSymmetricBandMatrix : public AWBandMatrix<elem_type>,
                                public AWBase<elem_type,
                                              AWSymmetricBandMatrix<elem_type>,
                                              arma::Mat<elem_type> > {
  public:

    // We have a diamond inheritance issue: AWBandMatrix inherits
    // from AWBase<elem_type, AWBandMatrix, arma::Mat>, and AWSymmetricMatrix
    // inherits from AWBase<elem_type, AWSymmetricBandMatrix, arma::Mat>.
    // So we have to explicitly use the SymmtricBandMatrix implementations.
    
    using armawrap::AWBase<elem_type,
                           AWSymmetricBandMatrix<elem_type>,
                           arma::Mat<elem_type> >::Nrows;
    using armawrap::AWBase<elem_type,
                           AWSymmetricBandMatrix<elem_type>,
                           arma::Mat<elem_type> >::Ncols; 
    using armawrap::AWBase<elem_type, 
                           AWSymmetricBandMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();
    using armawrap::AWBase<elem_type, 
                           AWSymmetricBandMatrix<elem_type>, 
                           arma::Mat<elem_type> >::get_at_ref;


    AWSymmetricBandMatrix() : AWBandMatrix<elem_type>() { }

    inline virtual ~AWSymmetricBandMatrix() { }

    AWSymmetricBandMatrix(int nelem, int bandwidth) :
      AWBandMatrix<elem_type>(nelem, bandwidth, bandwidth) { }

    template<typename T>
    inline AWSymmetricBandMatrix(const arma::Base<elem_type, T> &X) : arma::Mat<elem_type>() {
      operator=(X);
    }

    /*
     * AWSymmetricBandMatrix i(100, 10);
     * i = (some expression);
     */
    template<typename T> 
    inline const AWSymmetricBandMatrix & 
    operator=(const arma::Base<elem_type, T> &X) {

      arma::Mat<elem_type> exp = X.get_ref();

      // start at 0 so the main diagonal is copied
      for (int i = 0; i <= this->lower; i++) { arma::Mat<elem_type>::diag(-i) = exp.diag(-i); };
      for (int i = 1; i <= this->upper; i++) { arma::Mat<elem_type>::diag( i) = exp.diag( i); };

      return *this;
    }

    /*
     * AWSymmetricBandMatrix i(100, 10);
     * i = 99; // set all non-zero elements to 99.
     */
    inline const AWSymmetricBandMatrix &
    operator=(const elem_type k) {

      for (int i = 0; i <= this->lower; i++) { arma::Mat<elem_type>::diag(-i).fill(k); };
      for (int i = 1; i <= this->upper; i++) { arma::Mat<elem_type>::diag( i).fill(k); }; 

      return *this;
    }

    inline void ReSize(const int nelem, const int b) {
      AWBandMatrix<elem_type>::ReSize(nelem, b, b);
    } 
  };
}

#endif /* __MATRIX_BAND_HPP__ */
