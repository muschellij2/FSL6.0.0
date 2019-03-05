#ifndef __SUBVIEW_HPP__
#define __SUBVIEW_HPP__

/*
 * Sub-matrix views. Controlling assignent to subviews of special matrix
 * types (e.g. a subview of a diagonal matrix may only modify values along
 * the diagonal of the parent matrix).
 */

namespace armawrap {

  /*
   * This function is used (well, it should be used) by all pieces of code
   * which create AWSubView instances. Newmat supports zero-sized submatrices
   * for convenience. Armadillo does too, but it will throw errors on zero-
   * sized submatrices of different shapes .... yeah ... For example, you will
   * get an error if you try to assign a (0x5) subview to a (3x0) subview.
   *
   * This function makes sure that if a zero-sized subview has been requested,
   * its size is fixed to (0x0).
   */
  void FixSubviewIndices(int &first_row, int &last_row, int &first_col, int &last_col) {
    if      (last_row == first_row - 1) { last_col = first_col - 1; }
    else if (last_col == first_col - 1) { last_row = first_row - 1; }
  }


  /*
   * The diag_assign and sym_assign functions are used by some of the
   * subview_assign functions, for copying from the diagonal of a
   * matrix/expression to a the diagonal of a subview, or for assigning a scalar
   * value to a subview diagonal.
   */

  // assignment to a matrix/expression
  template<typename elem_type, typename wT, typename T>
  inline void diag_assign(AWSubView<elem_type, wT> &sv, const T &val, const int diag) {
    sv.diag(diag) = arma::diagvec(val, diag);
  }

  // assignment to a scalar
  template<typename elem_type, typename wT>
  inline void diag_assign(AWSubView<elem_type, wT> &sv, const elem_type val, const int diag) {
    sv.diag(diag).fill(val);
  }

  template<typename elem_type, typename wT, typename T>
  inline void sym_assign(AWSubView<elem_type, wT> &sv, const T &val) {
    sv.get_at_ref().operator=(val);
    sv.target.submat(sv.first_col - 1,
                     sv.first_row - 1,
                     sv.last_col  - 1,
                     sv.last_row  - 1) = val.t();
  }

  template<typename elem_type, typename wT>
  inline void sym_assign(AWSubView<elem_type, wT> &sv, const elem_type val) {
    sv.fill(val);
    sv.target.submat(sv.first_col - 1,
                     sv.first_row - 1,
                     sv.last_col  - 1,
                     sv.last_row  - 1).fill(val);
  }


  /*
   * The easy ones - AWMatrix/AWRowVector/AWColVector.
   */


  template<typename elem_type, typename T>
  inline void subview_assign(AWSubView<elem_type, AWMatrix<elem_type> > &sv, const arma::Base<elem_type, T> &val) {
    return sv.get_at_ref().operator=(val);
  }

  template<typename elem_type>
  inline void subview_assign(AWSubView<elem_type, AWMatrix<elem_type> > &sv, const elem_type val) {
    return sv.get_at_ref().fill(val);
  }

  template<typename elem_type, typename T>
  inline void subview_assign(AWSubView<elem_type, AWRowVector<elem_type> > &sv, const arma::Base<elem_type, T> &val) {
    return sv.get_at_ref().operator=(val);
  }

  template<typename elem_type>
  inline void subview_assign(AWSubView<elem_type, AWRowVector<elem_type> > &sv, const elem_type val) {
    return sv.get_at_ref().fill(val);
  }

  template<typename elem_type, typename T>
  inline void subview_assign(AWSubView<elem_type, AWColVector<elem_type> > &sv, const arma::Base<elem_type, T> &val) {
    return sv.get_at_ref().operator=(val);
  }

  template<typename elem_type>
  inline void subview_assign(AWSubView<elem_type, AWColVector<elem_type> > &sv, const elem_type val) {
    return sv.get_at_ref().fill(val);
  }


  /*
   * AWIdentityMatrix
   */


  template<typename elem_type, typename T>
  inline void subview_assign(AWSubView<elem_type, AWIdentityMatrix<elem_type> > &sv, const T &val)  {
    throw AWException("AWIdentityMatrix subviews are read-only");
  }


  /*
   * AWDiagonalMatrix
   */


  template<typename elem_type, typename T>
  inline void subview_assign(AWSubView<elem_type, AWDiagonalMatrix<elem_type> > &sv, const T &val)  {

    // complain if the subview does
    // not overlap with the main diagonal.
    if ((sv.last_row < sv.first_col) ||
        (sv.last_col < sv.first_row))
      throw AWException("Invalid subview on AWDiagonalMatrix");

    diag_assign(sv, val, sv.first_row - sv.first_col);
  }


  /*
   * AWUpperTriangularMatrix
   */


  template<typename elem_type, typename T>
  inline void subview_assign(AWSubView<elem_type, AWUpperTriangularMatrix<elem_type> > &sv, const T &val)  {

    // complain if the subview does
    // not overlap with the upper triangle
    if (sv.first_row > sv.last_col)
      throw AWException("Invalid subview on AWUpperTriangularMatrix");

    // if the subview is wholly within the
    // upper triangle, it's a straight copy
    if (sv.last_row <= sv.first_col)
      sv.get_at_ref().operator=(val);

    // otherwise copy one diagonal at a time,
    // ignoring everything below the main
    // diagonal of the original matrix.
    // There's probably a better way to do this.
    else {
      int first_diag = sv.first_row - sv.first_col;
      int last_diag  = sv.last_col  - sv.first_col;
      for (int diag = first_diag; diag <= last_diag; diag++) {
        diag_assign(sv, val, diag);
      }
    }
  }


  /*
   * AWLowerTriangularMatrix
   */


  template<typename elem_type, typename T>
  inline void subview_assign(AWSubView<elem_type, AWLowerTriangularMatrix<elem_type> > &sv, const T &val) {

    // See comments in the AWUpperTriangularMatrix
    // implementation above for an explanation
    if (sv.last_row < sv.first_col)
      throw AWException("Invalid subview on AWLowerTriangularMatrix");

    if (sv.last_col <= sv.first_row)
      sv.get_at_ref().operator=(val);

    else {
      int first_diag = sv.first_row - sv.last_row;
      int last_diag  = sv.first_row - sv.first_col;
      for (int diag = first_diag; diag <= last_diag; diag++) {
        diag_assign(sv, val, diag);
      }
    }
  }


  template<typename elem_type, typename T>
  inline void subview_assign(AWSubView<elem_type, AWSymmetricMatrix<elem_type> > &sv, const T &val) {

    // If the subivew lies wholly on one
    // side of the diagonal, assign to
    // both sides of the diagonal
    if ((sv.last_row <= sv.first_col) ||
        (sv.last_col <= sv.first_row)) {
      sym_assign(sv, val);
    }

    // If the subview overlaps with the main
    // diagonal, it's a bit more complicated
    else {

      // Newmat stores symmetric matrices as
      // lower triangular matrices, so we
      // give these values precedence. Copy
      // all values from the bottom left up
      // to the main diagonal of the original
      // matrix.
      int first_diag = sv.first_row  - sv.last_row;
      int main_diag  = sv.first_row - sv.first_col;

      for (int diag = first_diag; diag <= main_diag; diag++) {
        diag_assign(sv, val, diag);
      }

      // Then, using a new subview with a main
      // diagonal that overlaps with the main
      // diagonal of the original matrix, copy
      // values from the lower triangle to the
      // upper triangle of this new subview.
      arma::subview<elem_type> msv = sv.target.submat(sv.first_col - 1,
                                                      sv.first_col - 1,
                                                      sv.last_row  - 1,
                                                      sv.last_row  - 1);

      for (unsigned int diag = 1; diag <= (msv.n_rows - 1); diag++) {
        msv.diag(diag) = msv.diag(-diag);
      }
    }
  }

  /*
   * Functions which return the number of elements in a subview
   * of different matrix types.
   */

  /*
   * Default for simple types (e.g.
   * Matrix, RowVector, ColumnVector)
   */

  template<typename elem_type, typename wT>
  inline int subview_storage(const AWSubView<elem_type, wT> &sv) {
    return sv.AWBase<elem_type,
                     AWSubView<elem_type, wT>,
                     arma::subview<elem_type> >::Storage();
  }


  /*
   * Helper functions for calculating storage
   * of subviews on triangular matrices
   */

  inline int _lt_storage(int fr, int lr, int fc, int lc) {

    int storage = 0;
    int ncols;
    int row;

    // For each row, count the number of
    // columns which are before or on the
    // diagonal
    for (row = fr; row <= lr; row++) {

      ncols    = (lc - fc + 1)         -
                 std::max(0, lc - row) -
                 std::max(0, fc - row);
      storage += std::max(0, ncols);
    }
    return storage;
  }

  inline int _ut_storage(int sz, int fr, int lr, int fc, int lc) {

    int storage = 0;
    int ncols;
    int row;

    // For each row, count the number of
    // on or after the diagonal
    for (row = fr; row <= lr; row++) {

      ncols    = std::max(0, 1 + lc - row) -
                 std::max(0,     fc - row);
      storage += std::max(0, ncols);
    }

    return storage;
  }


  /*
   * Storage sizes of symmetric/triangular matrix sub-views.
   */


  template<typename elem_type>
  inline int subview_storage(
    const AWSubView<elem_type, AWSymmetricMatrix<elem_type> > &sv) {
    return _lt_storage(sv.first_row, sv.last_row, sv.first_col, sv.last_col);
  }


  template<typename elem_type>
  inline int subview_storage(
    const AWSubView<elem_type, AWLowerTriangularMatrix<elem_type> > &sv) {
    return _lt_storage(sv.first_row, sv.last_row, sv.first_col, sv.last_col);
  }


  template<typename elem_type>
  inline int subview_storage(
    const AWSubView<elem_type, AWUpperTriangularMatrix<elem_type> > &sv) {
    return _ut_storage(
      sv.target.Nrows(), sv.first_row, sv.last_row, sv.first_col, sv.last_col);
  }

  /*
   * Band matrices are unimplemented at the moment
   * (as there may not be any need for them).
   */

  /*
   * The AWSubView class
   */


  template<typename elem_type, typename wT>
  class AWSubView : public arma::subview<elem_type>,
                    public armawrap::AWBase<elem_type,
                                            AWSubView<elem_type, wT>,
                                            arma::subview<elem_type> > {

  public:

    using armawrap::AWBase<elem_type,
                           AWSubView<elem_type, wT>,
                           arma::subview<elem_type> >::operator<<;

    using armawrap::AWBase<elem_type,
                           AWSubView<elem_type, wT>,
                           arma::subview<elem_type> >::operator();

    using armawrap::AWBase<elem_type,
                           AWSubView<elem_type, wT>,
                           arma::subview<elem_type> >::t;

    wT &      target;
    const int first_row;
    const int last_row;
    const int first_col;
    const int last_col;

    /* If not NULL, this pointer is deleted in the destructor */
    const wT *ptr;

    /*
     * An AWSubView object is returned by all sub-matrix methods - see base.hpp.
     */
    inline AWSubView(const wT & target,
                     const int  first_row,
                     const int  last_row,
                     const int  first_col,
                     const int  last_col) :
      arma::subview<elem_type>(target,
                               first_row - 1,
                               first_col - 1,
                               last_row  - first_row + 1,
                               last_col  - first_col + 1),
      target(const_cast<wT&>(target)),
      first_row(first_row),
      last_row( last_row),
      first_col(first_col),
      last_col( last_col),
      ptr(NULL) { }

    /*
     * This constructor is used when a subview of an expression is
     * accessed - see exprs.hpp.
     */
    inline AWSubView(const wT & target,
                     const int  first_row,
                     const int  last_row,
                     const int  first_col,
                     const int  last_col,
                     const wT * ptr) :
      arma::subview<elem_type>(target,
                               first_row - 1,
                               first_col - 1,
                               last_row  - first_row + 1,
                               last_col  - first_col + 1),
      target(const_cast<wT&>(target)),
      first_row(first_row),
      last_row( last_row),
      first_col(first_col),
      last_col( last_col),
      ptr(ptr) { }

    inline ~AWSubView() {
      if (ptr != NULL) {
        delete ptr;
      }
    }

    /*
     * AWMatrix a(10, 10);
     * a.SubMatrix(4, 6, 4, 6) = (some expression);
     */
    template<typename T>
    inline const AWSubView &
    operator=(const arma::Base<elem_type, T> &X) {
      subview_assign(*this, X);
      return *this;
    }

    inline const AWSubView &
    operator=(const AWSubView &X) {
      subview_assign(*this, X);
      return *this;
    }

    /*
     * AWMatrix a(10, 10);
     * a.SubMatrix(4, 6, 4, 6) = (some scalar);
     */
    inline const AWSubView &
    operator=(const elem_type k) {
      subview_assign(*this, k);
      return *this;
    }

    inline void operator+=(const elem_type k) { subview_assign(*this, (*this) + k); }
    inline void operator-=(const elem_type k) { subview_assign(*this, (*this) - k); }
    inline void operator/=(const elem_type k) { subview_assign(*this, (*this) / k); }
    inline void operator*=(const elem_type k) { subview_assign(*this, (*this) * k); }

    template<typename T>
    inline void operator+=(const arma::Base<elem_type, T> &X) { subview_assign(*this, (*this) + X.get_ref()); }

    template<typename T>
    inline void operator-=(const arma::Base<elem_type, T> &X) { subview_assign(*this, (*this) - X.get_ref()); }

    template<typename T>
    inline void operator/=(const arma::Base<elem_type, T> &X) { subview_assign(*this, (*this) / X.get_ref()); }


    inline elem_type operator()(const int i) const {
      // This code does the following:
      //
      //  1. Take the 1D index, turn it into
      //     row/column indices into the
      //     subview
      //
      //  2. Offset the row/columns so they
      //     index the target matrix
      int svRow     = (i - 1) / this->Ncols() + 1;
      int svCol     = (i - 1) % this->Ncols() + 1;
      int targetRow = svRow + first_row - 1;
      int targetCol = svCol + first_col - 1;
      return target.operator()(targetRow, targetCol);
    }

    inline AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value> operator()(const int i) {
      int svRow     = (i - 1) / this->Ncols() + 1;
      int svCol     = (i - 1) % this->Ncols() + 1;
      int targetRow = svRow + first_row - 1;
      int targetCol = svCol + first_col - 1;
      return target.operator()(targetRow, targetCol);
    }

    inline elem_type operator()(const int row, const int col) const {
      return target.operator()(row + first_row, col + first_col);
    }

    inline AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value> operator()(const int row, const int col) {
      return target.operator()(row + first_row, col + first_col);
    }


    inline AWSubView<elem_type, wT> SubMatrix(const int first_row,
                                              const int last_row,
                                              const int first_col,
                                              const int last_col) {

      int fr = first_row + this->first_row - 1;
      int lr = last_row  + this->first_row - 1;;
      int fc = first_col + this->first_col - 1;
      int lc = last_col  + this->first_col - 1;;

      FixSubviewIndices(fr, lr, fc, lc);

      return AWSubView<elem_type, wT>(target, fr, lr, fc, lc);
    }

    inline AWSubView<elem_type, wT> Rows(const int first, const int last) {
      return SubMatrix(first, last, 1, this->last_col - this->first_col + 1);
    }


    inline AWSubView<elem_type, wT> SymSubMatrix(const int first, const int last) {
      return SubMatrix(first, last, first, last);
    }

    inline AWSubView<elem_type, wT> Columns(const int first, const int last) {
      return SubMatrix(1, this->last_row - this->first_row + 1, first, last);
    }

    inline AWSubView<elem_type, wT> Row(const int row) {
      return SubMatrix(row, row, 1, this->last_col - this->first_col + 1);
    }

    inline AWSubView<elem_type, wT> Column(const int col) {
      return SubMatrix(1, this->last_row - this->first_row + 1, col, col);
    }

    inline virtual int Storage() const {
      return subview_storage(*this);
    }
  };
}

#endif /* __SUBVIEW_HPP__ */
