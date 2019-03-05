/**
 * The AWBase class is a mixin of sorts, which provides a bunch of methods that
 * are common to the AWMatrix, AWRowVector, and AWColVector classes.
 */
#ifndef __BASE_HPP__
#define __BASE_HPP__

#include <math.h>

namespace armawrap {

  template<typename elem_type, /* data type (e.g. double)      */
           typename wT,        /* armawrap type                */
           typename aT>        /* corresponding armadillo type */

  class AWBase {

  private:

    /*
     * The horrible things below have the effect that calls to eval(*this)
     * on Matrix types will just return a reference to the matrix, but
     * calls to eval(*this) on expressions will force evaluation.
     */

    static const bool isMat = arma::is_Mat<aT>::value;

    template<typename this_type>
    inline typename enable_if< this_type::isMat, const aT &>::type
    eval(const this_type &self) const { return get_at_ref(); }

    template<typename this_type>
    inline typename enable_if< this_type::isMat, aT &>::type
    eval(this_type &self) { return get_at_ref(); }

    template<typename this_type>
    inline typename enable_if<!this_type::isMat, const arma::Mat<elem_type> >::type
    eval(const this_type &self) const { return get_at_ref().eval(); }

    template<typename this_type>
    inline typename enable_if<!this_type::isMat, arma::Mat<elem_type> >::type
    eval(this_type &self) { return get_at_ref().eval(); }

  public:

    inline const wT & get_wt_ref() const { return static_cast<const wT&>(*this);        }
    inline const aT & get_at_ref() const { return static_cast<const aT&>(get_wt_ref()); }
    inline wT & get_wt_ref()             { return static_cast<      wT&>(*this);        }
    inline aT & get_at_ref()             { return static_cast<      aT&>(get_wt_ref()); }

    /*
     * The Storage() method may be overridden by subclasses, as there
     * are type-specific differences in storage mechanisms between
     * armadillo and newmat.
     */
    inline int         Nrows()   const { return (int)eval(*this).n_rows; }
    inline int         Ncols()   const { return (int)eval(*this).n_cols; }
    virtual inline int Storage() const { return (int)eval(*this).n_elem; }


    // Code which uses Store() is going to break,
    // because Newmat stores matrices row-wise,
    // whereas Armadillo stores matrices column-wise.
    inline elem_type * Store()   const { return const_cast<elem_type *>(get_at_ref().memptr()); }


    /*
     * The real implementations for these functions are
     * in functions_base.hpp, as standalone functions.
     * These methods are defined as virtual so sub-classes
     * may provide type-specific implementations - see e.g.
     * the DiagonalMatrix class.
     */
    inline virtual elem_type               Minimum()              const { return armawrap::Minimum(             eval(*this)); }
    inline virtual elem_type               Maximum()              const { return armawrap::Maximum(             eval(*this)); }
    inline virtual elem_type               MinimumAbsoluteValue() const { return armawrap::MinimumAbsoluteValue(eval(*this)); }
    inline virtual elem_type               MaximumAbsoluteValue() const { return armawrap::MaximumAbsoluteValue(eval(*this)); }
    inline virtual elem_type               SumSquare()            const { return armawrap::SumSquare(           eval(*this)); }
    inline virtual elem_type               SumAbsoluteValue()     const { return armawrap::SumAbsoluteValue(    eval(*this)); }
    inline virtual elem_type               Sum()                  const { return armawrap::Sum(                 eval(*this)); }
    inline virtual elem_type               Trace()                const { return armawrap::Trace(               eval(*this)); }
    inline virtual elem_type               Norm1()                const { return armawrap::Norm1(               eval(*this)); }
    inline virtual elem_type               NormInfinity()         const { return armawrap::NormInfinity(        eval(*this)); }
    inline virtual elem_type               NormFrobenius()        const { return armawrap::NormFrobenius(       eval(*this)); }
    inline virtual elem_type               Determinant()          const { return armawrap::Determinant(         eval(*this)); }
    inline virtual AWLogAndSign<elem_type> LogDeterminant()       const { return armawrap::LogDeterminant(      eval(*this)); }
    inline virtual bool                    IsZero()               const { return armawrap::IsZero(              eval(*this)); }

    inline elem_type Minimum1(int &i) const {

      elem_type val;
      int row, col;
      val = Minimum2(row, col);
      i = (row - 1) * Ncols() + col;
      return val;
    }

    inline elem_type Maximum1(int &i) const {
      elem_type val;
      int row, col;
      val = Maximum2(row, col);
      i = (row - 1) * Ncols() + col;
      return val;
    }

    inline elem_type MinimumAbsoluteValue1(int &i) const {
      elem_type val;
      int row, col;
      val = MinimumAbsoluteValue2(row, col);
      i = (row - 1) * Ncols() + col;
      return val;
    }

    inline elem_type MaximumAbsoluteValue1(int &i) const {
      elem_type val;
      int row, col;
      val = MaximumAbsoluteValue2(row, col);
      i = (row - 1) * Ncols() + col;
      return val;
    }

    inline elem_type Minimum2(int &row, int &col) const {
      arma::uword urow, ucol;
      elem_type val = get_at_ref().eval().min(urow, ucol);
      row = urow + 1;
      col = ucol + 1;
      return val;
    }

    inline elem_type Maximum2(int &row, int &col) const {
      arma::uword urow, ucol;
      elem_type val = get_at_ref().eval().max(urow, ucol);
      row = urow + 1;
      col = ucol + 1;
      return val;
    }

    inline elem_type MinimumAbsoluteValue2(int &row, int &col) const {
      arma::uword urow, ucol;
      elem_type val = arma::abs<aT>(get_at_ref()).eval().min(urow, ucol);
      row = urow + 1;
      col = ucol + 1;
      return val;
    }

    inline elem_type MaximumAbsoluteValue2(int &row, int &col) const {
      arma::uword urow, ucol;
      elem_type val = arma::abs<aT>(get_at_ref()).eval().max(urow, ucol);
      row = urow + 1;
      col = ucol + 1;
      return val;
    }

    /*
     * Transpose
     */
    inline AWOp<typename armawrap_type_map<aT>::type, arma::op_htrans> t() const {
      return AWOp<typename armawrap_type_map<aT>::type, arma::op_htrans>(get_at_ref());
    }

    inline AWOp<typename armawrap_type_map<aT>::type, arma::op_htrans> Transpose() const {
      return AWOp<typename armawrap_type_map<aT>::type, arma::op_htrans>(get_at_ref());
    }


    /*
     * Horizontal/vertical concatenation (see also operator_ampersand.hpp and
     * operator_pipe.hpp).
     */

    template<typename T>
    inline void operator&=(const arma::Base<elem_type, T> &X) {
      get_at_ref().insert_rows(Nrows(), X);
    }

    template<typename T>
    inline void operator|=(const arma::Base<elem_type, T> &X) {
      get_at_ref().insert_cols(Ncols(), X);
    }


    /*
     * Reverse a matrix/vector
     */
    inline AWOp<arma::Op<aT, arma::op_flipud>, arma::op_fliplr> Reverse() const {

      // The inner operation is created on the heap
      // so it won't be deleted when we leave this
      // scope. It is deleted by the AWOp object.
      arma::Op<aT, arma::op_flipud> *inner =
        new arma::Op<aT, arma::op_flipud>(get_at_ref());

      return AWOp<arma::Op<aT, arma::op_flipud>, arma::op_fliplr>(true, *inner);
    }

    /*
     * Matrix value access/assignment - see operator_call.hpp.
     */

    inline elem_type operator()(const int i) const {
      return AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value>(get_wt_ref(), i);
    }

    inline elem_type operator()(const int row, const int col) const {
      return AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value>(get_wt_ref(), row, col);
    }

    inline AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value> operator()(const int i) {
      return AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value>(get_wt_ref(), i);
    }

    inline AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value> operator()(const int row, const int col) {
      return AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value>(get_wt_ref(), row, col);
    }

    /*
     * 0-indexed matrix value access/assignment - again, see operator_call.hpp.
     */

    inline elem_type element(const int i) const {
      return operator()(i + 1);
    }

    inline elem_type element(const int row, const int col) const {
      return operator()(row + 1, col + 1);
    }

    inline AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value> element(const int i) {
      return operator()(i + 1);
    }

    inline AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value> element(const int row, const int col) {
      return operator()(row + 1, col + 1);
    }

    /*
     * Initialising matrix/row/column values via the << operator - see
     * operator_insert.hpp
     */
    inline AWInsertManager<elem_type, wT> operator<<(const elem_type val) {
      return AWInsertManager<elem_type, wT>(get_wt_ref(), val);
    }

    inline AWInsertManager<elem_type, wT> operator<<(const int val) {
      return AWInsertManager<elem_type, wT>(get_wt_ref(), val);
    }

    /*
     * Inserting via a pointer or array
     */
    inline void operator<<(const elem_type* vals) {

      AWInsertManager<elem_type, wT> ins(get_wt_ref(), vals[0]);
      for (int i = 1; i < Storage(); i++) {
        ins << vals[i];
      }
    }

    template <size_t N>
    inline void operator<<(const elem_type (&vals)[N]) {

      AWInsertManager<elem_type, wT> ins(get_wt_ref(), vals[0]);
      for (int i = 1; i < Storage(); i++) {
        ins << vals[i];
      }
    }

    /*
     * Inserting via a matrix/expression.
     */
    template<typename T>
    inline void operator<<(const arma::Base<elem_type, T> &X) {

      get_wt_ref().operator=(X);
    }


    /*
     * Sub matrix views.
     */

    inline AWSubView<elem_type, wT> SubMatrix(const int first_row,
                                              const int last_row,
                                              const int first_col,
                                              const int last_col) const {

      int fr = first_row;
      int lr = last_row;
      int fc = first_col;
      int lc = last_col;

      FixSubviewIndices(fr, lr, fc, lc);

      return AWSubView<elem_type, wT>(get_wt_ref(), fr, lr, fc, lc);
    }


    inline AWSubView<elem_type, wT> SubMatrix(const int first_row,
                                              const int last_row,
                                              const int first_col,
                                              const int last_col) {

      int fr = first_row;
      int lr = last_row;
      int fc = first_col;
      int lc = last_col;

      FixSubviewIndices(fr, lr, fc, lc);

      return AWSubView<elem_type, wT>(get_wt_ref(), fr, lr, fc, lc);
    }

    inline AWSubView<elem_type, wT> SymSubMatrix(const int first, const int last) const {
      return SubMatrix(first, last, first, last);
    }

    inline AWSubView<elem_type, wT> SymSubMatrix(const int first, const int last) {
      return SubMatrix(first, last, first, last);
    }

    inline AWSubView<elem_type, wT> Rows(const int first, const int last) const {
      return SubMatrix(first, last, 1, Ncols());
    }

    inline AWSubView<elem_type, wT> Rows(const int first, const int last) {
      return SubMatrix(first, last, 1, Ncols());
    }

    inline AWSubView<elem_type, wT> Columns(const int first, const int last) const {
      return SubMatrix(1, Nrows(), first, last);
    }

    inline AWSubView<elem_type, wT> Columns(const int first, const int last) {
      return SubMatrix(1, Nrows(), first, last);
    }

    inline AWSubView<elem_type, wT> Row(const int row) const {
      return SubMatrix(row, row, 1, Ncols());
    }

    inline AWSubView<elem_type, wT> Row(const int row) {
      return SubMatrix(row, row, 1, Ncols());
    }

    inline AWSubView<elem_type, wT> Column(const int col) const {
      return SubMatrix(1, Nrows(), col, col);
    }

    inline AWSubView<elem_type, wT> Column(const int col) {
      return SubMatrix(1, Nrows(), col, col);
    }

    /*
     * Matrix resizing. Resizing of 1D vectors and diagonal matrices
     * is implemented in the respective subclasses (see vector.hpp and
     * matrix_diag.hpp).
     */

    inline void ReSize(const int nrows, const int ncols) {
      const_cast<aT&>(get_at_ref()).set_size(nrows, ncols);
    }

    template<typename elem_type2, typename wT2, typename aT2>
    inline void ReSize(const AWBase<elem_type2, wT2, aT2> &X) {
      const_cast<aT&>(get_at_ref()).set_size(X.Nrows(), X.Ncols());
    }

    inline void CleanUp() {
      const_cast<aT&>(get_at_ref()).set_size(0, 0);
    }

    /*
     * Conversion to other types
     */

    inline elem_type AsScalar() const {
      return arma::as_scalar(get_at_ref());
    }

    inline typename armawrap_asrow_type_map< wT>::type AsRow()      const { return armawrap::AsRow<     wT>(get_wt_ref()); }
    inline typename armawrap_ascol_type_map< wT>::type AsColumn()   const { return armawrap::AsColumn<  wT>(get_wt_ref()); }
    inline typename armawrap_asdiag_type_map<wT>::type AsDiagonal() const { return armawrap::AsDiagonal<wT>(get_wt_ref()); }

    inline typename armawrap_asmat_type_map< wT>::type AsMatrix(const int nrows, const int ncols) const {
      return armawrap::AsMatrix<wT>(get_wt_ref(), nrows, ncols);
    }

    /*
     * Not sure about this.
     */
    inline wT &  Evaluate()         { return get_wt_ref(); }
    inline void  Release()          {}
    inline void  ReleaseAndDelete() {}
    inline wT &  ForReturn()        { return get_wt_ref(); }
  };
}

#endif
