#ifndef __OPERATOR_CALL_HPP__
#define __OPERATOR_CALL_HPP__

#include <cmath>

/*
 * Call operator ( () ) for matrix element access and assignment.
 */


namespace armawrap {

  /*
   * Standalone functions for each matrix type are implemented to:
   *
   *   - Look up the value of an element (read only) - lookup*D
   *   - Assign a value to an element                - assign*D
   *   - Obtain a reference to an element            - ref*D
   *
   * The ref*D functions are only implemented for the 'simple' matrix types
   * (AWMatrix, AWColVector, AWRowVector).
   */

  /*
   * AWMatrix
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWMatrix<elem_type> &target, const int i) {
    arma::uword row = (i - 1) / target.Ncols();
    arma::uword col = (i - 1) % target.Ncols();
    return target.get_at_ref()(row, col);
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWMatrix<elem_type> &target, const int row, const int col) {
    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWMatrix<elem_type> &target, const int i, const elem_type val) {
    arma::uword row = (i - 1) / target.Ncols();
    arma::uword col = (i - 1) % target.Ncols();
    target.get_at_ref()(row, col) = val;
  }

  template<typename elem_type>
  inline void assign2D(AWMatrix<elem_type> &target, const int row, const int col, const elem_type val) {
    target.get_at_ref()(row - 1, col - 1) = val;
  }

  template<typename elem_type>
  inline elem_type & ref1D(AWMatrix<elem_type> &target, const int i) {
    arma::uword row = (i - 1) / target.Ncols();
    arma::uword col = (i - 1) % target.Ncols();
    return target.get_at_ref()(row, col);
  }

  template<typename elem_type>
  inline elem_type & ref2D(AWMatrix<elem_type> &target, const int row, const int col) {
    return target.get_at_ref()(row - 1, col - 1);
  }

  /*
   * AWRowVector
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWRowVector<elem_type> &target, const int i) {
    return target.get_at_ref()(i - 1);
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWRowVector<elem_type> &target, const int row, const int col) {
    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWRowVector<elem_type> &target, const int i, const elem_type val) {
    target.get_at_ref()(i - 1) = val;
  }

  template<typename elem_type>
  inline void assign2D(AWRowVector<elem_type> &target, const int row, const int col, const elem_type val) {
    target.get_at_ref()(row - 1, col - 1) = val;
  }

  template<typename elem_type>
  inline elem_type & ref1D(AWRowVector<elem_type> &target, const int i) {
    return target.get_at_ref()(i - 1);
  }

  template<typename elem_type>
  inline elem_type & ref2D(AWRowVector<elem_type> &target, const int row, const int col) {
    return target.get_at_ref()(row - 1, col - 1);
  }

  /*
   * AWColVector
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWColVector<elem_type> &target, const int i) {
    return target.get_at_ref()(i - 1);
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWColVector<elem_type> &target, const int row, const int col) {
    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWColVector<elem_type> &target, const int i, const elem_type val) {
    target.get_at_ref()(i - 1) = val;
  }

  template<typename elem_type>
  inline void assign2D(AWColVector<elem_type> &target, const int row, const int col, const elem_type val) {
    target.get_at_ref()(row - 1, col - 1) = val;
  }

  template<typename elem_type>
  inline elem_type & ref1D(AWColVector<elem_type> &target, const int i) {
    return target.get_at_ref()(i - 1);
  }

  template<typename elem_type>
  inline elem_type & ref2D(AWColVector<elem_type> &target, const int row, const int col) {
    return target.get_at_ref()(row - 1, col - 1);
  }

  /*
   * AWIdentityMatrix
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWIdentityMatrix<elem_type> &target, const int i) {
    return target.get_at_ref()(i - 1, i - 1);
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWIdentityMatrix<elem_type> &target, const int row, const int col) {
    static elem_type default_val = 0;
    if (row != col) return default_val;
    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWIdentityMatrix<elem_type> &target, const int i, const elem_type val) {
    throw AWException("Attempt to assign to element of IdentityMatrix");
  }

  template<typename elem_type>
  inline void assign2D(AWIdentityMatrix<elem_type> &target, const int row, const int col, const elem_type val) {
    throw AWException("Attempt to assign to element of IdentityMatrix");
  }

  /*
   * DiagonalMatrix
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWDiagonalMatrix<elem_type> &target, const int i) {
    return target.get_at_ref()(i - 1, i - 1);
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWDiagonalMatrix<elem_type> &target, const int row, const int col) {
    static elem_type default_val = 0;
    if (row != col) return default_val;
    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWDiagonalMatrix<elem_type> &target, const int i, const elem_type val) {
    target.get_at_ref()(i - 1, i - 1) = val;
  }

  template<typename elem_type>
  inline void assign2D(AWDiagonalMatrix<elem_type> &target, const int row, const int col, const elem_type val) {
    if (row != col) return;
    target.get_at_ref()(row - 1, col - 1) = val;
  }

  template<typename elem_type>
  inline elem_type & ref1D(AWDiagonalMatrix<elem_type> &target, const int i) {
    return target.get_at_ref()(i - 1, i - 1);
  }

  template<typename elem_type>
  inline elem_type & ref2D(AWDiagonalMatrix<elem_type> &target, const int row, const int col) {
    if (row != col) throw AWException("Attempt to obtain a reference to a non-diagonal element in an AWDiagonalMatrix");
    return target.get_at_ref()(row - 1, col - 1);
  }

  /*
   * The next few functions are used for 1D indexing of triangular matrices
   * (converting from a 1D sequential index to 2D row/column indices).
   */

  /*
   * Calculate the nth triangular number.
   */
  inline double _tri(double n) {
    return (n * (n + 1)) / 2.0;
  }

  /*
   * For an upper triangular matrix of size base*base, calculate the
   * index (starting from 0) of the diagonal element on the nth row.
   */
  inline double _ut_diag(double base, double n) {
    return _tri(base) - _tri(base - n);
  }

  /*
   * For an upper triangular matrix of size base*base, calculate the diagonal
   * number (the row number, starting from 0) of the element at the given 1D
   * index. Take the floor of the result.
   */
  inline double _ut_diag_root(double base, double n) {
    return ((2 * base + 1) - sqrt(pow(-2 * base - 1, 2) - 8*n)) / 2.0;
  }

  /*
   * For a lower triangular matrix of any size, calculate the index
   * (starting from 0) of the diagonal element on the nth row.
   */
  inline double _lt_diag(double n) {
    return _tri(n + 1) - 1;
  }

  /*
   * For a lower triangular matrix of any size, calculate the row number
   * of the element at the given 1D index. Take the ceiling of the result.
   */
  inline double _lt_diag_root(double n) {
    return -1.5 + sqrt(pow(1.5, 2) + 2 * n);
  }

  /*
   * UpperTriangularMatrix
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWUpperTriangularMatrix<elem_type> &target, const int i) {
    arma::uword row = floor(_ut_diag_root(target.Nrows(), i - 1));
    arma::uword col = row + i - 1 - _ut_diag(target.Nrows(), row);
    return target.get_at_ref()(row, col);
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWUpperTriangularMatrix<elem_type> &target, const int row, const int col) {
    static elem_type default_val = 0;
    if (row > col) return default_val;
    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWUpperTriangularMatrix<elem_type> &target, const int i, const elem_type val) {
    arma::uword row = floor(_ut_diag_root(target.Nrows(), i - 1));
    arma::uword col = row + i - 1 - _ut_diag(target.Nrows(), row);
    target.get_at_ref()(row, col) = val;
  }

  template<typename elem_type>
  inline void assign2D(AWUpperTriangularMatrix<elem_type> &target, const int row, const int col, const elem_type val) {
    if (row > col) throw AWException("Access to lower half of UpperTriangularMatrix");
    target.get_at_ref()(row - 1, col - 1) = val;
  }

  /*
   * LowerTriangularMatrix
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWLowerTriangularMatrix<elem_type> &target, const int i) {
    arma::uword row = ceil(_lt_diag_root(i - 1));
    arma::uword col = row - (_lt_diag(row) - (i - 1));
    return target.get_at_ref()(row, col);
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWLowerTriangularMatrix<elem_type> &target, const int row, const int col) {
    static elem_type default_val = 0;
    if (row < col) return default_val;
    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWLowerTriangularMatrix<elem_type> &target, const int i, const elem_type val) {
    arma::uword row = ceil(_lt_diag_root(i - 1));
    arma::uword col = row - (_lt_diag(row) - (i - 1));
    target.get_at_ref()(row, col) = val;
  }

  template<typename elem_type>
  inline void assign2D(AWLowerTriangularMatrix<elem_type> &target, const int row, const int col, const elem_type val) {
    if (row < col) throw AWException("Access to upper half of LowerTriangularMatrix");
    target.get_at_ref()(row - 1, col - 1) = val;
  }

  /*
   * SymmetricMatrix. Newmat stores symmetric matrix
   * in the same way as it stores lower triangular
   * matrices, so we use the same 1D/2D conversion.
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWSymmetricMatrix<elem_type> &target, const int i) {
    arma::uword row = ceil(_lt_diag_root(i - 1));
    arma::uword col = row - (_lt_diag(row) - (i - 1));
    return target.get_at_ref()(row, col);
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWSymmetricMatrix<elem_type> &target, const int row, const int col) {
    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWSymmetricMatrix<elem_type> &target, const int i, const elem_type val) {
    arma::uword row = ceil(_lt_diag_root(i - 1));
    arma::uword col = row - (_lt_diag(row) - (i - 1));
    target.get_at_ref()(row, col) = val;
    target.get_at_ref()(col, row) = val;
  }

  template<typename elem_type>
  inline void assign2D(AWSymmetricMatrix<elem_type> &target, const int row, const int col, const elem_type val) {
    target.get_at_ref()(row - 1, col - 1) = val;
    target.get_at_ref()(col - 1, row - 1) = val;
  }

  /*
   * BandMatrix
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWBandMatrix<elem_type> &target, const int i) {
    throw AWException("1D indexing is not yet implemented on BandMatrix types.");
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWBandMatrix<elem_type> &target, const int row, const int col) {

    static elem_type default_val = 0;

    AWMatrixBandWidth bw = target.BandWidth();

    int diag = col - row;

    if (diag < 0 && (-diag > bw.lower)) return default_val;
    if (diag > 0 && ( diag > bw.upper)) return default_val;

    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWBandMatrix<elem_type> &target, const int i, const elem_type val) {
    throw AWException("1D indexing is not yet implemented on BandMatrix types.");
  }

  template<typename elem_type>
  inline void assign2D(AWBandMatrix<elem_type> &target, const int row, const int col, const elem_type val) {

    AWMatrixBandWidth bw = target.BandWidth();

    int diag = col - row;

    if (diag < 0 && (-diag > bw.lower)) throw AWException("Access to off-band region of BandMatrix");
    if (diag > 0 && ( diag > bw.upper)) throw AWException("Access to off-band region of BandMatrix");

    target.get_at_ref()(row - 1, col - 1) = val;
  }

  /*
   * SymmetricBandMatrix
   */

  template<typename elem_type>
  inline const elem_type & lookup1D(AWSymmetricBandMatrix<elem_type> &target, const int i) {
    throw AWException("1D indexing is not yet implemented on BandMatrix types.");
  }

  template<typename elem_type>
  inline const elem_type & lookup2D(AWSymmetricBandMatrix<elem_type> &target, const int row, const int col) {

    static elem_type default_val = 0;

    AWMatrixBandWidth bw = target.BandWidth();

    int diag = col - row;

    if (diag < 0 && (-diag > bw.lower)) return default_val;
    if (diag > 0 && ( diag > bw.upper)) return default_val;

    return target.get_at_ref()(row - 1, col - 1);
  }

  template<typename elem_type>
  inline void assign1D(AWSymmetricBandMatrix<elem_type> &target, const int i, const elem_type val) {
    arma::uword row = (i - 1) / target.Ncols();
    arma::uword col = (i - 1) % target.Ncols();

    AWMatrixBandWidth bw = target.BandWidth();

    int diag = col - row;

    if (diag < 0 && (-diag > bw.lower)) throw AWException("Access to off-band region of SymmetricBandMatrix");
    if (diag > 0 && ( diag > bw.upper)) throw AWException("Access to off-band region of SymmetricBandMatrix");

    target.get_at_ref()(row, col) = val;
    target.get_at_ref()(col, row) = val;
  }

  template<typename elem_type>
  inline void assign2D(AWSymmetricBandMatrix<elem_type> &target, const int row, const int col, const elem_type val) {

    AWMatrixBandWidth bw = target.BandWidth();

    int diag = col - row;

    if (diag < 0 && (-diag > bw.lower)) throw AWException("Access to off-band region of SymmetricBandMatrix");
    if (diag > 0 && ( diag > bw.upper)) throw AWException("Access to off-band region of SymmetricBandMatrix");

    target.get_at_ref()(row - 1, col - 1) = val;
    target.get_at_ref()(col - 1, row - 1) = val;
  }


  /**
   * Element access on the armawrap matrix types triggers creation of an
   * AWCallManager object. This is because different matrix types require
   * different lookup strategies and, more importantly, assignment to elements
   * of some matrix types needs to result in side effects (the SymmetricMatrix
   * is a prime example here).
   *
   * Lookup and assignment is performed by the standalone lookup1D, lookup2D,
   * assign1D and assign2D functions defined above.
   *
   * The AWCallManagerBase class contains most of the AWCallManager
   * functionality. Two partial template specialisations, called
   * AWCallManager, follow the AWCallManagerBase definition.
   */
  template<typename elem_type, typename wT>
  class AWCallManagerBase {
  protected:

    wT &      target;
    const int i;
    const int row;
    const int col;

  public:

    inline AWCallManagerBase(const wT &target, const int i) :
      target(const_cast<wT &>(target)),
      i(i),
      row(-1),
      col(-1) { };

    inline AWCallManagerBase(const wT &target, const int row, const int col) :
      target(const_cast<wT &>(target)),
      i(-1),
      row(row),
      col(col) { };

    inline operator const elem_type & () const {
      if (i != -1) return lookup1D(target, i);
      else         return lookup2D(target, row, col);
    }
  };


  /**
   * If the matrix is of a 'simple' type, it is possible to obtain references
   * to individual elements (elem_type &). This is not possible for the other
   * matrix types.
   */
  template<typename elem_type, typename wT, bool isSimpleType> class AWCallManager;

  /*
   * Default AWCallManager, for non-simple matrix types.
   */
  template<typename elem_type, typename wT>
  class AWCallManager<elem_type, wT, false> :
    public AWCallManagerBase<elem_type, wT> {

  public:

    inline AWCallManager(const wT &target, const int i) :
      AWCallManagerBase<elem_type, wT>(target, i) { }

    inline AWCallManager(const wT &target, const int row, const int col) :
      AWCallManagerBase<elem_type, wT>(target, row, col) { }


    inline elem_type operator=(AWCallManager &other) {
      return operator=(static_cast<elem_type>(other));
    }

    inline elem_type operator=(elem_type val) {
      if (this->i != -1) assign1D(this->target, this->i, val);
      else               assign2D(this->target, this->row, this->col, val);
      return val;
    }

    inline void operator+=(elem_type val) {

      elem_type oldval;

      if (this->i != -1) oldval = lookup1D(this->target, this->i);
      else               oldval = lookup2D(this->target, this->row, this->col);

      if (this->i != -1) assign1D(this->target, this->i,              oldval + val);
      else               assign2D(this->target, this->row, this->col, oldval + val);
    }

    inline void operator-=(elem_type val) {

      elem_type oldval;

      if (this->i != -1) oldval = lookup1D(this->target, this->i);
      else               oldval = lookup2D(this->target, this->row, this->col);

      if (this->i != -1) assign1D(this->target, this->i,              oldval - val);
      else               assign2D(this->target, this->row, this->col, oldval - val);
    }


    inline elem_type operator++(int) {
      elem_type val = *this;
      operator+=(1);
      return val;
    }

    inline elem_type operator--(int) {
      elem_type val = *this;
      operator-=(1);
      return val;
    }


    inline void operator*=(elem_type val) {

      elem_type oldval;

      if (this->i != -1) oldval = lookup1D(this->target, this->i);
      else               oldval = lookup2D(this->target, this->row, this->col);

      if (this->i != -1) assign1D(this->target, this->i,              oldval * val);
      else               assign2D(this->target, this->row, this->col, oldval * val);
    }


    inline void operator/=(elem_type val) {

      elem_type oldval;

      if (this->i != -1) oldval = lookup1D(this->target, this->i);
      else               oldval = lookup2D(this->target, this->row, this->col);

      if (this->i != -1) assign1D(this->target, this->i,              oldval / val);
      else               assign2D(this->target, this->row, this->col, oldval / val);
    }

    template<typename eT2> inline bool operator==(eT2 val) { return static_cast<elem_type>(*this) == val; }
    template<typename eT2> inline bool operator!=(eT2 val) { return static_cast<elem_type>(*this) != val; }
    template<typename eT2> inline bool operator<( eT2 val) { return static_cast<elem_type>(*this) <  val; }
    template<typename eT2> inline bool operator>( eT2 val) { return static_cast<elem_type>(*this) >  val; }
    template<typename eT2> inline bool operator<=(eT2 val) { return static_cast<elem_type>(*this) <= val; }
    template<typename eT2> inline bool operator>=(eT2 val) { return static_cast<elem_type>(*this) >= val; }
  };

  /*
   * Partial template specialisation for 'simple' matrix types. The only
   * difference between this specialisation and the above is that this one
   * defines a conversion operator to (elem_type &), thus can be used
   * to obtain a reference to matrix elements.
   */
  template<typename elem_type, typename wT>
  class AWCallManager<elem_type, wT, true> :
    public AWCallManagerBase<elem_type, wT> {
  public:

    inline AWCallManager(const wT &target, const int i) :
      AWCallManagerBase<elem_type, wT>(target, i) { }

    inline AWCallManager(const wT &target, const int row, const int col) :
      AWCallManagerBase<elem_type, wT>(target, row, col) { }

    inline operator elem_type & () {
      if (this->i != -1) return ref1D(this->target, this->i);
      else               return ref2D(this->target, this->row, this->col);
    }


    elem_type * operator&() {
      return &(elem_type &)(*this);
    }

    inline elem_type operator=(const elem_type val) {
      (elem_type &)(*this) = (elem_type)val;
      return val;
    }

    inline elem_type operator=(const AWCallManager &other) {
      elem_type val = (elem_type)other;
      (elem_type &)(*this) = val;
      return val;
    }
  };
}


/*
 * These extensions to various standard namespace functions are necessary
 * so that expressions such as the following will compile:
 *
 *   Matrix a;
 *   ...
 *   cout << std::max(a(4,6), 5.0) << endl;
 *
 * And this is also why every armawrap file contains the 'namespace armawrap
 * {...}'  declaration!
 */
namespace std {

  template<typename elem_type, typename wT, typename T>
  inline T max(const T X, const armawrap::AWCallManager<elem_type, wT, armawrap::is_simple_armawrap_type<wT>::value> &Y) {
    return std::max<T>(X, static_cast<elem_type>(Y));
  }

  template<typename elem_type, typename wT, typename T>
  inline T max(const armawrap::AWCallManager<elem_type, wT, armawrap::is_simple_armawrap_type<wT>::value> &X, const T Y) {
    return std::max<T>(static_cast<elem_type>(X), Y);
  }

  template<typename elem_type, typename wT1, typename wT2>
  inline elem_type max(const armawrap::AWCallManager<elem_type, wT1, armawrap::is_simple_armawrap_type<wT1>::value> &X,
                       const armawrap::AWCallManager<elem_type, wT2, armawrap::is_simple_armawrap_type<wT2>::value> &Y) {
    return std::max<elem_type>(static_cast<elem_type>(X), Y);
  }

  template<typename elem_type, typename wT, typename T>
  inline T min(const T X, const armawrap::AWCallManager<elem_type, wT, armawrap::is_simple_armawrap_type<wT>::value> &Y) {
    return std::min<T>(X, static_cast<elem_type>(Y));
  }

  template<typename elem_type, typename wT, typename T>
  inline T min(const armawrap::AWCallManager<elem_type, wT, armawrap::is_simple_armawrap_type<wT>::value> &X, const T Y) {
    return std::min<T>(static_cast<elem_type>(X), Y);
  }

  template<typename elem_type, typename wT1, typename wT2>
  inline elem_type min(const armawrap::AWCallManager<elem_type, wT1, armawrap::is_simple_armawrap_type<wT1>::value> &X,
                       const armawrap::AWCallManager<elem_type, wT2, armawrap::is_simple_armawrap_type<wT2>::value> &Y) {
    return std::min<elem_type>(static_cast<elem_type>(X), Y);
  }

  template<typename elem_type, typename wT1, typename wT2>
  inline void swap(const armawrap::AWCallManager<elem_type, wT1, armawrap::is_simple_armawrap_type<wT1>::value> &a,
                   const armawrap::AWCallManager<elem_type, wT2, armawrap::is_simple_armawrap_type<wT2>::value> &b) {

    // the function arguments must be const,
    // so temporaries can be passed in
    armawrap::AWCallManager<elem_type, wT1, armawrap::is_simple_armawrap_type<wT1>::value> ca =
      const_cast<armawrap::AWCallManager<elem_type, wT1, armawrap::is_simple_armawrap_type<wT1>::value> &>(a);
    armawrap::AWCallManager<elem_type, wT2, armawrap::is_simple_armawrap_type<wT2>::value> cb =
      const_cast<armawrap::AWCallManager<elem_type, wT2, armawrap::is_simple_armawrap_type<wT2>::value> &>(b);

    elem_type va = a;
    ca = (elem_type)b;
    cb = va;
  }

  template<typename elem_type, typename wT>
  inline bool isfinite(const armawrap::AWCallManager<elem_type, wT, armawrap::is_simple_armawrap_type<wT>::value> &val) {
    return isfinite((elem_type)val);
  }

  template<typename elem_type, typename wT>
  inline bool isinf(const armawrap::AWCallManager<elem_type, wT, armawrap::is_simple_armawrap_type<wT>::value> &val) {
    return isinf((elem_type)val);
  }

  template<typename elem_type, typename wT>
  inline bool isnan(const armawrap::AWCallManager<elem_type, wT, armawrap::is_simple_armawrap_type<wT>::value> &val) {
    return isnan((elem_type)val);
  }
}

#endif /* __OPERATOR_CALL_HPP__ */
