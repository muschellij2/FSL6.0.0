#ifndef __OPERATOR_STREAM_HPP__
#define __OPERATOR_STREAM_HPP__

#include <iostream>


namespace armawrap {

  /*
   * This function has been written to produce identical output
   * to the NEWMAT::operator<< function (see newmat9.cpp).
   */
  template<typename wT>
  inline std::ostream & doPrint(std::ostream &s, const wT &X) {

    int  w = s.width();
    int  p = s.precision();

    s.setf(std::ios::fixed, std::ios::floatfield);

    for (int ri = 1; ri <= X.Nrows(); ri++) {
      for (int ci = 1; ci <= X.Ncols(); ci++) {
        s.width(w);
        s.precision(p);
        s << X(ri, ci) << " ";
      }
      s << "\n";
    }

    return s;
  }

  /*
   * Only the diagonal values are printed for identity/diagonal matrices.
   */
  template<typename wT>
  inline std::ostream & doPrintDiag(std::ostream &s, const wT &X) {
    int  w = s.width();
    int  p = s.precision();

    s.setf(std::ios::fixed, std::ios::floatfield);

    for (int ri = 1; ri <= X.Nrows(); ri++) {

      for (int space = 1; space < ri; space++) { s << " "; }

      s.width(w);
      s.precision(p);
      s << X(ri) << " ";
      s << "\n";
    }

    return s;
  }

  /*
   * values off the band for band matrices are not printed.
   */
  template<typename wT>
  inline std::ostream & doPrintBand(std::ostream &s, const wT &x) {

    int  w = s.width();
    int  p = s.precision();

    AWMatrixBandWidth bw = x.BandWidth();

    s.setf(std::ios::fixed, std::ios::floatfield);

    for (int ri = 1; ri <= x.Nrows(); ri++) {
      for (int ci = 1; ci <= x.Ncols(); ci++) {

        int diag = ci - ri;

        if (((diag < 0) && (-diag > bw.lower)) || 
            ((diag > 0) && ( diag > bw.upper))) {
          s << " ";
        }
        else {
          s.width(w);
          s.precision(p);
          s << x(ri, ci) << " ";
        }
      }
      s << "\n";
    }

    return s;
  }

  /*
   * Only the relevant half is printed for upper/lower triangular matrices.
   */
  template<typename wT>
  inline std::ostream & doPrintTri(std::ostream &s, const wT&X, bool upper) {

    int  w = s.width();
    int  p = s.precision();

    s.setf(std::ios::fixed, std::ios::floatfield);

    for (int ri = 1; ri <= X.Nrows(); ri++) {

      if (upper) { for (int space = 1; space < ri; space++) { s << " "; }; }

      int ci;
      int ce;

      if (upper) ci = ri;
      else       ci = 1;

      if (upper) ce = X.Nrows();
      else       ce = ri;

      for (; ci <= ce; ci++) {
        s.width(w);
        s.precision(p);
        s << X(ri, ci) << " ";
      }
      s << "\n";
    }

    return s;
  }


  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWMatrix<elem_type> &X) {
    return doPrint<AWMatrix<elem_type> >(s, X);
  }

  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWRowVector<elem_type> &X) {
    return doPrint<AWRowVector<elem_type> >(s, X);
  }

  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWColVector<elem_type> &X) {
    return doPrint<AWColVector<elem_type> >(s, X);
  }

  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWIdentityMatrix<elem_type> &X) {
    return doPrintDiag<AWIdentityMatrix<elem_type> >(s, X);
  }

  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWDiagonalMatrix<elem_type> &X) {
    return doPrintDiag<AWDiagonalMatrix<elem_type> >(s, X);
  }


  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWUpperTriangularMatrix<elem_type> &X) {
    return doPrintTri<AWUpperTriangularMatrix<elem_type> >(s, X, true);
  }

  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWLowerTriangularMatrix<elem_type> &X) {
    return doPrintTri<AWLowerTriangularMatrix<elem_type> >(s, X, false);
  }

  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWSymmetricMatrix<elem_type> &X) {
    return doPrint<AWSymmetricMatrix<elem_type> >(s, X);
  }

  template<typename elem_type>
  inline std::ostream &
  operator<<(std::ostream &s, const AWBandMatrix<elem_type> &X) {
    return doPrintBand<AWBandMatrix<elem_type> >(s, X);
  }

  template<typename elem_type, typename wT>
  inline std::ostream &
  operator<<(std::ostream &s, const AWSubView<elem_type, wT> &X) {
    return doPrint<AWMatrix<elem_type> >(s, X);
  }

  /*
   * Insertion from an input stream
   */

  template<typename elem_type, typename wT>
  inline std::istream &
  operator>>(std::istream &s, const AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value> &X) {
  
    elem_type v;
    s >> v;

    // the input parameter must be const, to allow references
    // to temporary objects, which the AWCallmanager always is
    const_cast<AWCallManager<elem_type, wT, is_simple_armawrap_type<wT>::value> &>(X) = v;
    return s;
  }
}

#endif /* __OPERATOR_STREAM_HPP__ */
