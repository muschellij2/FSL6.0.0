#ifndef __EXPRS_HPP__
#define __EXPRS_HPP__

/*
 * Subclasses of all the armadillo operation/expression types.
 * These are necessary so that Newmat methods may be called
 * on the results of expressions, e.g.
 *
 *   Matrix a, b;
 *   ...
 *   cout << (a + b).Sum() << endl;
 */

namespace armawrap {

  /*
   * The AWExpr class is the base class for all AW expression classes.
   * Its sole purpose for existence is so that submatrix views may be
   * called on the results of expressions, e.g.:
   *
   *   Matrix a, b;
   *   ...
   *   cout << (a + b).Row(4).Sum() << endl;
   */
  template<typename elem_type, typename aT, typename wT>
  class AWExpr : public armawrap::AWBase<elem_type, aT, wT> {
  public:
    
    inline AWSubView<elem_type, arma::Mat<elem_type> >
    SubMatrix(const int first_row,
              const int last_row,
              const int first_col,
              const int last_col) const {

      int fr = first_row;
      int lr = last_row;
      int fc = first_col;
      int lc = last_col;

      FixSubviewIndices(fr, lr, fc, lc);

      // Armadillo subview objects do not accept expression objects -
      // they are converted into matrices (i.e. evaluated) when
      // passed to the arma::subview constructor. There's a problem
      // here though - if a temporary expression object is passed
      // to the subview constructor, it is evaluated, but the
      // resulting matrix is also a temporary object which goes out
      // of scope when the arma::subview constructor completes.
      // This causes random segfaults with code like:
      //
      //   Matrix a;
      //   ...
      //   std::cout << a.Column(3) - 10.0).Row(1).AsScalar() << std::endl;
      //
      // because by the time the AsScalar method is called on the second
      // subview object, the temporary matrix created when that subview
      // was constructed might have already been destroyed.
      //
      // I'm circumventing this behaviour by putting evaluated
      // expression matrices on the heap, before passing them to the
      // subview constructor. The AWSubView object will delete the
      // matrix when it goes out of scope.
      arma::Mat<elem_type> *eval = new arma::Mat<elem_type>(this->get_at_ref());
      return AWSubView<elem_type, arma::Mat<elem_type> >(*eval,
                                                         fr, lr, fc, lc,
                                                         eval); 
    }


    inline AWSubView<elem_type, arma::Mat<elem_type> >
    SubMatrix(const int first_row, 
              const int last_row, 
              const int first_col, 
              const int last_col) {

      int fr = first_row;
      int lr = last_row;
      int fc = first_col;
      int lc = last_col;

      FixSubviewIndices(fr, lr, fc, lc);
      
      arma::Mat<elem_type> *eval = new arma::Mat<elem_type>(this->get_at_ref());
      return AWSubView<elem_type, arma::Mat<elem_type> >(*eval, 
                                                         fr, lr, fc, lc,
                                                         eval);
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    SymSubMatrix(const int first, const int last) const {
      return SubMatrix(first, last, first, last);
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    SymSubMatrix(const int first, const int last) {
      return SubMatrix(first, last, first, last);
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    Rows(const int first, const int last) const {
      return SubMatrix(first, last, 1, this->Ncols());
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    Rows(const int first, const int last) {
      return SubMatrix(first, last, 1, this->Ncols());
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    Columns(const int first, const int last) const {
      return SubMatrix(1, this->Nrows(), first, last);
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    Columns(const int first, const int last) {
      return SubMatrix(1, this->Nrows(), first, last);
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    Row(const int row) const {
      return SubMatrix(row, row, 1, this->Ncols());
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    Row(const int row) {
      return SubMatrix(row, row, 1, this->Ncols());
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    Column(const int col) const {
      return SubMatrix(1, this->Nrows(), col, col);
    }

    inline AWSubView<elem_type, arma::Mat<elem_type> >
    Column(const int col) {
      return SubMatrix(1, this->Nrows(), col, col);
    } 
  };
                                         

  template<typename T1, typename T2, typename eglue_type>
  class AWEGlue : public arma::eGlue<T1, T2, eglue_type>,
                  public AWExpr<typename T1::elem_type, 
                                AWEGlue<T1, T2, eglue_type>,
                                arma::eGlue<T1, T2, eglue_type> > {  
  public:

    using armawrap::AWBase<typename T1::elem_type, 
                           AWEGlue<T1, T2, eglue_type>, 
                           arma::eGlue<T1, T2, eglue_type> >::operator();
    
    using armawrap::AWBase<typename T1::elem_type, 
                           AWEGlue<T1, T2, eglue_type>, 
                           arma::eGlue<T1, T2, eglue_type> >::t;

    inline virtual ~AWEGlue() { }
    
    inline AWEGlue(const T1 &A, const T2 &B) : arma::eGlue<T1, T2, eglue_type>(A, B) { }
  };


  template<typename T1, typename T2, typename glue_type>
  class AWGlue : public arma::Glue<T1, T2, glue_type>, 
                 public AWExpr<typename T1::elem_type, 
                               AWGlue<T1, T2, glue_type>,
                               arma::Glue<T1, T2, glue_type> > {  
  public:

    using armawrap::AWBase<typename T1::elem_type, 
                           AWGlue<T1, T2, glue_type>, 
                           arma::Glue<T1, T2, glue_type> >::operator();
    
    using armawrap::AWBase<typename T1::elem_type, 
                           AWGlue<T1, T2, glue_type>, 
                           arma::Glue<T1, T2, glue_type> >::t;

    inline virtual ~AWGlue() { }
    
    inline AWGlue(const T1 &A, const T2 &B)                        : arma::Glue<T1, T2, glue_type>(A, B)      { }
    inline AWGlue(const T1 &A, const T2 &B, const arma::uword aux) : arma::Glue<T1, T2, glue_type>(A, B, aux) { }
  };


  template<typename T, typename op_type>
  class AWOp : public arma::Op<T, op_type>, 
               public AWExpr<typename T::elem_type, 
                             AWOp<T, op_type>,
                             arma::Op<T, op_type> > {

  private:

    /* 
     * If this field is set in the constructor, 
     * it is deleted in the destructor.
     */
    const T *ptr;

  public:

    using armawrap::AWBase<typename T::elem_type, 
                           AWOp<T, op_type>,
                           arma::Op<T, op_type> >::operator();

    using armawrap::AWBase<typename T::elem_type, 
                           AWOp<T, op_type>,
                           arma::Op<T, op_type> >::t;

    /*
     * If a variant of the constructor with the 'ignored' parameter is used,
     * it is assumed that the given T has been created on the heap. A pointer
     * to it is saved and, when the destructor for this AWOp object is called,
     * the pointer is deleted.
     *
     * This is basically a horrible hack, implemented for the purposes of
     * enabling nested expressions. If a nested Op object has been created as
     * a temporary, or goes out of scope before the expression is evaluated,
     * it will be deleted, which will result in all sorts of wonderful errors.
     * 
     * See the AWBase::Reverse method (in base.hpp) for an example of where 
     * this constructor is used.
     */

    inline AWOp(                    const T &A) : arma::Op<T, op_type>(A), ptr(NULL) { }
    inline AWOp(const bool ignored, const T &A) : arma::Op<T, op_type>(A), ptr(&A)   { }
    
    inline AWOp(                    const T &A, const typename T::elem_type in_aux) :
      arma::Op<T, op_type>(A, in_aux), ptr(NULL) { }
    inline AWOp(const bool ignored, const T &A, const typename T::elem_type in_aux) :
      arma::Op<T, op_type>(A, in_aux), ptr(&A)   { }
    
    inline AWOp(                    const T &A, const typename T::elem_type in_aux, const arma::uword in_aux_uword_a, const arma::uword in_aux_uword_b) :
      arma::Op<T, op_type>(A, in_aux, in_aux_uword_a, in_aux_uword_b), ptr(NULL) { }
    inline AWOp(const bool ignored, const T &A, const typename T::elem_type in_aux, const arma::uword in_aux_uword_a, const arma::uword in_aux_uword_b) :
      arma::Op<T, op_type>(A, in_aux, in_aux_uword_a, in_aux_uword_b), ptr(&A)   { } 
    
    inline AWOp(                    const T &A, const arma::uword in_aux_uword_a, const arma::uword in_aux_uword_b) :
      arma::Op<T, op_type>(A, in_aux_uword_a, in_aux_uword_b), ptr(NULL) { }
    inline AWOp(const bool ignored, const T &A, const arma::uword in_aux_uword_a, const arma::uword in_aux_uword_b) :
      arma::Op<T, op_type>(A, in_aux_uword_a, in_aux_uword_b), ptr(&A)   { } 
    
    inline AWOp(                    const T &A, const arma::uword in_aux_uword_a, const arma::uword in_aux_uword_b, const arma::uword in_aux_uword_c, const char junk) :
      arma::Op<T, op_type>(A, in_aux_uword_a, in_aux_uword_b, in_aux_uword_c, junk), ptr(NULL) { }
    inline AWOp(const bool ignored, const T &A, const arma::uword in_aux_uword_a, const arma::uword in_aux_uword_b, const arma::uword in_aux_uword_c, const char junk) :
      arma::Op<T, op_type>(A, in_aux_uword_a, in_aux_uword_b, in_aux_uword_c, junk), ptr(&A)   { } 
    

    inline virtual ~AWOp() {
      
      if (ptr != NULL) {
        delete ptr;
      }
    }
  };

  template<typename T, typename eop_type>
  class AWEOp : public arma::eOp<T, eop_type>, 
                public AWExpr<typename T::elem_type, 
                              AWEOp<T, eop_type>,
                              arma::eOp<T, eop_type> > {
    
  public:

    using armawrap::AWBase<typename T::elem_type, 
                           AWEOp<T, eop_type>,
                           arma::eOp<T, eop_type> >::operator();
    using armawrap::AWBase<typename T::elem_type, 
                           AWEOp<T, eop_type>,
                           arma::eOp<T, eop_type> >::t;

    inline virtual ~AWEOp() { }

    inline AWEOp(const T &A)                                : arma::eOp<T, eop_type>(A)    { }
    inline AWEOp(const T &A, const typename T::elem_type k) : arma::eOp<T, eop_type>(A, k) { }
  };
}

#endif /* __EXPRS_HPP__ */
