#ifndef __MATRIX_DIAG_HPP__
#define __MATRIX_DIAG_HPP__

/*
 * Identity and diagonal matrix types.
 */

#include <cmath>
#include <exception>

namespace armawrap {

  template<typename elem_type>
    class AWIdentityMatrix : public arma::Mat<elem_type>, 
                             public armawrap::AWBase<elem_type, 
                                                     AWIdentityMatrix<elem_type>,
                                                     arma::Mat<elem_type> >{ 
  public:

    using armawrap::AWBase<elem_type, 
                           AWIdentityMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator<<;

    using armawrap::AWBase<elem_type, 
                           AWIdentityMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();

    using armawrap::AWBase<elem_type, 
                           AWIdentityMatrix<elem_type>, 
                           arma::Mat<elem_type> >::t;

    using armawrap::AWBase<elem_type, 
                           AWIdentityMatrix<elem_type>, 
                           arma::Mat<elem_type> >::ReSize; 

    /* AWIdentityMatrix i; */
    inline AWIdentityMatrix() : arma::Mat<elem_type>() {}

    inline virtual ~AWIdentityMatrix() { }

    /*
     * AWIdentityMatrix i(100);
     */
    inline AWIdentityMatrix(const int nelems) : arma::Mat<elem_type>(nelems, nelems) { 
      operator=(1);
    }

    /* AWIdentityMatrix a(some expression); */
    template<typename T>
    inline AWIdentityMatrix(const arma::Base<elem_type, T> &X) : arma::Mat<elem_type>() {
      operator=(X);
    }

    /*
     * AWIdentityMatrix i;
     * i = (some expression);
     */
    template<typename T> 
    inline const AWIdentityMatrix & 
    operator=(const arma::Base<elem_type, T> &X) {
      arma::Mat<elem_type>::operator=(arma::diagmat(X.get_ref()));
      return *this;
    }

    /*
     * AWIdentityMatrix i(100);
     * i = 99; // set all diagonal elements to 99.
     */
    inline const AWIdentityMatrix &
    operator=(const elem_type k) {

      arma::Mat<elem_type>::zeros();
      arma::Col<elem_type> diag = arma::Col<elem_type>(arma::Mat<elem_type>::n_rows);
      diag.fill(k);
      arma::Mat<elem_type>::diag() = diag;
      return *this;
    }

    /* 
     * Newmat only stores a single value for identity 
     * matrices, but armadillo stores a full matrix.
     */
    inline virtual int Storage() const { return 1; }

    inline void ReSize(int nelem) {
      ReSize(nelem, nelem);
    }

    /* 
     * Many of the common methods defined in 
     * AWBase are overridden here, as they have 
     * trivial results for identity matrices.
     */
    inline virtual elem_type Minimum()              const { return operator()(1); }
    inline virtual elem_type Maximum()              const { return operator()(1); }
    inline virtual elem_type MinimumAbsoluteValue() const { return std::abs(operator()(1)); }
    inline virtual elem_type MaximumAbsoluteValue() const { return std::abs(operator()(1)); }
    inline virtual elem_type SumSquare()            const { return this->Nrows() * std::pow(Minimum(), 2); }
    inline virtual elem_type SumAbsoluteValue()     const { return this->Nrows() * MinimumAbsoluteValue(); }
    inline virtual elem_type Sum()                  const { return this->Nrows() * Minimum(); }
    inline virtual elem_type Trace()                const { return Sum(); }
    inline virtual elem_type Norm1()                const { return MinimumAbsoluteValue(); }
    inline virtual elem_type NormInfinity()         const { return MinimumAbsoluteValue(); }
    inline virtual elem_type NormFrobenius()        const { return std::pow(SumSquare(), 0.5); }
    inline virtual bool      IsZero()               const { return Minimum() == 0; } 
  };


  template<typename elem_type> 
    class AWDiagonalMatrix : public arma::Mat<elem_type>, 
                             public armawrap::AWBase<elem_type, 
                                                     AWDiagonalMatrix<elem_type>,
                                                     arma::Mat<elem_type> >{ 

  public:

    using armawrap::AWBase<elem_type, 
                           AWDiagonalMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator<<;

    using armawrap::AWBase<elem_type, 
                           AWDiagonalMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();

    using armawrap::AWBase<elem_type, 
                           AWDiagonalMatrix<elem_type>, 
                           arma::Mat<elem_type> >::t;

    using armawrap::AWBase<elem_type, 
                           AWDiagonalMatrix<elem_type>, 
                           arma::Mat<elem_type> >::ReSize; 

    /* AWDiagonalMatrix i; */
    inline AWDiagonalMatrix() : arma::Mat<elem_type>() {}

    inline virtual ~AWDiagonalMatrix() { }

    /*
     * AWDiagonalMatrix i(100);
     */
    inline AWDiagonalMatrix(const int nelems) : arma::Mat<elem_type>(nelems, nelems, arma::fill::zeros) { }

    /*
     * AWDiagonalMatrix i(some expression);
     */ 
    template<typename T>
    inline AWDiagonalMatrix(const arma::Base<elem_type, T> &X) : arma::Mat<elem_type>() {
      operator=(X);
    }

    /*
     * AWDiagonalMatrix i;
     * i = (some expression);
     */
    template<typename T> 
    inline const AWDiagonalMatrix & 
    operator=(const arma::Base<elem_type, T> &X) {
      arma::Mat<elem_type>::operator=(arma::diagmat(X.get_ref()));
      return *this;
    }

    /*
     * AWDiagonalMatrix i(100);
     * i = 99; // set all diagonal elements to 99.
     */
    inline const AWDiagonalMatrix &
    operator=(const elem_type k) {
      arma::Mat<elem_type>::zeros();
      arma::Col<elem_type> diag = arma::Col<elem_type>(arma::Mat<elem_type>::n_rows);
      diag.fill(k);

      arma::Mat<elem_type>::diag() = diag;
      return *this;
    }

    /* 
     * Newmat only stores diagonal values for diagonal
     * matrices, but armadillo stores a full matrix.
     */
    inline virtual int Storage() const { return arma::Mat<elem_type>::n_rows; }

    inline void ReSize(int nelem) {
      ReSize(nelem, nelem);
    }

    /* 
     * Many of the common methods defined in AWBase are 
     * overridden here for performance reasons. The 
     * AWDiagonalMatrix is stored as a full square matrix, 
     * but we can improve performance by implementing 
     * these functions to work only on the diagonal.
     */
    inline virtual elem_type               Minimum()              const { return armawrap::Minimum(             this->diag()); }
    inline virtual elem_type               Maximum()              const { return armawrap::Maximum(             this->diag()); }
    inline virtual elem_type               MinimumAbsoluteValue() const { return armawrap::MinimumAbsoluteValue(this->diag()); }
    inline virtual elem_type               MaximumAbsoluteValue() const { return armawrap::MaximumAbsoluteValue(this->diag()); }
    inline virtual elem_type               SumSquare()            const { return armawrap::SumSquare(           this->diag()); }
    inline virtual elem_type               SumAbsoluteValue()     const { return armawrap::SumAbsoluteValue(    this->diag()); }
    inline virtual elem_type               Sum()                  const { return armawrap::Sum(                 this->diag()); }
    inline virtual elem_type               Trace()                const { return Sum(); }
    inline virtual elem_type               Norm1()                const { return MaximumAbsoluteValue(); }
    inline virtual elem_type               NormInfinity()         const { return MaximumAbsoluteValue(); }
    inline virtual elem_type               NormFrobenius()        const { return std::pow(SumSquare(), 0.5); }
    inline virtual bool                    IsZero()               const { return armawrap::IsZero(this->diag()); } 
  };
}

#endif /* __MATRIX_DIAG_HPP__ */
