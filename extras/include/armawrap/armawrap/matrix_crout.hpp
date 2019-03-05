#ifndef __MATRIX_CROUT_HPP__
#define __MATRIX_CROUT_HPP__

/*
 * CroutMatrix - a proxy object which encapsulates the LUP decomposition of a
 * matrix.
 */

namespace armawrap {

  /*
   * An instance of this class is returned by the AWCroutMatrix.i() method.
   * It supports two operations:
   *
   *  -  calculation of the inverse of the matrix originally passed to the 
   *     AWCroutMatrix
   *
   *  - Solving of a linear equation, when multiplied with another matrix.
   *
   * Both of these operations use the LUP decomposition of the original 
   * matrix, as stored in the AWCroutMatrix.
   *
   * See http://en.wikipedia.org/wiki/LU_decomposition
   */

  template<typename elem_type>
  class AWInvertedCroutMatrix : public arma::Mat<elem_type>,
                                public armawrap::AWBase<elem_type,
                                                        AWInvertedCroutMatrix<elem_type>,
                                                        arma::Mat<elem_type> > {

  private:
    friend class AWCroutMatrix<elem_type>;

    const AWCroutMatrix<elem_type> &A;

    /*
     * Calculate the inverse of a matrix from its LUP
     * decomposition.  Given a matrix A, we solve the following 
     * equation for X:
     *
     *    AX = I
     *
     * Given the LUP decomposition of A, this is equivalent to:
     *
     *    LUX = PI
     * 
     * Which is now in the same form as that solved by the approach
     * implemented in the operator* method. So we use the same 
     * method to here solve for X.
     */ 
    void solve() {
      arma::Mat<elem_type> Y  = arma::solve(A.L, A.P * arma::eye(A.L.n_rows, A.L.n_rows));
      arma::Mat<elem_type> Ai = arma::solve(A.U, Y);
      arma::Mat<elem_type>::operator=(Ai.get_ref());
    }

  public:
    using armawrap::AWBase<elem_type, 
                           AWInvertedCroutMatrix<elem_type>, 
                           arma::Mat<elem_type> >::operator();

    using armawrap::AWBase<elem_type, 
                           AWInvertedCroutMatrix<elem_type>, 
                           arma::Mat<elem_type> >::t;


    inline          AWInvertedCroutMatrix(const AWCroutMatrix<elem_type> &A) : A(A) { }
    virtual inline ~AWInvertedCroutMatrix() { }

    
    inline AWInvertedCroutMatrix operator=(const AWCroutMatrix<elem_type> &A) {
      arma::Mat<elem_type> Y  = arma::solve(A.L, A.P * arma::eye(A.L.n_rows, A.L.n_rows));
      arma::Mat<elem_type> Ai = arma::solve(A.U, Y);

      arma::Mat<elem_type>::operator=(Ai.get_ref());
      return *this;
    } 
    

    /*
     * When the inverse of an AWCroutMatrix is multiplied by another matrix,
     * we treat it as a system of linear equations to be solved:
     *
     *   X = A.i()*B
     *
     *   AX = B
     *
     * Given the LUP decomposition of A,  this is equivalent to:
     *   
     *   LUX = PB
     *
     * Which can be solved in two stages. First, we solve the following for Y:
     *
     *   LY = PB
     *
     * Then we solve for X:
     *
     *   UX = Y
     * 
     */
    template<typename T> inline arma::Mat<elem_type> operator*(const T &B) {

      arma::Mat<elem_type> Y = arma::solve(A.L, A.P * B);
      arma::Mat<elem_type> X = arma::solve(A.U, Y);

      return X;
    }
  };

  /*
   * An AWCroutMatrix encapsulates the LUP decomposition of a matrix A.
   * It is convenient for a situation where a system of linear equations
   * needs to be solved multiple times using A.
   */
  template<typename elem_type> class AWCroutMatrix {

  private:

    bool                             singular;
    arma::Mat<            elem_type> L;
    arma::Mat<            elem_type> U;
    arma::Mat<            elem_type> P;
    AWInvertedCroutMatrix<elem_type> I;

    friend class AWInvertedCroutMatrix<elem_type>;
    
  public:
    
    inline virtual ~AWCroutMatrix() { }

    template<typename T>
    inline AWCroutMatrix(const arma::Base<elem_type, T> &A) : I(*this) {
      singular = arma::lu(L, U, P, A);
      I.solve();
    }

    inline bool IsSingular() {
      return singular;
    }

    inline const AWInvertedCroutMatrix<elem_type> & i() const {
      return I;
    }
  };
}

#endif /* __MATRIX_CROUT_HPP__ */
