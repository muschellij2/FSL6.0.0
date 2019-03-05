#ifndef __TESTS_HPP__
#define __TESTS_HPP__

#ifndef SIZE
#define SIZE 250
#endif

#ifndef IDX
#define IDX 200
#endif

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdint.h>

#ifdef  USING_NEWMAT_API
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
using namespace NEWMAT;


/*
 * consistent cross-platform RNG
 *
 * https://en.wikipedia.org/wiki/Linear_congruential_generator
 */
uint64_t rand_seed = 1234567890;
uint64_t awrand() {
  rand_seed = (1103515245 * rand_seed + 12345) % 2147483648;
  return rand_seed;
}


inline void randu(Matrix &X) {
  for (int ci = 1; ci <= X.Ncols(); ci++) {
    for (int ri = 1; ri <= X.Nrows(); ri++) {
      X(ri, ci) = (float)awrand() / RAND_MAX;
    }
  }
}

inline void randu(RowVector &X) {
  for (int ci = 1; ci <= X.Ncols(); ci++) {
    X(ci) = (float)awrand() / RAND_MAX;
  }
}

inline void randu(ColumnVector &X) {
  for (int ri = 1; ri <= X.Nrows(); ri++) {
    X(ri) = (float)awrand() / RAND_MAX;
  }
}

inline void randu(DiagonalMatrix &X) {
  for (int ri = 1; ri <= X.Nrows(); ri++) {
    X(ri) = (float)awrand() / RAND_MAX;
  }
}

inline void randu(IdentityMatrix &X) {
  X = (float)awrand() / RAND_MAX;
}

inline void randu(UpperTriangularMatrix &X) {
  for (int ri = 1; ri <= X.Nrows(); ri++) {
    for (int ci = ri; ci <= X.Nrows(); ci++) {
      X(ri, ci) = (float)awrand() / RAND_MAX;
    }
  }
}

inline void randu(LowerTriangularMatrix &X) {
  for (int ri = 1; ri <= X.Nrows(); ri++) {
    for (int ci = 1; ci <= ri; ci++) {
      X(ri, ci) = (float)awrand() / RAND_MAX;
    }
  }
}

inline void randu(SymmetricMatrix &X) {
  for (int ri = 1; ri <= X.Nrows(); ri++) {
    for (int ci = 1; ci <= ri; ci++) {
      X(ri, ci) = (float)awrand() / RAND_MAX;
    }
  }
}

inline void randu(BandMatrix &X) {

  MatrixBandWidth bw = X.BandWidth();

  for (  int ri = 1; ri <= X.Nrows(); ri++) {
    for (int ci = 1; ci <= X.Nrows(); ci++) {

      int diag = ci - ri;

      if ((diag < 0) && (-diag > bw.lower)) continue;
      if ((diag > 0) && ( diag > bw.upper)) continue;

      X(ri, ci) = (float)awrand() / RAND_MAX;
    }
  }
}

inline void randu(SymmetricBandMatrix &X) {

  MatrixBandWidth bw = X.BandWidth();

  for (  int ri = 1; ri <= X.Nrows(); ri++) {
    for (int ci = 1; ci <= X.Nrows(); ci++) {

      if (ri < ci)  continue;

      int diag = ci - ri;

      if ((diag < 0) && (-diag > bw.lower)) continue;
      if ((diag > 0) && ( diag > bw.upper)) continue;

      X(ri, ci) = (float)awrand() / RAND_MAX;
    }
  }
}


inline void abs(Matrix &X) {
  for (int ri = 1; ri <= X.Nrows(); ri++) {
  for (int ci = 1; ci <= X.Ncols(); ci++) {
    X(ri, ci) = fabs(X(ri, ci));
  }}
}

inline void abs(UpperTriangularMatrix &X) {
  for (int ri = 1; ri <= X.Nrows(); ri++) {
  for (int ci = ri; ci <= X.Ncols(); ci++) {
    X(ri, ci) = fabs(X(ri, ci));
  }}
}



template<typename T>
inline void printmat(const char *name, const T &a) {
  std::cout << name << SumSquare(a) << ", (" << a.Nrows() << ", " << a.Ncols() << " [" << a.Storage() << "])" << std::endl;
}

#else
#include "armadillo"
template<typename T>
inline void randu(T &a) {
  for (int i = 0; i < a.n_elem; i++) {
    a(i) = awrand();
  }
}

template<typename T>
inline void printmat(const char *name, const T &a) {
  std::cout << name << arma::accu(a) << ", (" << a.n_rows << ", " << a.n_cols << " [" << a.n_elem << "])" << std::endl;
}
#endif /* USING_NEWMAT_API */

#endif /* __TESTS_HPP__ */
