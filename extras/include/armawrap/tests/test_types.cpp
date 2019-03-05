#include "tests.hpp"

#if LIB == NEWMAT_API
void printdims(const char *name, GeneralMatrix  &b) {
  cout << name << b.Nrows() << " x " << b.Ncols() << ": " << b.Storage() << endl;
}
#elif LIB == ARMAWRAP_API
template<typename eT, typename wT, typename aT>
void printdims(const char *name, armawrap::AWBase<eT, wT, aT> &b) {
  cout << name << b.Nrows() << " x " << b.Ncols() << ": " << b.Storage() << endl;
}
#endif


Matrix doThing(Matrix &m) {

  m = 1345;

  return m.t() * 34;
}

int main(int argv, char *argc[]) {

  Matrix a = Matrix(SIZE, SIZE);
  IdentityMatrix i(SIZE);

  i = 90;

  printmat("a: ", a);
  printmat("i: ", i);


  Matrix                 m(10, 10);
  RowVector             rv(10);
  ColumnVector          cv(10);
  IdentityMatrix        im(10);
  DiagonalMatrix         d(10);
  UpperTriangularMatrix ut(10);
  LowerTriangularMatrix lt(10);
  BandMatrix             b(5, 3, 2);
  UpperBandMatrix       ub(5, 3);
  LowerBandMatrix       lb(5, 2);

  Matrix mi = IdentityMatrix(10);

  mi(4,4) = 99;
  mi(2,3) = 123;

  printdims("m:  ",  m);
  printdims("rv: ", rv);
  printdims("cv: ", cv);
  printdims("i:  ", im);
  printdims("d:  ",  d);
  printdims("ut: ", ut);
  printdims("lt: ", lt);
  printdims("b:  ",  b);
  printdims("ub: ", ub);
  printdims("lb  ", lb);
  
  cout << "mi: " << endl;
  cout <<  mi    << endl;


  //test that a RowVector can be passed to
  // a function which accepts a Matrix refr
  Matrix mr = doThing(rv);

  cout << "mr: " << endl;
  cout << mr     << endl;

  return 0;
}

