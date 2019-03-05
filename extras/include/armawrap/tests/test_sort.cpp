#include "tests.hpp"

#include <typeinfo>

template<typename T>
void doSort(const char *s, T &X) {

  cout << endl << "Sort test (" << typeid(X).name() << ")" << endl << endl;
  printmat(s, X);
  cout << X.SubMatrix(1, X.Nrows(), 1, X.Ncols()) << endl;

  T Xo = X;

  cout << "SortAscending("  << s << ")" << endl;
  SortAscending(Xo);
  printmat(" -> ", Xo);
  cout << Xo.SubMatrix(1, X.Nrows(), 1, X.Ncols()) << endl;

  Xo = X;

  cout << "SortDescending(" << s << ")" << endl;
  SortDescending(Xo);
  printmat(" -> ", Xo);
  cout << Xo.SubMatrix(1, X.Nrows(), 1, X.Ncols()) << endl; 
}

int main (int argc, char *argv[]) {

  Matrix                 m(5, 6);
  DiagonalMatrix         d(5);
  ColumnVector           c(5);
  RowVector              r(5);

  randu(m);
  randu(d);
  randu(c);
  randu(r);

  doSort(" m: ", m);
  doSort(" d: ", d);
  doSort(" c: ", c);
  doSort(" r: ", r);

  return 0;
}
