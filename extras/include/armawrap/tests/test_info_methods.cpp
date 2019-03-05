#include "tests.hpp"


template<typename T>
void testthings(T &X) {
  cout << "X.Nrows():" << X.Nrows()  << endl;
  cout << "X.Ncols():" << X.Ncols()  << endl;

  int i, j;
  cout << "X.Minimum:               " << X.Minimum()                   << endl;
  cout << "Minimum(X):              " << Minimum(X)                    << endl;
  cout << "X.Maximum:               " << X.Maximum()                   << endl;
  cout << "Maximum(X):              " << Maximum(X)                    << endl;
  cout << "X.MinimumAbsoluteValue:  " << X.MinimumAbsoluteValue()      << endl;
  cout << "MinimumAbsoluteValue(X): " << MinimumAbsoluteValue(X)       << endl;
  cout << "X.MaximumAbsoluteValue:  " << X.MaximumAbsoluteValue()      << endl;
  cout << "MaximumAbsoluteValue(X): " << MaximumAbsoluteValue(X)       << endl;
  cout << "X.Minimum1:              " << X.Minimum1(i)                 << endl;
  cout << "     index:              " << i                             << endl;
  cout << "X.Minimum2:              " << X.Minimum2(i, j)              << endl;
  cout << "     index:              " << i << "," << j                 << endl;
  cout << "X.Maximum1:              " << X.Maximum1(i)                 << endl;
  cout << "     index:              " << i                             << endl;
  cout << "X.Maximum2:              " << X.Maximum2(i, j)              << endl;
  cout << "     index:              " << i << "," << j                 << endl;
  cout << "X.MinimumAbsoluteValue1: " << X.MinimumAbsoluteValue1(i)    << endl;
  cout << "     index:              " << i                             << endl;
  cout << "X.MinimumAbsoluteValue2: " << X.MinimumAbsoluteValue2(i, j) << endl;
  cout << "     index:              " << i << "," << j                 << endl;
  cout << "X.MaximumAbsoluteValue1: " << X.MaximumAbsoluteValue1(i)    << endl;
  cout << "     index:              " << i                             << endl;
  cout << "X.MaximumAbsoluteValue2: " << X.MaximumAbsoluteValue2(i, j) << endl;
  cout << "     index:              " << i << "," << j                 << endl;
}

int main (int argc, char *argv[]) {
  Matrix       a(SIZE-5, SIZE+5);
  RowVector    r(SIZE);
  ColumnVector c(SIZE);


  randu(a);
  randu(r);
  randu(c);

  cout << "Testing matrix methods ..."     << endl << endl;
  testthings<Matrix>(a);
  cout << "Testing row vector methods ..." << endl << endl;
  testthings<RowVector>(r);
  cout << "Testing col vector methods ..." << endl << endl;
  testthings<ColumnVector>(c);
}
