#include "tests.hpp"

int main(int argc, char *argv[]) {

  Matrix m(8, 8);
  randu(m);
  
  Matrix b = m.SubMatrix(3, 5, 4, 7);
  Matrix c = m.SubMatrix(3, 5, 4, 7).Reverse();

  Matrix d = m.SubMatrix(5, 4, 1, 8);

  m.Rows(8, 7) = d;

  ColumnVector v(0);

  m.SubMatrix(6, 5, 4, 4) += v;

  printmat("m: ", m);
  printmat("b: ", b);
  printmat("c: ", c);

  cout << "m: " << endl;
  cout <<  m    << endl;

  cout << "b: " << endl;
  cout <<  b    << endl;

  cout << "c: " << endl;
  cout <<  c    << endl;
  
  return 0;
}
