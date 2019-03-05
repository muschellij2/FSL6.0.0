#include "tests.hpp"

int main (int argc, char *argv[]) {
  Matrix       a(SIZE, SIZE);
  ColumnVector x(SIZE);
  RowVector    r(SIZE);

  Matrix ra(SIZE, SIZE*2);

  randu(a);
  randu(r);
  randu(ra);

  x = 23;

  Matrix cpa = a;
  Matrix ax;
  ax = a + 4.0;
  ax = 4  + ax;
  ax = ax / 3.0;
  ax = ax - 423;
  ax = ax - a;
  ax = ax - a + 8;
  ax = 423.0 - ax;
  ax = ax * 431.0;
  ax = 431.0 * ax;
  ax = -ax;

  ax = ax * x;

  Matrix axm = ax.AsDiagonal();

  axm = (a + axm) + a;
  axm = a + (a + axm);
  axm = (a * axm) * a;
  axm = a * (a * axm);
  axm = (a - axm) - a;
  axm = a - (a - axm);

  axm = (a + axm) + 100.0;
 

  Matrix c  = (a * axm) * x / 100.0 - (x * 4.0);
  Matrix cc = c | c | c;
  Matrix cv = c & c & c;


  cc |= c;
  cv &= c;

  cv += 3214;

  Matrix co(cv.Nrows(), cv.Ncols());
  randu(co);

  cv -= co;

  printmat("a: ", a);
  printmat("x: ", x);
  printmat("cpa: ", cpa);
  printmat("ax: ", ax);
  printmat("axm: ", axm);
  printmat("c: ", c);
  printmat("cc: ", cc);
  printmat("cv: ", cv);
  printmat("ra: ", ra);

  cout << "(a - axm).Maximum():   " << (a - axm).Maximum() << endl;
  cout << "(a / axm).Maximum():   " << (a / 888).Maximum() << endl;
  cout << "(a * axm).Maximum(): " << (a * axm).Maximum() << endl;
  cout << "(a + axm).Nrows():   " << (a + axm).Minimum()   << endl;
  cout << "(a & axm).Minimum():   " << (a & axm).Minimum() << endl;
  cout << "(a | axm).Sum()    :   " << (a | axm).Sum()     << endl;


  return 0;
}
