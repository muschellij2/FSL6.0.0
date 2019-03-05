#include "tests.hpp"

int main (int argc, char *argv[]) {
  Matrix       a(SIZE, SIZE);
  ColumnVector x(SIZE);
  RowVector    r(SIZE);

  randu(a);
  randu(r);

  x = 23;

  Matrix cpa = a;
  Matrix ax;
  ax = a + 4.0;
  ax = 4  + ax;
  ax = ax / 3.0;
  ax = ax - 423;
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
 
  Matrix c = (a * axm) * x / 100.0 - (x * 4.0);

  ColumnVector m(SIZE);
  randu(m);

  cout << "axm != a: " << (axm != a) << endl;
  cout << "cpa == a: " << (cpa == a) << endl;
  cout << "m == m:   " << (m == m)   << endl;
  cout << "m != m:   " << (m != m)   << endl;

  return 0;
}
