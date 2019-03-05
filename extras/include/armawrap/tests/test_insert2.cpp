#include "tests.hpp"


template <typename T> void insert(T type) {

  double          x;
  double          y;
  double          z;
  double          *u;
  ColumnVector    v( 4);
  ColumnVector    v2(4);
  Matrix          m( 5, 5);
  SymmetricMatrix s( 5);
  DiagonalMatrix  d( 5);

  randu(m);
  randu(v2);

  x = 1;
  y = 2;
  z = 3;

  v << x << y << z << 1.0;
  cout << v << endl;

  v << 0 << 1 << 2 << 3;
  cout << v << endl;

  v.Row(3) << m(3, 4);
  cout << v << endl;

  v << m.Column(4).Rows(1, 4);
  cout << v << endl;

  v << m(1, 1) << m(2, 2) << m(3, 3) << m(4, 4);
  cout << v << endl;

  v << 1 << m(2, 2) << m(3, 3) << m(4, 4);
  cout << v << endl;

  v << 1.0 << m(2, 2) << m(3, 3) << m(4, 4);
  cout << v << endl;

  v << m(1, 1) << m(2, 2) << m(3, 3) << 9;
  cout << v << endl;

  v << m(1, 1) << m(2, 2) << m(3, 3) << 9.0;
  cout << v << endl;

  v << v2(1) << v2(2) << v2(3) << v2(4);
  cout << v << endl;

  s <<  1 <<  2 <<  3 << 4 << 5
    <<  6 <<  7 <<  8 << 9
    << 10 << 11 << 12
    << 13 << 14
    << 15;
  cout << s << endl;

  d << 10 << 11 << 12 << 13 << 14;
  cout << d << endl;

  u = (double *)malloc(4 * sizeof(Real));

  for (int i = 0; i < 4; i++)
      u[i] = i;

  v << u;

  free(u);

  cout << v << endl;
}


int main(int argc, char *argv[]) {

  char type;

  type = 'a';

  insert(type);

  return 0;
}
