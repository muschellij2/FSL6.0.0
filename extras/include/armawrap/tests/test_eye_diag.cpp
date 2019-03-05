#include "tests.hpp"

int main(int argc, char *argv[]) {

  IdentityMatrix a(10);
  IdentityMatrix b(5);
  IdentityMatrix c(10);
  DiagonalMatrix d(10);
  DiagonalMatrix e(5);
  DiagonalMatrix f(10);
  IdentityMatrix  g(10);
  DiagonalMatrix h(10);

  a = 1;
  b = 2;
  c = 3;
  d = 4;
  e = 5;
  f = 6;

  a *= 2;
  b *= 3;
  c *= 4;
  d *= 5;
  e *= 6;
  f *= 7;

  d(5) = 25;
  e(5) = 25;
  f(5) = 25;

  g = SP(a, c);
  h = SP(d, f);

  printmat("a: ", a);
  printmat("b: ", b);
  printmat("c: ", c);
  printmat("d: ", d);
  printmat("e: ", e);
  printmat("f: ", f);
  printmat("g: ", g);
  printmat("h: ", h);

  cout << a << endl;
  cout << b << endl;
  cout << c << endl;
  cout << d << endl;
  cout << e << endl;
  cout << f << endl;
  cout << g << endl;
  cout << h << endl;

  cout << d(5) << endl;
  cout << e(5) << endl;
  cout << f(5) << endl;
}
