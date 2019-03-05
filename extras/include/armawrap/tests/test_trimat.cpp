#include "tests.hpp"

int main(int argc, char *argv[]) {

  UpperTriangularMatrix a(5);
  LowerTriangularMatrix b(5);

  a = 15;
  b = 25;

  printmat("a: ", a);
  printmat("b: ", b);

  cout << a << endl;
  cout << b << endl;

  cout << "a(3, 3): " << a(3, 3) << endl;
  cout << "a(3, 4): " << a(3, 4) << endl;
  cout << "b(3, 2): " << b(3, 2) << endl;
  cout << "b(3, 3): " << b(3, 3) << endl;

}
