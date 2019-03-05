#include <stdio.h>

#include "tests.hpp"

void mod(double &val) {
  val += 99999;
}

void modptr(double *val) {
  *val *= 100;
}

int main(int argc, char *argv[]) {

  Matrix a(6, 6);

  DiagonalMatrix d(4);

  randu(a);
  randu(d);

  cout << "a(3, 4): " << a(3, 4) << endl;
  cout << "d(3, 3): " << d(3, 3) << endl;

  a(3, 4) = 12345.0;
  d(3, 3) = 54321.0;
  a(3, 3) += 1234;

  d(4, 4) += 1234;

  a(4, 2) = d(2, 2);
  a(4, 5) = a(1, 3);
  d(2, 2) = a(5, 3);

  mod(a(2, 6));
  mod(d(3, 3));

  modptr(&a(2, 6));
  modptr(&d(4, 4));

  cout << "a: " << endl;
  cout <<  a    << endl;
  cout << "d: " << endl;
  cout <<  d    << endl;


  // cast is necessary
  printf("Printf: %e\n", (double)a(3, 5));

  return 0;
}
