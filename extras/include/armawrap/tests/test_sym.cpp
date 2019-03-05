#include "tests.hpp"


int main(int argc, char *argv[]) {

  SymmetricMatrix m(3);

  m(1, 1) = 10;
  m(1, 2) = 11;
  m(1, 3) = 12;
  m(2, 1) = 13;
  m(2, 2) = 14;
  m(2, 3) = 15;
  m(3, 1) = 16;
  m(3, 2) = 17;
  m(3, 3) = 18;
  cout << m << endl;

  m << 20 << 21 << 22 << 23 << 24 << 25;
  cout << m << endl;

  m.Row(1) << 1;
  m.Row(2) << 2 << 3;
  m.Row(3) << 4 << 5 << 6;
  cout << m << endl;
}
